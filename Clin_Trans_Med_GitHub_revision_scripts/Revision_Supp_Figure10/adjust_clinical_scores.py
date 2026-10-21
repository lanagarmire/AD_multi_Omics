"""Reproduce the clinical notebook's residual boxplots, additionally adjusting education.

Final clinical residual model with age, sex and education.
The original complete-case cohort and pooled fitting/scaling are retained.
Run this file directly. Required inputs are documented in README.md.
"""
from pathlib import Path
import itertools
import json
import hashlib
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'inputs/cn_emci_lmci_ad_clin_pseudotime.csv'
OUT = ROOT / 'results'
FIG = ROOT / 'figures'
VARIABLES = [
    'CDRSB', 'ADAS11', 'ADAS13', 'ADASQ4', 'MMSE', 'RAVLT_immediate',
    'RAVLT_learning', 'RAVLT_forgetting', 'RAVLT_perc_forgetting',
    'LDELTOTAL', 'TRABSCOR', 'FAQ',
]
GROUPS = ['CN', 'EMCI1', 'EMCI2', 'LMCI1', 'LMCI2', 'AD']
LABELS = dict(zip(['CN', 'emcisubtype1', 'emcisubtype2',
                  'lmcisubtype1', 'lmcisubtype2', 'AD'], GROUPS))
COLORS = dict(zip(GROUPS, ['#4E79A7', '#F28E2B', '#E15759', '#59A14F', '#B07AA1', '#EDC948']))
STEM = 'Supplementary Figure10_adjustedclinical'


def stars(p):
    return '****' if p <= 1e-4 else '***' if p <= .001 else '**' if p <= .01 else '*' if p <= .05 else 'ns'


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(INPUT)
    assert source['PID'].is_unique
    # Cell 0 calls complete.cases(clinical) across the supplied table.
    # Preserve this filter, including its exclusions for non-modeled columns,
    # so adding education is the only change to the adjustment analysis.
    keep = source.notna().all(axis=1)
    clinical = source.loc[keep].copy()
    assert clinical['SNF'].isin(LABELS).all()
    clinical['Group'] = clinical['SNF'].map(LABELS)
    assert set(clinical.PTGENDER.unique()) <= {0, 1}
    selection = source[['PID', 'SNF']].copy()
    selection['included'] = keep
    selection['missing_columns'] = source.apply(lambda r: ';'.join(r.index[r.isna()]), axis=1)
    selection.to_csv(OUT / 'participant_selection.csv', index=False)
    clinical.to_csv(OUT / 'analysis_cohort.csv', index=False)
    counts = clinical.groupby('Group').size().reindex(GROUPS).rename('n').reset_index()
    counts.to_csv(OUT / 'group_sample_sizes.csv', index=False)

    long, coefficients, audit, baseline = [], [], [], []
    for variable in VARIABLES:
        x = sm.add_constant(clinical[['AGE', 'PTGENDER', 'PTEDUCAT']].astype(float))
        model = sm.OLS(clinical[variable].astype(float), x).fit()
        residual = model.resid
        # R scale(resid(model)) uses the sample SD (n - 1).
        z = (residual - residual.mean()) / residual.std(ddof=1)
        for pid, group, raw, pred, resid, value in zip(
                clinical.PID, clinical.Group, clinical[variable], model.fittedvalues, residual, z):
            long.append(dict(PID=pid, Group=group, Variable=variable,
                             raw_score=raw, fitted_score=pred, Residual=resid, Zadj=value))
        ci = model.conf_int()
        for term in x:
            coefficients.append(dict(Variable=variable, term=term, coefficient=model.params[term],
                                     SE=model.bse[term], CI95_lower=ci.loc[term, 0],
                                     CI95_upper=ci.loc[term, 1], p=model.pvalues[term],
                                     n=len(clinical), R2=model.rsquared))
        old_x = sm.add_constant(clinical[['AGE', 'PTGENDER']].astype(float))
        old_model = sm.OLS(clinical[variable].astype(float), old_x).fit()
        old_z = (old_model.resid - old_model.resid.mean()) / old_model.resid.std(ddof=1)
        for group in GROUPS:
            mask = clinical.Group.eq(group)
            baseline.append(dict(Variable=variable, Group=group, n=int(mask.sum()),
                                 age_sex_mean=old_z[mask].mean(),
                                 age_sex_education_mean=z[mask].mean(),
                                 age_sex_education_median=z[mask].median()))
        covariate_correlation = max(abs(np.corrcoef(z, x[c])[0, 1]) for c in ['AGE', 'PTGENDER', 'PTEDUCAT'])
        audit.append(dict(Variable=variable, n=len(z), z_mean=z.mean(),
                          z_sample_sd=z.std(ddof=1),
                          max_abs_covariate_correlation=covariate_correlation,
                          passed=bool(abs(z.mean()) < 1e-10 and abs(z.std(ddof=1)-1) < 1e-10
                                      and covariate_correlation < 1e-10)))
    long = pd.DataFrame(long)
    long.to_csv(OUT / 'age_sex_education_adjusted_scores.csv', index=False)
    long.pivot(index=['PID', 'Group'], columns='Variable', values='Zadj').to_csv(OUT / 'adjusted_scores_wide.csv')
    pd.DataFrame(coefficients).to_csv(OUT / 'OLS_coefficients.csv', index=False)
    pd.DataFrame(baseline).to_csv(OUT / 'age_sex_vs_age_sex_education_group_summary.csv', index=False)
    pd.DataFrame(audit).to_csv(OUT / 'residual_verification.csv', index=False)
    assert all(row['passed'] for row in audit)

    tests = []
    for variable in VARIABLES:
        sub = long.loc[long.Variable.eq(variable)]
        rows = []
        for g1, g2 in itertools.combinations(GROUPS, 2):
            a = sub.loc[sub.Group.eq(g1), 'Zadj']
            b = sub.loc[sub.Group.eq(g2), 'Zadj']
            # Same large-sample, tie-corrected, continuity-corrected two-sided
            # rank-sum test used by R wilcox.test for these group sizes (>50).
            result = mannwhitneyu(a, b, alternative='two-sided', method='asymptotic', use_continuity=True)
            rows.append(dict(Variable=variable, group1=g1, group2=g2,
                             n1=len(a), n2=len(b), U=result.statistic, p=result.pvalue))
        q = multipletests([r['p'] for r in rows], method='fdr_bh')[1]
        for row, value in zip(rows, q):
            row.update(p_BH=float(value), significance=stars(value))
        tests.extend(rows)
    tests = pd.DataFrame(tests)
    tests.to_csv(OUT / 'pairwise_Wilcoxon_BH.csv', index=False)
    assert len(tests) == 12 * 15

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(2, 6, figsize=(22, 10))
    # ggplot facet_wrap alphabetically orders these character-valued facets.
    for ax, variable in zip(axes.flat, sorted(VARIABLES)):
        sub = long.loc[long.Variable.eq(variable)]
        arrays = [sub.loc[sub.Group.eq(g), 'Zadj'].to_numpy() for g in GROUPS]
        boxes = ax.boxplot(arrays, positions=np.arange(1, 7), widths=.73, patch_artist=True,
                           medianprops={'color': '#333333', 'linewidth': 1.5},
                           boxprops={'edgecolor': '#444444', 'linewidth': 1},
                           whiskerprops={'color': '#444444', 'linewidth': 1},
                           capprops={'color': '#444444', 'linewidth': 1},
                           flierprops={'marker': 'o', 'markersize': 2.5, 'markerfacecolor': '#333333',
                                       'markeredgecolor': '#333333', 'alpha': .7})
        for patch, group in zip(boxes['boxes'], GROUPS):
            patch.set_facecolor(COLORS[group])
        lo, hi = sub.Zadj.min(), sub.Zadj.max()
        for j, row in enumerate(tests.loc[tests.Variable.eq(variable)].itertuples()):
            x1, x2 = GROUPS.index(row.group1)+1, GROUPS.index(row.group2)+1
            y = hi + .5 + j * .4
            ax.plot([x1, x1, x2, x2], [y-.1, y, y, y-.1], color='#333333', linewidth=.7)
            ax.text((x1+x2)/2, y+.015, row.significance, ha='center', va='bottom', fontsize=8)
        ax.set_ylim(lo-.35, hi+.5+14*.4+.55)
        ax.set_xlim(.4, 6.6)
        ax.set_xticks(range(1, 7), GROUPS, rotation=45, ha='right')
        ax.set_title(variable, fontsize=10, pad=8,
                     bbox={'facecolor': '#D9D9D9', 'edgecolor': '#777777', 'boxstyle': 'square,pad=.45'})
        ax.grid(True, color='#e9e9e9', linewidth=.8)
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_color('#777777')
    fig.supylabel('Age-, sex- and education-adjusted z-score residual', x=.012, fontsize=14)
    fig.supxlabel('Subtype', y=.015, fontsize=13)
    fig.legend(handles=[Patch(facecolor=COLORS[g], edgecolor='#444444', label=g) for g in GROUPS],
               title='Group', loc='center right', bbox_to_anchor=(.999, .5), frameon=False)
    fig.subplots_adjust(left=.045, right=.94, top=.95, bottom=.09, hspace=.23, wspace=.25)
    for suffix in ['png', 'pdf', 'svg']:
        fig.savefig(FIG / f'{STEM}.{suffix}', dpi=300, facecolor='white')
    plt.close(fig)

    metadata = dict(
        source_notebook='Revision_Code_to_Github/5_ADNI_clinical_scores.ipynb',
        source_cells=[0, 24, 25, 29], input_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        model='raw_score ~ intercept + AGE + PTGENDER + PTEDUCAT',
        fit_population='All six groups pooled, retaining original notebook complete-case filtering',
        scaling='Residual mean centered and divided by pooled sample SD (ddof=1)',
        test='Two-sided Wilcoxon rank-sum; asymptotic, tie and continuity corrected',
        multiplicity='BH adjustment over 15 pairwise comparisons separately for each of 12 measures',
        n_input=len(source), n_analyzed=len(clinical), n_excluded=int((~keep).sum()),
        n_excluded_for_added_education=int(source.loc[keep, 'PTEDUCAT'].isna().sum()),
        group_counts=clinical.groupby('Group').size().reindex(GROUPS).to_dict(),
        all_residual_checks_pass=True,
        versions=dict(python=platform.python_version(), numpy=np.__version__,
                      pandas=pd.__version__, scipy=scipy.__version__, statsmodels=statsmodels.__version__),
    )
    (OUT / 'analysis_metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    caption = ('**Supplementary Figure 10. Age-, sex- and education-adjusted cognitive and functional '
               'measures across ADNI subtypes.** Boxplots show standardized residuals for 12 measures '
               'across CN, EMCI1, EMCI2, LMCI1, LMCI2 and AD. Each raw score was regressed on age, sex '
               'and years of education in the pooled complete-case cohort, and residuals were centered '
               'and scaled by their sample standard deviation. Boxes show the median and interquartile '
               'range; whiskers extend to the most extreme values within 1.5 times the interquartile '
               'range, and points show outliers. Brackets indicate two-sided pairwise Wilcoxon rank-sum '
               'tests with BH adjustment across 15 comparisons per measure '
               '(*adjusted P ≤0.05, **≤0.01, ***≤0.001, ****≤0.0001; ns, adjusted P >0.05).')
    (ROOT / 'legend.md').write_text(caption+'\n')
    print(json.dumps(metadata, indent=2))
    return metadata


if __name__ == '__main__':
    run()
