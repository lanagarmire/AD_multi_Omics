#!/usr/bin/env python3
"""Reproduce visit-level dementia endpoints and within-stage survival comparisons.

Final 90-month analysis. Inputs are resolved relative to this script.
"""
import argparse
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import lifelines
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test, logrank_test, proportional_hazard_test

GROUPS = ['EMCI1', 'EMCI2', 'LMCI1', 'LMCI2']
COLORS = dict(zip(GROUPS, ['red', '#FFA600', '#0000FF', 'black']))
HERE = Path(__file__).resolve().parent
DEFAULT_PROJECT = HERE / 'inputs'
DEFAULT_VISITS = HERE / 'inputs/ADNIMERGE.csv'

def display(value):
    print(value.to_string())

def participant_outcome(records, expected_baseline, horizon=96.0):
    """Return an endpoint or an explicit exclusion reason; no diagnosis imputation."""
    if records.empty:
        return {'exclusion': 'no_visit_records'}
    dated = records.dropna(subset=['EXAMDATE']).sort_values('EXAMDATE')
    baseline = dated[dated['VISCODE'].eq('bl')].drop_duplicates(
        subset=['EXAMDATE', 'DX', 'DX_bl'])
    if len(baseline) != 1:
        return {'exclusion': 'missing_or_ambiguous_baseline'}
    b = baseline.iloc[0]
    if pd.notna(b['DX']) and b['DX'] == 'Dementia':
        return {'exclusion': 'baseline_dementia'}
    if pd.isna(b['DX']) or b['DX'] != 'MCI':
        return {'exclusion': 'baseline_not_observed_MCI'}
    if pd.isna(b['DX_bl']) or b['DX_bl'] != expected_baseline:
        return {'exclusion': 'baseline_group_mismatch'}
    valid = dated[dated['DX'].isin(['CN', 'MCI', 'Dementia'])].copy()
    if (valid.loc[valid['EXAMDATE'] < b['EXAMDATE'], 'DX'] == 'Dementia').any():
        return {'exclusion': 'dementia_before_baseline'}
    valid = valid[valid['EXAMDATE'] >= b['EXAMDATE']]
    if (valid.groupby('EXAMDATE')['DX'].nunique() > 1).any():
        return {'exclusion': 'conflicting_same_date_diagnoses'}
    valid = valid.drop_duplicates(subset=['EXAMDATE', 'DX'])
    valid['month'] = (valid['EXAMDATE'] - b['EXAMDATE']).dt.days / (365.25 / 12)
    if not (valid['month'] > 0).any():
        return {'exclusion': 'no_postbaseline_diagnosis'}
    dementia = valid.loc[valid['DX'].eq('Dementia'), 'month']
    first_event = dementia.min() if len(dementia) else np.nan
    last_month = valid['month'].max()
    event = int(pd.notna(first_event) and first_event <= horizon)
    duration = float(first_event if event else min(last_month, horizon))
    return {'exclusion': '', 'duration_months': duration, 'event': event,
            'baseline_date': b['EXAMDATE'], 'first_dementia_month': first_event,
            'last_observed_month': last_month}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, default=DEFAULT_PROJECT)
    parser.add_argument('--visits', type=Path, default=DEFAULT_VISITS)
    parser.add_argument('--output', type=Path, default=HERE)
    args = parser.parse_args()
    SOURCE_PROJECT = args.project
    PROJECT = args.output
    ADNIMERGE_CSV = args.visits
    MAX_MONTHS = 90.0
    WINDOW = f'0_{MAX_MONTHS:g}'
    TICKS = np.unique(np.append(np.arange(0, MAX_MONTHS + 0.01, 12), MAX_MONTHS))
    FIGS = PROJECT / 'figs'
    FIGS.mkdir(parents=True, exist_ok=True)
    (PROJECT / 'extdata').mkdir(parents=True, exist_ok=True)
    visits = pd.read_csv(ADNIMERGE_CSV, low_memory=False)
    visits = visits.rename(columns={'DX.bl': 'DX_bl'})
    required = {'PTID', 'EXAMDATE', 'VISCODE', 'DX_bl', 'DX'}
    if not required.issubset(visits.columns):
        raise ValueError(f'Missing columns: {required - set(visits.columns)}')
    for column in ['PTID', 'VISCODE', 'DX_bl', 'DX']:
        visits[column] = visits[column].astype('string').str.strip()
    visits['EXAMDATE'] = pd.to_datetime(visits['EXAMDATE'], errors='coerce')

    members = []
    for group in GROUPS:
        names = pd.read_csv(SOURCE_PROJECT / 'extdata' / f'{group[:-1]}_{group[-1]}_names.csv', dtype=str)
        if 'x' not in names:
            raise ValueError(f'Missing participant column x in {group} names')
        ids = names['x'].dropna().str.strip()
        ids = ids[ids.ne('')].drop_duplicates()
        members.extend({'PTID': ptid, 'group': group} for ptid in ids)
    membership = pd.DataFrame(members)
    if membership['PTID'].duplicated().any():
        raise ValueError('Participant belongs to multiple groups; resolve membership before analysis.')
    visits = visits[visits['PTID'].isin(membership['PTID'])].copy()
    print('Visits with invalid/missing dates:', int(visits['EXAMDATE'].isna().sum()))
    print('Observed diagnosis labels:', sorted(visits['DX'].dropna().unique()))
    unknown = set(visits['DX'].dropna()) - {'CN', 'MCI', 'Dementia', ''}
    if unknown:
        raise ValueError(f'Unrecognized diagnosis labels; explicitly map before analysis: {unknown}')


    records_by_id = dict(tuple(visits.groupby('PTID')))
    audit = pd.DataFrame([{**member, **participant_outcome(
        records_by_id.get(member['PTID'], visits.iloc[:0]), member['group'][:-1], MAX_MONTHS)}
        for member in members])
    cohort = audit[audit['exclusion'].eq('')].copy()
    if cohort.empty:
        raise ValueError('No eligible participants; inspect audit and source diagnosis coding.')
    cohort['event'] = cohort['event'].astype(int)
    assert cohort['PTID'].is_unique
    assert cohort['duration_months'].between(0, MAX_MONTHS, inclusive='both').all()
    summary = cohort.groupby('group').agg(n=('PTID', 'size'), events=('event', 'sum'))
    summary['censored'] = summary['n'] - summary['events']
    display(summary.reindex(GROUPS))
    display(audit.groupby(['group', 'exclusion']).size().rename('n').to_frame())
    audit.to_csv(PROJECT / f'extdata/KM_dementia_{WINDOW}_participant_audit.csv', index=False)
    cohort.to_csv(PROJECT / f'extdata/KM_dementia_{WINDOW}_participants.csv', index=False)
    summary.to_csv(PROJECT / f'extdata/KM_dementia_{WINDOW}_group_summary.csv')


    # Endpoint checks: event, reversion, missing DX, and administrative censoring.
    def example(dates, diagnoses):
        return pd.DataFrame({'EXAMDATE': pd.to_datetime(dates), 'DX': diagnoses,
                             'DX_bl': 'EMCI', 'VISCODE': ['bl'] + ['followup']*(len(dates)-1)})
    r = participant_outcome(example(['2020-01-01', '2021-01-01', '2022-01-01'],
                                     ['MCI', 'Dementia', 'MCI']), 'EMCI')
    assert r['event'] == 1 and 11 < r['duration_months'] < 13
    r = participant_outcome(example(['2020-01-01', '2021-01-01', '2022-01-01'],
                                     ['MCI', 'MCI', None]), 'EMCI')
    assert r['event'] == 0 and 11 < r['duration_months'] < 13
    r = participant_outcome(example(['2020-01-01', '2025-01-01'], ['MCI', 'Dementia']), 'EMCI', 48)
    assert r['event'] == 0 and r['duration_months'] == 48
    assert participant_outcome(example(['2020-01-01'], ['MCI']), 'EMCI')['exclusion']
    assert participant_outcome(example(['2020-01-01'], ['Dementia']), 'EMCI')['exclusion'] == 'baseline_dementia'
    print('Endpoint checks passed.')

    r = participant_outcome(example(['2020-01-01', '2025-01-01'], ['MCI', 'Dementia']), 'EMCI', 96)
    assert r['event'] == 1 and 59 < r['duration_months'] < 61
    r = participant_outcome(example(['2020-01-01', '2029-01-01'], ['MCI', 'Dementia']), 'EMCI', 96)
    assert r['event'] == 0 and r['duration_months'] == 96
    print('96-month endpoint checks passed.')


    fitters = {}
    for group in GROUPS:
        data = cohort[cohort['group'].eq(group)]
        if data.empty:
            raise ValueError(f'No eligible participants in {group}')
        fitters[group] = KaplanMeierFitter(alpha=0.05).fit(
            data['duration_months'], event_observed=data['event'], label=group)
        assert (fitters[group].survival_function_.diff().dropna() <= 1e-12).all().all()

    # At risk immediately before each tick: events/censors at the tick are included.
    at_risk = pd.DataFrame({g: [int((cohort.loc[cohort['group'].eq(g), 'duration_months'] >= t).sum())
                               for t in TICKS] for g in GROUPS}, index=TICKS).T
    at_risk.columns.name = 'Month'
    display(at_risk)
    at_risk.to_csv(PROJECT / f'extdata/KM_dementia_{WINDOW}_numbers_at_risk.csv')

    # Effect sizes: subtype 2 relative to subtype 1, separately within each stage.
    cox_rows = []
    for stage in ['EMCI', 'LMCI']:
        data = cohort[cohort['group'].str.startswith(stage)].copy()
        data['subtype2'] = data['group'].eq(stage + '2').astype(int)
        model_data = data[['duration_months', 'event', 'subtype2']]
        model = CoxPHFitter().fit(model_data, 'duration_months', 'event')
        row = model.summary.loc['subtype2']
        ph = proportional_hazard_test(model, model_data, time_transform='rank')
        cox_rows.append({'stage': stage, 'contrast': f'{stage}2 / {stage}1',
            'hazard_ratio': row['exp(coef)'],
            'CI95_lower': row['exp(coef) lower 95%'],
            'CI95_upper': row['exp(coef) upper 95%'],
            'cox_p_unadjusted': row['p'],
            'PH_rank_test_p': ph.summary.loc['subtype2', 'p']})
    cox = pd.DataFrame(cox_rows)
    cox.to_csv(PROJECT / f'extdata/KM_dementia_{WINDOW}_within_stage_cox.csv', index=False)
    print(cox.to_string(index=False, formatters={
        column: '{:.2f}'.format for column in ['hazard_ratio', 'CI95_lower', 'CI95_upper']}))
    # Refit at 48 months on the same eligible cohort for plot annotations.
    annotation_rows = []
    for horizon in sorted({min(48.0, MAX_MONTHS), MAX_MONTHS}):
        for stage in ['EMCI', 'LMCI']:
            data = cohort[cohort['group'].str.startswith(stage)].copy()
            data['event'] = ((data['event'] == 1) &
                             (data['duration_months'] <= horizon)).astype(int)
            data['duration_months'] = data['duration_months'].clip(upper=horizon)
            data['subtype2'] = data['group'].eq(stage + '2').astype(int)
            model = CoxPHFitter().fit(data[['duration_months', 'event', 'subtype2']],
                                      'duration_months', 'event')
            row = model.summary.loc['subtype2']
            annotation_rows.append({'stage': stage, 'horizon_months': horizon,
                'hazard_ratio': row['exp(coef)'],
                'CI95_lower': row['exp(coef) lower 95%'],
                'CI95_upper': row['exp(coef) upper 95%']})
    annotation_table = pd.DataFrame(annotation_rows)
    annotation_table.to_csv(PROJECT / f'extdata/KM_dementia_{WINDOW}_plot_HRs.csv', index=False)
    hr_labels = {}
    for stage in ['EMCI', 'LMCI']:
        labels = [f'{row.horizon_months:g} months: HR={row.hazard_ratio:.2f} '
                  f'(95% CI {row.CI95_lower:.2f}–{row.CI95_upper:.2f})'
                  for row in annotation_table[annotation_table.stage.eq(stage)].itertuples()]
        hr_labels[stage] = f'{stage}2 / {stage}1: ' + '; '.join(labels)

    for cumulative, stem in [(False, f'Figure_4E_KM_dementia_free_survival_{WINDOW}'),
                             (True, f'Figure_4E_KM_dementia_probability_{WINDOW}')]:
        fig, ax = plt.subplots(figsize=(14, 8))
        for group, km in fitters.items():
            plotter = km.plot_cumulative_density if cumulative else km.plot_survival_function
            plotter(ax=ax, color=COLORS[group], linewidth=2, ci_show=True,
                    ci_alpha=0.12, show_censors=True, censor_styles={'marker': '+', 'ms': 6})
        ax.set(xlim=(0, MAX_MONTHS), ylim=(0, 1.02), xticks=TICKS,
               xlabel='Months from baseline', ylabel=('KM-estimated dementia probability (1 − S(t))'
               if cumulative else 'Dementia-free survival probability S(t)'))
        ax.set_title('\n'.join(hr_labels.values()), fontsize=12)
        ax.tick_params(labelsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', title='MCI subgroup')
        table = ax.table(cellText=at_risk.values, rowLabels=GROUPS,
                         colLabels=[f'{t:g}' for t in TICKS], cellLoc='center',
                         bbox=[0, -0.42, 1, 0.27])
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        ax.text(0, -0.12, 'Number at risk immediately before month', transform=ax.transAxes)
        fig.subplots_adjust(bottom=0.34)
        fig.savefig(FIGS / f'{stem}.png', dpi=300, bbox_inches='tight')
        fig.savefig(FIGS / f'{stem}.pdf', bbox_inches='tight')
        plt.close(fig)



    omnibus = multivariate_logrank_test(cohort['duration_months'], cohort['group'], cohort['event'])
    omnibus_table = omnibus.summary.copy()
    omnibus_table['comparison'] = 'All four groups'
    display(omnibus_table)
    omnibus_table.to_csv(PROJECT / f'extdata/KM_dementia_{WINDOW}_global_logrank.csv', index=False)

    comparisons = []
    for a, b in [('EMCI1', 'EMCI2'), ('LMCI1', 'LMCI2')]:
        x = cohort[cohort['group'].eq(a)]
        y = cohort[cohort['group'].eq(b)]
        test = logrank_test(x['duration_months'], y['duration_months'], x['event'], y['event'])
        comparisons.append({'comparison': f'{a} vs {b}', 'chi_square': test.test_statistic,
                            'df': 1, 'p_unadjusted': test.p_value})
    comparisons = pd.DataFrame(comparisons)
    # Holm step-down adjustment for exactly the two declared contrasts.
    raw = comparisons['p_unadjusted'].to_numpy()
    order = np.argsort(raw)
    adjusted = np.empty(len(raw))
    adjusted[order] = np.minimum(1, np.maximum.accumulate(raw[order] * np.arange(len(raw), 0, -1)))
    comparisons['p_Holm'] = adjusted
    display(comparisons)
    comparisons.to_csv(PROJECT / f'extdata/KM_dementia_{WINDOW}_within_stage_logrank.csv', index=False)

    report = summary.reindex(GROUPS).copy()
    for group in GROUPS:
        data = cohort[cohort['group'].eq(group)]
        report.loc[group, 'censored_before_horizon'] = ((data['event'] == 0) & (data['duration_months'] < MAX_MONTHS)).sum()
        report.loc[group, 'administratively_censored'] = ((data['event'] == 0) & (data['duration_months'] == MAX_MONTHS)).sum()
        reverse = KaplanMeierFitter().fit(data['duration_months'], 1-data['event'])
        report.loc[group, 'reverse_KM_median_followup_months'] = reverse.median_survival_time_
        report.loc[group, 'at_risk_at_horizon'] = (data['duration_months'] >= MAX_MONTHS).sum()
        report.loc[group, 'fraction_at_risk_at_horizon'] = (data['duration_months'] >= MAX_MONTHS).mean()
    display(report)
    report.to_csv(PROJECT / f'extdata/KM_dementia_{WINDOW}_followup_summary.csv')

    estimates = []
    for group, km in fitters.items():
        data = cohort[cohort['group'].eq(group)]
        for month in [48.0, MAX_MONTHS]:
            risk = int((data['duration_months'] >= month).sum())
            supported = risk > 0
            survival = float(km.predict(month)) if supported else np.nan
            ci = km.confidence_interval_survival_function_.loc[:month].iloc[-1]
            estimates.append({'group': group, 'month': month, 'n_at_risk': risk,
                'survival': survival, 'KM_dementia_probability': 1-survival,
                'dementia_probability_CI95_lower': 1-float(ci.iloc[1]) if supported else np.nan,
                'dementia_probability_CI95_upper': 1-float(ci.iloc[0]) if supported else np.nan})
    estimates = pd.DataFrame(estimates).drop_duplicates(['group', 'month'])
    display(estimates)
    estimates.to_csv(PROJECT / f'extdata/KM_dementia_{WINDOW}_estimates.csv', index=False)



    # Dedicated pairwise plots, annotated with the two prespecified log-rank tests.
    for stage in ['EMCI', 'LMCI']:
        groups = [stage + '1', stage + '2']
        result = comparisons[comparisons['comparison'].eq(f'{groups[0]} vs {groups[1]}')].iloc[0]
        fig, ax = plt.subplots(figsize=(10, 7))
        for group in groups:
            fitters[group].plot_cumulative_density(ax=ax, color=COLORS[group], ci_show=True,
                ci_alpha=0.15, show_censors=True, censor_styles={'marker': '+', 'ms': 5})
        ax.set(xlim=(0, MAX_MONTHS), ylim=(0, 1), xticks=TICKS,
            xlabel='Months from baseline', ylabel='KM-estimated dementia probability (1 - S(t))',
            title=(f'{hr_labels[stage].replace("; ", chr(10))}\n'
                   f'Log-rank p={result.p_unadjusted:.3g}; Holm p={result.p_Holm:.3g}'))
        ax.grid(alpha=0.25)
        ax.table(cellText=at_risk.loc[groups].values, rowLabels=groups,
            colLabels=[f'{t:g}' for t in TICKS], cellLoc='center', bbox=[0, -0.32, 1, 0.18])
        ax.text(0, -0.12, 'Number at risk immediately before month', transform=ax.transAxes)
        fig.subplots_adjust(bottom=0.28)
        for suffix in ['png', 'pdf']:
            fig.savefig(FIGS / f'{stage}_within_stage_KM_dementia_probability_{WINDOW}.{suffix}',
                        dpi=300, bbox_inches='tight')
        plt.close(fig)

    metadata = {'python': platform.python_version(), 'lifelines': lifelines.__version__,
        'pandas': pd.__version__, 'numpy': np.__version__, 'horizon_months': MAX_MONTHS,
        'visit_source': str(ADNIMERGE_CSV.resolve()), 'membership_source': str(SOURCE_PROJECT.resolve()),
        'n_eligible': len(cohort), 'n_events': int(cohort.event.sum()), 'n_excluded': len(audit)-len(cohort)}
    (PROJECT / 'run_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')

if __name__ == "__main__":
    main()
