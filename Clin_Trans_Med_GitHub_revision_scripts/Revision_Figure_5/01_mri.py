#!/usr/bin/env python3
"""Baseline regional-volume comparisons with OLS and HC3 robust inference."""
from pathlib import Path
from itertools import combinations
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from patsy import build_design_matrices

HERE = Path(__file__).resolve().parent
REVIEW = HERE.parent
SOURCE = HERE / 'inputs/cn_emci_lmci_ad_clin_pseudotime.csv'
VISITS = HERE / 'inputs/ADNIMERGE.csv'
GROUPS = ['CN', 'EMCI1', 'EMCI2', 'LMCI1', 'LMCI2', 'AD']
REGIONS = ['MidTemp_bl', 'Fusiform_bl', 'Entorhinal_bl']
COLORS = ['#4E79A7','#F28E2B','#E15759','#59A14F','#B07AA1','#EDC948']
for directory in ['extdata','figs']: (HERE/directory).mkdir(exist_ok=True)
labels = pd.read_csv(SOURCE)
labels['subtype'] = labels.SNF.replace({'emcisubtype1':'EMCI1','emcisubtype2':'EMCI2','lmcisubtype1':'LMCI1','lmcisubtype2':'LMCI2'})
assert labels.PID.is_unique
visits = pd.read_csv(VISITS, low_memory=False)
baseline = visits.loc[visits.VISCODE.eq('bl')].copy()
assert baseline.PTID.is_unique, 'Resolve duplicate baseline visits before fitting.'
cols = ['PTID','AGE','PTGENDER','ICV_bl','FLDSTRENG_bl','FSVERSION_bl'] + REGIONS
cohort = labels[['PID','subtype']].merge(baseline[cols],left_on='PID',right_on='PTID',how='left',validate='one_to_one')
for col in ['AGE','ICV_bl'] + REGIONS: cohort[col] = pd.to_numeric(cohort[col],errors='coerce')
for col in ['ICV_bl'] + REGIONS: cohort.loc[cohort[col]<=0,col] = np.nan
cohort['subtype'] = pd.Categorical(cohort.subtype,categories=GROUPS,ordered=True)
cohort.to_csv(HERE/'extdata/merged_baseline_cohort.csv',index=False)
scanner = []
for col in ['FLDSTRENG_bl','FSVERSION_bl']:
    tab = pd.crosstab(cohort.subtype,cohort[col].fillna('Missing'))
    tab.to_csv(HERE/f'extdata/{col}_by_subtype.csv')
    scanner.append({'variable':col,'levels':int(cohort[col].nunique()),'missing':int(cohort[col].isna().sum())})
rows, sizes, adjusted = [], [], []
FORMULA = 'C(subtype) + AGE + C(PTGENDER) + ICV_bl + C(FSVERSION_bl)'

def contrasts(model, ref, region, kind):
    grid = pd.DataFrame([dict(ref,subtype=g) for g in GROUPS])
    design = np.asarray(build_design_matrices([model.model.data.design_info],grid)[0])
    results=[]
    for a,b in combinations(range(6),2):
        test = model.t_test(design[b]-design[a])
        ci=np.asarray(test.conf_int()).ravel()
        results.append(dict(region=region,group1=GROUPS[a],group2=GROUPS[b],contrast=f'{GROUPS[b]} - {GROUPS[a]}',
            difference_mm3=float(np.asarray(test.effect).item()),CI95_lower=float(ci[0]),CI95_upper=float(ci[1]),
            p=float(np.asarray(test.pvalue).item()),analysis=kind))
    return results

for region in REGIONS:
    # Preserve the reference age, sex and ICV from the original analysis sample.
    base=cohort.dropna(subset=[region,'subtype','AGE','PTGENDER','ICV_bl']).copy()
    ref={'AGE':float(base.AGE.mean()),'ICV_bl':float(base.ICV_bl.mean()),
         'PTGENDER':base.PTGENDER.mode().iloc[0]}
    data=base.dropna(subset=['FSVERSION_bl']).copy()
    ref['FSVERSION_bl']=data.FSVERSION_bl.mode().iloc[0]
    raw=smf.ols(f'{region} ~ {FORMULA}',data=data).fit()
    assert np.linalg.matrix_rank(raw.model.exog)==raw.model.exog.shape[1]
    assert raw.df_resid>20
    model=raw.get_robustcov_results(cov_type='HC3',use_t=True)
    rows.extend(contrasts(model,ref,region,'age_sex_ICV_FSversion_HC3'))
    sizes.extend(dict(region=region,subtype=g,n=int((data.subtype==g).sum())) for g in GROUPS)
    same=data.copy()
    for column,value in ref.items(): same[column]=value
    data['adjusted_volume_mm3']=data[region]-raw.predict(data)+raw.predict(same)
    data['region']=region
    adjusted.append(data[['PID','subtype','region','adjusted_volume_mm3']])

results=pd.DataFrame(rows)
results['p_Holm_45']=multipletests(results.p,method='holm')[1]
primary=results[((results.group1=='EMCI1')&(results.group2=='EMCI2'))|((results.group1=='LMCI1')&(results.group2=='LMCI2'))].copy()
primary['p_Holm_6']=multipletests(primary.p,method='holm')[1]
results.to_csv(HERE/'extdata/all_pairwise_adjusted_comparisons.csv',index=False)
primary.to_csv(HERE/'extdata/within_stage_adjusted_comparisons.csv',index=False)
pd.DataFrame(sizes).to_csv(HERE/'extdata/sample_sizes.csv',index=False)
adjusted=pd.concat(adjusted,ignore_index=True)
adjusted.to_csv(HERE/'extdata/adjusted_plot_values.csv',index=False)

def stars(p):
    return '****' if p<=.0001 else '***' if p<=.001 else '**' if p<=.01 else '*' if p<=.05 else 'ns'

for all_pairs in [True]:
    fig,axes=plt.subplots(1,3,figsize=(19,7))
    for ax,region in zip(axes,REGIONS):
        d=adjusted[adjusted.region.eq(region)]
        vals=[d.loc[d.subtype.eq(g),'adjusted_volume_mm3'].values for g in GROUPS]
        boxes=ax.boxplot(vals,positions=np.arange(6),patch_artist=True,widths=.7,medianprops={'color':'#333333','linewidth':2})
        for patch,color in zip(boxes['boxes'],COLORS):patch.set_facecolor(color)
        ax.set_xticks(range(6),GROUPS,rotation=45,ha='right');ax.set_axisbelow(True);ax.grid(color='#e8e8e8')
        lo,hi=d.adjusted_volume_mm3.min(),d.adjusted_volume_mm3.max();step=(hi-lo)*.065
        tests=results[results.region.eq(region)] if all_pairs else primary[primary.region.eq(region)]
        tests=tests.sort_values(['group1','group2'])
        for k,r in enumerate(tests.itertuples()):
            a,b=GROUPS.index(r.group1),GROUPS.index(r.group2);y=hi+step*(k+1)
            ax.plot([a,a,b,b],[y-step*.15,y,y,y-step*.15],color='#444444',lw=1)
            p=r.p_Holm_45 if all_pairs else r.p_Holm_6
            ax.text((a+b)/2,y+step*.05,stars(p),ha='center',fontsize=10)
        ax.set_ylim(lo-step,hi+step*(len(tests)+3));ax.set_ylabel('Covariate-adjusted regional volume (mm³)')
        ax.add_patch(Rectangle((0,1),1,.085,transform=ax.transAxes,clip_on=False,facecolor='#d9d9d9',edgecolor='#555555'))
        ax.text(.5,1.04,region,transform=ax.transAxes,ha='center',va='center',fontsize=15)
    fig.text(.005,.96,'A',fontsize=30,fontweight='bold')
    fig.legend(handles=[Patch(facecolor=c,label=g) for c,g in zip(COLORS,GROUPS)],loc='center right',frameon=False,title='Subtype')
    family='45 pairwise comparisons' if all_pairs else '6 within-stage comparisons'
    adjustment='age, sex, baseline ICV and FreeSurfer version'
    fig.text(.5,.02,f'Adjusted for {adjustment}; OLS contrasts with HC3 standard errors. Holm correction across {family}.',ha='center',fontsize=11)
    fig.subplots_adjust(left=.065,right=.91,top=.87,bottom=.18,wspace=.25)
    stem='MRI_adjusted_all_groups' if all_pairs else 'MRI_adjusted_within_stage'
    stem += '_FSversion'
    for ext in ['png','pdf']:fig.savefig(HERE/f'figs/{stem}.{ext}',dpi=300,bbox_inches='tight')
    plt.close(fig)
(HERE/'run_metadata.json').write_text(json.dumps({'membership_source':str(SOURCE),'MRI_source':str(VISITS),'formula':FORMULA,'inference':'HC3, t reference','scanner_fields':scanner,'software':{'pandas':pd.__version__,'statsmodels':sm.__version__}},indent=2))
print(primary.to_string(index=False))
