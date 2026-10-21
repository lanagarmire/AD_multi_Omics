"""Reproduce the user-selected final fit and test fixed-trajectory LOWESS R².

Run with python311 environment. Permutation inference is conditional on the
selected trajectory and settings, not corrected for the prior parameter search.
"""
from pathlib import Path
import importlib.util, json, random
import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('analysis',ROOT/'clinical_trajectory.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)

def main():
    for folder in ['extdata', 'figs']: (ROOT/folder).mkdir(parents=True, exist_ok=True)
    np.random.seed(0); random.seed(0)
    clinical=pd.read_csv(a.INPUT)
    adjusted,cn=a.fit_cn_normative_scores(clinical,include_education=True)
    x,_,tree,root,trajectories=a.fit_clintrajan(adjusted)
    pseudotime=a.collapse_pseudotime(adjusted,trajectories)
    models=a.fit_subtype_models(pseudotime,include_education=True)
    curve=a.group_curve_data(x,adjusted,tree,trajectories)
    observed=float(curve['descriptive_r2'])
    print('REPRODUCED R2',observed,flush=True)
    if round(observed,2)!=0.78: raise ValueError('Selected settings did not reproduce 0.78; do not relabel.')
    for name,data in [('final_adjusted_scores',adjusted),('final_CN_models',cn),('final_pseudotime',pseudotime),('final_subtype_models',models)]:
        data.to_csv(ROOT/'extdata'/f'{name}.csv',index=False)
    order=np.argsort(curve['pseudotime']);xx=curve['pseudotime'][order]; yy=curve['group_code'][order]
    pd.DataFrame({'participant_index':curve['points'][order],'pseudotime':xx,'group_code':yy}).to_csv(ROOT/'extdata/final_curve_points.csv',index=False)
    rng=np.random.default_rng(20261019); B=9999; null=np.empty(B)
    for i in range(B):
        perm=rng.permutation(yy)
        smooth=lowess(perm,xx,frac=.30,it=1,return_sorted=True)
        fitted=np.interp(xx,smooth[:,0],smooth[:,1])
        null[i]=1-np.sum((perm-fitted)**2)/np.sum((perm-perm.mean())**2)
        if (i+1)%1000==0: print('Permutations',i+1,flush=True)
    k=int(np.sum(null>=observed));p=(k+1)/(B+1)
    result=dict(r2=observed,n_participants=len(yy),permutations=B,exceedances=k,p_value=p,seed=20261019,parameters=a.PRIMARY_TREE_PARAMETERS,root=int(root),terminal=int(curve['terminal_node']),test='One-sided permutation of ordinal clinical-group codes on fixed selected trajectory; LOWESS refitted each permutation',limitation='Conditional on the selected settings and trajectory; not corrected for parameter or trajectory selection, or residual covariate effects.')
    (ROOT/'extdata/final_permutation_result.json').write_text(json.dumps(result,indent=2)+'\n')
    pd.DataFrame({'permutation_r2':null}).to_csv(ROOT/'extdata/final_permutation_null.csv',index=False)
    a.plot_requested_panels(x,tree,root,adjusted,trajectories,models,file_tag='final',figure_title='Final selected ClinTrajan trajectory: 35 nodes; alpha=0.01; Mu=0.20; Lambda=0.010; leaf extension=0.25')
    print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__': main()
