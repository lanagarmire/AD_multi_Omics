"""Rebuttal visualization of the reproduced ADNI p180 transformation workflow."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import skew
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

H=Path(__file__).resolve().parent
stages=['imputed','log2','standardized','clipped']
data={mode:{stage:pd.read_csv(H/'inputs'/f'{mode}_{stage}.csv').set_index('RID')
            for stage in stages} for mode in ['fia','uplc']}
rows=[]
for mode,ds in data.items():
    for stage in stages:
        assert ds[stage].columns.equals(ds['imputed'].columns)
        assert ds[stage].index.equals(ds['imputed'].index)
    raw,log,z,clipped=[ds[s].to_numpy() for s in stages]
    assert np.allclose(log,np.log2(raw))
    assert np.allclose(z,(log-log.mean(0))/log.std(0,ddof=1))
    assert np.allclose(clipped,np.clip(z,-3,3))
    for j,name in enumerate(ds['imputed'].columns):
        rows.append(dict(platform=mode.upper(),metabolite=name,n=len(raw),
            skew_prelog=float(skew(raw[:,j],bias=False)),
            skew_log2=float(skew(log[:,j],bias=False)),
            lower_clipped=int((z[:,j]<-3).sum()),upper_clipped=int((z[:,j]>3).sum()),
            clipped_percent=100*float((abs(z[:,j])>3).mean())))
stats=pd.DataFrame(rows)
stats.to_csv(H/'feature_distribution_summary.csv',index=False)
zall=np.concatenate([data[m]['standardized'].to_numpy().ravel() for m in data])
call=np.clip(zall,-3,3)
summary=dict(features=len(stats),features_right_skewed_prelog=int((stats.skew_prelog>0).sum()),
    percent_right_skewed_prelog=100*float((stats.skew_prelog>0).mean()),
    median_skew_prelog=float(stats.skew_prelog.median()),median_skew_log2=float(stats.skew_log2.median()),
    median_absolute_skew_prelog=float(stats.skew_prelog.abs().median()),
    median_absolute_skew_log2=float(stats.skew_log2.abs().median()),
    values=len(zall),values_clipped=int((abs(zall)>3).sum()),
    percent_values_clipped=100*float((abs(zall)>3).mean()),
    n_fia=len(data['fia']['imputed']),n_uplc=len(data['uplc']['imputed']))
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'pdf.fonttype':42})
fig,axes=plt.subplots(2,2,figsize=(11,8))
ax=axes[0,0]
ax.boxplot([stats.skew_prelog,stats.skew_log2],
           showfliers=True,flierprops=dict(markersize=3,alpha=.4))
ax.set_xticks([1,2],['Before log₂','After log₂'])
ax.axhline(0,color='#888888',lw=.8)
ax.set(ylabel='Feature-wise skewness',title='A  Feature-wise skewness')
ax.set_yscale('symlog',linthresh=1)
ax.text(.97,.96,f'Median |skewness|\n{summary["median_absolute_skew_prelog"]:.3f} → {summary["median_absolute_skew_log2"]:.3f}',
        transform=ax.transAxes,ha='right',va='top',fontsize=10)
ax.grid(axis='y',alpha=.15)
representatives={}
for mode in ['fia']:
    s=stats[stats.platform.eq(mode.upper())]
    target=s.skew_prelog.median()
    rep=s.loc[(s.skew_prelog-target).abs().idxmin()]
    name=rep.metabolite
    representatives[mode]=name
    for col,stage in [(0,'imputed'),(1,'log2')]:
        ax=axes[0,1] if col==0 else axes[1,0]
        values=data[mode][stage][name].to_numpy()
        ax.hist(values,bins=40,density=True,color='#518ab5' if col==0 else '#db9450',alpha=.85)
        letter=['C','D'][col]
        label='Before log₂' if col==0 else 'After log₂'
        ax.set(title=f'{letter}  {name} — {label}',
               xlabel='Concentration' if col==0 else 'Log₂ concentration',ylabel='Density')
        sk=rep.skew_prelog if col==0 else rep.skew_log2
        ax.text(.97,.95,f'Skewness = {sk:.2f}',transform=ax.transAxes,ha='right',va='top')
        ax.grid(axis='y',alpha=.15)
ax=axes[1,1]
# Common z-score scale permits pooling standardized entries across metabolites.
bins=np.linspace(-5,5,151)
# Zoom the view, keeping density normalized to all entries, including tails.
weights=np.full(len(zall),1/(len(zall)*(bins[1]-bins[0])))
ax.hist(zall,bins=bins,weights=weights,histtype='step',lw=1.5,color='#518ab5',label='Before clipping')
ax.hist(call,bins=bins,weights=weights,histtype='step',lw=1.5,color='#db9450',label='After clipping')
for cutoff in [-3,3]:ax.axvline(cutoff,linestyle='--',lw=.8,color='#777777')
ax.set(title='B  Standardization and clipping',xlabel='Z score',ylabel='Density',xlim=(-5,5))
ax.legend(fontsize=8,loc='upper left')
ax.grid(axis='y',alpha=.15)
fig.suptitle('ADNI p180 metabolomics preprocessing',fontsize=16)
fig.tight_layout(rect=[0,0,1,.95],h_pad=2,w_pad=2)
for ext in ['png','pdf']:fig.savefig(H/f'metabolomics_preprocessing_rebuttal.{ext}',dpi=300,bbox_inches='tight')
plt.close(fig)
summary['representative_features']=representatives
(H/'distribution_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
