from pathlib import Path
import pickle,json
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from docx import Document
ROOT=Path(__file__).resolve().parent;R=ROOT/'results';F=ROOT/'figures'
def observed(stage):
 return pd.read_csv(R/f'{stage}_observed_validity.csv')
eobs=observed('EMCI');enull=pd.read_csv(R/'EMCI_permutation_replicates.csv');e={'k_range':eobs.k.to_numpy(),'score_sil_list':eobs.silhouette_affinity.to_numpy(),'score_chz_list':eobs.calinski_harabasz.to_numpy(),**{'null_'+short:{int(k):g[metric].to_numpy() for k,g in enull.groupby('k')} for short,metric in [('sil','silhouette_affinity'),('ch','calinski_harabasz')]}}
lobs=observed('LMCI');lnull=pd.read_csv(R/'LMCI_permutation_replicates.csv');b=pd.read_csv(R/'all_modality_baselines.csv')
plt.rcParams.update({'font.size':12,'axes.titlesize':13,'axes.labelsize':12,'pdf.fonttype':42})
fig,axes=plt.subplots(2,3,figsize=(24,13));exports=[]
methods=['SNF + spectral clustering','PCA + k-means','PCA + Ward hierarchical clustering'];colors=['#b6487e','#d39c26','#365b7b']
for row,stage in enumerate(['EMCI','LMCI']):
 ks=np.array(e['k_range'] if stage=='EMCI' else lobs.k)
 for col,(metric,short,obskey,label) in enumerate([('silhouette_affinity','sil','score_sil_list','Affinity silhouette'),('calinski_harabasz','ch','score_chz_list','Calinski–Harabasz index')]):
  observed=np.array(e[obskey] if stage=='EMCI' else lobs[metric]);nulls=[np.array(e['null_'+short][int(k)]) if stage=='EMCI' else lnull[lnull.k==k][metric].to_numpy() for k in ks]
  means=np.array([v.mean() for v in nulls]);lo=np.array([np.quantile(v,.025) for v in nulls]);hi=np.array([np.quantile(v,.975) for v in nulls]);ax=axes[row,col]
  ax.fill_between(ks,lo,hi,color='#d2dfeb',label='Null 95% interval');ax.plot(ks,means,'--',color='#55738e',label='Null mean');ax.plot(ks,observed,'o-',color=colors[0],label='Observed');ax.axvline(2,color='grey',linestyle=':',lw=1)
  k2_index=int(np.flatnonzero(ks==2)[0])
  p=(1+np.sum(nulls[k2_index]>=observed[k2_index]))/(1+len(nulls[k2_index]));ax.set(title=f'{stage}: maximum at k={int(ks[observed.argmax()])}',ylabel=label,xlabel='Number of clusters',xticks=ks)
  p_label='p < 0.001' if p<.001 else f'p = {p:.3f}'
  ax.text(.97,.96,f'k=2: observed = {observed[k2_index]:.3f}\n{p_label}',transform=ax.transAxes,ha='right',va='top',fontsize=11)
  if row==0 and col==0:ax.legend(fontsize=10,loc='lower left')
  exports.extend(dict(stage=stage,k=int(k),metric=metric,observed=float(o),null_mean=float(m),null_025=float(low),null_975=float(high)) for k,o,m,low,high in zip(ks,observed,means,lo,hi))
 ax=axes[row,2]
 order=[(method,modality) for method in methods for modality in ['Fused','Gene expression','Metabolomics']]
 pivot=b[b.stage==stage].pivot(index=['method','modality'],columns='k',values='native_silhouette').loc[order]
 sns.heatmap(pivot,ax=ax,annot=True,fmt='.3f',cmap='YlOrRd',vmin=0,vmax=.8,linewidths=.4,annot_kws={'fontsize':8},cbar_kws={'label':'Native silhouette','shrink':.8,'pad':.03})
 ax.set(title=f'{stage}: native silhouette scores',ylabel='',xlabel='Number of clusters')
 names={'SNF + spectral clustering':'Network + spectral','PCA + k-means':'PCA + k-means','PCA + Ward hierarchical clustering':'PCA + hierarchical (Ward)'}
 ax.set_yticklabels([names[method]+'\n'+modality for method,modality in order],rotation=0,fontsize=9)
 for boundary in [3,6]:ax.axhline(boundary,color='white',lw=3)
 ax.tick_params(axis='y',length=0)
 ax.collections[0].colorbar.ax.tick_params(labelsize=9)
 ax.collections[0].colorbar.ax.set_ylabel('Native silhouette',fontsize=10)
 b[b.stage==stage][['stage','method','modality','k','native_silhouette','native_metric']].to_csv(R/f'{stage}_integrated_DH_native_heatmap.csv',index=False)
 # The integrated silhouette curve and fused-network heatmap must match at every k.
 curve=e['score_sil_list'] if stage=='EMCI' else lobs.silhouette_affinity.to_numpy()
 assert np.allclose(pivot.loc[(methods[0],'Fused'),ks].to_numpy(),curve,rtol=0,atol=1e-12)
 for col in range(3):
  axes[row,col].text(-.18,1.07,(['A','B','C'] if row==0 else ['D','E','F'])[col],transform=axes[row,col].transAxes,fontweight='bold',fontsize=17);
  if col!=2:axes[row,col].grid(alpha=.15)
fig.suptitle('EMCI and LMCI: null comparisons and native-silhouette benchmarks across modalities',fontsize=18);fig.tight_layout(rect=[0,0,1,.96],w_pad=3,h_pad=3)
for ext in ['png','pdf']:fig.savefig(F/f'EMCI_LMCI_integrated_cluster_validity.{ext}',dpi=300,bbox_inches='tight')
plt.close(fig);pd.DataFrame(exports).to_csv(R/'EMCI_LMCI_integrated_null_plot_data.csv',index=False)
(ROOT/'integrated_figure_legend.md').write_text('Observed and 1000 complete-pipeline feature-permutation scores use identical original fusion defaults: affinity K=10, standalone mu=0.5, fusion input mu=0.3, fusion K=20, t=20, alpha=1. All spectral clustering uses seed=42 and n_init=10. Gene inputs are the original measured gene-expression profiles (20,032 features). Independent modality clustering; native affinity and PCA Euclidean silhouettes are different evaluation spaces. Archived downstream labels unchanged.\n')
