from pathlib import Path
import json
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from matplotlib.collections import LineCollection
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
A=P/'inputs/ADNI'
E=P/'inputs/EFIGA'
adni=pd.read_csv(A/'ADNI_umap.csv').merge(pd.read_csv(A/'ADNI_metadata.csv'),on='ID',validate='one_to_one')
efiga=pd.read_csv(E/'SELECTED_full_cohort_UMAP_coordinates.csv')
archived_labels=pd.read_csv(E/'SELECTED_labels.csv')[['ID','label']]
# Restrict EFIGA MCI to the final biomarker-available 149 participants.
biomarker_path=P/'inputs/patients.csv'
biomarker_cohort=pd.read_csv(biomarker_path)
assert biomarker_cohort.ID.is_unique and len(biomarker_cohort)==149
mci_ids=set(biomarker_cohort.ID)
assert biomarker_cohort.label.value_counts().to_dict()=={'EMCI1':51,'LMCI1':44,'EMCI2':34,'LMCI2':20}
comparison=archived_labels.merge(biomarker_cohort[['ID','label']],on='ID',suffixes=('_archived','_final'),validate='one_to_one')
assert len(comparison)==149
comparison['broad_archived']=comparison.label_archived.str[:4]
comparison['broad_final']=comparison.label_final.str[:4]
comparison['subtype_changed']=comparison.label_archived.ne(comparison.label_final)
comparison['broad_group_changed']=comparison.broad_archived.ne(comparison.broad_final)
comparison.to_csv(P/'EFIGA_archived_vs_final_labels.csv',index=False)
efiga=efiga.merge(biomarker_cohort[['ID','label']],on='ID',how='left',validate='one_to_one')
assert mci_ids <= set(efiga.loc[efiga.clinical_stage.eq('MCI'),'ID'])
excluded=efiga[efiga.clinical_stage.eq('MCI') & ~efiga.ID.isin(mci_ids)].copy()
excluded.to_csv(P/'EFIGA_excluded_MCI_without_final_biomarkers.csv',index=False)
efiga=efiga[efiga.clinical_stage.ne('MCI') | efiga.ID.isin(mci_ids)].copy()
assert efiga.clinical_stage.eq('MCI').sum()==149
assert efiga.clinical_stage.value_counts().to_dict()=={'CN':388,'MCI':149,'AD':130}
adni['highlight_group']=adni.clinical_stage.where(adni.clinical_stage.ne('MCI'),adni.big_label)
efiga['highlight_group']=efiga.clinical_stage.where(efiga.clinical_stage.ne('MCI'),efiga.label.str.replace(r'[12]$','',regex=True))
for name,d in [('ADNI',adni),('EFIGA',efiga)]:
 assert d.ID.is_unique and d[['UMAP1','UMAP2','clinical_stage','highlight_group']].notna().all().all()
 assert set(d.highlight_group)<= {'CN','AD','EMCI','LMCI'}
 d.to_csv(P/f'{name}_plot_coordinates.csv',index=False)
 print(name,d.groupby(['clinical_stage','highlight_group']).size().to_string())
colors={'CN':'#4E79A7','MCI':'#F28E2B','AD':'#E15759','EMCI':'#E69F00','LMCI':'#009E73'}
plt.rcParams.update({'font.size':14,'axes.titlesize':19,'axes.labelsize':16,'legend.fontsize':13,'pdf.fonttype':42,'ps.fonttype':42})
adni_nodes=pd.read_csv(A/'ADNI_tree_nodes_UMAP.csv').set_index('node')
adni_edges=pd.read_csv(A/'ADNI_tree_edges.csv').to_numpy(dtype=int)
adni_pca_nodes=pd.read_csv(A/'ADNI_tree_nodes_PCA.csv',index_col='node').to_numpy()
adni_root=int(json.loads((A/'ADNI_trajectory_config.json').read_text())['root'])
lengths=np.linalg.norm(adni_pca_nodes[adni_edges[:,1]]-adni_pca_nodes[adni_edges[:,0]],axis=1)
graph=csr_matrix((np.r_[lengths,lengths],(np.r_[adni_edges[:,0],adni_edges[:,1]],np.r_[adni_edges[:,1],adni_edges[:,0]])),shape=(len(adni_nodes),len(adni_nodes)))
node_pseudotime=dijkstra(graph,directed=False,indices=adni_root)
assert np.isfinite(node_pseudotime).all()
adni_nodes['pseudotime']=node_pseudotime
adni_nodes.to_csv(P/'ADNI_tree_nodes_UMAP_pseudotime.csv',index_label='node')
pt_norm=plt.Normalize(0,node_pseudotime.max())
fig,axes=plt.subplots(2,2,figsize=(16,12))
for row,(name,d) in enumerate([('ADNI',adni),('EFIGA',efiga)]):
 for col in [0,1]:
  ax=axes[row,col]; field='clinical_stage' if col==0 else 'highlight_group'
  groups=['CN','MCI','AD'] if col==0 else ['CN','AD','EMCI','LMCI']
  for group in groups:
   x=d[d[field].eq(group)]
   alpha=.3 if col==1 and group in ['CN','AD'] else .85
   ax.scatter(x.UMAP1,x.UMAP2,s=24 if col==0 or group in ['CN','AD'] else 38,color=colors[group],alpha=alpha,edgecolors='none',label=f'{group} (n={len(x)})',zorder=3 if group in ['EMCI','LMCI'] else 2)
  if row==0:
   segments=[adni_nodes.loc[[i,j],['UMAP1','UMAP2']].to_numpy() for i,j in adni_edges]
   path=LineCollection(segments,cmap='viridis',norm=pt_norm,linewidths=2,zorder=4)
   path.set_array(node_pseudotime[adni_edges].mean(axis=1));ax.add_collection(path)
   root=adni_nodes.loc[adni_root]
   ax.scatter(root.UMAP1,root.UMAP2,marker='*',c='black',s=125,zorder=5,label='_nolegend_')
   fig.colorbar(path,ax=ax,fraction=.045,pad=.025,label='Pseudotime (tree distance)')
  ax.set(title=f'{name}: '+('CN, MCI and AD' if col==0 else 'EMCI and LMCI highlighted'),xlabel='UMAP 1',ylabel='UMAP 2')
  ax.text(-.12,1.05,'ABCD'[row*2+col],transform=ax.transAxes,fontweight='bold',fontsize=28)
  ax.legend(loc='lower left',frameon=True,facecolor='white',framealpha=.92,markerscale=1.5)
  ax.spines[['top','right']].set_visible(False)
  # Identical coordinate windows within each cohort.
  for coord, setter in [('UMAP1',ax.set_xlim),('UMAP2',ax.set_ylim)]:
   lo,hi=d[coord].min(),d[coord].max();margin=(hi-lo)*.06;setter(lo-margin,hi+margin)
fig.tight_layout(pad=2.3,w_pad=3,h_pad=3)
for ext in ['png','pdf','svg']:fig.savefig(P/f'Metabolomics_UMAP_ADNI_EFIGA_2x2.{ext}',dpi=300,bbox_inches='tight')
(P/'input_manifest.json').write_text(json.dumps({'ADNI_coordinates':str(A/'ADNI_umap.csv'),'ADNI_metadata':str(A/'ADNI_metadata.csv'),'EFIGA_coordinates':str(E/'SELECTED_full_cohort_UMAP_coordinates.csv'),'EFIGA_subtypes':str(biomarker_path),'EFIGA_archived_subtypes':str(E/'SELECTED_labels.csv'),'label_subtype_changes':int(comparison.subtype_changed.sum()),'label_broad_group_changes':int(comparison.broad_group_changed.sum()),'EFIGA_final_broad_counts':comparison.broad_final.value_counts().to_dict(),'EFIGA_archived_broad_counts':comparison.broad_archived.value_counts().to_dict(),'tree_refitted':False,'ADNI_tree_PCA':str(A/'ADNI_tree_nodes_PCA.csv'),'ADNI_tree_edges':str(A/'ADNI_tree_edges.csv'),'ADNI_tree_display':str(P/'ADNI_tree_nodes_UMAP.csv'),'ADNI_trajectory_config':str(A/'ADNI_trajectory_config.json'),'trajectory_panels':['A','B'],'EFIGA_tree_displayed':False,'CN_AD_highlight_alpha':0.3,'EFIGA_MCI_biomarker_cohort':str(biomarker_path),'EFIGA_MCI_retained':149,'EFIGA_MCI_excluded':len(excluded),'coordinates_refitted':False},indent=2))

