"""Assemble Figure6 A–I from the selected 149-patient native Monocle3 trajectory and recomputed DE/GAM results."""
from pathlib import Path
import json,textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.collections import LineCollection
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'
plt.rcParams.update({'font.size':12,'axes.titlesize':15,'axes.labelsize':12,'xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':10,'pdf.fonttype':42,'ps.fonttype':42})
ORDER=['EMCI1','EMCI2','LMCI1','LMCI2'];COLORS=['#F8766D','#7CAE00','#00BFC4','#C77CFF']
d=pd.read_csv(ROOT/'inputs/patients.csv',index_col='ID')
nodes=pd.read_csv(ROOT/'inputs/selected_trajectory_nodes.csv',index_col='node')
edges=pd.read_csv(ROOT/'inputs/selected_trajectory_edges.csv')
# Draw the longest root-to-terminal path; full native graph remains in inputs.
import networkx as nx
graph=nx.Graph()
for _,v in edges.iterrows():graph.add_edge(v['from'],v['to'],weight=float(np.linalg.norm(nodes.loc[v['from'],['UMAP1','UMAP2']]-nodes.loc[v['to'],['UMAP1','UMAP2']])))
root_node=pd.read_csv(ROOT/'inputs/selected_trajectory_configuration.csv').iloc[0].root_node
distances,paths=nx.single_source_dijkstra(graph,root_node)
terminal=max([v for v in graph if graph.degree[v]==1],key=lambda v:distances[v])
main_path=paths[terminal];main_pairs={frozenset([a,b])for a,b in zip(main_path[:-1],main_path[1:])}
edges=edges.loc[[frozenset([a,b])in main_pairs for a,b in edges[['from','to']].to_numpy()]]
nodes=nodes.loc[main_path]
config=pd.read_csv(ROOT/'inputs/selected_trajectory_configuration.csv').iloc[0]
mapping=pd.read_csv(ROOT/'inputs/metabolite_mapping.csv').set_index('target_hmdb')
names=mapping.Name.to_dict()
full_labels=pd.read_csv(OUT/'full_metabolite_labels.csv').set_index('target_hmdb').figure_label.to_dict()
def display(h):
 return full_labels[h]
de=pd.read_csv(OUT/'EFIGA_broad_DE.csv');ad=pd.read_csv(OUT/'ADNI_broad_DE.csv')
comp=de[['HMDB','standardized_effect','P.Value','adj.P.Val']].merge(ad[['HMDB','standardized_effect','P.Value','adj.P.Val']],on='HMDB',suffixes=('_EFIGA','_ADNI'))
comp['same_direction']=comp.standardized_effect_EFIGA*comp.standardized_effect_ADNI>0
comp['selected_raw']=comp.same_direction&comp['P.Value_EFIGA'].lt(.05)&comp['P.Value_ADNI'].lt(.05)
comp['selected_FDR']=comp.same_direction&comp['adj.P.Val_EFIGA'].lt(.05)&comp['adj.P.Val_ADNI'].lt(.05)
comp['Name']=comp.HMDB.map(names);comp.to_csv(OUT/'panel_C_cross_cohort_comparison.csv',index=False)
vC=comp[comp.selected_raw].sort_values('standardized_effect_EFIGA');vC.to_csv(OUT/'panel_C_selected_raw_p005.csv',index=False)
gam=pd.read_csv(OUT/'shared54_original_style_pair_GAM.csv')
GCOLS=['EMCI1->EMCI2','EMCI1->LMCI1','EMCI1->LMCI2']
r2=gam.pivot(index='feature',columns='comparison',values='R2_adjusted')[GCOLS]
gq=gam.pivot(index='feature',columns='comparison',values='q_BH')[GCOLS]
r2_differences=r2.diff(axis=1).iloc[:,1:]
increasing_r2=r2_differences.gt(0).all(axis=1);decreasing_r2=r2_differences.lt(0).all(axis=1)
audit=r2.copy();audit['all_R2_gt_015']=r2.gt(.15).all(axis=1);audit['monotonic_R2']=increasing_r2|decreasing_r2;audit['selected']=audit.all_R2_gt_015 & audit.monotonic_R2
audit['direction']=np.where(increasing_r2,'increasing R2',np.where(decreasing_r2,'decreasing R2','nonmonotonic R2'));audit.to_csv(OUT/'panel_E_monotonic_R2_selection_audit.csv')
selected=r2.loc[audit.selected]
R2_VMAX=float(np.ceil(selected.to_numpy().max()*10)/10) if len(selected) else 1.0
features=selected.sort_values('EMCI1->LMCI2',ascending=False).index.tolist()
selected.loc[features].to_csv(OUT/'panel_E_R2.csv');gq.loc[features].to_csv(OUT/'panel_E_q_BH.csv')
pathrows=[]
for group in ['EMCI','LMCI']:
 tab=pd.read_csv(OUT/'metaboanalyst'/group/'pathway_results.csv',index_col=0).reset_index()
 tab=tab.rename(columns={tab.columns[0]:'pathway','Raw p':'p_raw','Hits':'hits','Impact':'impact'})
 tab['group']=group;pathrows.append(tab)
pathways=pd.concat(pathrows,ignore_index=True)
pathways['selected']=pathways.hits.gt(1)&pathways.impact.gt(.001)
pathways.to_csv(OUT/'panel_F_all_pathways.csv',index=False)
pathselected=pathways[pathways.selected].copy();pathselected['minus_log10_p']=-np.log10(pathselected.p_raw)
pathselected.to_csv(OUT/'panel_F_displayed_pathways.csv',index=False)
stat=pd.read_csv(ROOT/'inputs/biomarker_statistics.csv')
impactnorm=plt.Normalize(0,max(float(pathselected.impact.max()),1e-6))

def volcano(ax):
 x=de.logFC.to_numpy();q=de['adj.P.Val'].to_numpy();y=-np.log10(np.clip(q,1e-300,None));sig=q<.05
 for mask,color,label in [(sig&(x<0),'#00BFC4','Lower in LMCI'),(~sig,'#bbbbbb','FDR ≥ 0.05'),(sig&(x>0),'#F8766D','Higher in LMCI')]:ax.scatter(x[mask],y[mask],s=18,c=color,label=label)
 ax.axhline(-np.log10(.05),c='#777777',ls='--',lw=.8);ax.axvline(0,c='#aaaaaa',lw=.6)
 for positive in [False,True]:
  top=de[(de['adj.P.Val']<.05)&((de.logFC>0) if positive else (de.logFC<0))].nsmallest(5,'adj.P.Val').sort_values('adj.P.Val',ascending=False)
  for i,r in enumerate(top.itertuples()):
   ax.annotate(textwrap.fill(display(r.HMDB),31),(r.logFC,-np.log10(de.loc[de.HMDB.eq(r.HMDB),'adj.P.Val'].iloc[0])),xytext=(.98 if positive else .02,.25+.57*i/max(len(top)-1,1)),textcoords='axes fraction',ha='right' if positive else 'left',va='center',fontsize=9,arrowprops=dict(arrowstyle='-',color='#888888',lw=.5))
 ax.set_xlabel('logFC');ax.set_ylabel('−log₁₀(FDR)');ax.set_title('EFIGA differential metabolites: LMCI vs EMCI');ax.legend(fontsize=9,loc='upper left',bbox_to_anchor=(1.02,1),borderaxespad=0);ax.margins(x=.35,y=.22);ax.set_box_aspect(1)

def ellipse(ax,xy,color):
 vals,vec=np.linalg.eigh(np.cov(xy.T));ix=np.argsort(vals)[::-1];vals=vals[ix];vec=vec[:,ix]
 dims=2*np.sqrt(vals*chi2.ppf(.85,2));angle=np.degrees(np.arctan2(vec[1,0],vec[0,0]))
 ax.add_patch(Ellipse(xy.mean(axis=0),*dims,angle=angle,facecolor=color,edgecolor='none',alpha=.12))
 ax.add_patch(Ellipse(xy.mean(axis=0),*dims,angle=angle,fill=False,edgecolor=color,linestyle='--',linewidth=1))

def umap(ax,big=False):
 labels=d.label.str[:4] if big else d.label
 for g,color in zip(['EMCI','LMCI'] if big else ORDER,['#E41A1C','#377EB8'] if big else COLORS):
  xy=d.loc[labels.eq(g),['UMAP1','UMAP2']].to_numpy();ellipse(ax,xy,color);ax.scatter(xy[:,0],xy[:,1],color=color,s=15,alpha=.85,label=f'{g} (n={len(xy)})',zorder=3)
 coords=nodes[['UMAP1','UMAP2']];segments=[coords.loc[[r['from'],r['to']]].to_numpy() for _,r in edges.iterrows()]
 values=[nodes.loc[[r['from'],r['to']],'pseudotime'].mean() for _,r in edges.iterrows()]
 norm=plt.Normalize(float(d.pseudotime.min()),float(d.pseudotime.max()))
 ax.add_collection(LineCollection(segments,colors='black',linewidths=2.7,zorder=4))
 lc=LineCollection(segments,cmap='viridis',norm=norm,linewidths=1.6,zorder=5);lc.set_array(np.array(values));ax.add_collection(lc)
 ax.scatter(nodes.UMAP1,nodes.UMAP2,c=nodes.pseudotime,cmap='viridis',norm=norm,s=18,edgecolors='black',linewidths=.4,zorder=6)
 rt=nodes.loc[config.root_node];ax.scatter(rt.UMAP1,rt.UMAP2,marker='*',c='black',s=105,zorder=7)
 ax.annotate('Root',(rt.UMAP1,rt.UMAP2),xytext=(6,-12),textcoords='offset points',fontsize=10)
 ax.figure.colorbar(lc,ax=ax,fraction=.035,pad=.025,label='Monocle3 pseudotime')
 ax.legend(fontsize=10,loc='best');ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');ax.set_title('149 MCI: broad groups' if big else '149 MCI: four subtypes')

def direction(ax):
 yy=np.arange(len(vC));ax.barh(yy-.18,vC.standardized_effect_EFIGA,height=.34,color='#F8766D',label='EFIGA');ax.barh(yy+.18,vC.standardized_effect_ADNI,height=.34,color='#00BFC4',label='ADNI')
 ax.set_yticks(yy,[textwrap.fill(display(h),40) for h in vC.HMDB],fontsize=11);ax.invert_yaxis();ax.axvline(0,c='#888888',lw=.7)
 ax.legend(fontsize=11);ax.set_xlabel('Standardized LMCI − EMCI difference');ax.set_title('Concordant differential metabolites\nRaw p < 0.05 in EFIGA and ADNI')

def heatmap(fig,slot):
 if not features:
  ax=fig.add_subplot(slot);ax.text(.5,.5,'No metabolites meet adjusted R² > 0.15\nin all three EMCI1-referenced subsets',ha='center',va='center');ax.axis('off');return ax
 gs=slot.subgridspec(1,2,width_ratios=[1,.055],wspace=.15)
 ax=fig.add_subplot(gs[0,0]);cax=fig.add_subplot(gs[0,1])
 mat=selected.loc[features];im=ax.imshow(mat.to_numpy(),cmap='viridis',vmin=0,vmax=R2_VMAX,aspect='auto')
 ax.set_yticks(range(len(features)),[textwrap.fill(display(f),34) for f in features],fontsize=10)
 ax.set_xticks(range(3),['EMCI2','LMCI1','LMCI2'],fontsize=10)
 for i,f in enumerate(features):
  for j,c in enumerate(GCOLS):
   val=mat.loc[f,c];ax.text(j,i,f'{val:.2f}'+('*' if gq.loc[f,c]<.05 else ''),ha='center',va='center',fontsize=10,color='black' if val/R2_VMAX>.8 else 'white')
 ax.set_xticks(np.arange(-.5,3,1),minor=True);ax.set_yticks(np.arange(-.5,len(features),1),minor=True);ax.grid(which='minor',color='white',linewidth=.7);ax.tick_params(which='minor',bottom=False,left=False)
 ax.set_title('Metabolite–pseudotime associations',fontsize=14)
 ax.set_xlabel('Comparison group (reference: EMCI1)',fontsize=11)
 fig.colorbar(im,cax=cax,label='GAM adjusted R²');return ax

def pathway(ax,group):
 rows=pathselected[pathselected.group.eq(group)].sort_values('minus_log10_p',ascending=False)
 scores=rows.minus_log10_p.to_numpy();yy=np.arange(len(rows))
 ax.barh(yy,scores,color=plt.get_cmap('viridis')(impactnorm(rows.impact.to_numpy())),height=.55 * (.7 if group == 'LMCI' else 1))
 ax.set_yticks(yy,[textwrap.fill(v,27) for v in rows.pathway],fontsize=10);ax.invert_yaxis()
 # Keep the original vertical range so thinner LMCI bars remain visibly thinner.
 if group=='LMCI' and len(rows):ax.set_ylim(len(rows)-1+.55*.55,-.55*.55)
 ax.axvline(-np.log10(.05),c='#333333',ls='--',lw=1,label='p = 0.05')
 ax.set_xlim(0,max(1.55,float(scores.max())*1.15 if len(scores) else 1.55));ax.set_xlabel('−log₁₀(p)',fontsize=12)
 ax.set_title(f'{group}2 vs {group}1: KEGG pathways',fontsize=14);ax.tick_params(axis='x',labelsize=11)
 ax.spines[['top','right']].set_visible(False)

def split_pathways(fig,slot):
 gs=slot.subgridspec(2,2,width_ratios=[1,.055],height_ratios=[2,1],hspace=.60,wspace=.10)
 axes=[fig.add_subplot(gs[i,0]) for i in range(2)]
 for ax,g in zip(axes,['EMCI','LMCI']):pathway(ax,g)
 cax=fig.add_subplot(gs[:,1]);fig.colorbar(plt.cm.ScalarMappable(norm=impactnorm,cmap='viridis'),cax=cax,label='Pathway impact')
 return axes[0]

def biomarker(ax,marker,title):
 values=[d.loc[d.label.eq(g),marker].dropna().to_numpy() for g in ORDER];medians=[np.median(v) for v in values]
 for i,(v,c) in enumerate(zip(values,COLORS)):
  ax.boxplot(v,positions=[i],widths=.50,showfliers=False,medianprops={'color':c,'linewidth':1.8},boxprops={'color':c},whiskerprops={'color':c},capprops={'color':c})
  ax.scatter(np.random.default_rng(666+i).normal(i,.055,len(v)),v,s=12,c=c,alpha=.55)
 ax.plot(range(4),medians,'s--',color='#555555',lw=1,ms=3,label='Median');ax.legend(fontsize=10,loc='upper right')
 top=max(v.max() for v in values);bottom=min(v.min() for v in values);step=max(top-bottom,top*.1)*.08
 for k,r in enumerate(stat[stat.marker.eq(marker)].sort_values(['j','i']).itertuples()):
  height=top+(k+1)*step
  ax.plot([r.i,r.i,r.j,r.j],[height,height+step*.15,height+step*.15,height],c='#555555',lw=.6)
  ax.text((r.i+r.j)/2,height+step*.22,f'p={r.p_raw:.3g}',ha='center',fontsize=9)
 ax.set_ylim(bottom-step,top+8.5*step);ax.set_xticks(range(4),[f'{g}\nn={len(v)}' for g,v in zip(ORDER,values)],fontsize=11)
 ax.set_title(title+'\nOne-sided Wilcoxon; raw p',fontsize=15);ax.set_ylabel(title)

def proportion(ax):
 counts=[]
 for g in ORDER:
  v=d.loc[d.label.eq(g),'pTau181_positive'].dropna();counts.append((len(v),int(v.eq(1).sum()),int(v.eq(0).sum())))
 positive=np.array([p/n for n,p,z in counts]);negative=1-positive;xx=np.arange(4)
 ax.bar(xx,positive,color='#F8766D',label='Positive',width=.70);ax.bar(xx,negative,bottom=positive,color='#00BFC4',label='Negative',width=.70)
 for i,(p,n) in enumerate(zip(positive,negative)):
  ax.text(i,p/2,f'{100*p:.1f}%',ha='center',va='center',fontsize=12);ax.text(i,p+n/2,f'{100*n:.1f}%',ha='center',va='center',fontsize=12)
 ax.set_xticks(xx,[f'{g}\nn={c[0]}' for g,c in zip(ORDER,counts)],fontsize=11)
 ax.set_yticks([0,.25,.5,.75,1],['0%','25%','50%','75%','100%']);ax.set_ylim(0,1.05)
 ax.set_title('pTau181 positivity (≥2.63 pg/mL)');ax.set_ylabel('Proportion');ax.legend(fontsize=10,loc='upper left',bbox_to_anchor=(1.02,1),borderaxespad=0);ax.set_box_aspect(1)
 pd.DataFrame([dict(label=g,n=n,positive=p,negative=z,proportion=p/n) for g,(n,p,z) in zip(ORDER,counts)]).to_csv(OUT/'panel_G_positivity_counts.csv',index=False)

VARIANTS={'tau217_abratio': [('P_PT217','pTau217'),('P_AB42AB40ratio','Aβ42/Aβ40')]}
def render_panel(fig,slot,letter,variant):
 if letter=='E':ax=heatmap(fig,slot)
 elif letter=='F':ax=split_pathways(fig,slot)
 else:
  ax=fig.add_subplot(slot)
  if letter=='A':volcano(ax)
  elif letter=='B':umap(ax,True)
  elif letter=='C':direction(ax)
  elif letter=='D':umap(ax,False)
  elif letter=='G':proportion(ax)
  else:
   marker,title=VARIANTS[variant][0 if letter=='H' else 1];biomarker(ax,marker,title)
  ax.spines[['top','right']].set_visible(False)
 ax.text(-.14,1.11,letter,transform=ax.transAxes,fontsize=26,fontweight='bold')
 return ax

def enlarge_subplot_fonts(fig, factor=1.3):
 # Include annotations, tick labels, legends, and colorbar text.
 from matplotlib.text import Text
 fig.canvas.draw()
 for text in fig.findobj(match=Text):
  text.set_fontsize(text.get_fontsize()*factor)

def main():
 subplot_dir=ROOT/'figures'/'subplot'
 subplot_dir.mkdir(parents=True,exist_ok=True)
 for variant in VARIANTS:
  fig=plt.figure(figsize=(28,23))
  gs=fig.add_gridspec(3,3,width_ratios=[1.05,1.2,1.05],left=.06,right=.98,bottom=.065,top=.93,hspace=.50,wspace=.90)
  for i,letter in enumerate('ABCDEFGHI'):render_panel(fig,gs[i//3,i%3],letter,variant)
  for ext in ['png','pdf']:fig.savefig(ROOT/'figures'/f'Figure6_complete_{variant}.{ext}',dpi=300,bbox_inches='tight')
  plt.close(fig)
 for letter in 'ABCDEFG':
  if letter=='F':
   fig=plt.figure(figsize=(10,8))
   slot=fig.add_gridspec(1,1,left=.32,right=.90,bottom=.12,top=.86)[0,0]
  else:
   fig=plt.figure(figsize=(10,8),layout='constrained')
   slot=fig.add_gridspec(1,1)[0,0]
  render_panel(fig,slot,letter,'tau217_abratio')
  enlarge_subplot_fonts(fig)
  for ext in ['png','pdf']:
   fig.savefig(ROOT/'panels'/f'Figure6{letter}.{ext}',dpi=300,bbox_inches='tight')
   fig.savefig(subplot_dir/f'Figure6{letter}.{ext}',dpi=300,bbox_inches='tight')
  plt.close(fig)
 for variant in VARIANTS:
  for letter in ['H','I']:
   fig=plt.figure(figsize=(9,8),layout='constrained');render_panel(fig,fig.add_gridspec(1,1)[0,0],letter,variant)
   enlarge_subplot_fonts(fig)
   for ext in ['png','pdf']:
    fig.savefig(ROOT/'panels'/f'Figure6{letter}_{variant}.{ext}',dpi=300,bbox_inches='tight')
    fig.savefig(subplot_dir/f'Figure6{letter}_{variant}.{ext}',dpi=300,bbox_inches='tight')
    if variant=='tau217_abratio':fig.savefig(subplot_dir/f'Figure6{letter}.{ext}',dpi=300,bbox_inches='tight')
   plt.close(fig)
 for g in ['EMCI','LMCI']:
  fig,ax=plt.subplots(figsize=(9,5),layout='constrained');pathway(ax,g);fig.colorbar(plt.cm.ScalarMappable(norm=impactnorm,cmap='viridis'),ax=ax,label='Pathway impact')
  enlarge_subplot_fonts(fig)
  for ext in ['png','pdf']:
   fig.savefig(ROOT/'panels'/f'Figure6F_{g}.{ext}',dpi=300,bbox_inches='tight')
   fig.savefig(subplot_dir/f'Figure6F_{g}.{ext}',dpi=300,bbox_inches='tight')
  plt.close(fig)
 pd.DataFrame([dict(label=g,n_metabolites=int(d.label.eq(g).sum()),n_biomarkers=int(d.loc[d.label.eq(g),'P_PT181V2'].notna().sum()),pseudotime_mean=float(d.loc[d.label.eq(g),'pseudotime'].mean()),pseudotime_median=float(d.loc[d.label.eq(g),'pseudotime'].median()),**{m+'_median':float(d.loc[d.label.eq(g),m].median()) for m in ['P_PT181V2','P_PT217','P_AB42AB40ratio']}) for g in ORDER]).to_csv(OUT/'patient_group_summary.csv',index=False)
 summary=dict(n_metabolite_patients=len(d),n_biomarker_patients=int(d.P_PT181V2.notna().sum()),n_A_FDR_DE=int(de['adj.P.Val'].lt(.05).sum()),n_A_higher_LMCI=int((de['adj.P.Val'].lt(.05)&de.logFC.gt(0)).sum()),n_A_lower_LMCI=int((de['adj.P.Val'].lt(.05)&de.logFC.lt(0)).sum()),n_C_raw_concordant=int(comp.selected_raw.sum()),n_C_FDR_concordant=int(comp.selected_FDR.sum()),panel_C_names=[display(h) for h in vC.HMDB],panel_E_color_max=R2_VMAX,n_E_selected=len(features),panel_E_features=features,pathways=pathselected[['group','pathway','p_raw','FDR','impact','hits']].to_dict('records'),variant_names=list(VARIANTS),trajectory_root=str(config.root_node))
 (OUT/'figure_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
