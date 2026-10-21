"""Plot archived final modality scores with shared fused labels; no fitting or tuning."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
FIG=Path(__file__).resolve().parent
DIMS=[50,172,512,1024,5000,10000,20032]
def plot(df,rule):
 fig,axes=plt.subplots(1,2,figsize=(12,4.6),sharey=True)
 colors={'Integrated (fused)':'#1b9e77','Metabolomics':'#d95f02','Gene expression':'#7570b3'}
 for ax,stage in zip(axes,['EMCI','LMCI']):
  for name,color in colors.items():
   d=df[(df.stage==stage)&(df.label_rule==rule)&(df.modality==name)].set_index('transcriptomic_dimension').loc[DIMS]
   ax.plot(range(len(DIMS)),d.silhouette,'o-',label=name,color=color,linewidth=1.8,markersize=5)
  ax.set_title(stage);ax.set_xticks(range(len(DIMS)),[str(d) if d!=20032 else '20032\n(raw)' for d in DIMS],rotation=40)
  ax.set_xlabel('Transcriptomic dimension');ax.set_ylim(-.05,1);ax.grid(axis='y',alpha=.2);ax.axvspan(.75,1.25,color='#ffd166',alpha=.25)
 axes[0].set_ylabel('Affinity-based silhouette');axes[1].legend(loc='upper right',fontsize=9)
 fig.suptitle('Legacy silhouette: independently fitted clusters' if rule=='own_partition' else 'Legacy silhouette: shared fused labels within each dimension')
 fig.tight_layout();stem='bottleneck_legacy_silhouette_raw20032' if rule=='own_partition' else 'bottleneck_legacy_shared_fused_labels_raw20032'
 fig.savefig(FIG/f'{stem}.png',dpi=240,bbox_inches='tight');fig.savefig(FIG/f'{stem}.pdf',bbox_inches='tight');plt.close(fig)


if __name__=='__main__':
 plot(pd.read_csv(FIG/'inputs/silhouette_scores.csv'),'shared_fused_partition')
