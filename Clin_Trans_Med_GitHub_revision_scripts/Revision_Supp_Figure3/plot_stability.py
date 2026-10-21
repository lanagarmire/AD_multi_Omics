"""Plot final archived dimension/subsampling metrics without model fitting or selection."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent
INPUT=OUT/'inputs'
LABELS=['50','172','512','1024','5000','10000','20032']
def order(version):
    return [f'cached_{d}' for d in [50,172,512,1024,5000,10000]]+['uncompressed']
def make_composite_raw(summary):
    """Full-data A-C above transposed EMCI/LMCI subsampling rows."""
    full=pd.read_csv(INPUT/'figure_metrics_raw20032.csv')
    raw_order=order('raw20032')
    dims=[50,172,512,1024,5000,10000,20032]
    x=np.arange(7)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    fig,axes=plt.subplots(3,3,figsize=(18,13),sharex='col')
    top_specs=[
        ('ARI_mean_including_self','Adjusted Rand Index'),
        ('Dunn_legacy_graph_rows','Dunn Index'),
        ('conductance','Graph Conductance')]
    for col,(metric,title) in enumerate(top_specs):
        ax=axes[0,col]
        stage_values=[]
        for stage,color in [('EMCI','#377eb8'),('LMCI','#e41a1c')]:
            d=full[full.stage.eq(stage)].set_index('representation').loc[raw_order]
            values=d[metric].to_numpy()
            stage_values.append(values)
            ax.plot(x,values,'o-',color=color,label=stage)
        ax.plot(x,np.mean(stage_values,axis=0),'--',color='#333333',label='Mean across stages')
        ax.set_title(title,fontsize=16,fontweight='bold',pad=12)
        ax.text(-.08,1.06,chr(65+col),transform=ax.transAxes,fontsize=17,fontweight='bold')
        ax.set_ylabel('Cross-representation ARI' if col==0 else ('Dunn index' if col==1 else 'Graph conductance'))
        ax.legend(fontsize=8)
        ax.grid(axis='y',alpha=.2)
        # Reproduce the review table: equal-weight means across the two stages.
        # Its ARI excludes self-comparisons, unlike the plotted ARI series.
        table_metric='ARI_mean_excluding_self' if col==0 else metric
        means=full.groupby('representation')[table_metric].mean().loc[raw_order]
        table_label='Mean ARI excluding self' if col==0 else f'Mean {title.lower()}'
        ax.text(.5,-.085,table_label,transform=ax.transAxes,
                ha='center',va='top',fontsize=9)
        table=ax.table(cellText=[[f'{v:.4f}' for v in means]],
                       colLabels=[f'{d:,}' for d in dims],
                       cellLoc='center',bbox=[0,-.28,1,.17])
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        for (table_row,table_col),cell in table.get_celld().items():
            cell.set_edgecolor('#cccccc')
            cell.set_linewidth(.5)
            if table_row==0:
                cell.set_facecolor('#f0f0f0')
            if table_col==1:
                cell.set_facecolor('#fff0cc')
                cell.get_text().set_weight('bold')
    lower_specs=[
        ('ARI_fused_stability','ARI with full-data fused partition','Subsample partition stability'),
        ('Dunn_fused','Dunn index','Subsample cluster separation'),
        ('conductance_fused','Graph conductance','Subsample graph separation')]
    panel=3
    for row,stage in enumerate(['EMCI','LMCI'],start=1):
        for col,(metric,ylabel,title) in enumerate(lower_specs):
            ax=axes[row,col]
            d=summary[(summary.stage.eq(stage))&(summary.metric.eq(metric))].set_index('representation').loc[raw_order]
            y=d['median'].to_numpy()
            err=np.vstack([y-d.low.to_numpy(),d.high.to_numpy()-y])
            ax.errorbar(x,y,yerr=np.maximum(err,0),fmt='o-',capsize=3,color='#377eb8' if stage=='EMCI' else '#e41a1c')
            ax.set_title(f'{stage}: {title.title()}',fontsize=14,fontweight='bold',pad=12)
            ax.text(-.08,1.06,chr(65+panel),transform=ax.transAxes,fontsize=17,fontweight='bold')
            panel+=1
            ax.set_ylabel(ylabel)
            ax.grid(axis='y',alpha=.2)
            ax.set_xlabel('Transcriptomic dimension')
    for ax in axes.flat:
        ax.axvspan(.75,1.25,color='#ffd166',alpha=.3)
        ax.set_xticks(x,LABELS,rotation=45)
    fig.subplots_adjust(left=.07,right=.98,bottom=.085,top=.94,wspace=.24,hspace=.60)
    for ext,kwargs in [('png',dict(dpi=300)),('pdf',{})]:
        fig.savefig(OUT/f'Supplementary_Figure3_A_I.{ext}',bbox_inches='tight',**kwargs)
    plt.close(fig)


if __name__=='__main__':
    make_composite_raw(pd.read_csv(INPUT/'subsampling80_summary.csv'))
