"""Exact UpSet intersections of DE sets; export PNG/PDF and membership tables."""
from pathlib import Path
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
OUT=HERE/'results/overlap_plots'; OUT.mkdir(parents=True,exist_ok=True)

def upset(ax_top,ax_matrix,sets,title):
    names=list(sets);universe=set.union(*sets.values())
    memberships={f:tuple(int(f in sets[n]) for n in names) for f in sorted(universe)}
    counts=Counter(memberships.values())
    patterns=sorted(counts,key=lambda p:(-counts[p],p))
    xs=np.arange(len(patterns))
    colors=['#527f9c' if sum(p)==len(names) else '#8e9ba5' for p in patterns]
    ax_top.bar(xs,[counts[p] for p in patterns],color=colors)
    for x,p in zip(xs,patterns):ax_top.text(x,counts[p],str(counts[p]),ha='center',va='bottom',fontsize=8)
    ax_top.set_title(title,fontsize=15,fontweight='bold',pad=12);ax_top.set_ylabel('Features in exact intersection')
    ax_top.set_ylim(0,max([counts[p] for p in patterns]+[1])*1.2);ax_top.tick_params(axis='x',bottom=False,labelbottom=False)
    for x,p in zip(xs,patterns):
        ax_matrix.scatter([x]*len(names),range(len(names)),c=['#243746' if v else '#dddddd' for v in p],s=30)
        active=[i for i,v in enumerate(p) if v]
        if len(active)>1:ax_matrix.plot([x,x],[min(active),max(active)],color='#243746',zorder=0)
    ax_matrix.set_yticks(range(len(names)),[f'{n} (n={len(sets[n])})' for n in names],fontsize=8)
    ax_matrix.set_ylim(len(names)-.5,-.5);ax_matrix.set_xticks([])
    for ax in [ax_top,ax_matrix]:
        ax.set_xlim(-.6,max(len(patterns)-.4,.6))
        ax.spines[['top','right']].set_visible(False)
    rows=[dict(zip(names,p),intersection_n=counts[p]) for p in patterns]
    return memberships,rows

# Final compilation: four labeled panels, with retention shown inside each panel.
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
fig=plt.figure(figsize=(14,10.5))
outer=fig.add_gridspec(2,2,wspace=.42,hspace=.34)
models={'Baseline':'baseline','+ APOE ε4':'APOE4','+ Education':'education',
        '+ APOE ε4 + Education':'APOE4_education'}
folder=HERE/'results/additional_confounders'
for i,(stage,omics) in enumerate([('EMCI','MRNA'),('LMCI','MRNA'),('EMCI','META'),('LMCI','META')]):
    prefix=f'{stage}_{omics}'
    # All four models are required; missing results must not yield blank panels.
    sets={label:set(pd.read_csv(folder/f'{prefix}_{model}_DE.csv').feature)
          for label,model in models.items()}
    sub=outer[i//2,i%2].subgridspec(2,1,height_ratios=[3,1.4],hspace=.04)
    ax=fig.add_subplot(sub[0]);mat=fig.add_subplot(sub[1],sharex=ax)
    title=f'{stage} '+('Genes' if omics=='MRNA' else 'Metabolite')
    members,counts=upset(ax,mat,sets,title)
    ax.text(-.43,1.04,'ABCD'[i],transform=ax.transAxes,
            fontsize=21,fontweight='bold',ha='left',va='bottom')
    baseline=sets['Baseline'];joint=sets['+ APOE ε4 + Education']
    retained=len(baseline & joint)
    rate=100*retained/len(baseline) if baseline else float('nan')
    rate_label=f'{rate:.2f}'.rstrip('0').rstrip('.')
    ax.text(.79,.83,'Original DE retained after\nAPOE + Education:\n'
            f'{retained}/{len(baseline)} ({rate_label}%)',transform=ax.transAxes,
            fontsize=12,color='#2e86c1',ha='center',va='top',linespacing=1.35)
    pd.DataFrame(counts).to_csv(OUT/f'{prefix}_additional_confounders_intersections.csv',index=False)
    pd.DataFrame([dict(feature=f,**dict(zip(sets,p))) for f,p in members.items()]).to_csv(
        OUT/f'{prefix}_additional_confounders_membership.csv',index=False)
fig.subplots_adjust(left=.16,right=.98,bottom=.055,top=.95)
for ext in ['png','pdf']:
    fig.savefig(OUT/f'Supplementary_Figure8_APOE4_education_adj_upset.{ext}',dpi=220)
plt.close(fig)
print('Final Supplementary Figure 8 saved.')
