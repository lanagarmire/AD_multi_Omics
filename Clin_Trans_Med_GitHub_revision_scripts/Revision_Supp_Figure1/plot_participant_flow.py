"""One integrated participant flow; FIA/UPLC use their shared sample roster."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

H=Path(__file__).resolve().parent
counts=pd.read_csv(H/'inputs/preprocessing_counts.csv')
def n(cohort,step):
 row=counts[(counts.cohort==cohort)&(counts.step==step)&(counts.phase=='All')]
 assert len(row)==1
 return int(row.iloc[0].participants)
assert [n('FIA','Sample missingness'),n('FIA','Fasting'),n('Metabolomics','Bile acid overlap'),n('Paired omics','Omics overlap'),n('Paired omics','MCI selection')]==[1680,1573,1562,681,401]
external=pd.read_csv(H/'inputs/patient_group_summary.csv').set_index('label')
assert external.loc[['EMCI1','EMCI2','LMCI1','LMCI2'],'n_metabolites'].tolist()==[51,34,44,20]
assert external.n_biomarkers.sum()==external.n_metabolites.sum()==149
plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
fig,ax=plt.subplots(figsize=(17,11))
ax.set(xlim=(0,17),ylim=(0,11));ax.axis('off')
ink='#243b53';blue='#eaf2fa';green='#eaf5ee'
def subtypes(a,b,c,d):
 return f'EMCI-1 {a} · EMCI-2 {b}\nLMCI-1 {c} · LMCI-2 {d}'
def box(x,y,w,h,title,value,detail='',color=blue,fontsize=10.5):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.035,rounding_size=.07',facecolor=color,edgecolor='#91a4b7',linewidth=1))
 lines=[title,f'{value:,}' if isinstance(value,int) else value]
 if detail:lines.append(detail)
 ax.text(x+w/2,y+h/2,'\n'.join(lines),ha='center',va='center',fontsize=fontsize,color=ink,linespacing=1.3)
def arrow(x,y,xx,yy,label=''):
 ax.annotate('',xy=(xx,yy),xytext=(x,y),arrowprops={'arrowstyle':'->','color':'#6b7f91','lw':1.1})
 if label:ax.text((x+xx)/2+.10,(y+yy)/2,label,fontsize=9,color='#885649',ha='left',va='center',bbox={'facecolor':'white','edgecolor':'none','pad':1})
ax.text(.2,10.65,'Participant flow',fontsize=19,fontweight='bold',color=ink)
ax.text(.2,10.2,'ADNI',fontsize=13,fontweight='bold',color=ink)
ax.text(12.0,10.2,'EFIGA',fontsize=13,fontweight='bold',color=ink)

box(.25,8.85,4.15,1.05,'P180 records',1790,'ADNI1 849 · ADNI2/GO 941\nFIA + UPLC')
box(5.55,8.85,4.15,1.05,'Gene samples',744,'GO 295 · ADNI2 449',green)
box(.25,7.25,4.15,.95,'QC',1680,'157 P180 features')
arrow(2.325,8.85,2.325,8.2,'QC −58 · missing −1\nRepeats −51 records')
box(.25,5.8,4.15,.9,'Fasting',1573)
arrow(2.325,7.25,2.325,6.7,'Nonfasting −107')
box(.25,4.15,4.15,1.15,'Metabolomics',1562,'CN 374 · MCI 790 · AD 299\nSMC 99 · 172 metabolites')
arrow(2.325,5.8,2.325,5.3,'Unmatched −11')
box(5.55,7.25,4.15,.95,'Gene profiles',744,'20,032 genes',green)
arrow(7.625,8.85,7.625,8.2)
box(5.55,4.15,4.15,1.15,'Paired omics',681,'CN 239 · MCI 401 · AD 41',green)
arrow(4.4,4.725,5.55,4.725)
ax.text(4.975,5.1,'−881',ha='center',fontsize=9,color='#885649')
arrow(7.625,7.25,7.625,5.3,'Unpaired −63')
box(5.55,2.55,4.15,1.1,'MCI',401,subtypes(108,93,85,115),green)
arrow(7.625,4.15,7.625,3.65,'CN/AD −280')

# One branch point links the actual MCI output directly to all downstream analyses.
arrow(7.625,2.55,7.625,2.15)
ax.plot([1.6,7.7],[2.15,2.15],color='#6b7f91',lw=1.1)
outputs=[(.25,2.7,'Differential',401,subtypes(108,93,85,115)+'\nGene + metabolite'),
         (3.3,2.7,'Trajectory',1274,subtypes(108,93,85,115)+'\nCN 534 · AD 339'),
         (6.35,2.7,'Survival',387,subtypes(104,87,84,112))]
for x,w,title,value,detail in outputs:
 box(x,.55,w,1.4,title,value,detail,fontsize=9.5)
 arrow(x+w/2,2.15,x+w/2,1.95)

box(12.25,8.85,4.35,1.05,'Available cohort',692,'CN 388 · MCI 156 · AD 130\nUnclassified 18',green)
box(12.25,7.25,4.35,.95,'MCI',156,color=green)
arrow(14.425,8.85,14.425,8.2,'Non-MCI −536')
box(12.25,5.8,4.35,.9,'Complete biomarkers',149,'QC / scaling',green)
arrow(14.425,7.25,14.425,6.7,'Missing −7')
box(12.25,4.15,4.35,1.15,'projectR transfer',149,'54 shared metabolites',green)
arrow(9.7,3.1,12.25,4.725,'ADNI reference')
arrow(14.425,5.8,14.425,5.3,'Excluded 0')
box(12.55,2.35,3.8,1.25,'Differential / trajectory',149,subtypes(51,34,44,20),green)
arrow(14.425,4.15,14.45,3.6)
box(12.55,.55,3.8,1.4,'Biomarkers',149,subtypes(51,34,44,20),green)
arrow(14.45,2.35,14.45,1.95,'Excluded 0')
fig.savefig(H/'participant_flow.pdf',bbox_inches='tight')
fig.savefig(H/'participant_flow.png',dpi=240,bbox_inches='tight')

