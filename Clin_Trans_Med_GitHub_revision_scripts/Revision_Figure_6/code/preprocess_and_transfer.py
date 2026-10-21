"""Replay the selected EFIGA-only preprocessing and fixed hierarchical transfer rule."""
from pathlib import Path
import ast,json,itertools
import numpy as np,pandas as pd
from scipy.stats import rankdata,mannwhitneyu,chi2_contingency
from sklearn.preprocessing import StandardScaler
ROOT=Path(__file__).resolve().parents[1];I=ROOT/'inputs';O=ROOT/'results';ORDER=np.array(['EMCI1','EMCI2','LMCI1','LMCI2'])
# projectR::cluster2pattern and least-squares projection numerical equivalents.
def cluster2pattern_exact(X_samples_features,labels,classes):
 data=np.asarray(X_samples_features,float).T;labels=np.asarray(labels,str);p=np.zeros((data.shape[0],len(classes)))
 for j,c in enumerate(classes):
  sub=data[:,labels==c];mean_cluster=sub.mean(axis=0)
  for i,x in enumerate(sub):
   if np.std(x)>0 and np.std(mean_cluster)>0:
    v=np.corrcoef(x,mean_cluster)[0,1];p[i,j]=v if np.isfinite(v) else 0
 return p

def projectR_exact(target,patterns):return np.linalg.lstsq(patterns,np.asarray(target,float).T,rcond=None)[0].T

def main():
 ref=pd.read_csv(ROOT/'reference/EFIGA.csv',index_col=0);raw=pd.read_csv(I/'EFIGA_MCI_positive2000_QC_imputed.csv',index_col=0).loc[ref.index];mapping=pd.read_csv(I/'metabolite_mapping.csv').set_index('target_hmdb')
 assert raw.shape==(149,2000) and np.isfinite(raw).all().all()
 profile=np.sort(raw.to_numpy(),axis=1).mean(axis=0);q=np.vstack([np.interp(rankdata(v)-1,np.arange(len(v)),profile)for v in raw.to_numpy()]);z=StandardScaler().fit_transform(q);z=pd.DataFrame(z,index=raw.index,columns=raw.columns)
 shared=z.loc[:,mapping.loc[ref.columns,'feature_name']];shared.columns=ref.columns;err=float(np.max(np.abs(shared.to_numpy()-ref.to_numpy())));assert err<1e-10
 # Preserve the original CSV precision after independently verifying normalization.
 ref.to_csv(I/'EFIGA.csv',index_label='ID');ref.to_csv(I/'shared54.csv',index_label='ID')
 a=pd.read_csv(I/'ADNI.csv',index_col=0).loc[:,ref.columns];y=pd.read_csv(I/'ADNI_labels.csv',index_col=0).loc[a.index,'label'].to_numpy();big=np.array([v[:4]for v in y]);scores=[]
 for mask,labels,classes in [(np.ones(len(y),bool),big,['EMCI','LMCI']),(big=='EMCI',y,['EMCI1','EMCI2']),(big=='LMCI',y,['LMCI1','LMCI2'])]:
  patterns=cluster2pattern_exact(a.to_numpy()[mask],labels[mask],classes);projection=projectR_exact(ref.to_numpy(),patterns);scores.append(projection[:,1]-projection[:,0]);pd.DataFrame(patterns,index=ref.columns,columns=classes).to_csv(O/('projection_patterns_'+classes[0]+'.csv'))
 scores=np.array(scores).T;expected=pd.read_csv(ROOT/'reference/transfer_scores.csv',index_col=0).loc[ref.index];np.testing.assert_allclose(scores,expected,atol=1e-10,rtol=1e-10)
 rules=json.loads((I/'protocol.json').read_text())['selected_rank_cuts'];sort=np.argsort(scores[:,0],kind='stable');labels=np.empty(149,int)
 for ids,col,cut,offset in [(sort[:rules['broad_EMCI']],1,rules['EMCI1_within_EMCI'],0),(sort[rules['broad_EMCI']:],2,rules['LMCI1_within_LMCI'],2)]:
  ids=ids[np.argsort(scores[ids,col],kind='stable')];labels[ids[:cut]]=offset;labels[ids[cut:]]=offset+1
 d=pd.read_csv(I/'selected149_patients.csv',index_col=0).loc[ref.index];assert np.array_equal(ORDER[labels],d.label);d['label']=ORDER[labels];d.drop(columns=['UMAP1','UMAP2','pseudotime'],errors='ignore').to_csv(I/'patients.csv',index_label='ID')
 pd.DataFrame(scores,index=ref.index,columns=expected.columns).to_csv(O/'transfer_scores.csv',index_label='ID');d[['label']].assign(broad_label=d.label.str[:4]).to_csv(O/'transferred_labels.csv',index_label='ID')
 rows=[]
 for marker in ['P_PT181V2','P_PT217','P_AB42AB40ratio','pTau181_positive']:
  for i,j in itertools.combinations(range(4),2):
   g1,g2=ORDER[i],ORDER[j];v=d.loc[d.label==g1,marker].dropna().to_numpy();w=d.loc[d.label==g2,marker].dropna().to_numpy();row=dict(marker=marker,group1=g1,group2=g2,i=i,j=j,n1=len(v),n2=len(w))
   if marker=='pTau181_positive':
    table=np.array([[sum(v==0),sum(v==1)],[sum(w==0),sum(w==1)]]);stat,p,df,expect=chi2_contingency(table,correction=False);row.update(alternative='two-sided',test='Pearson chi-square (no Yates correction)',minimum_expected_count=expect.min(),negative1=table[0,0],positive1=table[0,1],negative2=table[1,0],positive2=table[1,1])
   else:
    alt='greater'if marker=='P_AB42AB40ratio'else'less';stat,p=mannwhitneyu(v,w,alternative=alt,method='asymptotic',use_continuity=True);row.update(alternative=alt,test='Wilcoxon rank-sum (tie and continuity corrected)')
   row.update(statistic=stat,p_raw=p);rows.append(row)
 stats=pd.DataFrame(rows);expected=pd.read_csv(ROOT/'reference/biomarker_statistics.csv');match=stats.merge(expected,on=['marker','group1','group2'],suffixes=('_new','_reference'));assert len(match)==24;np.testing.assert_allclose(match.p_raw_new,match.p_raw_reference,atol=1e-12,rtol=1e-12);stats.to_csv(I/'biomarker_statistics.csv',index=False);stats.to_csv(O/'biomarker_pairwise_statistics.csv',index=False)
 (O/'transfer_verification.json').write_text(json.dumps({'n_EFIGA':149,'n_ADNI':len(a),'counts':d.label.value_counts().to_dict(),'normalization_max_error':err,'selected_labels_reproduced':True,'transfer_scores_reproduced':True,'biomarker_p_values_reproduced':True,'ADNI_additional_normalization':False},indent=2));print('Preprocessing, label transfer and biomarker statistics verified.')
if __name__=='__main__':main()
