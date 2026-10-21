from pathlib import Path
import analyze as m
import numpy as np,pandas as pd,json,time,multiprocessing
from concurrent.futures import ProcessPoolExecutor
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans,AgglomerativeClustering
from sklearn.metrics import silhouette_score,calinski_harabasz_score
from scipy import stats
R=m.R

def fast_affinity(a,b):
 # Same normalized Euclidean distances and snf affinity kernel as make_affinity.
 # Gram evaluation avoids repeated O(n^2 * 20032) scalar distance loops.
 aff=[]
 for x in [a,b]:
  z=np.nan_to_num(stats.zscore(x,ddof=1));v=np.square(z).sum(1)
  d=np.sqrt(np.maximum(v[:,None]+v[None,:]-2*z@z.T,0));np.fill_diagonal(d,0)
  aff.append(m.snf.compute.affinity_matrix(d,K=10,mu=.3))
 return m.snf.snf(aff,K=20,t=20,alpha=1.)

def permute(i):
 rng=np.random.default_rng(np.random.SeedSequence([42,i]))
 a=np.column_stack([rng.permutation(m.G_M[:,j]) for j in range(m.G_M.shape[1])]);b=np.column_stack([rng.permutation(m.G_T[:,j]) for j in range(m.G_T.shape[1])]);x=np.concatenate([a,b],1);af=fast_affinity(a,b)
 return [dict(replicate=i,k=k,silhouette_affinity=float(m.snf.metrics.silhouette_score(af,lab)),calinski_harabasz=float(calinski_harabasz_score(x,lab))) for k in m.KS for lab in [m.labels(af,k)]]

def main():
 baseline=[]
 for stage in ['EMCI','LMCI']:
  a,b,ids,_=m.load(stage);raw=np.concatenate([a,b],1);af=m.affinity(a,b);obs=[]
  z=PCA(n_components=.8,svd_solver='full').fit_transform(StandardScaler().fit_transform(raw))
  for k in m.KS:
   lab=m.labels(af,k);obs.append(dict(k=k,silhouette_affinity=float(m.snf.metrics.silhouette_score(af,lab)),calinski_harabasz=float(calinski_harabasz_score(raw,lab))))
   for method in ['SNF + spectral clustering','PCA + k-means','PCA + Ward hierarchical clustering']:
    l=lab if method.startswith('SNF') else (KMeans(n_clusters=k,n_init=100,random_state=42).fit_predict(z) if 'k-means' in method else AgglomerativeClustering(n_clusters=k,linkage='ward').fit_predict(z))
    baseline.append(dict(stage=stage,method=method,k=k,native_silhouette=float(m.snf.metrics.silhouette_score(af,l) if method.startswith('SNF') else silhouette_score(z,l)),native_metric='Affinity' if method.startswith('SNF') else 'PCA Euclidean'))
  pd.DataFrame(obs).to_csv(R/f'{stage}_observed_validity.csv',index=False);print(stage,'observed',obs[0],flush=True)
  path=R/f'{stage}_original_defaults_permutations.csv';rows=pd.read_csv(path).to_dict('records') if path.exists() else [];done={r['replicate'] for r in rows};start=time.time()
  with ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn'),initializer=m.initialize,initargs=(a,b)) as pool:
   for j,batch in enumerate(pool.map(permute,[i for i in range(1000) if i not in done],chunksize=2)):
    rows.extend(batch)
    if (j+1)%25==0:pd.DataFrame(rows).to_csv(path,index=False);print(stage,'permutations',len(done)+j+1,'elapsed',round(time.time()-start),flush=True)
  d=pd.DataFrame(rows);d.to_csv(path,index=False);d.to_csv(R/f'{stage}_permutation_replicates.csv',index=False)
  summary=[]
  for o in obs:
   for key in ['silhouette_affinity','calinski_harabasz']:
    v=d[d.k==o['k']][key];summary.append(dict(k=o['k'],metric=key,observed=o[key],null_mean=v.mean(),p=(1+(v>=o[key]).sum())/1001))
  pd.DataFrame(summary).to_csv(R/f'{stage}_null_summary.csv',index=False)
 pd.DataFrame(baseline).to_csv(R/'combined_modality_baselines.csv',index=False)
 (m.ROOT/'configuration.json').write_text(json.dumps(dict(affinity_K=10,standalone_mu=.5,fusion_input_mu=.3,fusion_K=20,t=20,alpha=1,spectral_seed=42,n_init=10,permutations_per_cohort=1000,gene_representation='Original measured gene expression',metric='snf.metrics.silhouette_score'),indent=2))
 print('DONE',flush=True)
if __name__=='__main__':main()
