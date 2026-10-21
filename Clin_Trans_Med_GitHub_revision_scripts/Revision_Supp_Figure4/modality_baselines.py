from pathlib import Path
import importlib.util
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans,AgglomerativeClustering
from sklearn.metrics import silhouette_score
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('analysis',ROOT/'analyze.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
combined=pd.read_csv(ROOT/'results/combined_modality_baselines.csv');rows=combined[['stage','method','k','native_silhouette','native_metric']].assign(modality='Fused').to_dict('records')
for stage in ['EMCI','LMCI']:
 gene,metab,ids,_=m.load(stage)
 for modality,x in [('Gene expression',gene),('Metabolomics',metab)]:
  aff=m.snf.make_affinity(x,metric='euclidean',K=10,mu=.5);scaled=StandardScaler().fit_transform(x);pca=PCA(n_components=.8,svd_solver='full',random_state=42);z=pca.fit_transform(scaled)
  for method in ['SNF + spectral clustering','PCA + k-means','PCA + Ward hierarchical clustering']:
   for k in range(2,11):
    lab=m.labels(aff,k) if method.startswith('SNF') else (KMeans(n_clusters=k,n_init=10,random_state=42).fit_predict(z) if 'k-means' in method else AgglomerativeClustering(n_clusters=k,linkage='ward').fit_predict(z))
    score=m.snf.metrics.silhouette_score(aff,lab) if method.startswith('SNF') else silhouette_score(z,lab)
    rows.append(dict(stage=stage,method=method,modality=modality,k=k,native_silhouette=score,native_metric='Affinity' if method.startswith('SNF') else 'PCA Euclidean'))
  print(stage,modality,'done',flush=True)
d=pd.DataFrame(rows);d['affinity_K']=10;d['standalone_mu']=.5;d['fusion_input_mu']=.3;d['fusion_K']=20;d['fusion_iterations']=20;d['spectral_seed']=42;d['spectral_n_init']=10;d.to_csv(ROOT/'results/all_modality_baselines.csv',index=False);print(d.query('k==2').to_string(index=False))
