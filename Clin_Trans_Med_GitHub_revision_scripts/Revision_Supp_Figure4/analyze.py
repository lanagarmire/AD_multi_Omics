"""Reproduce the supplied EMCI procedure for LMCI and compare combined-modality models."""
import os
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import json,pickle,hashlib,sys,time,multiprocessing
import numpy as np,pandas as pd
import snf
from sklearn.cluster import spectral_clustering,KMeans,AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score,calinski_harabasz_score,davies_bouldin_score,adjusted_rand_score
from scipy.spatial.distance import pdist,squareform
from concurrent.futures import ProcessPoolExecutor
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parent; EXT=ROOT/'inputs'
R=ROOT/'results';FIG=ROOT/'figures'
for p in [R,FIG]:p.mkdir(exist_ok=True)
KS=list(range(2,11)); PARAMS=dict(affinity_K=10,mu=.3,fusion_K=20,t=20,alpha=1.)
G_M=None;G_T=None

def load(stage):
 original=EXT
 m=pd.read_csv(original/f'mrna20k_{stage.lower()}.csv');t=pd.read_csv(original/f'all_{stage.lower()}_metab.csv')
 assert m['index'].is_unique and t.PID.is_unique and set(m['index'])==set(t.PID)
 assert np.array_equal(m['index'].to_numpy(),t.PID.to_numpy()), 'Original row orders must agree'
 a=m.drop(columns=['Unnamed: 0','index']).to_numpy(float);b=t.drop(columns=['PID','Label']).to_numpy(float)
 assert np.isfinite(a).all() and np.isfinite(b).all()
 return a,b,t.PID.to_numpy(),t.Label.to_numpy()

def affinity(m,t):return snf.snf(snf.make_affinity([m,t],metric='euclidean',K=10,mu=.3),K=20,t=20,alpha=1.)
def labels(aff,k):return spectral_clustering(aff,n_clusters=k,random_state=42,n_init=10)
def initialize(m,t):
 global G_M,G_T;G_M=m;G_T=t;threadpool_limits(limits=1)

