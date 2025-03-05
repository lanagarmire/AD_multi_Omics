import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier
from sklearn.model_selection import GridSearchCV
from unioncom import UnionCom

intersect_adni = pd.read_csv("extdata/intersect_adni_metab.csv",index_col=0) 
intersect_EFIGA_mci = pd.read_csv("extdata/intersect_EFIGA_mci_metab.csv", index_col=0)  


adni_mrna_merge = pd.read_csv("extdata/adni_mrna_merge.csv",index_col=0)


conditions = [
    (adni_mrna_merge['subtype'] == "EMCI1"),
    (adni_mrna_merge['subtype'] == "EMCI2") ,
    (adni_mrna_merge['subtype'] == "LMCI1") ,
    (adni_mrna_merge['subtype'] == "LMCI2")
    ]

values = ['EMCI', 'EMCI', 'LMCI', 'LMCI']

adni_mrna_merge['subtype_big'] = np.select(conditions, values)

adni_mrna_merge = adni_mrna_merge.sort_values(by=['Unnamed: 0'])
labels = adni_mrna_merge['subtype_big']

intersect_adni = intersect_adni.sort_index()
data1 = intersect_adni.to_numpy()
data2 = intersect_EFIGA_mci.to_numpy()

#normalization
from sklearn import preprocessing
min_max_scaler = preprocessing.MinMaxScaler()
ss_scaler = preprocessing.StandardScaler()


norm_data1 = min_max_scaler.fit_transform(data1) ### adni
norm_data2 = min_max_scaler.fit_transform(data2) # EFIGA

# Initialize UnionCom

uc = UnionCom.UnionCom(output_dim=11,epoch_pd=2000,epoch_DNN=0,distance_mode='geodesic')


#computer and extracts the match results from the cor_pairs.
integrated_data = uc.fit_transform(dataset=[norm_data1,norm_data2]) 
part_1 = integrated_data[0]
part_2 = integrated_data[1]


cor_pairs = uc.match(dataset=[norm_data1,norm_data2])
norm_pat_match = cor_pairs[0]

def get_indices(lst, targets):
    indices = []
    for target in targets:
        if target in lst:
            indices.append(lst.index(target))
    return indices


EMCI_ID = list(np.where(np.isin(labels,'EMCI'))[0])
LMCI_ID = list(np.where(np.isin(labels,'LMCI'))[0])
clus_in = EMCI_ID


matched = []   


thre = 0.995


for index in clus_in:

    match_vec = norm_pat_match[index]
# similarity scores in match_vec are greater than a threshold value
    matched_index = np.where(match_vec > thre*max(match_vec))
    matched.append(matched_index[0])

    flat_list = [item for sublist in matched for item in sublist]

new_flat_list = list(set(sorted(flat_list)))






### HERE THE NUMBER SHOULD BE THE NUMBER OF SAMPLES
UnionCom_predict =  ['LMCI'] * 10
for index in new_flat_list:
    UnionCom_predict[index] = 'EMCI'
    
intersect_EFIGA_mci['labels_matching_UC_big'] = UnionCom_predict
intersect_EFIGA_mci.to_csv('extdata/labels_matching_uc.csv')