import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from unioncom import UnionCom
import os


intersect_adni = pd.read_csv("extdata/intersect_adni_metab.csv", index_col=0)
intersect_adni['Label'] = intersect_adni['Label'].replace({
    'subtype1': 'LMCI1',
    'subtype2': 'LMCI2'
 })
LMCI_conditions = intersect_adni['Label'].isin(['LMCI1', 'LMCI2'])
intersect_adni_LMCI = intersect_adni[LMCI_conditions]
labels_ADNI = intersect_adni_LMCI['Label']
intersect_adni_LMCI.drop(columns=['Label'], inplace=True)

intersect_EFIGA_mci = pd.read_csv("extdata/intersect_EFIGA_mci_metab.csv", index_col=0)

EFIGA_LMCI_samples = intersect_EFIGA_mci.loc[
    intersect_EFIGA_mci['big label'] == 'LMCI'
]
data2_LMCI = EFIGA_LMCI_samples.to_numpy()
data2_LMCI = data2_LMCI[:, 1:]

min_max_scaler = MinMaxScaler()
norm_data1 = min_max_scaler.fit_transform(intersect_adni_LMCI.to_numpy())
norm_data2 = min_max_scaler.fit_transform(data2_LMCI)

uc = UnionCom.UnionCom(output_dim=11, epoch_pd=2000, epoch_DNN=0, distance_mode='l2')
integrated_data = uc.fit_transform(dataset=[norm_data1, norm_data2])
cor_pairs = uc.match(dataset=[norm_data1, norm_data2])
norm_pat_match = cor_pairs[0]


LMCI1_ID = list(np.where(np.isin(labels_ADNI, 'LMCI1'))[0])
LMCI2_ID = list(np.where(np.isin(labels_ADNI, 'LMCI2'))[0])


matched_LMCI1, matched_LMCI2 = [], []
thre = 0.995
ties = 0

for index in range(norm_pat_match.shape[0]):
    match_vec = norm_pat_match[index]
    if max(LMCI1_ID) >= len(match_vec) or max(LMCI2_ID) >= len(match_vec):
        continue

    max_score_LMCI1 = max(match_vec[LMCI1_ID])
    max_score_LMCI2 = max(match_vec[LMCI2_ID])

    if abs(max_score_LMCI1 - max_score_LMCI2) < 1e-6:
        ties += 1
        print(f"Tie detected for sample index {index}.")

    if max_score_LMCI1 >= max_score_LMCI2:
        matched_LMCI1.append(index)
    else:
        matched_LMCI2.append(index)


flat_list_LMCI1 = list(set(sorted(matched_LMCI1)))
flat_list_LMCI2 = list(set(sorted(matched_LMCI2)))

intersect_EFIGA_mci['labels_matching_LMCI'] = 'Unmatched'
efiga_lmci_indices = intersect_EFIGA_mci.loc[intersect_EFIGA_mci['big label'] == 'LMCI'].index

for i, idx in enumerate(efiga_lmci_indices):
    if i in flat_list_LMCI1:
        intersect_EFIGA_mci.at[idx, 'labels_matching_LMCI'] = 'LMCI1'
    elif i in flat_list_LMCI2:
        intersect_EFIGA_mci.at[idx, 'labels_matching_LMCI'] = 'LMCI2'

unmatched_indices = intersect_EFIGA_mci.loc[
    (intersect_EFIGA_mci['big label'] == 'LMCI') & 
    (intersect_EFIGA_mci['labels_matching_LMCI'] == 'Unmatched')
].index

for idx in unmatched_indices:
    col_idx = list(efiga_lmci_indices).index(idx)
    similarity_scores = norm_pat_match[:, col_idx]
    highest_similarity_index = np.argmax(similarity_scores)
    highest_similarity_label = labels_ADNI.iloc[highest_similarity_index]
    intersect_EFIGA_mci.at[idx, 'labels_matching_LMCI'] = highest_similarity_label

final_unmatched_count = intersect_EFIGA_mci.loc[
    (intersect_EFIGA_mci['big label'] == 'LMCI') & 
    (intersect_EFIGA_mci['labels_matching_LMCI'] == 'Unmatched')
].shape[0]


intersect_EFIGA_mci.to_csv('extdata/labels_matching_uc_refined.csv')
