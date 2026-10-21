# Revision_Supp_Figure11

Final UMAP comparison using archived coordinates and final Figure 6 labels. No UMAP or trajectory settings are optimized.

Run from this directory with Python 3.11 or R as applicable:

```bash
python plot_umap.py
```

Required local inputs (not distributed): inputs/patients.csv (final Figure 6 cohort); inputs/ADNI/{ADNI_umap.csv,ADNI_metadata.csv,ADNI_tree_nodes_UMAP.csv,ADNI_tree_nodes_PCA.csv,ADNI_tree_edges.csv,ADNI_trajectory_config.json}; inputs/EFIGA/{SELECTED_full_cohort_UMAP_coordinates.csv,SELECTED_labels.csv}.

Only scripts are distributed. Input and generated output files must remain outside version control. Fixed analysis settings are retained; no tuning or parameter search is run.
