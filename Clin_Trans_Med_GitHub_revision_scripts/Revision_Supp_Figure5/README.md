# Revision_Supp_Figure5

Final modality-separation plot from fixed final scores only. Encoder fitting, dimension search and parameter selection are excluded.

Run from this directory with Python 3.11 or R as applicable:

```bash
python plot_modality_contribution.py
```

Required local inputs (not distributed): inputs/silhouette_scores.csv, the final bottleneck sensitivity score table with stage, label_rule, modality, transcriptomic_dimension and silhouette columns.

Only scripts are distributed. Input and generated output files must remain outside version control. Fixed analysis settings are retained; no tuning or parameter search is run.
