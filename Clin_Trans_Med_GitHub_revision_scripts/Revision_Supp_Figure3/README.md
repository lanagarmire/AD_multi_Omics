# Revision_Supp_Figure3

Final full-data and subsampling panels A–I only, using fixed archived metrics. No encoder training, dimension ranking, search or parameter selection. Panels J–K belong to the previous revision and are intentionally excluded.

Run from this directory with Python 3.11 or R as applicable:

```bash
python plot_stability.py
```

Required local inputs (not distributed): inputs/figure_metrics_raw20032.csv and inputs/subsampling80_summary.csv. Previous-revision J–K heatmap panels are intentionally excluded.

Only scripts are distributed. Input and generated output files must remain outside version control. Fixed analysis settings are retained; no tuning or parameter search is run.

The layout follows compiled S3 A–I: three columns, full-data metrics above EMCI and LMCI subsampling rows. G–I are labeled LMCI; the compiled raster incorrectly repeats EMCI. Output stem: Supplementary_Figure3_A_I.
