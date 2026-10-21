# Revision_Supp_Figure8

Final differential-feature overlap across baseline, APOE, education and joint adjustment models, retaining the published covariates and thresholds.

Run from this directory with Python 3.11 or R as applicable:

```bash
Rscript 01_baseline_APOE.R
python 02_prepare_education.py
Rscript 03_adjust_confounders.R
python 04_plot_overlap.py
```

Required local inputs (not distributed): inputs/ADNIMERGE.csv, cn_emci_lmci_ad_clin_pseudotime.csv, mrna20k_{emci,lmci}.csv, all_{emci,lmci}_metab_hmdb_lilikoi.csv and mrna_{emci,lmci}_deconv.tsv. R requires limma.

Only scripts are distributed. Input and generated output files must remain outside version control. Fixed analysis settings are retained; no tuning or parameter search is run.

Final design: four panels A–D (EMCI genes, LMCI genes, EMCI metabolites, LMCI metabolites), bold panel titles, and blue joint-adjustment retention annotations inside each plot, matching the final compilation. Output stem: Supplementary_Figure8_APOE4_education_adj_upset.
