# Revision_Supp_Figure10

Final age-, sex- and education-adjusted cognitive/functional distributions.

Run from this directory with Python 3.11 or R as applicable:

```bash
python adjust_clinical_scores.py
```

Required local inputs (not distributed): inputs/cn_emci_lmci_ad_clin_pseudotime.csv. Python requires numpy, pandas, scipy, statsmodels and matplotlib.

Only scripts are distributed. Input and generated output files must remain outside version control. Fixed analysis settings are retained; no tuning or parameter search is run.
