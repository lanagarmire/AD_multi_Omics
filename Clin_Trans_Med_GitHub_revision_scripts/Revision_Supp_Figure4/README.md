# Revision_Supp_Figure4

Final 1,000-replicate feature-permutation tests and modality benchmarks. The 2–10-cluster curves are retained as the published validity analysis; they do not tune or replace the final subtype model.

Run from this directory with Python 3.11 or R as applicable:

```bash
bash run.sh
```

Required local inputs (not distributed): inputs/mrna20k_{emci,lmci}.csv and inputs/all_{emci,lmci}_metab.csv. Python requires snfpy, scikit-learn, scipy, numpy, pandas, threadpoolctl, matplotlib, seaborn and python-docx.

Only scripts are distributed. Input and generated output files must remain outside version control. Fixed analysis settings are retained; no tuning or parameter search is run.
