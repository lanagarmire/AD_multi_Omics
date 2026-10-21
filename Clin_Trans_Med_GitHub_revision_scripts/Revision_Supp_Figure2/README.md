# Revision_Supp_Figure2

Final metabolite transformation visualization. No raw data, transformation snapshots or output plots are included.

Run from this directory with Python 3.11 or R as applicable:

```bash
python plot_preprocessing.py
```

Required local inputs (not distributed): inputs/{fia,uplc}_{imputed,log2,standardized,clipped}.csv with RID and metabolite columns.

Only scripts are distributed. Input and generated output files must remain outside version control. Fixed analysis settings are retained; no tuning or parameter search is run.
