# Revision_Table_1

Final characteristics of the paired-omics cohort.

Run from this directory with Python 3.11 or R as applicable:

```bash
Rscript generate_table.R
python export_word.py
```

Required local inputs (not distributed): inputs/cn_emci_lmci_ad_clin_pseudotime.csv, inputs/ADNIMERGE.csv and inputs/681p_172f.csv. Python Word export uses the standard library.

Only scripts are distributed. Input and generated output files must remain outside version control. Fixed analysis settings are retained; no tuning or parameter search is run.
