# Revision_Figure_5

Panels A–E: final MRI model, education-adjusted fixed clinical trajectory and conditional permutation test, and 90-month dementia-free survival. The permutation inference remains conditional on the selected trajectory; earlier selection is not accounted for.

Run from this directory with Python 3.11 or R as applicable:

```bash
python 01_mri.py
python 02_clinical_trajectory.py
python 03_survival.py
```

Required local inputs (not distributed): inputs/cn_emci_lmci_ad_clin_pseudotime.csv; inputs/ADNIMERGE.csv; inputs/extdata/ subtype membership files used by survival; ClinTrajan helper modules clintraj_eltree.py and clintraj_qi.py and their dependencies installed or placed in clintra/.

Only scripts are distributed. Input and generated output files must remain outside version control. Fixed analysis settings are retained; no tuning or parameter search is run.
