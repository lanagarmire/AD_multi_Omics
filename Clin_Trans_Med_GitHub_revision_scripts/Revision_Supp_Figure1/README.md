# Revision_Supp_Figure1

Final integrated participant-flow layout. Reads audited cohort counts without rerunning raw-data audits.

Run from this directory with Python 3.11 or R as applicable:

```bash
python plot_participant_flow.py
```

Required local inputs (not distributed): inputs/preprocessing_counts.csv and inputs/patient_group_summary.csv from the finalized cohort audit.

Only scripts are distributed. Input and generated output files must remain outside version control. Fixed analysis settings are retained; no tuning or parameter search is run.
