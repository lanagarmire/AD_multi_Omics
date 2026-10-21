# Revision_Supp_Figure6

Final joint-clustering heatmap, extracted from the final plotting cell of the linked analysis notebook. Uses the frozen final contingency table; SNF parameter comparisons and trial fits are excluded.

```bash
Rscript plot_joint_clustering.R
```

Required private input: inputs/contingency_table.csv, a numeric 4-by-4 matrix with row and column names EMCI1, EMCI2, LMCI1, LMCI2 in that order. Rows are the final aligned unsupervised groups and columns are the inferred subtypes. Total n=401 and diagonal agreement=287. No table or plot is included. This is a plotting-only export; it does not rerun joint clustering. Requires ggplot2, reshape2 and dplyr.
