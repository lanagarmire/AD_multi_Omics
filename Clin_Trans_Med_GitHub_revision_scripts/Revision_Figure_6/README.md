# Revision_Figure_6

Final 149-participant EFIGA preprocessing, fixed label transfer, native Monocle3 trajectory, differential metabolites, GAMs, pathways and Figure 6. Historical transfer calibration and graph-setting searches are excluded.

Run from this directory with Python 3.11 or R as applicable:

```bash
bash run.sh
```

Required local inputs (not distributed): Privately provide inputs/ matrices, patient roster, mapping and fixed protocol; reference/ numerical reference files used by reproducibility assertions; resources/metaboanalyst/{EMCI,LMCI}/{compound_db.qs,hsa.qs}. These files are not shipped. The scripts retain the frozen verification contracts of the final pipeline. Use PYTHON_BIN, RSCRIPT_BIN, TRAJECTORY_RSCRIPT and PATHWAY_RSCRIPT to select compatible environments. R packages include monocle3, uwot, limma, mgcv and MetaboAnalystR; see code imports for the full list.

Only scripts are distributed. Input and generated output files must remain outside version control. Fixed analysis settings are retained; no tuning or parameter search is run.
