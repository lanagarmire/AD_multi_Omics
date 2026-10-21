# Revision scripts for GitHub

Script-only submission export. No participant data, input matrices, cached results, figures, notebook outputs, execution logs or environments are included.

Numbering follows the final compiled Supplementary Figure 1–11 legends, including education-adjusted Supplementary Figure 10 and UMAP Supplementary Figure 11. Run commands and required private inputs are listed in each folder. Run each analysis from its folder. Python 3.11 is recommended. Install the packages imported by the relevant scripts in compatible environments; this export does not provide a tested environment lockfile.

## Order

- `Revision_Figure_5`: Panels A–E: final MRI model, education-adjusted fixed clinical trajectory and conditional permutation test, and 90-month dementia-free survival. The permutation inference remains conditional on the selected trajectory; earlier selection is not accounted for.
- `Revision_Figure_6`: Final 149-participant EFIGA preprocessing, fixed label transfer, native Monocle3 trajectory, differential metabolites, GAMs, pathways and Figure 6. Historical transfer calibration and graph-setting searches are excluded.
- `Revision_Supp_Figure1`: Final integrated participant-flow layout. Reads audited cohort counts without rerunning raw-data audits.
- `Revision_Supp_Figure2`: Final metabolite transformation visualization. No raw data, transformation snapshots or output plots are included.
- `Revision_Supp_Figure3`: Final full-data and subsampling panels A–I only, using fixed archived metrics. No encoder training, dimension ranking, search or parameter selection. Panels J–K belong to the previous revision and are intentionally excluded.
- `Revision_Supp_Figure4`: Final 1,000-replicate feature-permutation tests and modality benchmarks. The 2–10-cluster curves are retained as the published validity analysis; they do not tune or replace the final subtype model.
- `Revision_Supp_Figure5`: Final modality-separation plot from fixed final scores only. Encoder fitting, dimension search and parameter selection are excluded.
- `Revision_Supp_Figure6`: Final joint-clustering heatmap from the fixed final contingency table; parameter comparison cells excluded.
- `Revision_Supp_Figure8`: Final differential-feature overlap across baseline, APOE, education and joint adjustment models, retaining the published covariates and thresholds.
- `Revision_Supp_Figure10`: Final age-, sex- and education-adjusted cognitive/functional distributions.
- `Revision_Supp_Figure11`: Final UMAP comparison using archived coordinates and final Figure 6 labels. No UMAP or trajectory settings are optimized.
- `Revision_Table_1`: Final characteristics of the paired-omics cohort.

## Coverage and exclusions

The supplied analysis folder contains final revision code for main Figures 5–6, Supplementary Figures 1–6, 8, 10–11 and Table 1. Supplementary Figure 3 panels J–K and Supplementary Figures 7 and 9 belong to the previous revision and are intentionally excluded from this revision-only export. Other unchanged main figures/tables are outside this export. Supplementary Figures 3, 5 and 6 are plot-only exports of the final fixed metrics; their encoder/parameter comparison fitting workflows are excluded.

Excluded: autoencoder training/retraining, transfer-threshold searches, graph/hyperparameter sweeps, cutoff searches, exploratory MRI models, superseded pipelines, provenance archives, execution/display notebooks, report assembly, and verification-only utilities. Published fixed settings and statistical validity tests remain. Removing search code does not make the historical choices prespecified or change the inferential limitations in the source analysis.

## Inputs and outputs

Place required data privately at the portable relative paths documented in each folder. Input files are deliberately not included. Figure 6 additionally requires frozen numerical references and pathway resources described in its README. Output-writing statements remain so that users with authorized inputs can rerun the analyses. Scientific pipelines were not rerun end to end. Selected plotting/final-model scripts were smoke-checked in temporary directories; no generated artifacts are included.

Validation covers Python syntax, R syntax, shell syntax, privacy scans and final-script dependency review; it does not establish numerical reproduction without the private inputs and compatible scientific environments.

## Final supplementary design reference

The design reference is the package's `figures/Supp_fig_final_compiles/` collection. Numbering is preserved, with no renumbering of the intentional exclusions. S1/S2/S4/S5/S10/S11 retain their final source layouts. S3 exports the A–I three-column layout with the final metrics, stage colors and 172-dimensional highlight; the third-row stage label is correctly LMCI (the compiled raster repeats EMCI there). S6 uses the final four-by-four heatmap layout. S8 uses the compiled four-panel UpSet design with bold A–D labels and blue retention annotations. S3/S8 styling is reconstructed in plotting code from the compiled design; it is not a claim of pixel-identical raster export across environments.
