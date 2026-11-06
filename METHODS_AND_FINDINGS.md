# Methods & preliminary findings

## Research question
How are plasma protein expression measurements associated with time since the first positive COVID-19 swab among dialysis patients? The objective is *exploratory association*, not causal inference or clinical prediction.

## Design and sources
Secondary analysis of real, published plasma Olink NPX measurements and linked clinical metadata. Patients contribute repeated samples; their measurements are not independent. `SampleID` maps samples to metadata, `Individual_ID` identifies participants, `Panel`+`Assay` identifies panel-specific assay targets.

## Quality control findings
- 307 unique sample identifiers join to the metadata records without an observed mismatch.
- 436 assay definitions (90 CM, 87 CVD2, 87 CVD3, 87 IR, 85 In).
- 275 absent sample/assay records relative to a reference set built from observed assays, mostly from 3 missing whole panels (two missing CVD2, one missing IR).
- Null-value inspection alone would not have revealed structurally absent assay records.
- The published NPX files are already processed and normalized. We did not impute or renormalize.

## Exploratory sample-level PCA
Assays complete across all 307 samples are selected to avoid imputation. Standardize assay features (mean zero, unit SD); run 2-component PCA. Plate and case-status groupings are overlaid as descriptive covariates. PCA does not isolate technical batch effects or prove causal effects.

## Full assay screen
The pipeline fits per-assay `NPX ~ time_days + C(Plate_ID)` using OLS with patient-clustered robust standard errors. It uses COVID-positive records with known time since first positive swab. BH FDR is calculated across all fitted assays. This is a high-throughput prioritization analysis, not the final longitudinal model. Patient-level fixed/mixed effects, time nonlinearity, follow-up selection, and robustness of hits need closer study.

## Initial output summary
See `outputs/run_summary.csv` and `outputs/all_assay_longitudinal_cluster_robust.csv`. In this exploratory screen, 64 of 436 assay time coefficients had BH-adjusted p<0.05. This is a statistical screen, not proof that 64 proteins are biological biomarkers.

## Why model choice matters: IL6
A previous patient random-intercept model estimated ~-0.066 NPX/day (p<0.001). A separate patient-cluster-robust OLS model in the pipeline estimated ~-0.0435 NPX/day (nominal p=0.045, BH FDR=0.166). These are **different models** and differ in how they use between-person and within-person information and uncertainty. Their disagreement is a reason to perform sensitivity analysis and not overstate an IL6 result.

## Open scientific tasks
- Nonlinear and within-patient time modeling for shortlisted assays.
- Clinical severity associations with participant-level and repeated-sample design.
- Protein functional annotations/pathway enrichment against measured assay background.
- Figure captions, biological interpretation, and formal limitations for a paper-style report.
- If desired, a separate real raw-data/genomics exercise to support further RTI competencies.
