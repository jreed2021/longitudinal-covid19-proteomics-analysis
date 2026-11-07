# Methods and Preliminary Findings

**Analysis author:** Janessa Reed, MS  
**Project:** Independent, AI-assisted secondary analysis of published COVID-19 proteomics data  
**Original research:** Gisby et al. (2021), *Longitudinal proteomic profiling of dialysis patients with COVID-19 reveals markers of severity and predictors of death*. [eLife article](https://doi.org/10.7554/eLife.64827) · [Original research code](https://github.com/jackgisby/longitudinal_olink_proteomics)

## Research question and scope

Among COVID-19-positive dialysis patients, how is relative plasma protein expression (Olink NPX) associated with days since the first positive COVID-19 swab? This is an exploratory secondary analysis, not a reproduction of the original investigators' statistical models, a causal analysis, or a clinical biomarker validation study.

The authoritative workflow is [`notebooks/01_Guided_COVID19_Proteomics_Analysis.ipynb`](notebooks/01_Guided_COVID19_Proteomics_Analysis.ipynb). The notebook is intended to be run top to bottom using the two original plasma CSV files, obtained separately from the source study. Tables and figures shown in the executed notebook constitute the published project outputs; any local `outputs/` and `figures/` files are regenerated when it runs.

## Study data and observational units

The supplied processed datasets contain **133,577 assay measurements**, **307 sample identifiers**, and **105 unique individual identifiers**. The study's plasma cohort includes **256 COVID-positive samples from 55 patients** and **51 control samples**. There are **436 observed panel–assay definitions** across five Olink panels.

Each measurement is associated with a **sample** (`SampleID`) and a **person** (`Individual_ID`). Repeated observations from the same patient are correlated and must not be treated as independent. Assay identity is handled using `Panel` and `Assay` together.

One patient, **C105**, has negative and positive sample-level status labels. The original research documentation explains sampling before and after the first positive swab; these labels are preserved, not automatically recoded or discarded. `WHO_Severity_Peak` remains constant within each patient's records, as expected for an overall disease-course variable.

## Data integrity and missingness

- All **307** sample identifiers matched their metadata in the notebook's sample reconciliation.
- **275** structurally absent sample–panel–assay combinations were identified relative to the panel–assay definitions observed elsewhere in the file.
- Three sample–panel combinations are entirely absent: two CVD2 panels and one IR panel, accounting for **261** absent assay combinations; the other **14** absent combinations are isolated assay omissions.
- These are **absent records**, not `NaN` values in existing measurements. The notebook does not automatically impute them.
- The source NPX data are already normalized relative log2 protein-expression measurements. The notebook does not perform raw instrument-level QC or normalization.

## Exploratory quality control and PCA

Sample-level median NPX distributions are visualized by laboratory plate within each Olink panel. Because NPX values are assay-relative, differences between panels should not be interpreted as differences in absolute concentrations.

For descriptive principal component analysis (PCA), the notebook constructs a sample-by-assay matrix and retains **261 assay columns observed in all 307 samples**. Assay columns are standardized prior to two-component PCA. **PC1 explains 18.82%** and **PC2 explains 8.59%** of variance (combined **27.41%**). Plate groups overlap considerably in the displayed first two components, but this does **not** establish that plate effects are absent.

This PCA intentionally differs from the original paper's approach, which used missing-value imputation; results are not expected to be identical.

## Exploratory assay-wise longitudinal screening

The notebook screens the COVID-positive observations with known time since first positive swab using an ordinary least squares model for each panel–assay definition:

`NPX ~ time_days + C(Plate_ID)`

Standard errors are clustered by `Individual_ID` to account for within-patient dependence in the uncertainty estimates. Benjamini–Hochberg false-discovery-rate (BH-FDR) adjustment is applied across the **436** fitted assay tests. These regressions do not estimate patient-specific trajectories or fully remove between-patient confounding.

All **436** screens completed; **64** had BH-adjusted p-values below **0.05**. These are **exploratory time associations**, not independently validated biomarkers.

For **IL6**, the assay-wise screening model estimated approximately **−0.0435 NPX units/day** (nominal *p* ≈ 0.045; BH-adjusted *p* ≈ 0.166); therefore IL6 did **not** meet the screening FDR threshold. An earlier, separate random-intercept mixed-effects model estimated approximately **−0.066 NPX units/day** (*p* < 0.001). The models use different assumptions, so their estimates should not be described as interchangeable or as proof of a robust IL6 time trend.

## Relationship to the original publication

Gisby et al. used longitudinal mixed-effects models, time from **symptom onset**, additional clinical covariates, nonlinear time functions, and severity-by-time interactions for their longitudinal research questions. This notebook uses **time from first positive swab**, linear exploratory association screens, and patient-clustered standard errors. Consequently, its count of 64 exploratory associations should **not** be directly compared with the article's counts as if they were the same analysis.

The original publication includes separate severity, survival, and machine-learning analyses that this portfolio notebook does not reproduce.

## Limitations and next scientific work

- The cohort consists of patients with end-stage kidney disease receiving dialysis, so generalizability is limited.
- Samples are unevenly distributed across follow-up time; later time points are sparse.
- Plate allocation and patient composition are related because repeated samples from each participant were placed on the same plate.
- Cluster-robust standard errors do not by themselves account for informative follow-up, nonlinear time trends, or confounding in the estimated coefficient.
- Restricting PCA to complete assay columns excludes partially observed assays and differs from the published preprocessing method.
- Before drawing biological conclusions about shortlisted proteins, examine within-patient models, nonlinear time, severity, sensitivity analyses, and relevant technical QC information.
- No raw-mass-spectrometry processing, GWAS, genotype QC, or multi-omics integration was conducted.

## Attribution and AI disclosure

The original investigators designed the study, collected patient samples, generated proteomic measurements, and published the research. This repository contains an **independent secondary analysis**, not their original study or analysis pipeline. ChatGPT assisted with code development, debugging, methodology explanations, and documentation. Responsibility for describing and checking the published analysis remains with the repository author.
