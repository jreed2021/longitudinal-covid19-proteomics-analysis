# Longitudinal COVID-19 Proteomics: Reproducible Python Secondary Analysis

**Project Author:** Janessa Reed, MS  
**Project Type:** Independent Secondary Analysis of Published Biomedical Research  
**Tools:** Python, Pandas, NumPy, Statsmodels, Scikit-learn, Matplotlib, Seaborn  
**Research Data Source:** Gisby et al. (2021)

## AI Assistance Disclosure

This project was completed with assistance from ChatGPT for Python code development, troubleshooting, statistical methodology explanations, workflow organization, and documentation.

The analysis uses authentic, publicly available proteomics data originally generated and published by Gisby and colleagues. AI assistance did not generate or replace the original biological measurements.

All findings presented are exploratory secondary analyses and should not be interpreted as independently validated clinical discoveries.


## Original data and research attribution

**Data:** Published plasma Olink NPX proteomics and sample metadata, from *Longitudinal proteomic profiling of dialysis patients with COVID-19 reveals markers of severity and predictors of death* (Gisby et al., eLife, 2021). Original research: https://elifesciences.org/articles/64827 . Study authors' repository: https://github.com/jackgisby/longitudinal_olink_proteomics . Cite the original authors and obey the dataset's reuse license.

## Project overview

This repository documents a reproducible, independent educational secondary analysis of published longitudinal COVID-19 proteomics data using the Jupyter Notebook as the single supported analytical workflow. The original investigators conducted the clinical study and produced the Olink measurements; this portfolio project focuses on subsequent computational analysis.

## What this project demonstrates
- Intake, schema and patient/sample ID reconciliation, and structural assay completeness checking.
- Panel/sample QC and observed-expression summaries (no automatically deleted or imputed measurements).
- Exploratory PCA on **assays observed in every sample** (no missing values imputed).
- 436 separate exploratory longitudinal association models with patient-clustered robust standard errors, plate adjustment, and Benjamini–Hochberg FDR correction.
- Explicit documentation of design limitations and the distinction between sample-level and patient-level inference.

## Reproducing the analysis (primary Jupyter Notebook workflow)

The authoritative analysis is `notebooks/01_Guided_COVID19_Proteomics_Analysis.ipynb`. The notebook generates CSV outputs and plots locally when executed; these generated files are not required to be stored in the repository. The analysis should be regenerated from the notebook, not from a separate Python script.

1. Download this repository and open the project folder in VS Code or Jupyter.
2. Use Python 3.10+ and install dependencies: `python -m pip install -r requirements.txt`. Select the corresponding Jupyter kernel.
3. Download the original published plasma data from the authors' source and put **both files locally** in `data/` using the exact filenames `plasma_npx_level.csv` and `plasma_sample_level.csv`. They are not bundled in this repository.
4. Open `notebooks/01_Guided_COVID19_Proteomics_Analysis.ipynb` and select **Run All** (or run its cells sequentially). The working directory should be the project root or the `notebooks/` folder.
5. Review the analysis tables, PCA visualizations, and final summary in the notebook. Running it also generates CSV tables in `outputs/` and figures in `figures/` locally. These generated files do not need to be committed to GitHub. Save the executed notebook if you want the results visible in GitHub’s notebook preview.

The notebook checks source-data structure and patient/sample mappings and does not alter the original input CSV files. Reproducible execution of this exploratory workflow is distinct from independent biological or clinical validation.

## Preliminary results from supplied real files

- 133,577 rows, 307 samples, 105 distinct patients, and 436 observed panel-assay definitions.
- 275 structurally absent sample–assay combinations; three sample-panel combinations entirely absent.
- 261 assay columns are complete across all 307 samples and used in PCA. PC1 explains 18.82% and PC2 explains 8.59% of variance in this standardized subset.
- 436 assay-specific exploratory time models were fitted on the COVID-positive group, using 55 distinct patients. 64 time associations meet exploratory BH-FDR < 0.05.
- **IL6 caution:** the earlier exploratory mixed-effects random-intercept model gave a slope ~-0.066 NPX/day with p<0.001, but the independent patient-clustered OLS screening model in the notebook gives a slope ~-0.0435, nominal p~0.045 and FDR~0.166. They estimate time associations under different models. This difference makes model specification and sensitivity analysis essential; do not claim IL6 is a robust FDR-significant longitudinal discovery.

## Statistical limitations and next validation work

The supplied NPX values are already normalized and cannot demonstrate raw instrument processing. This notebook **does not** do GWAS, genotype QC, genome imputation, QTL mapping, or multi-omics integration. Those require different source datasets and tools. The study population is a specific dialysis cohort; generalization is limited.

The per-assay models are `NPX ~ time_days + C(Plate_ID)` estimated by OLS with patient-clustered robust inference. This is a fast exploratory association screen; it is **not** a substitute for a fully validated longitudinal mixed-effects model. Patient-clustered standard errors help account for within-patient dependence in uncertainty estimates but do not solve informative follow-up, time-varying confounding, nonlinear trends, or between-patient confounding in the estimated slope. A patient fixed-effect or mixed-effects sensitivity model and nonlinear-time comparison should be done for shortlisted assays before research claims.

Plate effects are confounded with patient composition and must not automatically be interpreted as technical bias. PCA is exploratory and is conducted on 261 complete-assay columns only, so the subset may not represent all 436 assays. This project does not infer laboratory-provided limits of detection that were not supplied. No biological causal claims should be made from these associations.

For a fuller study, assess QC metadata from primary authors, time and outcome missingness by patient, assay-specific outliers, severity association in repeated-measures design, and pathway enrichment against suitable background. The already completed per-sample `WHO_Severity_Peak` field is patient-consistent within the positive subset; severity analyses still need design-aware modeling and careful interpretation.

## Folder layout

- `data/`: instructions for obtaining the original public CSV inputs; source CSVs are excluded from this GitHub-ready package.
- `notebooks/01_Guided_COVID19_Proteomics_Analysis.ipynb`: primary, executable analysis with explanations and results.
- `outputs/`: locally generated integrity, QC, PCA, and model tables; not committed to this repository.
- `figures/`: selected quality-control visualizations derived from the published research data. Additional figures may be generated locally when the notebook runs.
- `METHODS_AND_FINDINGS.md`: concise methods and interpretation notes.

## Resume-ready description (accurate scope)

*Developed a reproducible Jupyter Notebook workflow for secondary analysis of published longitudinal Olink proteomics data (133,577 measurements, 307 samples), including metadata reconciliation, panel-level completeness audits, sample QC, PCA and exploratory assay-wise longitudinal association screening with patient-clustered uncertainty estimates and FDR correction.*

## Attribution and scope

This is an **independent educational secondary analysis**, not the original investigators' pipeline and not affiliated with RTI International. The scientific dataset and original study are by Gisby and colleagues. Links to the original repository and DOI are above; credit them when reusing or presenting the dataset. The notebook workflow and summaries here were developed for this portfolio exercise with AI-assisted code development and analytical review. Verify outputs before treating them as scientific conclusions.

The statistical results presented in the Jupyter Notebook are exploratory research findings, not clinical advice or independently validated biomarkers. Patient C105 contributes a pre-positive-swab NEGATIVE sample and subsequent POSITIVE samples; those original sample classifications are preserved, not forcibly made patient-constant. Peak WHO severity is a different patient-level variable.
