# ASO Model Metabolomics

## Overview

This repository contains the data and code used for the analysis presented in our study: https://www.biorxiv.org/content/10.1101/2024.06.07.597975v1.abstract on the metabolomics of ASO mouse models with SPF and germ-free microbiomes.

## Data and Code Availability

All raw data are available on Zenodo: [Zenodo Link]([https://zenodo.org/record/...](https://zenodo.org/records/10841477))

### Data Availability

The Zenodo deposit includes:
- Raw metabolomics data
- Processed data files
- Analysis scripts
- Supplementary materials


## Repository Structure

```md
The repository is organized as follows:

ASO_model_metabolomics/
│
#################################################################################
#######   Code for processing raw and performing core statistical tests   #######

├── Sia - github/
│   ├── 5370_Gut_brain_axis_study_Sub1_Pilot_study_Report_2020-04-27.pdf
│   ├── Cal_Tech_PD_Mice/
│   │   ├── 5370_Gut_brain_axis_study_Sub1_Pilot_study_Data_2020-04-27.xlsx
│   │   ├── 5370_Gut_brain_axis_study_Sub1_Pilot_study_Report_2020-04-27.pdf
│   │   ├── MxP-Quant-500-kit_MetaboINDICATOR™_V01-2019.xlsx
│   │   └── Sample_reconciliation_messages/
│   │       ├── RE_Sample_reception_(BC_No._5370_Duke).msg
│   │       ├── Sample_information_From_PD_Mice_2019_LHM_AKa_DSo.xlsx
│   │       ├── Sample_information_From_PD_Mice_2019_LHM_original.xlsx
│   │       └── Sample_overview_complete.xlsx
│   └── ...
│
├── data/
│   ├── 5370_Cal_Tech_PD_Mouse_Model_Gut_brain_axis_study_09112020.xlsx
│   ├── Metadata_Metabolomics_2019_checked_2020_updated.xlsx
│   ├── Metadata_Metabolomics_2019.xlsx
│   ├── normalized/
│   │   ├── 5370_Gut_brain_axis_study_Sub1_Pilot_study_Data_2020-04-27.xlsx
│   └── unnormalized/
│       ├── 5370_Gut_brain_axis_study_Sub1_not_normalized_2020-11-04.xlsx
│       ├── 5370_Gut_brain_axis_study_Sub1_not_normalized_2020-11-05.xlsx
│   └── ...
│
├── figures/
│   ├── August_2021/ ...
│
├── results/
│   ├── caltech.PD_mice.122622.xlsx
│
├── scripts/
│   ├── Cal_Tech_PD_Mice.q500.analysis.R
│   ├── Cal_Tech_PD_Mice.q500.figures.R
│   ├── Cal_Tech_PD_Mice.q500.preprocess.R
│   ├── Cal_Tech_PD_Mice.q500.regression.R
│   └── impute.knn.obs.sel.R


#################################################################################
##############   Code for interpretation and figure generation   ##############

├── ASO_model_metabolomics.Rproj
│
├── src/
│   ├── analyte_class_correlations.R
│   ├── analyte_class_enrichment.R
│   ├── body-weight-correlation-analysis.html
│   ├── body-weight-correlation-analysis.qmd
│   ├── body-weight-correlation-analysis.R
│   ├── boxplots.R
│   ├── correlations_vis.R
│   ├── correlations.R
│   ├── functions.R
│   ├── interaction_effect_plots.R
│   ├── load_packages.R
│   └── tissue_shared_analytes.R
│
├── figures/
│   ├── Analyte_class_estimates/
│   ├── Correlations/
│   ├── Heatmaps/
│   ├── Interaction_plots/
│   ├── Venn_Diagrams/
│   └── ...
│
├── files/
│   ├── 2023-01-22_body-weight-correlations.rds
│   ├── 2023-01-24_lm-comparisons-BWQuant.rds
│   ├── 5370_Cal_Tech_PD_Mouse_Model_Gut_brain_axis_study_09112020.xlsx
│   ├── analyte_class_colors_grouped.rds
│   ├── analyte_class_colors.RData
│   ├── caltech.PD_mice_G_by_M_results.xlsx
│   ├── caltech.PD_mice.022223.xlsx
│   ├── caltechdata.csv
│   ├── Metadata_Metabolomics_2019_checked_2020_updated.xlsx
│   ├── plasma_correlation_data.rds
│   └── ...
│
├── renv/ ... # for reproducibility 
├── README.md  # This file
└── LICENSE
```

## Data Dictionary

### Raw Data

- **5370_Gut brain axis study_Sub1_Pilot study_gut content data_not normalized_2020-11-05.xlsx**:
  - Raw data for gut samples.
- **5370_Gut brain axis study_Sub1_Pilot study_Data_not normalized_2020-11-04.xlsx**:
  - Raw data for all other samples (Plasma & brain)

### Processed Data

- **5370_Gut_brain_axis_study_Sub1_Pilot_study_Data_2020-04-27.xlsx**:
  - Contains aggregated data for the gut-brain axis study, sub1 pilot study conducted on 2020-04-27.
- **caltechdata.csv**:
  - Imputated and unimputed abundances for all metabolites and samples.

### Results and Statistical Reports

- **caltech.PD_mice.022223.xlsx**:
  - Contains statistical results for linear models, ANOVAs, t-tests, and conditional testing.
- **plasma_correlation_data.rds**:
  - Correlation results for each metabolite between plasma and all other tissue levels.
- **2023-01-22_body-weight-correlations.rds**:
  - Mouse body weight correlation results.
- **2023-01-24_lm-comparisons-BWQuant.rds**:
  - Contains linear model results for body weight quantification.


### Supplementary Data

- **Metadata_Metabolomics_2019_checked_2020_updated.xlsx**:
  - Metadata for samples and mice.

### How to cite this material
DOI:10.5281/zenodo.10841427
