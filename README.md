This repository contains the analysis code for Muller et al, "A meta-analysis study of the robustness of gut microbiome-metabolome associations".

Citation: Muller, E., Algavi, Y.M. & Borenstein, E. A meta-analysis study of the robustness and universality of gut microbiome-metabolome associations. Microbiome 9, 203 (2021). https://doi.org/10.1186/s40168-021-01149-z

### Abstract  

*Background:* Microbiome-metabolome studies of the human gut have been gaining popularity in recent years, mostly due to accumulating evidence of the interplay between gut microbes, metabolites, and host health. Statistical and machine learning-based methods have been widely applied to analyze such paired microbiome-metabolome data, in the hope of identifying metabolites that are governed by the composition of the microbiome. Such metabolites can be likely modulated by microbiome-based interventions, offering a route for promoting gut metabolic health. Yet, to date, it remains unclear whether findings of microbially-associated metabolites in any single study carry over to other studies or cohorts, and how universal are microbiome-metabolites links.     

*Results:* In this study, we addressed this challenge by performing a comprehensive meta-analysis to identify human gut metabolites that are ‘universally’ well-predicted by the gut microbiota. To this end, we processed data from 1733 samples from 10 independent human gut microbiome-metabolome studies, focusing initially on healthy subjects, and implemented a machine learning pipeline to predict metabolite levels in each dataset based on the composition of the microbiome. Comparing the predictability of each metabolite across datasets, we found 97 universally well-predicted metabolites. These include metabolites involved in important microbial pathways such as bile acid transformations and polyamines metabolism. Importantly, however, other metabolites exhibited large variation in predictability across datasets, suggesting a cohort- or study-specific relationship between the microbiome and the metabolite. Comparing taxonomic contributors to different models, we found that some universally well-predicted metabolites were predicted by markedly different sets of taxa across datasets, suggesting that some microbially associated metabolites may be governed by different members of the microbiome in different cohorts. We finally examined whether models trained on a control group of a given study successfully predicted the metabolite’s level in the disease group of the same study, identifying several metabolites where the model was not transferable, indicating a shift in microbial metabolism in disease-associated dysbiosis.  

*Conclusions:* Combined, our findings provide a better understanding of the link between the microbiome and metabolites and allow researchers to put identified microbially-associated metabolites within the context of other studies.

### Scripts and files overview

The scripts and files included in this repository are detailed below.

```
root
|-- code                                # All the scripts used for the analysis
|   |-- Analysis_Notebook.Rmd           # Main analysis script, includes all manuscript results and plots
|   |-- Custom_Formatting.css           # Notebook formatting
|   |-- Utilities.R                     # Utility functions for miscellaneous tasks
|   |-- Machine_Learning_Pipeline_For_Condor.R  # Machine learning pipeline (meant to be used with HTCondor or another distributed job management system)
|   |-- Subgroup_Meta-Analysis_Mixed_Effects.R  # Functions for subgroup meta-analysis
|-- data                                # All the data used for the anlaysis including intermediate data
|   |-- Processed_Data_29092020.RData   # Processed microbiome-metabolome datasets from 10 studies, loaded to R
|   |-- Combined_Metabolome_Predictability_Results.zip                             # Intermediate results from various machine learning pipelines, per dataset per metabolite (should be unzipped before running)
|   |-- REM_Results.tsv                 # A processed table with all REM model results
|   |-- Full_REM_Results.RData          # Raw REM model results
|   |-- MelonnPan_supp_table_3_withHMDB.xlsx   # Table related to Mallick et al. 2019
|   |-- Zierer_et_al_2018               # Files related to Zierer et al. 2018 
|   |-- Subgroup_REM_Results.tsv        # A processed table with all subgroup-REM model results
|   |-- Full_Subgroup_REM_Results.RData # Raw subgroup-REM model results
|   |-- Permutation_Feature_Importances.tsv # Raw feature-importance analysis results
|   |-- Pairwise_Contrib_Results.tsv    # Raw feature-importance analysis results for pairwise studies comparisons
|   |-- Pairwise_Cross_Pred_Results.tsv # Raw cross-predictability results for pairwise studies comparisons
|   |-- Pairwise_Cross_Pred_Results_Healthy_Disease.tsv # Raw cross-predictability results for healthy vs. disease analysis
|-- README.md
```

