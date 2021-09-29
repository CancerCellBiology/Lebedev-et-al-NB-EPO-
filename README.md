# Lebedev-et-al-NB-EPO-
Codes from "Growth factor signaling predicts therapy resistance mechanisms and defines neuroblastoma subtypes" Oncogene 2021
https://www.nature.com/articles/s41388-021-02018-7
This repository contains codes used for UMAP-based clustering and differential genes analysis. Exapmle files for each code are also provided.
Requirements:
Python 3.7
Codes were written and optimized for use with Spider 4.1
Python libraries:
pandas
numpy
statsmodels
scipy
umap
hdbscan
matplotlib
seaborn

The Logistic Regression contains code and data files to replicate the results.
Requirements:
pandas
sklearn
numpy
matplotlib

Download codes and excel files to a working folder. Set working folder for each code.
1. “Logistic Regression” code uses dataset (Cangelosi) with known outcomes and gene expression to create predictive models based on expression of certain genes (Diff genes for example). Models for MYCN amplified and MYCN non-amplified are created separately (use status parameter). Code selects top 20% genes for a final model, this parameter can be adjusted in the code. The code will create file with gene coefficients in a predictive model, which can be used for further analysis. This code also yields precession recall graph for created model using test part of the dataset.
2. “Predict outcome” code uses models created during the first step to predict outcomes for other datasets. Gene expression in datasets should be normalized by dividing gene expression for each sample by mean gene expression in the dataset (provided datasets are already normalized). Initial “score” parameter should be set equal to model intercept parameter from the first code. As a result this code yields data files with calculated event (death) probabilities.
3. “Plot ROC curves” code plots ROC curves and calculates AUC values for each datasets. Data files created in previous step should be used as input. To calculate ROC and AUC you need datasets with known outcomes.


