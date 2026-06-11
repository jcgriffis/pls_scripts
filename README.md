
This repository contains generic versions of the PLSR/PLSC scripts used for the main analyses in Griffis et al., 2019 (Cell Reports - DOI: 10.1016/j.celrep.2019.07.100)

1. PLSC_Code

The code contained in this folder implements partial least squared correlation (PLSC) analyses. PLSC analyses aim to link two different data tables, and are in this sense similar to canonical correlation analyses (CCA). Specifically, PLSC attempts to identify linear combinations of the variables in each table (i.e. latent variables – LVs) that maximally account for the covariance between the two data tables. A detailed description of the PLSC method can be found in Griffis e t al., (2019 – Cell Reports) and in the papers referenced there.

This implementation of PLSC performs permutation testing to determine the significance of the LVs, and bootstrap ratios are computed to determine the significance of the variable loadings on each LV. These procedures are also described in detail by Griffis et al., (2019 – Cell Reports).

plsc_script.m – this is an input script for running PLSC analyses. Usage is explained in detail in the script header.

2. PLSR_Code

The code contained in this folder implements partial least squared regression (PLSR) analyses. PLSR is a multivariate regression technique that is similar to principal components regression (PCR), but differs in that the predictor matrix is decomposed to identify variance components that maximally account for variance in the dependent variable (i.e. rather than variance in the predictor matrix). A detailed description of the PLSR method can be found in Griffis et al., (2019 – Cell Reports) and in the papers referenced there.

plsr_script_fe.m - This implementation of PLSR fits an explanatory model that utilizes the full dataset (i.e. rather than a predictive model that makes predictions on held-out data). A single leave-one-out cross-validation loop is used to determine the optimal number of components for the final model based on the change in cross-validation model error associated with the addition of each component, as in Griffis et al., (2019 – Cell Reports). Model R-squareds correspond to the variance explained in the dependent variable, computed using the sums-of-squares formulation. Beta weights are returned for individual variables, and bootstrapped confidence intervals on the beta estimates are used to determine the significance of individual variables as in Griffis et al., (2019 – Cell Reports). Usage is explained in detail in the script header.

plsr_script_oos.m – This implementation of PLSR fits predictive models that attempt to generalize to unseen data. It performs nested leave-one-out cross-validation in which the optimal number of components for the training model fit in the outer loops is determined within a second inner loop. Model R-squareds correspond to the cross-validation R-squared computed using the sums-of-squares formulation (e.g. see Poldrack et al., 2019 – JAMA Psychiatry). Usage is explained in detail in the script header. Note: out-of-sample prediction may be less optimal than PCRR due to both (1) optimistic bias of PCRR due to PCA performed prior to train-test splits/arbitrary component selection, and (2) lack of a regularization term in PLSR (meaning betas are unbiased, but may be less stable than penalized betas obtained from ridge).
