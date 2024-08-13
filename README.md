# causalblb_test

## bias_std_err

**Goal**: Understand why increasing sample size with a fixed gamma produces worse coverage in previous iterations. We want to understand the bias in our results, so we are saving the entire bootstrap replication for each subset for each replication.

**Data Dictionary**:

- **n**: The full sample size
- **gamma**: The exponent for the subset size
- **subsets**: The number of subsets for the BLB
- **prop_form**: Whether the propensity score is correctly specified or misspecified
- **estim_subset**: The estimate of the ATE in each subset, calculated using a DR estimator. In this case, both propensity and outcome model are correctly specified.
- **sd_subset**: An estimate of the SE of the ATE in each subset. Calculated using the EIF.
- **boot_reps**: A vector of bootstrap replicates of the ATE. Here, the nuisance functions are re-estimated for each bootstrap resample.
- **boot_rep_num**: An index of the bootstrap resample
- **subset_num**: An index of the subset number
- **replication**: An index of the replication number
- **estim_full**: The estimate of the ATE on the full dataset
- **sd_full**: An estimated of the *variance* (misnamed column) on the full dataset

## high_n

**Goal**: Check if coverage stays approximately the same by increasing the number of subsets while increasing the total sample size *n*.

**Data Dictionary**:

- **n**: The full sample size
- **gamma**: The exponent for the subset size
- **subsets**: The number of subsets for the BLB
- **prop_form**: Whether the propensity score is correctly specified or misspecified
- **estim_subset**: The estimate of the ATE in each subset, calculated using a DR estimator. In this case, both propensity and outcome model are correctly specified.
- **sd_subset**: An estimate of the SE of the ATE in each subset. Calculated using the EIF.
- **boot_reps**: A vector of bootstrap replicates of the ATE. 
- **boot_rep_num**: An index of the bootstrap resample
- **subset_num**: An index of the subset number
- **replication**: An index of the replication number
- **estim_full**: The estimate of the ATE on the full dataset
- **sd_full**: An estimated of the *variance* (misnamed column) on the full dataset


## svm

**Goal**: Check the coverage compared to the full bootstrap using SVM to estimate the nuisance functions. 

**Data Dictionary**:


## weighted_bootstrap

**Goal**: Do separate treatment and control group bootstrappings. 

**Data Dictionary**:





