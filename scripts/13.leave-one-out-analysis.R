# 13. leave one out analysis

## 1) carry out the exponential model, but leaving one of the 4 sera out of the analysis.
##.   - only worth doing for Kendra map really
## 2) repeat for each sera
## 3) compare the attack rates for each analysis to the full dataset

## Modified by AI assistant on 2026-07-24 using GPT-5. Removed the hard-coded
## local project path so this script runs from the repository root.
project_dir <- normalizePath(getwd())

## Kendra and Debbink maps
## There are 6 rows total in the input table, the first 2 rows use the IC50 data for the two maps
for(run_number in 2:5){ 
  ## Obs ver 2 and 3 -- fit with each version
  for(j in 2){
    obs_vers <- c(2,3)
    obs_ver <- obs_vers[j]
    fp_rate_set <- 0.001
    fp_rate_fixed <- 1
    rerun_mcmc <- TRUE
    source("code/main_exponential_model_loo.R")
    gc()
  }
}
## Use obs ver 3 for the main results
# obs_ver <- 3
# fp_rate_set <- 0.001
# fp_rate_fixed <- 1
# source("code/generate_tables_exponential.R")