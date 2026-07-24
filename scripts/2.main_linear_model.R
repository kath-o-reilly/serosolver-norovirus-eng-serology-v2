## Modified by AI assistant on 2026-07-24 using GPT-5. Set the default linear
## model workflow to obs_ver 3 and removed the hard-coded local project path.
project_dir <- normalizePath(getwd())


## Kendra and Debbink maps
## There are 6 rows total in the input table, the first 2 rows use the IC50 data for the two maps
for(run_number in 2){ 
  ## Obs ver 2 and 3 -- fit with each version
  for(j in 2){
    obs_vers <- c(2,3)
    obs_ver <- obs_vers[j]
    fp_rate_set <- 0.001
    fp_rate_fixed <- 1
    rerun_mcmc <- FALSE
    source("code/main_linear_model_aux.R")
    gc()
  }
}
## Use obs ver 3 for the main results
obs_ver <- 3
fp_rate_set <- 0.001
fp_rate_fixed <- 1
source("code/generate_tables_linear.R")
