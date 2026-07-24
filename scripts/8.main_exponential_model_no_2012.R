## Created and modified by AI assistant through 2026-07-24 using GPT-5.
## Runs the obs-ver-3 no-2012 sensitivity from the repository root.
project_dir <- normalizePath(getwd())
obs_ver <- 3
fp_rate_set <- 0.001
fp_rate_fixed <- 1

for(i in 1:1){
  run_number <- i
  rerun_mcmc <- TRUE
  
  source("code/main_exponential_model_no_2012_aux.R")
  rm(list=ls())
  gc()
}

obs_ver <- 3
fp_rate_set <- 0.001
fp_rate_fixed <- 1
source("code/generate_tables_exponential_no_2012.R")
