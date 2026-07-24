## Created by AI assistant on 2026-07-22 using GPT-5. Copied from the main
## exponential workflow to run the obs-ver-3 no-2012 sensitivity analysis.
setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2")
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
