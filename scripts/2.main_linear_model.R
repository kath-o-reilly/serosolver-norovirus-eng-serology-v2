## Modified by AI assistant on 2026-07-22 using GPT-5. Set the default linear
## model workflow to obs_ver 3 for spike-and-slab observation-model outputs.
setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2")
obs_ver <- 3
fp_rate_set <- 0.001
fp_rate_fixed <- 1

for(i in 2:2 ){ #6){
  run_number <- i
  rerun_mcmc <- TRUE
  source("code/main_linear_model_aux.R")
  
  rm(list=ls())
  gc()
}
obs_ver <- 3
fp_rate_set <- 0.001
fp_rate_fixed <- 1
source("code/generate_tables_linear.R")
