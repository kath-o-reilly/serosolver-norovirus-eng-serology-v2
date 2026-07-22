setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2")
obs_ver <- 3
fp_rate_set <- 0.001
fp_rate_fixed <- 1

for(i in 1:1){
  run_number <- i
  rerun_mcmc <- TRUE
  
  source("code/main_exponential_model_aux.R")
  
  rm(list=ls())
  gc()
}
source("code/generate_tables_exponential.R")
