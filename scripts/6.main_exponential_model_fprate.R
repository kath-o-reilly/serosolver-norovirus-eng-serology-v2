setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2")

for(i in 1:3){
  fp_rates <- c(0.001,0.01,0.02)
  run_number <- 1
  rerun_mcmc <- TRUE
  obs_ver <- 3
  fp_rate_set <- fp_rates[i]
  fp_rate_fixed <- 1
  
  source("code/main_exponential_model_fprate_aux.R")
  
  rm(list=ls())
  gc()
}
