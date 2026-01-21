setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2")
for(i in 1:1){
  run_number <- i
  rerun_mcmc <- FALSE
  
  source("code/main_exponential_model_aux.R")
  
  rm(list=ls())
  gc()
}
source("code/generate_tables_exponential.R")
