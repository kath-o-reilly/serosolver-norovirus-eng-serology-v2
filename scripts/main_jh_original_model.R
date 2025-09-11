setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2")

for(i in 1:6){
  run_number <- i
  rerun_mcmc <- TRUE
  source("code/main_script_v4_original_model.R")
  
  rm(list=ls())
  gc()
}
source("code/generate_tables_linear.R")
