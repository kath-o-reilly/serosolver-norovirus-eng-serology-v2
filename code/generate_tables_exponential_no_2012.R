## Created by AI assistant on 2026-07-22 using GPT-5. Copied from the
## exponential table script to summarize the no-2012 sensitivity analysis.
library(tidyverse)
source("code/scenario_paths.R")

results_scenario_dir <- file.path(
  "results",
  "sensitivity_analyses",
  build_scenario_id(
    "sensitivity_analyses",
    "exponential",
    obs_ver,
    sensitivity_name = "exclude_2012"
  )
)

## Combine and save parameter estimates
theta_files <- list.files(results_scenario_dir, pattern="theta_chain", full.names=TRUE, recursive=TRUE)
theta_ests_all <- NULL
for(file in theta_files){
  theta_ests <- read_csv(file)
  theta_ests_tmp <- theta_ests %>% dplyr::group_by(name_key, biomarker_group_label,run_name,map) %>% dplyr::summarize(post_mean=mean(est),
                                                                                           lower95=quantile(est,0.025),
                                                                                           lower50=quantile(est,0.25),
                                                                                           upper50=quantile(est,0.975),
                                                                                           upper95=quantile(est,0.975)) %>%
  dplyr::mutate(est=paste0(signif(post_mean,3), " (", signif(lower95,3), ", ", signif(upper95,3), ")")) %>%
    dplyr::select(name_key,biomarker_group_label, run_name,map,est) %>%
    dplyr::mutate(map=if_else(map=="debbink","Debbink","Kendra")) %>% pivot_wider(names_from=biomarker_group_label,values_from=est) %>%
    dplyr::rename("Parameter"=name_key, "Data"=run_name, "Map"=map)
theta_ests_all <- bind_rows(theta_ests_all, theta_ests_tmp)
}
write_csv(theta_ests_all, file.path(results_scenario_dir, "parameter_estimates_all.csv"))

## Combine and save attach rate estimates
ar_files <- list.files(results_scenario_dir, pattern="attack_rates.csv", full.names=TRUE, recursive=TRUE)
ar_ests_all <- NULL
for(file in ar_files){
  ar_ests <- read_csv(file)
  ar_ests_summary <- ar_ests %>%dplyr::group_by(time,run_name,map) %>% dplyr::summarize(post_med=median(V1),
                                                                                           lower95=quantile(V1,0.025),
                                                                                           lower50=quantile(V1,0.25),
                                                                                           upper50=quantile(V1,0.975),
                                                                                           upper95=quantile(V1,0.975)) %>%
    dplyr::mutate(est=paste0(signif(post_med,3), " (", signif(lower95,3), ", ", signif(upper95,3), ")")) %>%
    dplyr::select(time,run_name,map, est) %>%
    dplyr::mutate(map=if_else(map=="debbink","Debbink","Kendra")) %>% 
    dplyr::rename("Data"=run_name, "Map"=map,"Year"=time,"Estimate"=est)
  ar_ests_all <- bind_rows(ar_ests_all, ar_ests_summary)
}
write_csv(ar_ests_all, file.path(results_scenario_dir, "ar_estimates_all.csv"))
