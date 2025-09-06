library(tidyverse)

## Combine and save parameter estimates
theta_files <- list.files("results_exponential",pattern="theta_chain",full.names=TRUE)
theta_ests_all <- NULL
for(file in theta_files){
  theta_ests <- read_csv(file)
  theta_ests_tmp <- theta_ests %>% group_by(name_key, biomarker_group_label,run_name,map) %>% dplyr::summarize(post_mean=mean(est),
                                                                                           lower95=quantile(est,0.025),
                                                                                           lower50=quantile(est,0.25),
                                                                                           upper50=quantile(est,0.975),
                                                                                           upper95=quantile(est,0.975)) %>%
  mutate(est=paste0(signif(post_mean,3), " (", signif(lower95,3), ", ", signif(upper95,3), ")")) %>%
  select(name_key,biomarker_group_label, run_name,map,est) %>%
  mutate(map=if_else(map=="debbink","Debbink","Kendra")) %>% pivot_wider(names_from=biomarker_group_label,values_from=est) %>%
  rename("Parameter"=name_key, "Data"=run_name, "Map"=map)
theta_ests_all <- bind_rows(theta_ests_all, theta_ests_tmp)
}
write_csv(theta_ests_all, "results_exponential/parameter_estimates_all.csv")

## Combine and save attach rate estimates
ar_files <- list.files("results_exponential",pattern="attack_rates.csv",full.names=TRUE)
ar_ests_all <- NULL
for(file in ar_files){
  ar_ests <- read_csv(file)
  ar_ests_summary <- ar_ests %>% group_by(time,run_name,map) %>% dplyr::summarize(post_med=median(V1),
                                                                                           lower95=quantile(V1,0.025),
                                                                                           lower50=quantile(V1,0.25),
                                                                                           upper50=quantile(V1,0.975),
                                                                                           upper95=quantile(V1,0.975)) %>%
    mutate(est=paste0(signif(post_med,3), " (", signif(lower95,3), ", ", signif(upper95,3), ")")) %>%
    select(time,run_name,map, est) %>%
    mutate(map=if_else(map=="debbink","Debbink","Kendra")) %>% 
    rename("Data"=run_name, "Map"=map,"Year"=time,"Estimate"=est)
  ar_ests_all <- bind_rows(ar_ests_all, ar_ests_summary)
}
write_csv(ar_ests_all, "results_exponential/ar_estimates_all.csv")