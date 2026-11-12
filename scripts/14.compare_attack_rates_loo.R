## Created by AI assistant on 2026-07-23 using GPT-5. Compares attack rate

## sensitivity analysis.
## compare the annula attack rates from the full anlaysis
## to the analysis generated where each biomarker has been left out

library(tidyverse)

main_ar <- read_csv(
  "results/model_versions/exponential_obs_ver_3/ic50_real/ic50_real_attack_rates.csv",
  show_col_types = FALSE
) %>%
  mutate(scenario = "All sera")

loo_2002_ar <- read_csv(
  "results/sensitivity_analyses/exponential_obs_ver_3_block_loo/ic50_real_2002/ic50_real_2002_attack_rates.csv",
  show_col_types = FALSE
) %>%
  mutate(scenario = "No 2002 sera")

loo_2006_ar <- read_csv(
  "results/sensitivity_analyses/exponential_obs_ver_3_block_loo/ic50_real_2006/ic50_real_2006_attack_rates.csv",
  show_col_types = FALSE
) %>%
  mutate(scenario = "No 2006 sera")

loo_2009_ar <- read_csv(
  "results/sensitivity_analyses/exponential_obs_ver_3_block_loo/ic50_real_2009/ic50_real_2009_attack_rates.csv",
  show_col_types = FALSE
) %>%
  mutate(scenario = "No 2009 sera")

loo_2012_ar <- read_csv(
  "results/sensitivity_analyses/exponential_obs_ver_3_block_loo/ic50_real_2012/ic50_real_2012_attack_rates.csv",
  show_col_types = FALSE
) %>%
  mutate(scenario = "No 2012 sera")

attack_rates <- bind_rows(main_ar, loo_2002_ar, loo_2006_ar, loo_2009_ar, loo_2012_ar) %>%
  group_by(scenario, time) %>%
  dplyr::summarise(
    median = median(V1),
    lower = quantile(V1, 0.025),
    upper = quantile(V1, 0.975),
    .groups = "drop"
  )

dir.create("figures/comparisons/sensitivity_loo", recursive = TRUE, showWarnings = FALSE)

p <- ggplot(attack_rates, aes(x = time, y = median, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(breaks = seq(min(attack_rates$time), max(attack_rates$time), by = 1)) +
  labs(x = "Year", y = "Attack rate", color = NULL, fill = NULL) +
  theme_classic()

ggsave("figures/comparisons/sensitivity_loo/attack_rates_main_vs_loo.pdf", p, width = 7, height = 4, units = "in", dpi = 300)
ggsave("figures/comparisons/sensitivity_loo/attack_rates_main_vs_loo.png", p, width = 7, height = 4, units = "in", dpi = 300)

# we want to report the % difference compared to the main model


tmp <- attack_rates %>% filter(scenario=="All sera"|scenario=="No 2002 sera")
chisq.test(tmp$median[tmp$scenario=="No 2002 sera"],tmp$median[tmp$scenario=="Main model"]) 

diff_tab <- (tmp$median[tmp$scenario=="No 2002 sera"]/tmp$median[tmp$scenario=="Main model"]) 
diff_tab <- rbind(diff_tab,attack_rates$median[attack_rates$scenario=="No 2006 sera"]/attack_rates$median[attack_rates$scenario=="All sera"]) 
diff_tab <- rbind(diff_tab,attack_rates$median[attack_rates$scenario=="No 2009 sera"]/attack_rates$median[attack_rates$scenario=="All sera"]) 
diff_tab <- rbind(diff_tab,attack_rates$median[attack_rates$scenario=="No 2012 sera"]/attack_rates$median[attack_rates$scenario=="All sera"]) 

write.csv(diff_tab,"figures/comparisons/sensitivity_loo/attack_rates_compare.csv")

ar_years <- attack_rates$median[attack_rates$scenario=="All sera"]
ar_years <- rbind(ar_years,attack_rates$median[attack_rates$scenario=="No 2002 sera"])
ar_years <- rbind(ar_years,attack_rates$median[attack_rates$scenario=="No 2006 sera"])
ar_years <- rbind(ar_years,attack_rates$median[attack_rates$scenario=="No 2009 sera"])
ar_years <- rbind(ar_years,attack_rates$median[attack_rates$scenario=="No 2012 sera"])

write.csv(ar_years,"figures/comparisons/sensitivity_loo/attack_rates_median.csv")


# end