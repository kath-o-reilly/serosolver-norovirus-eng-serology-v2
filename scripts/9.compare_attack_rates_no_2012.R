## Created by AI assistant on 2026-07-23 using GPT-5. Compares attack rate
## estimates from the main exponential obs-ver-3 scenario and the no-2012
## sensitivity analysis.

library(tidyverse)

main_ar <- read_csv(
  "results/model_versions/exponential_obs_ver_3/ic50_real/ic50_real_attack_rates.csv",
  show_col_types = FALSE
) %>%
  mutate(scenario = "Main model")

no_2012_ar <- read_csv(
  "results/sensitivity_analyses/exponential_obs_ver_3_exclude_2012/ic50_real/ic50_real_attack_rates.csv",
  show_col_types = FALSE
) %>%
  mutate(scenario = "No 2012 data")

attack_rates <- bind_rows(main_ar, no_2012_ar) %>%
  group_by(scenario, time) %>%
  dplyr::summarise(
    median = median(V1),
    lower = quantile(V1, 0.025),
    upper = quantile(V1, 0.975),
    .groups = "drop"
  )

dir.create("figures/comparisons/no_2012", recursive = TRUE, showWarnings = FALSE)

p <- ggplot(attack_rates, aes(x = time, y = median, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(breaks = seq(min(attack_rates$time), max(attack_rates$time), by = 1)) +
  labs(x = "Year", y = "Attack rate", color = NULL, fill = NULL) +
  theme_classic()

ggsave("figures/comparisons/no_2012/attack_rates_main_vs_no_2012.pdf", p, width = 7, height = 4, units = "in", dpi = 300)
ggsave("figures/comparisons/no_2012/attack_rates_main_vs_no_2012.png", p, width = 7, height = 4, units = "in", dpi = 300)
