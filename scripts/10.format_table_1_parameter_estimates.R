## Created by AI assistant on 2026-07-23 using GPT-5. Formats the
## generated IC50 parameter estimate tables into the manuscript Table 1 layout.
library(tidyverse)

exp_file <- file.path("results", "model_versions", "exponential_obs_ver_3", "parameter_estimates_all.csv")
lin_file <- file.path("results", "model_versions", "linear_obs_ver_3", "parameter_estimates_all.csv")
out_file <- file.path("results", "model_versions", "table_1_parameter_estimates_formatted.csv")

param_key <- tribble(
  ~Parameter, ~Symbol, ~Meaning,
  "Long-term boost", "μ_long", "Long-term antibody boosting",
  "Transient boost", "μ_short", "Short-term antibody boosting",
  "Long-term cross-reactivity gradient", "σ_long", "Cross-reaction (retained)",
  "Transient cross-reactivity gradient", "σ_short", "Cross reaction (acute)",
  "Antigenic seniority parameter", "τ", "Immune imprinting, ie. proportion decrease in boost per successive infection",
  "Transient waning rate", "ω", "Acute waning",
  "Observation error (standard deviation)", "ε", "Standard deviation of the observation error distribution"
)

table_1 <- bind_rows(
  read_csv(exp_file, show_col_types = FALSE) %>% mutate(model = "exponential model"),
  read_csv(lin_file, show_col_types = FALSE) %>% mutate(model = "linear model")
) %>%
  filter(Data %in% c("ic50_real", "ic50_kendra")) %>%
  inner_join(param_key, by = "Parameter") %>%
  mutate(column = paste("Parameter estimate", model, Map, sep = " - ")) %>%
  select(Symbol, Meaning, Parameter, column, IC50) %>%
  pivot_wider(names_from = column, values_from = IC50) %>%
  arrange(match(Parameter, param_key$Parameter)) %>%
  select(
    Symbol,
    Meaning,
    `Parameter estimate - exponential model - Debbink`,
    `Parameter estimate - exponential model - Kendra`,
    `Parameter estimate - linear model - Debbink`,
    `Parameter estimate - linear model - Kendra`
  )

write_csv(table_1, out_file)
