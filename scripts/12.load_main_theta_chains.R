## Created by AI assistant on 2026-07-24 using GPT-5. Short helper script to
## load theta MCMC chains for the four main obs-ver-3 IC50 model scenarios.

library(tidyverse)
library(serosolver)

setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/")
source("code/scenario_paths.R")
fig_incidence_file <- function(...) file.path(incidence_wd, paste0(...))

main_theta_scenarios <- tribble(
  ~scenario, ~waning_model, ~obs_ver, ~run_name, ~burnin,
  "Exponential / Debbink", "exponential", 3, "ic50_real", 300000,
  "Exponential / Kendra", "exponential", 3, "ic50_kendra", 300000,
  "Linear / Debbink", "linear", 3, "ic50_real", 300000,
  "Linear / Kendra", "linear", 3, "ic50_kendra", 300000
) %>%
  mutate(
    paths = pmap(
      list(waning_model, obs_ver, run_name),
      ~ build_scenario_paths("model_versions", ..1, ..2, ..3)
    ),
    chain_wd = map_chr(paths, "chain_wd")
  )

missing_chain_dirs <- main_theta_scenarios %>%
  filter(!dir.exists(chain_wd))
if (nrow(missing_chain_dirs) > 0) {
  stop("Missing chain directories: ", paste(missing_chain_dirs$chain_wd, collapse = "; "))
}

theta_chains <- main_theta_scenarios %>%
  mutate(
    theta_chain = map2(
      chain_wd,
      burnin,
      ~ load_theta_chains(
        location = .x,
        burnin = .y,
        estimated_only = TRUE,
        convert_mcmc = FALSE,
        verbose = FALSE
      )$chain %>% as_tibble()
    )
  )

theta_chains_long <- theta_chains %>%
  select(scenario, waning_model, obs_ver, run_name, theta_chain) %>%
  unnest(theta_chain)

theta_chains %>%
  transmute(scenario, run_name, n_draws = map_int(theta_chain, nrow)) %>%
  print(n = Inf)

p_compare_likelihoods <- theta_chains_long %>% select(scenario, likelihood) %>% ggplot() + geom_density(aes(x=likelihood,fill=scenario),alpha=0.5) + theme_use + scale_y_continuous(expand=c(0,0)) +
  ylab("Posterior density") + xlab("Log likelihood")
ggsave(p_compare_likelihoods,filename="~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/final_figures/compare_likelihoods.pdf",width=5,height=4,dpi=300)
