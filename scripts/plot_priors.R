## Modified by AI assistant on 2026-07-24 using GPT-5. Moved the prior plotting
## workflow into scripts/, added posterior overlays, and removed hard-coded
## local project paths.
library(tidyverse)
library(ggplot2)
library(patchwork)
library(serosolver)

project_dir <- normalizePath(getwd())
source("code/simulate_priors.R")
source("code/scenario_paths.R")

theme_use <- theme_classic() + theme(
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 6),
  strip.text = element_text(size = 8),
  plot.title = element_text(size = 10)
)

dir.create("figures/priors", recursive = TRUE, showWarnings = FALSE)

plot_posteriors <- TRUE
posterior_n_draws <- 100
prior_fill <- "grey75"
prior_line <- "grey35"
posterior_cols <- c("exponential" = "#0072B2", "linear" = "#D55E00")
source_cols <- c(
  "Prior" = prior_line,
  "Exponential posterior" = posterior_cols[["exponential"]],
  "Linear posterior" = posterior_cols[["linear"]]
)
source_fills <- c(
  "Prior" = prior_fill,
  "Exponential posterior" = posterior_cols[["exponential"]],
  "Linear posterior" = posterior_cols[["linear"]]
)
source_legend_data <- tibble(
  source = factor(rep(names(source_cols), each = 2), levels = names(source_cols)),
  name = rep(c(0, 0.01), times = length(source_cols)),
  mean = -1,
  lower = -1,
  upper = -1
)
source_legend_layers <- function() {
  list(
    geom_ribbon(
      data = source_legend_data,
      aes(x = name, ymin = lower, ymax = upper, fill = source, group = source),
      inherit.aes = FALSE,
      alpha = 0.25,
      show.legend = TRUE
    ),
    geom_line(
      data = source_legend_data,
      aes(x = name, y = mean, col = source, group = source),
      inherit.aes = FALSE,
      show.legend = TRUE
    )
  )
}
source_scales <- function(include_fill = TRUE, show_guide = TRUE) {
  source_guide <- if (show_guide) {
    guide_legend(
      override.aes = list(
        fill = unname(source_fills),
        alpha = 1
      )
    )
  } else {
    "none"
  }
  scales <- list(
    scale_colour_manual(name = NULL, values = source_cols, limits = names(source_cols), drop = FALSE, guide = source_guide)
  )
  if (include_fill) {
    scales <- c(scales, list(
      scale_fill_manual(name = NULL, values = source_fills, limits = names(source_fills), drop = FALSE, guide = "none")
    ))
  }
  scales
}
collect_guides <- function(p) {
  p + plot_layout(guides = "collect") & theme(legend.position = "bottom")
}
posterior_scenarios <- tibble(
  model = c("exponential", "linear"),
  waning_model = c("exponential", "linear"),
  obs_ver = c(3, 3),
  run_name = c("ic50_real", "ic50_real"),
  burnin = c(300000, 300000)
)

simulate_and_plot_prior_traj <- function(N, priors, times, model = c("exponential", "linear")) {
  model <- match.arg(model)
  traj <- matrix(0, ncol = length(times), nrow = N)
  for (i in 1:N) {
    if (model == "exponential") {
      traj[i, ] <- priors$boost_long[i] + priors$boost_short[i] * exp(-times * priors$wane_short[i])
    } else {
      traj[i, ] <- priors$boost_long[i] + pmax(0, priors$boost_short[i] * (1 - times * priors$wane_short[i]))
    }
  }
  traj <- as.data.frame(traj)
  colnames(traj) <- times
  traj$samp <- 1:N
  draws <- traj %>% pivot_longer(-samp)
  draws$name <- as.numeric(draws$name)
  traj <- traj %>%
    pivot_longer(-samp) %>%
    group_by(name) %>%
    dplyr::summarize(
      mean = median(value),
      lower50 = quantile(value, 0.25),
      upper50 = quantile(value, 0.75),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    )
  traj$name <- as.numeric(traj$name)
  list(traj = traj, draws = draws)
}

load_posterior_theta <- function(use_model) {
  scenario <- posterior_scenarios %>% filter(model == use_model) %>% slice(1)
  chain_wd <- build_scenario_paths(
    analysis_family = "model_versions",
    waning_model = scenario$waning_model,
    obs_ver = scenario$obs_ver,
    run_name = scenario$run_name
  )$chain_wd
  theta <- load_theta_chains(
    location = chain_wd,
    burnin = scenario$burnin,
    estimated_only = TRUE,
    convert_mcmc = FALSE,
    verbose = FALSE
  )$chain %>%
    as_tibble() %>%
    select(samp_no, chain_no, boost_long, boost_short, wane_short)
  theta %>% slice_sample(n = min(posterior_n_draws, nrow(theta)))
}

add_posterior_overlay <- function(p_lhs, p_rhs, posterior_draws, model, title_prefix = "") {
  posterior_label <- ifelse(model == "exponential", "Exponential posterior", "Linear posterior")
  posterior_traj <- simulate_and_plot_prior_traj(nrow(posterior_draws), posterior_draws, times, model)
  posterior_traj$traj$Assay <- "IC50"
  posterior_traj$draws$Assay <- "IC50"
  p_lhs <- p_lhs +
    geom_ribbon(data = posterior_traj$traj, aes(x = name, ymin = lower, ymax = upper, fill = posterior_label), alpha = 0.18) +
    geom_line(data = posterior_traj$traj, aes(x = name, y = mean, col = posterior_label))
  p_rhs <- p_rhs +
    geom_line(
      data = posterior_traj$draws,
      aes(x = name, y = value, group = samp, col = posterior_label),
      alpha = 0.25,
      show.legend = FALSE
    )
  list(lhs = p_lhs + ggtitle(paste0(title_prefix, "prior with posterior overlay")),
       rhs = p_rhs + ggtitle("Prior draws with posterior overlay"))
}

N <- 10000
times <- seq(0, 10, by = 0.1)

priors <- simulate_priors(N) %>% pivot_wider(names_from = name, values_from = value)
priors_ic50 <- priors[priors$biomarker_group == 1, ]
# priors_avidity <- priors[priors$biomarker_group == 2, ]

traj_exponential_ic50 <- simulate_and_plot_prior_traj(N, priors_ic50, times, "exponential")
# traj_exponential_avidity <- simulate_and_plot_prior_traj(N, priors_avidity, times, "exponential")
traj_exponential <- bind_rows(
  traj_exponential_ic50$traj %>% mutate(Assay = "IC50")
  # traj_exponential_avidity$traj %>% mutate(Assay = "Avidity")
)
draws_exponential <- bind_rows(
  traj_exponential_ic50$draws %>% mutate(Assay = "IC50")
  # traj_exponential_avidity$draws %>% mutate(Assay = "Avidity")
)
posterior_exponential <- if (plot_posteriors) load_posterior_theta("exponential") else NULL

p_lhs <- ggplot(traj_exponential) +
  geom_ribbon(aes(x = name, ymin = lower, ymax = upper, fill = "Prior"), alpha = 0.25) +
  geom_ribbon(aes(x = name, ymin = lower50, ymax = upper50, fill = "Prior"), alpha = 0.5) +
  geom_line(aes(x = name, y = mean, col = "Prior")) +
  source_legend_layers() +
  source_scales() +
  labs(x = "Time since boost (years)", y = "Antibody level") +
  theme_use +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 16, by = 2)) +
  coord_cartesian(ylim = c(0, 16)) +
  facet_wrap(~Assay, ncol = 1) +
  ggtitle("Prior median/50% quantiles/95% quantiles")

p_rhs <- ggplot(traj_exponential) +
  geom_line(
    data = draws_exponential %>% filter(samp %in% 1:100),
    aes(x = name, y = value, group = samp, col = "Prior"),
    alpha = 0.25,
    show.legend = FALSE
  ) +
  source_scales(include_fill = FALSE, show_guide = FALSE) +
  labs(x = "Time since boost (years)", y = "Antibody level") +
  theme_use +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 16, by = 2)) +
  coord_cartesian(ylim = c(0, 16)) +
  facet_wrap(~Assay, ncol = 1) +
  ggtitle("Prior draws")

if (plot_posteriors) {
  overlay <- add_posterior_overlay(p_lhs, p_rhs, posterior_exponential, "exponential", "Exponential ")
  p_lhs <- overlay$lhs
  p_rhs <- overlay$rhs
}

ggsave("figures/priors/prior_exponential_traj.png", collect_guides(p_lhs | p_rhs), height = 5, width = 7)
ggsave("figures/priors/prior_exponential_traj.pdf", collect_guides(p_lhs | p_rhs), height = 5, width = 7)

p_lhs_ic50 <- ggplot(traj_exponential %>% filter(Assay == "IC50")) +
  geom_ribbon(aes(x = name, ymin = lower, ymax = upper, fill = "Prior"), alpha = 0.25) +
  geom_ribbon(aes(x = name, ymin = lower50, ymax = upper50, fill = "Prior"), alpha = 0.5) +
  geom_line(aes(x = name, y = mean, col = "Prior")) +
  source_legend_layers() +
  source_scales() +
  labs(x = "Time since boost (years)", y = "Antibody level") +
  theme_use +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 16, by = 2)) +
  coord_cartesian(ylim = c(0, 16)) +
  ggtitle("Prior median/50% quantiles/95% quantiles")

p_rhs_ic50 <- ggplot(traj_exponential %>% filter(Assay == "IC50")) +
  geom_line(
    data = draws_exponential %>% filter(samp %in% 1:100, Assay == "IC50"),
    aes(x = name, y = value, group = samp, col = "Prior"),
    alpha = 0.25,
    show.legend = FALSE
  ) +
  source_scales(include_fill = FALSE, show_guide = FALSE) +
  labs(x = "Time since boost (years)", y = "Antibody level") +
  theme_use +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 16, by = 2)) +
  coord_cartesian(ylim = c(0, 16)) +
  ggtitle("Prior draws")

if (plot_posteriors) {
  overlay <- add_posterior_overlay(p_lhs_ic50, p_rhs_ic50, posterior_exponential, "exponential", "Exponential ")
  p_lhs_ic50 <- overlay$lhs
  p_rhs_ic50 <- overlay$rhs
}

ggsave("figures/priors/prior_exponential_traj_ic50.png", collect_guides(p_lhs_ic50 | p_rhs_ic50), height = 3, width = 7)
ggsave("figures/priors/prior_exponential_traj_ic50.pdf", collect_guides(p_lhs_ic50 | p_rhs_ic50), height = 3, width = 7)

p_exponential_lhs_ic50 <- p_lhs_ic50
p_exponential_rhs_ic50 <- p_rhs_ic50

priors <- simulate_priors_original_model(N) %>% pivot_wider(names_from = name, values_from = value)
priors_ic50 <- priors[priors$biomarker_group == 1, ]
# priors_avidity <- priors[priors$biomarker_group == 2, ]

traj_linear_ic50 <- simulate_and_plot_prior_traj(N, priors_ic50, times, "linear")
# traj_linear_avidity <- simulate_and_plot_prior_traj(N, priors_avidity, times, "linear")
traj_linear <- bind_rows(
  traj_linear_ic50$traj %>% mutate(Assay = "IC50")
  # traj_linear_avidity$traj %>% mutate(Assay = "Avidity")
)
draws_linear <- bind_rows(
  traj_linear_ic50$draws %>% mutate(Assay = "IC50")
  # traj_linear_avidity$draws %>% mutate(Assay = "Avidity")
)
posterior_linear <- if (plot_posteriors) load_posterior_theta("linear") else NULL

p_lhs <- ggplot(traj_linear) +
  geom_ribbon(aes(x = name, ymin = lower, ymax = upper, fill = "Prior"), alpha = 0.25) +
  geom_ribbon(aes(x = name, ymin = lower50, ymax = upper50, fill = "Prior"), alpha = 0.5) +
  geom_line(aes(x = name, y = mean, col = "Prior")) +
  source_legend_layers() +
  source_scales() +
  labs(x = "Time since boost (years)", y = "Antibody level") +
  theme_use +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 16, by = 2)) +
  coord_cartesian(ylim = c(0, 16)) +
  facet_wrap(~Assay, ncol = 1) +
  ggtitle("Prior median/50% quantiles/95% quantiles")

p_rhs <- ggplot(traj_linear) +
  geom_line(
    data = draws_linear %>% filter(samp %in% 1:100),
    aes(x = name, y = value, group = samp, col = "Prior"),
    alpha = 0.25,
    show.legend = FALSE
  ) +
  source_scales(include_fill = FALSE, show_guide = FALSE) +
  labs(x = "Time since boost (years)", y = "Antibody level") +
  theme_use +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 16, by = 2)) +
  coord_cartesian(ylim = c(0, 16)) +
  facet_wrap(~Assay, ncol = 1) +
  ggtitle("Prior draws")

if (plot_posteriors) {
  overlay <- add_posterior_overlay(p_lhs, p_rhs, posterior_linear, "linear", "Linear ")
  p_lhs <- overlay$lhs
  p_rhs <- overlay$rhs
}

ggsave("figures/priors/prior_linear_traj.png", collect_guides(p_lhs | p_rhs), height = 5, width = 7)
ggsave("figures/priors/prior_linear_traj.pdf", collect_guides(p_lhs | p_rhs), height = 5, width = 7)

p_lhs_ic50 <- ggplot(traj_linear %>% filter(Assay == "IC50")) +
  geom_ribbon(aes(x = name, ymin = lower, ymax = upper, fill = "Prior"), alpha = 0.25) +
  geom_ribbon(aes(x = name, ymin = lower50, ymax = upper50, fill = "Prior"), alpha = 0.5) +
  geom_line(aes(x = name, y = mean, col = "Prior")) +
  source_legend_layers() +
  source_scales() +
  labs(x = "Time since boost (years)", y = "Antibody level") +
  theme_use +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 16, by = 2)) +
  coord_cartesian(ylim = c(0, 16)) +
  ggtitle("Prior median/50% quantiles/95% quantiles")

p_rhs_ic50 <- ggplot(traj_linear %>% filter(Assay == "IC50")) +
  geom_line(
    data = draws_linear %>% filter(samp %in% 1:100, Assay == "IC50"),
    aes(x = name, y = value, group = samp, col = "Prior"),
    alpha = 0.25,
    show.legend = FALSE
  ) +
  source_scales(include_fill = FALSE, show_guide = FALSE) +
  labs(x = "Time since boost (years)", y = "Antibody level") +
  theme_use +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_y_continuous(breaks = seq(0, 16, by = 2)) +
  coord_cartesian(ylim = c(0, 16)) +
  theme() +
  ggtitle("Prior draws")

if (plot_posteriors) {
  overlay <- add_posterior_overlay(p_lhs_ic50, p_rhs_ic50, posterior_linear, "linear", "Linear ")
  p_lhs_ic50 <- overlay$lhs
  p_rhs_ic50 <- overlay$rhs
}

ggsave("figures/priors/prior_linear_traj_ic50.png", collect_guides(p_lhs_ic50 | p_rhs_ic50), height = 3, width = 7)
ggsave("figures/priors/prior_linear_traj_ic50.pdf", collect_guides(p_lhs_ic50 | p_rhs_ic50), height = 3, width = 7)

p_comb_ic50 <- wrap_plots(
  p_exponential_lhs_ic50,
  p_exponential_rhs_ic50 + guides(colour = "none", fill = "none"),
  p_lhs_ic50 + guides(colour = "none", fill = "none"),
  p_rhs_ic50 + guides(colour = "none", fill = "none"),
  ncol = 2,
  guides = "collect"
) & theme(legend.position = "bottom")

ggsave("figures/priors/prior_comb_traj_ic50.pdf", p_comb_ic50, height = 6, width = 7)
ggsave("figures/priors/prior_comb_traj_ic50.png", p_comb_ic50, height = 6, width = 7)
