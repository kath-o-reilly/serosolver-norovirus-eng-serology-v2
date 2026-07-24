## Modified by AI assistant on 2026-07-24 using GPT-5. This script mirrors the
## prior plotting workflow, compares exponential versus linear cross-reactivity
## decline over antigenic distance, and optionally overlays posterior draws.
library(tidyverse)
library(ggplot2)
library(patchwork)
library(serosolver)

setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2/")
source("code/simulate_priors.R")
source("code/scenario_paths.R")

theme_use <- theme_classic() + theme(
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  legend.title = element_text(size = 6),
  legend.text = element_text(size = 6),
  strip.text = element_text(size = 8),
  plot.title = element_text(size = 8),
  legend.position = "bottom"
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

simulate_cr_distance_traj <- function(N, priors, distances, arm = c("short", "long"),
                                      distance_model = c("exponential", "linear")) {
  arm <- match.arg(arm)
  distance_model <- match.arg(distance_model)
  cr <- matrix(0, ncol = length(distances), nrow = N)
  for (i in 1:N) {
    cr_gradient <- if (arm == "short") priors$cr_short[i] else priors$cr_long[i]
    if (distance_model == "exponential") {
      cr[i, ] <- exp(-distances * cr_gradient)
    } else {
      cr[i, ] <- pmax(0, 1 - distances * cr_gradient)
    }
  }
  cr <- as.data.frame(cr)
  colnames(cr) <- distances
  cr$samp <- 1:N
  draws <- cr %>% pivot_longer(-samp)
  draws$name <- as.numeric(draws$name)
  cr <- cr %>%
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
  cr$name <- as.numeric(cr$name)
  list(traj = cr, draws = draws)
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
    select(samp_no, chain_no, cr_short, cr_long)
  theta %>% slice_sample(n = min(posterior_n_draws, nrow(theta)))
}

N <- 10000
distances <- seq(0, 10, by = 0.1)

plot_cross_reactivity_prior <- function(priors, title_prefix,
                                        distance_model = c("exponential", "linear"),
                                        posterior_draws = NULL) {
  distance_model <- match.arg(distance_model)
  model_label <- ifelse(distance_model == "exponential", "Exponential", "Linear")
  priors_ic50 <- priors[priors$biomarker_group == 1, ]
  # priors_avidity <- priors[priors$biomarker_group == 2, ]

  cr_short_ic50 <- simulate_cr_distance_traj(N, priors_ic50, distances, "short", distance_model)
  # cr_short_avidity <- simulate_cr_distance_traj(N, priors_avidity, distances, "short", distance_model)
  cr_long_ic50 <- simulate_cr_distance_traj(N, priors_ic50, distances, "long", distance_model)
  # cr_long_avidity <- simulate_cr_distance_traj(N, priors_avidity, distances, "long", distance_model)

  cr_traj <- bind_rows(
    cr_short_ic50$traj %>% mutate(Arm = "Short-term"),
    # cr_short_avidity$traj %>% mutate(Arm = "Short-term", Assay = "Avidity"),
    cr_long_ic50$traj %>% mutate(Arm = "Long-term")
    # cr_long_avidity$traj %>% mutate(Arm = "Long-term", Assay = "Avidity")
  )
  cr_draws <- bind_rows(
    cr_short_ic50$draws %>% mutate(Arm = "Short-term"),
    # cr_short_avidity$draws %>% mutate(Arm = "Short-term", Assay = "Avidity"),
    cr_long_ic50$draws %>% mutate(Arm = "Long-term")
    # cr_long_avidity$draws %>% mutate(Arm = "Long-term", Assay = "Avidity")
  )
  cr_traj$Arm <- factor(cr_traj$Arm, levels = c("Short-term", "Long-term"))
  cr_draws$Arm <- factor(cr_draws$Arm, levels = c("Short-term", "Long-term"))
  cr_traj$`Distance model` <- model_label
  cr_draws$`Distance model` <- model_label

  p_lhs <- ggplot(cr_traj) +
    geom_ribbon(aes(x = name, ymin = lower, ymax = upper, fill = "Prior"), alpha = 0.25) +
    geom_ribbon(aes(x = name, ymin = lower50, ymax = upper50, fill = "Prior"), alpha = 0.45) +
    geom_line(aes(x = name, y = mean, col = "Prior")) +
    source_legend_layers() +
    source_scales() +
    labs(x = "Antigenic distance", y = "Relative cross-reactivity") +
    theme_use +
    scale_x_continuous(breaks = seq(0, 10, by = 2)) +
    coord_cartesian(ylim = c(0, 1)) +
    facet_wrap(~Arm, ncol = 1) +
    ggtitle(paste0(title_prefix, "\n prior median/50% quantiles/95% quantiles"))

  p_rhs <- ggplot(cr_traj) +
    geom_line(
      data = cr_draws %>% filter(samp %in% 1:100),
      aes(x = name, y = value, group = samp, col = "Prior"),
      alpha = 0.25,
      show.legend = FALSE
    ) +
    source_scales(include_fill = FALSE, show_guide = FALSE) +
    labs(x = "Antigenic distance", y = "Relative cross-reactivity") +
    theme_use +
    scale_x_continuous(breaks = seq(0, 10, by = 2)) +
    coord_cartesian(ylim = c(0, 1)) +
    facet_wrap(~Arm, ncol = 1) +
    ggtitle(paste0(title_prefix, " draws"))

  if (!is.null(posterior_draws)) {
    posterior_label <- ifelse(distance_model == "exponential", "Exponential posterior", "Linear posterior")
    cr_short_posterior <- simulate_cr_distance_traj(nrow(posterior_draws), posterior_draws, distances, "short", distance_model)
    cr_long_posterior <- simulate_cr_distance_traj(nrow(posterior_draws), posterior_draws, distances, "long", distance_model)
    cr_traj_posterior <- bind_rows(
      cr_short_posterior$traj %>% mutate(Arm = "Short-term"),
      cr_long_posterior$traj %>% mutate(Arm = "Long-term")
    )
    cr_draws_posterior <- bind_rows(
      cr_short_posterior$draws %>% mutate(Arm = "Short-term"),
      cr_long_posterior$draws %>% mutate(Arm = "Long-term")
    )
    cr_traj_posterior$Arm <- factor(cr_traj_posterior$Arm, levels = c("Short-term", "Long-term"))
    cr_draws_posterior$Arm <- factor(cr_draws_posterior$Arm, levels = c("Short-term", "Long-term"))

    p_lhs <- p_lhs +
      geom_ribbon(data = cr_traj_posterior, aes(x = name, ymin = lower, ymax = upper, fill = posterior_label), alpha = 0.18) +
      geom_line(data = cr_traj_posterior, aes(x = name, y = mean, col = posterior_label)) +
      ggtitle(paste0(title_prefix, "\n prior with posterior overlay"))
    p_rhs <- p_rhs +
      geom_line(
        data = cr_draws_posterior,
        aes(x = name, y = value, group = samp, col = posterior_label),
        alpha = 0.25,
        show.legend = FALSE
      ) +
      ggtitle(paste0(title_prefix, "\n draws with posterior overlay"))
  }

  list(lhs = p_lhs, rhs = p_rhs, plot = p_lhs | p_rhs)
}

priors_exponential <- simulate_priors(N) %>% pivot_wider(names_from = name, values_from = value)
posterior_exponential <- if (plot_posteriors) load_posterior_theta("exponential") else NULL
p_exponential <- plot_cross_reactivity_prior(
  priors_exponential,
  "Exponential cross-reactivity distance model",
  "exponential",
  posterior_exponential
)
ggsave(
  "figures/priors/prior_cross_reactivity_distance_exponential_ic50.png",
  collect_guides(p_exponential$plot),
  height = 5,
  width = 7
)
ggsave(
  "figures/priors/prior_cross_reactivity_distance_exponential_ic50.pdf",
  collect_guides(p_exponential$plot),
  height = 5,
  width = 7
)

priors_linear <- simulate_priors_original_model(N) %>% pivot_wider(names_from = name, values_from = value)
posterior_linear <- if (plot_posteriors) load_posterior_theta("linear") else NULL
p_linear <- plot_cross_reactivity_prior(
  priors_linear,
  "Linear cross-reactivity distance model",
  "linear",
  posterior_linear
)
ggsave(
  "figures/priors/prior_cross_reactivity_distance_linear_ic50.png",
  collect_guides(p_linear$plot),
  height = 5,
  width = 7
)
ggsave(
  "figures/priors/prior_cross_reactivity_distance_linear_ic50.pdf",
  collect_guides(p_linear$plot),
  height = 5,
  width = 7
)

p_combined <- wrap_plots(
  p_exponential$lhs,
  p_exponential$rhs + guides(colour = "none", fill = "none"),
  p_linear$lhs + guides(colour = "none", fill = "none"),
  p_linear$rhs + guides(colour = "none", fill = "none"),
  ncol = 2,
  guides = "collect"
) & theme(legend.position = "bottom")

ggsave(
  "figures/priors/prior_cross_reactivity_distance_combined_ic50.png",
  p_combined,
  height = 10,
  width = 7
)
ggsave(
  "figures/priors/prior_cross_reactivity_distance_combined_ic50.pdf",
  p_combined,
  height = 10,
  width = 7
)
