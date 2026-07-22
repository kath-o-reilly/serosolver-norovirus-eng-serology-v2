# Created and modified by AI assistant on 2026-07-22 using GPT-5. Compares
# available ic50_real false-positive-rate sensitivity outputs using fprate_XX
# folders without sim_type as a directory layer.

library(tidyverse)
library(data.table)
library(serosolver)
library(ggplot2)
source("code/scenario_paths.R")

analysis_name <- "ic50_real"
output_dir <- file.path("figures", "comparisons", "fprate_scenarios", analysis_name)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

model_registry <- tibble::tribble(
  ~model_id, ~model_label, ~waning_model, ~obs_ver, ~observation_model, ~fp_rate, ~fp_rate_label, ~fallback_burnin, ~model_col,
  "fprate_0.001", "False-positive rate = 0.001", "exponential", 3, "spike-and-slab observation model", 0.001, "0.001", 300000, "#0072B2",
  "fprate_0.01", "False-positive rate = 0.01", "exponential", 3, "spike-and-slab observation model", 0.01, "0.01", 300000, "#D55E00",
  "fprate_0.02", "False-positive rate = 0.02", "exponential", 3, "spike-and-slab observation model", 0.02, "0.02", 300000, "#009E73"
) %>%
  mutate(
    scenario_id = purrr::pmap_chr(
      list(waning_model, obs_ver, fp_rate),
      ~build_scenario_id("fprate_scenarios", ..1, ..2, ..3)
    ),
    results_dir = file.path("results", "fprate_scenarios", scenario_id),
    local_dir = file.path("local_data", "fprate_scenarios", scenario_id),
    chain_dir = file.path(local_dir, "chains", analysis_name),
    figures_dir = file.path(local_dir, "figures", analysis_name),
    settings_file = file.path(chain_dir, paste0(analysis_name, "_serosolver_settings.RData")),
    theta_file = file.path(results_dir, analysis_name, paste0(analysis_name, "_theta_chains.csv")),
    out_res_file = file.path(figures_dir, "out_res.rda"),
    kinetics_file = NA_character_
  )

has_infection_chains <- function(chain_dir) {
  length(list.files(
    chain_dir,
    pattern = paste0("^", analysis_name, "_[0-9]+_infection_histories[.]csv$"),
    full.names = TRUE
  )) > 0
}

infer_burnin <- function(theta_file, fallback_burnin) {
  if (!file.exists(theta_file)) {
    return(fallback_burnin)
  }

  theta_chain <- suppressMessages(readr::read_csv(theta_file, show_col_types = FALSE))
  if (!"samp_no" %in% names(theta_chain) || nrow(theta_chain) == 0) {
    return(fallback_burnin)
  }

  max(0, min(theta_chain$samp_no, na.rm = TRUE) - 1)
}

load_settings <- function(settings_file) {
  env <- new.env(parent = emptyenv())
  load(settings_file, envir = env)
  env$serosolver_settings
}

load_out_res <- function(out_res_file) {
  if (!file.exists(out_res_file)) {
    return(NULL)
  }

  env <- new.env(parent = emptyenv())
  load(out_res_file, envir = env)
  env$res
}

included_models <- model_registry %>%
  rowwise() %>%
  mutate(
    has_settings = file.exists(settings_file),
    has_theta = file.exists(theta_file) || file.exists(out_res_file),
    has_inf = has_infection_chains(chain_dir),
    include = has_settings && has_theta && has_inf
  ) %>%
  ungroup()

message("Included models: ", paste(included_models$model_label[included_models$include], collapse = "; "))
message("Skipped models: ", paste(included_models$model_label[!included_models$include], collapse = "; "))

model_registry <- included_models %>%
  filter(include) %>%
  mutate(
    burnin = purrr::map2_dbl(theta_file, fallback_burnin, infer_burnin),
    model_label = as.character(model_label)
  )

if (nrow(model_registry) == 0) {
  stop("No complete model outputs found for ", analysis_name, ".")
}

model_cols <- setNames(model_registry$model_col, model_registry$model_label)
model_labels_wrapped <- setNames(stringr::str_wrap(names(model_cols), width = 28), names(model_cols))
noro_key <- c("2002" = "FH-2002", "2006" = "DH-2006", "2009" = "NO-2009", "2012" = "SY-2012")
noro_levels <- unname(noro_key)
parameter_labels <- c(
  "boost_long" = "Long-term boost",
  "boost_short" = "Transient boost",
  "wane_short" = "Transient waning rate",
  "cr_long" = "Long-term cross-reactivity gradient",
  "cr_short" = "Transient cross-reactivity gradient",
  "obs_sd" = "Observation error (standard deviation)",
  "antigenic_seniority" = "Antigenic seniority parameter"
)
theme_use <- theme_classic() +
  theme(
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )

load_model_chains <- function(row, estimated_only = FALSE) {
  settings <- load_settings(row$settings_file)
  chains <- load_mcmc_chains(
    row$chain_dir,
    settings$par_tab,
    burnin = row$burnin,
    estimated_only = estimated_only
  )
  list(settings = settings, chains = chains)
}

make_attack_rates <- function(row) {
  settings <- load_settings(row$settings_file)
  inf_chain <- load_infection_chains(location = row$chain_dir, burnin = row$burnin)$chain
  inf_chain <- pad_inf_chain(inf_chain, indivs = unique(settings$antibody_data$individual))
  colnames(inf_chain) <- c("individual", "time", "infected", "samp_no", "chain_no")

  norovirus_years <- 2001:2012
  inf_chain$time <- norovirus_years[inf_chain$time - min(inf_chain$time, na.rm = TRUE) + 1]

  ages <- settings$antibody_data %>%
    select(individual, birth) %>%
    distinct() %>%
    tidyr::expand_grid(time = norovirus_years) %>%
    mutate(age = time - birth)

  sample_times <- settings$antibody_data %>%
    select(individual, sample_time) %>%
    group_by(individual) %>%
    filter(sample_time == max(sample_time)) %>%
    ungroup() %>%
    distinct()

  inf_chain %>%
    as_tibble() %>%
    filter(time >= 2001, time <= 2012) %>%
    left_join(ages, by = c("individual", "time")) %>%
    left_join(sample_times, by = "individual") %>%
    filter(age >= 0, time <= sample_time) %>%
    group_by(samp_no, chain_no, time) %>%
    summarize(attack_rate = sum(infected) / n(), n_alive = n(), .groups = "drop") %>%
    group_by(time) %>%
    summarize(
      median = median(attack_rate),
      lower = quantile(attack_rate, 0.025),
      upper = quantile(attack_rate, 0.975),
      n_alive = median(n_alive),
      .groups = "drop"
    ) %>%
    mutate(
      model_id = row$model_id,
      model_label = row$model_label,
      waning_model = row$waning_model,
      observation_model = row$observation_model
    )
}

make_parameter_summaries <- function(row) {
  if (file.exists(row$theta_file)) {
    return(suppressMessages(readr::read_csv(row$theta_file, show_col_types = FALSE)) %>%
      filter(biomarker_group_label == "IC50") %>%
      group_by(label, name_key, names_new) %>%
      summarize(
        median = median(est),
        lower = quantile(est, 0.025),
        upper = quantile(est, 0.975),
        .groups = "drop"
      ) %>%
      mutate(
        parameter = if_else(is.na(name_key), names_new, name_key),
        model_id = row$model_id,
        model_label = row$model_label,
        waning_model = row$waning_model,
        observation_model = row$observation_model
      ))
  }

  out_res <- load_out_res(row$out_res_file)
  if (is.null(out_res$all_diagnostics$theta_estimates)) {
    stop("No theta-chain CSV or out_res theta estimates found for ", row$model_label)
  }

  as_tibble(as.data.frame(out_res$all_diagnostics$theta_estimates)) %>%
    filter(names %in% names(parameter_labels)) %>%
    transmute(
      label = parameter_labels[names],
      name_key = parameter_labels[names],
      names_new = names,
      median = median,
      lower = lower95_CrI,
      upper = upper95_CrI,
      parameter = parameter_labels[names],
      model_id = row$model_id,
      model_label = row$model_label,
      waning_model = row$waning_model,
      observation_model = row$observation_model
    )
}

make_kinetics_curves <- function(row) {
  if (!is.na(row$kinetics_file) && file.exists(row$kinetics_file)) {
    env <- new.env(parent = emptyenv())
    load(row$kinetics_file, envir = env)
    kinetics_plot <- get(ls(env)[1], envir = env)

    return(kinetics_plot$data %>%
      as_tibble() %>%
      mutate(
        strain = factor(strain, levels = noro_levels),
        model_id = row$model_id,
        model_label = row$model_label,
        waning_model = row$waning_model,
        observation_model = row$observation_model,
        time_since_infection = sample_time - min(sample_time, na.rm = TRUE)
      ))
  }

  regenerated_curves <- tryCatch({
    model <- load_model_chains(row, estimated_only = FALSE)
    antibody_data <- tidyr::expand_grid(
      individual = 1,
      sample_time = seq(2002, 2012, by = 0.1),
      biomarker_id = as.numeric(names(noro_key)),
      biomarker_group = 1
    )

    plot_estimated_antibody_model(
      model$chains$theta_chain,
      antibody_data = antibody_data,
      demographics = NULL,
      settings = model$settings,
      by_group = FALSE,
      nsamp = 100,
      antigenic_map = model$settings$antigenic_map %>% filter(inf_times >= 2002),
      possible_exposure_times = 2002:2012,
      par_tab = model$settings$par_tab %>% mutate(stratification = NA),
      solve_times = seq(2002, 2012, by = 0.1)
    )$data %>%
      mutate(
        strain = factor(noro_key[as.character(biomarker_id)], levels = noro_levels),
        model_id = row$model_id,
        model_label = row$model_label,
        waning_model = row$waning_model,
        observation_model = row$observation_model,
        time_since_infection = sample_time - 2002
      )
  }, error = function(err) {
    message("Falling back to saved antibody model data for ", row$model_label, ": ", conditionMessage(err))
    NULL
  })

  if (!is.null(regenerated_curves)) {
    return(regenerated_curves)
  }

  out_res <- load_out_res(row$out_res_file)
  if (!is.null(out_res$plot_antibody_model$data)) {
    return(out_res$plot_antibody_model$data %>%
      as_tibble() %>%
      filter(biomarker_group == 1) %>%
      mutate(
        strain = factor(noro_key[as.character(biomarker_id)], levels = noro_levels),
        model_id = row$model_id,
        model_label = row$model_label,
        waning_model = row$waning_model,
        observation_model = row$observation_model,
        time_since_infection = sample_time - min(sample_time, na.rm = TRUE)
      ))
  }

  stop("Could not build kinetics curves for ", row$model_label)
}

make_observation_histograms <- function(row) {
  out_res <- load_out_res(row$out_res_file)
  predictions <- out_res$antibody_predictions

  if (is.null(predictions$all_predictions_obs_draws)) {
    model <- load_model_chains(row, estimated_only = FALSE)
    predictions <- plot_antibody_predictions(
      model$chains$theta_chain,
      model$chains$inf_chain,
      settings = model$settings
    )
  }

  predicted_observations <- predictions$all_predictions_obs_draws %>%
    as_tibble() %>%
    filter(biomarker_group == 1) %>%
    transmute(
      value = prediction,
      distribution = "Posterior predictive observations"
    )

  observed_data <- predictions$all_predictions %>%
    as_tibble() %>%
    filter(biomarker_group == 1) %>%
    distinct(individual, sample_time, repeat_number, biomarker_id, biomarker_group, measurement) %>%
    transmute(
      value = measurement,
      distribution = "Observed data"
    )

  bind_rows(predicted_observations, observed_data) %>%
    mutate(
      distribution = factor(distribution, levels = c("Observed data", "Posterior predictive observations")),
      model_id = row$model_id,
      model_label = row$model_label,
      waning_model = row$waning_model,
      observation_model = row$observation_model
    )
}

attack_rates <- purrr::pmap_dfr(model_registry, function(...) make_attack_rates(tibble::tibble(...)))
parameter_summaries <- purrr::pmap_dfr(model_registry, function(...) make_parameter_summaries(tibble::tibble(...)))
kinetics_curves <- purrr::pmap_dfr(model_registry, function(...) make_kinetics_curves(tibble::tibble(...)))
observation_histograms <- purrr::pmap_dfr(model_registry, function(...) make_observation_histograms(tibble::tibble(...)))

required_models <- model_registry$model_label
attack_rates <- attack_rates %>% mutate(model_label = factor(model_label, levels = required_models))
parameter_summaries <- parameter_summaries %>% mutate(model_label = factor(model_label, levels = required_models))
kinetics_curves <- kinetics_curves %>% mutate(model_label = factor(model_label, levels = required_models))
observation_histograms <- observation_histograms %>% mutate(model_label = factor(model_label, levels = required_models))

stopifnot(
  all(required_models %in% as.character(unique(attack_rates$model_label))),
  all(required_models %in% as.character(unique(parameter_summaries$model_label))),
  all(required_models %in% as.character(unique(kinetics_curves$model_label))),
  all(required_models %in% as.character(unique(observation_histograms$model_label)))
)

p_attack_rates <- ggplot(attack_rates) +
  geom_ribbon(aes(x = time, ymin = lower, ymax = upper, fill = model_label), alpha = 0.18) +
  geom_line(aes(x = time, y = median, color = model_label), linewidth = 0.65) +
  geom_point(aes(x = time, y = median, color = model_label), size = 1.1) +
  scale_color_manual(name = "False-positive rate", values = model_cols, labels = model_labels_wrapped) +
  scale_fill_manual(name = "False-positive rate", values = model_cols, labels = model_labels_wrapped) +
  scale_x_continuous(breaks = 2001:2012) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  xlab("Norovirus year") +
  ylab("Attack rate (infections per capita)") +
  theme_use +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))

p_parameter_summaries <- ggplot(parameter_summaries) +
  geom_linerange(
    aes(x = model_label, ymin = lower, ymax = upper, color = model_label),
    position = position_dodge(width = 0.5),
    linewidth = 0.35
  ) +
  geom_point(
    aes(x = model_label, y = median, color = model_label, fill = model_label),
    shape = 21,
    size = 2.2,
    stroke = 0.35
  ) +
  facet_wrap(~parameter, scales = "free_y", ncol = 3) +
  scale_color_manual(name = "False-positive rate", values = model_cols, labels = model_labels_wrapped) +
  scale_fill_manual(name = "False-positive rate", values = model_cols, labels = model_labels_wrapped) +
  xlab(NULL) +
  ylab("Posterior median and 95% CrI") +
  theme_use +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))

p_kinetics_curves <- ggplot(kinetics_curves) +
  geom_ribbon(aes(x = time_since_infection, ymin = lower, ymax = upper, fill = model_label), alpha = 0.16) +
  geom_line(aes(x = time_since_infection, y = mean, color = model_label), linewidth = 0.55) +
  facet_wrap(~strain, ncol = 2) +
  scale_color_manual(name = "False-positive rate", values = model_cols, labels = model_labels_wrapped) +
  scale_fill_manual(name = "False-positive rate", values = model_cols, labels = model_labels_wrapped) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  xlab("Time since infection (years)") +
  ylab("Antibody level (log2 IC50)") +
  theme_use +
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))

p_observation_histograms <- ggplot(observation_histograms) +
  geom_histogram(
    aes(x = value, y = after_stat(density), fill = distribution),
    binwidth = 0.5,
    boundary = 0,
    position = "identity",
    alpha = 0.45,
    color = "grey30",
    linewidth = 0.15
  ) +
  facet_wrap(~model_label, labeller = label_wrap_gen(width = 28)) +
  scale_fill_manual(values = c("Observed data" = "grey25", "Posterior predictive observations" = "#56B4E9")) +
  xlab("Observed or predicted observation (log2 IC50)") +
  ylab("Density") +
  theme_use

ggsave(file.path(output_dir, "attack_rates_fprate_scenarios_ic50_real.pdf"), p_attack_rates, width = 9, height = 4.5, units = "in", dpi = 300)
ggsave(file.path(output_dir, "attack_rates_fprate_scenarios_ic50_real.png"), p_attack_rates, width = 9, height = 4.5, units = "in", dpi = 300)
ggsave(file.path(output_dir, "kinetics_posteriors_fprate_scenarios_ic50_real.pdf"), p_parameter_summaries, width = 9, height = 5.5, units = "in", dpi = 300)
ggsave(file.path(output_dir, "kinetics_posteriors_fprate_scenarios_ic50_real.png"), p_parameter_summaries, width = 9, height = 5.5, units = "in", dpi = 300)
ggsave(file.path(output_dir, "kinetics_curves_fprate_scenarios_ic50_real.pdf"), p_kinetics_curves, width = 9, height = 5.5, units = "in", dpi = 300)
ggsave(file.path(output_dir, "kinetics_curves_fprate_scenarios_ic50_real.png"), p_kinetics_curves, width = 9, height = 5.5, units = "in", dpi = 300)
ggsave(file.path(output_dir, "residuals_fprate_scenarios_ic50_real.pdf"), p_observation_histograms, width = 8, height = 5, units = "in", dpi = 300)
ggsave(file.path(output_dir, "residuals_fprate_scenarios_ic50_real.png"), p_observation_histograms, width = 8, height = 5, units = "in", dpi = 300)
ggsave(file.path(output_dir, "predicted_observation_histograms_fprate_scenarios_ic50_real.pdf"), p_observation_histograms, width = 8, height = 5, units = "in", dpi = 300)
ggsave(file.path(output_dir, "predicted_observation_histograms_fprate_scenarios_ic50_real.png"), p_observation_histograms, width = 8, height = 5, units = "in", dpi = 300)

save(
  model_registry,
  attack_rates,
  parameter_summaries,
  kinetics_curves,
  observation_histograms,
  p_attack_rates,
  p_parameter_summaries,
  p_kinetics_curves,
  p_observation_histograms,
  file = file.path(output_dir, "fprate_scenario_comparison_data_ic50_real.rda")
)

message("Saved false-positive-rate comparison figures to ", output_dir)
