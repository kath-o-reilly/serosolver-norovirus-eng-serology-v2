## Created by AI assistant on 2026-07-23 using GPT-5. Reads saved posterior
## predictive histogram data and compares exponential obs-ver-2 versus obs-ver-3
## IC50 predictions in a small faceted figure.
library(tidyverse)
library(ggplot2)

setwd("~/Documents/GitHub/serosolver-norovirus-eng-serology-v2")
source("code/scenario_paths.R")

run_name <- "ic50_real"
n_draws_to_plot <- 3
draws_to_plot <- NULL

output_dir <- file.path("figures", "comparisons", "model_versions", run_name)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

model_registry <- tibble(
  obs_ver = c(2, 3),
  model_label = c("Original observation model", "Updated observation model")
) %>%
  mutate(
    save_wd = map_chr(
      obs_ver,
      ~ build_scenario_paths(
        analysis_family = "model_versions",
        waning_model = "exponential",
        obs_ver = .x,
        run_name = run_name
      )$save_wd
    ),
    hist_file = file.path(save_wd, "antibody_prediction_histogram_draws_data.rds")
  )

missing_files <- model_registry$hist_file[!file.exists(model_registry$hist_file)]
if (length(missing_files) > 0) {
  stop(
    "Missing histogram data files. Re-run the relevant aux scripts first:\n",
    paste(missing_files, collapse = "\n")
  )
}

read_histogram_data <- function(hist_file, model_label) {
  dat <- readRDS(hist_file) %>%
    as_tibble() %>%
    filter(biomarker_group == 1)

  selected_draws <- if (is.null(draws_to_plot)) {
    unique(dat$name)[seq_len(min(n_draws_to_plot, n_distinct(dat$name)))]
  } else {
    draws_to_plot
  }

  dat <- dat %>%
    filter(name %in% selected_draws) %>%
    mutate(
      model_label = model_label,
      draw_panel = paste0("Posterior draw ", match(name, selected_draws))
    )

  bind_rows(
    dat %>%
      transmute(
        model_label,
        draw_panel,
        value = value,
        distribution = "Posterior predictive"
      ),
    dat %>%
      distinct(model_label, draw_panel, i, measurement) %>%
      transmute(
        model_label,
        draw_panel,
        value = measurement,
        distribution = "Observed"
      )
  )
}

plot_data <- pmap_dfr(
  model_registry %>% select(hist_file, model_label),
  read_histogram_data
) %>%
  mutate(
    model_label = factor(model_label, levels = model_registry$model_label),
    draw_panel = factor(draw_panel, levels = paste0("Posterior draw ", seq_len(n_draws_to_plot))),
    distribution = factor(distribution, levels = c("Observed", "Posterior predictive"))
  )

p_hist_compare <- ggplot(plot_data) +
  geom_histogram(
    aes(x = value, y = after_stat(density), fill = distribution),
    binwidth = 0.5,
    boundary = 0,
    position = "identity",
    alpha = 0.5,
    color = "grey30",
    linewidth = 0.15
  ) +
  facet_grid(draw_panel ~ model_label) +
  scale_fill_manual(values = c("Observed" = "grey25", "Posterior predictive" = "#56B4E9")) +
  xlab("Observed or predicted observation (log2 IC50)") +
  ylab("Density") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_dir, "prediction_histogram_draws_exponential_obs_ver_2_vs_3_ic50_real.pdf"),
  p_hist_compare,
  width = 6,
  height = 7,
  units = "in",
  dpi = 300
)
ggsave(
  file.path(output_dir, "prediction_histogram_draws_exponential_obs_ver_2_vs_3_ic50_real.png"),
  p_hist_compare,
  width = 6,
  height = 7,
  units = "in",
  dpi = 300
)
