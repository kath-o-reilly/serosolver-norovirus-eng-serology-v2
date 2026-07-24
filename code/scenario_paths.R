## Modified by AI assistant on 2026-07-23 using GPT-5. Centralizes output
## paths for model-version, false-positive-rate, and sensitivity analyses, with
## MCMC figures written to the top-level figures directory.

format_fp_rate_for_path <- function(fp_rate_set) {
  format(fp_rate_set, scientific = FALSE, trim = TRUE)
}

build_scenario_id <- function(analysis_family,
                              waning_model,
                              obs_ver,
                              fp_rate_set = NULL,
                              sensitivity_name = NULL) {
  if (!analysis_family %in% c("model_versions", "fprate_scenarios", "sensitivity_analyses")) {
    stop("analysis_family must be 'model_versions', 'fprate_scenarios', or 'sensitivity_analyses'.")
  }

  if (analysis_family == "fprate_scenarios") {
    if (is.null(fp_rate_set)) {
      stop("fp_rate_set is required for fprate_scenarios.")
    }
    return(paste0(waning_model, "_obs_ver_", obs_ver, "_fprate_", format_fp_rate_for_path(fp_rate_set)))
  }

  if (analysis_family == "sensitivity_analyses") {
    if (is.null(sensitivity_name)) {
      stop("sensitivity_name is required for sensitivity_analyses.")
    }
    return(paste0(waning_model, "_obs_ver_", obs_ver, "_", sensitivity_name))
  }

  paste0(waning_model, "_obs_ver_", obs_ver)
}

build_scenario_paths <- function(analysis_family,
                                 waning_model,
                                 obs_ver,
                                 run_name,
                                 fp_rate_set = NULL,
                                 project_dir = getwd(),
                                 sim_type = NULL,
                                 sensitivity_name = NULL) {
  scenario_id <- build_scenario_id(
    analysis_family = analysis_family,
    waning_model = waning_model,
    obs_ver = obs_ver,
    fp_rate_set = fp_rate_set,
    sensitivity_name = sensitivity_name
  )

  local_sim_dir <- file.path(project_dir, "local_data", analysis_family, scenario_id)
  results_run_dir <- file.path(project_dir, "results", analysis_family, scenario_id, run_name)
  figures_run_root <- file.path(project_dir, "figures", analysis_family, scenario_id)

  list(
    scenario_id = scenario_id,
    local_sim_dir = local_sim_dir,
    chain_wd = file.path(local_sim_dir, "chains", run_name),
    save_wd = file.path(figures_run_root, run_name),
    results_run_dir = results_run_dir,
    results_prefix = file.path(results_run_dir, run_name),
    figures_run_root = figures_run_root,
    kinetics_wd = file.path(figures_run_root, "kinetics", run_name),
    model_schematic_wd = file.path(figures_run_root, "model_schematic", run_name),
    incidence_wd = file.path(figures_run_root, "incidence", run_name),
    infection_trends_wd = file.path(figures_run_root, "infection_trends", run_name)
  )
}
