## Part 2. Alamar CNS biomarker aging and trajectory clustering

source("packages.R")
source("config.R")
source("functions.R")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)

alamar <- readr::read_csv(alamar_file, show_col_types = FALSE)
covariates <- readr::read_csv(covariate_file, show_col_types = FALSE)
domains <- readr::read_csv(domain_file, show_col_types = FALSE)

check_columns(alamar, c(id_var, exam_var, kit_var), "alamar_file")
check_columns(covariates, c(id_var, exam_var, age_var, covariates_part2_part3), "covariate_file")
check_columns(domains, c("protein_name", "domain"), "domain_file")

alamar <- alamar |>
  dplyr::rename(idno = dplyr::all_of(id_var), exam = dplyr::all_of(exam_var), kit_sheet = dplyr::all_of(kit_var))

if (outlier_var %in% names(alamar)) {
  alamar <- alamar |> dplyr::filter(is.na(.data[[outlier_var]]) | .data[[outlier_var]] == FALSE)
}
if (all(c(abeta42_var, abeta40_var) %in% names(alamar))) {
  alamar[[abeta_ratio_var]] <- alamar[[abeta42_var]] - alamar[[abeta40_var]]
}

nonprotein <- c("idno", "exam", "kit_sheet", outlier_var, "targetName", "sample_median", "sample_iqr")
protein_cols <- setdiff(names(alamar), nonprotein)
protein_cols <- protein_cols[vapply(alamar[protein_cols], is.numeric, logical(1))]
protein_cols <- setdiff(protein_cols, c(abeta40_var, abeta42_var, apoe4_protein_var))

bio <- alamar |>
  dplyr::select(idno, exam, kit_sheet, dplyr::all_of(protein_cols)) |>
  dplyr::left_join(
    covariates |> dplyr::rename(idno = dplyr::all_of(id_var), exam = dplyr::all_of(exam_var), age = dplyr::all_of(age_var)),
    by = c("idno", "exam")
  ) |>
  dplyr::filter(!is.na(age)) |>
  dplyr::mutate(idno = factor(idno))

for (v in covariates_part2_part3) bio[[v]] <- clean_factor(bio[[v]])
bio <- add_age_spline(bio, age_center = age_center, df = spline_df)

protein_long <- bio |>
  tidyr::pivot_longer(dplyr::all_of(protein_cols), names_to = "protein_name", values_to = "protein_value") |>
  dplyr::group_by(protein_name) |>
  dplyr::mutate(protein_z = as.numeric(scale(protein_value))) |>
  dplyr::ungroup() |>
  dplyr::filter(is.finite(protein_z))

age_grid <- 55:95
model_rows <- list()
trajectory_rows <- list()
gap_rows <- list()

for (protein in protein_cols) {
  dat <- protein_long |>
    dplyr::filter(protein_name == protein) |>
    dplyr::select(idno, exam, age, age_centered_70, age_ns1, age_ns2, age_ns3, protein_name, protein_z, dplyr::all_of(covariates_part2_part3)) |>
    dplyr::filter(stats::complete.cases(.))

  if (nrow(dat) < 50) next

  linear_formula <- as.formula(paste("protein_z ~ age_centered_70 +", paste(covariates_part2_part3, collapse = " + "), "+ (1 | idno)"))
  spline_formula <- as.formula(paste("protein_z ~ age_ns1 + age_ns2 + age_ns3 +", paste(covariates_part2_part3, collapse = " + "), "+ (1 | idno)"))

  fit_linear <- fit_lmer(linear_formula, dat, reml = FALSE)
  fit_spline <- fit_lmer(spline_formula, dat, reml = FALSE)
  if (is.null(fit_linear) || is.null(fit_spline)) next

  lrt <- anova(fit_linear, fit_spline)
  nonlinearity_p <- lrt$`Pr(>Chisq)`[2]
  selected_model <- if (!is.na(nonlinearity_p) && nonlinearity_p < 0.05) fit_spline else fit_linear
  selected_age_form <- if (identical(selected_model, fit_spline)) "natural spline" else "linear"

  age_only_formula <- if (selected_age_form == "natural spline") {
    protein_z ~ age_ns1 + age_ns2 + age_ns3 + (1 | idno)
  } else {
    protein_z ~ age_centered_70 + (1 | idno)
  }
  age_only_model <- fit_lmer(age_only_formula, dat, reml = FALSE)
  expected <- predict(age_only_model, newdata = dat, re.form = NA)

  model_rows[[protein]] <- data.frame(
    protein_name = protein,
    display_label = label_biomarker(protein),
    selected_age_form = selected_age_form,
    nonlinearity_p = nonlinearity_p,
    n_participants = dplyr::n_distinct(dat$idno),
    n_observations = nrow(dat)
  )
  trajectory_rows[[protein]] <- marginal_age_prediction(selected_model, dat, age_grid, age_center) |>
    dplyr::mutate(protein_name = protein, display_label = label_biomarker(protein))
  gap_rows[[protein]] <- dat |>
    dplyr::mutate(biomarker_expected_for_age = expected, biomarker_age_gap = protein_z - expected) |>
    dplyr::select(idno, exam, protein_name, protein_z, biomarker_expected_for_age, biomarker_age_gap)
}

model_selection <- dplyr::bind_rows(model_rows) |>
  dplyr::mutate(nonlinearity_q = p.adjust(nonlinearity_p, method = "fdr"))
trajectories <- dplyr::bind_rows(trajectory_rows)
biomarker_gaps <- dplyr::bind_rows(gap_rows)

trajectory_matrix <- trajectories |>
  dplyr::select(protein_name, age, estimate) |>
  tidyr::pivot_wider(names_from = age, values_from = estimate) |>
  dplyr::arrange(protein_name)

trajectory_ids <- trajectory_matrix$protein_name
trajectory_values <- as.matrix(trajectory_matrix[, -1])
trajectory_values <- t(apply(trajectory_values, 1, scale))
rownames(trajectory_values) <- trajectory_ids

hc <- hclust(dist(trajectory_values, method = "euclidean"), method = "complete")
clusters <- data.frame(
  protein_name = trajectory_ids,
  cluster_4 = as.integer(cutree(hc, k = 4)[trajectory_ids])
) |>
  dplyr::left_join(domains, by = "protein_name")

plot_data <- trajectories |> dplyr::left_join(clusters, by = "protein_name")
p <- ggplot2::ggplot(plot_data, ggplot2::aes(age, estimate, group = protein_name, color = factor(cluster_4))) +
  ggplot2::geom_line(alpha = 0.45) +
  ggplot2::facet_wrap(~ cluster_4, ncol = 2, labeller = ggplot2::label_both) +
  ggplot2::labs(
    title = "The mean standardized levels of biomarkers by age within each of the four identified clusters",
    x = "Age, years",
    y = "Predicted standardized biomarker level",
    color = "Cluster"
  ) +
  ggplot2::theme_classic(base_size = 12)
ggplot2::ggsave(file.path(output_dir, "figures", "part2_biomarker_trajectory_clusters.png"), p, width = 10, height = 8, dpi = 300)

write_xlsx(
  list(model_selection = model_selection, trajectory_clusters = clusters, trajectory_predictions = trajectories),
  file.path(output_dir, "part2_biomarker_aging_results.xlsx")
)

saveRDS(
  list(protein_long = protein_long, model_selection = model_selection, trajectories = trajectories, clusters = clusters, biomarker_gaps = biomarker_gaps, protein_cols = protein_cols),
  file.path(output_dir, "part2_biomarker_aging_results.rds")
)
