## Part 3. Association between cognitive aging and biomarker aging

source("packages.R")
source("config.R")
source("functions.R")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)

part1 <- readRDS(file.path(output_dir, "part1_cognitive_aging_results.rds"))
part2 <- readRDS(file.path(output_dir, "part2_biomarker_aging_results.rds"))

cog_gaps <- part1$data |>
  dplyr::select(idno, exam, GCC_age_gap)

protein_gaps <- part2$biomarker_gaps |>
  dplyr::left_join(
    part2$protein_long |> dplyr::select(idno, exam, protein_name, dplyr::all_of(covariates_part2_part3)),
    by = c("idno", "exam", "protein_name")
  ) |>
  dplyr::distinct()

analysis_data <- protein_gaps |>
  dplyr::left_join(cog_gaps, by = c("idno", "exam")) |>
  dplyr::filter(!is.na(GCC_age_gap), !is.na(biomarker_age_gap)) |>
  dplyr::group_by(protein_name) |>
  dplyr::mutate(biomarker_age_gap_z = as.numeric(scale(biomarker_age_gap))) |>
  dplyr::ungroup()

## Attrition IPW.
## Baseline is defined as Exam 5 among participants contributing an analytic row.
## The denominator model estimates the probability of observation at each follow-up
## window using baseline covariates. Weights are truncated at the configured
## 1st and 99th percentiles.
baseline <- analysis_data |>
  dplyr::filter(exam == 5) |>
  dplyr::distinct(idno, .keep_all = TRUE) |>
  dplyr::select(idno, dplyr::all_of(covariates_part2_part3))

observed_windows <- analysis_data |>
  dplyr::distinct(idno, exam) |>
  dplyr::mutate(observed = 1)

ipw_data <- tidyr::expand_grid(idno = baseline$idno, exam = sort(unique(analysis_data$exam))) |>
  dplyr::left_join(observed_windows, by = c("idno", "exam")) |>
  dplyr::mutate(observed = ifelse(is.na(observed), 0, observed)) |>
  dplyr::left_join(baseline, by = "idno")

denominator_formula <- as.formula(paste("observed ~", paste(covariates_part2_part3, collapse = " + ")))
numerator_formula <- observed ~ 1

ipw_data$prob_denom <- NA_real_
ipw_data$prob_num <- NA_real_
for (visit in sort(unique(ipw_data$exam))) {
  dat <- ipw_data |> dplyr::filter(exam == visit)
  den <- glm(denominator_formula, data = dat, family = binomial())
  num <- glm(numerator_formula, data = dat, family = binomial())
  ipw_data$prob_denom[ipw_data$exam == visit] <- predict(den, type = "response")
  ipw_data$prob_num[ipw_data$exam == visit] <- predict(num, type = "response")
}

ipw_data <- ipw_data |>
  dplyr::mutate(
    ipw = prob_num / prob_denom,
    ipw_truncated = pmin(
      pmax(ipw, stats::quantile(ipw, ipw_lower_quantile, na.rm = TRUE)),
      stats::quantile(ipw, ipw_upper_quantile, na.rm = TRUE)
    )
  )

analysis_data <- analysis_data |>
  dplyr::left_join(ipw_data |> dplyr::select(idno, exam, ipw_truncated), by = c("idno", "exam"))

model_rows <- list()
for (protein in unique(analysis_data$protein_name)) {
  dat <- analysis_data |>
    dplyr::filter(protein_name == protein) |>
    dplyr::select(idno, GCC_age_gap, biomarker_age_gap_z, ipw_truncated, dplyr::all_of(covariates_part2_part3)) |>
    dplyr::filter(stats::complete.cases(.))

  if (nrow(dat) < 50) next

  formula <- as.formula(paste("GCC_age_gap ~ biomarker_age_gap_z +", paste(covariates_part2_part3, collapse = " + "), "+ (1 | idno)"))
  fit <- fit_lmer(formula, dat, weights = dat$ipw_truncated, reml = FALSE)
  if (is.null(fit)) next

  model_rows[[protein]] <- extract_fixed_terms(fit) |>
    dplyr::filter(term == "biomarker_age_gap_z") |>
    dplyr::mutate(
      protein_name = protein,
      display_label = label_biomarker(protein),
      n_participants = dplyr::n_distinct(dat$idno),
      n_observations = nrow(dat)
    )
}

results <- dplyr::bind_rows(model_rows) |>
  dplyr::mutate(q_value = p.adjust(p.value, method = "fdr"))

make_volcano(
  results,
  beta_col = "estimate",
  p_col = "p.value",
  q_col = "q_value",
  title = "Associations between cognition age gap and CNS biomarker age gap",
  xlab = "Beta for 1-SD higher biomarker age gap",
  file = file.path(output_dir, "figures", "part3_cognition_biomarker_age_gap_volcano.png")
)

ipw_summary <- data.frame(
  n_rows = nrow(ipw_data),
  min = min(ipw_data$ipw_truncated, na.rm = TRUE),
  p1 = stats::quantile(ipw_data$ipw_truncated, 0.01, na.rm = TRUE),
  median = stats::median(ipw_data$ipw_truncated, na.rm = TRUE),
  p99 = stats::quantile(ipw_data$ipw_truncated, 0.99, na.rm = TRUE),
  max = max(ipw_data$ipw_truncated, na.rm = TRUE)
)

write_xlsx(
  list(gap_associations = results, ipw_summary = ipw_summary),
  file.path(output_dir, "part3_cognition_biomarker_age_gap_results.xlsx")
)
