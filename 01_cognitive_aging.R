## Part 1. Cognitive aging

source("packages.R")
source("config.R")
source("functions.R")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "figures"), showWarnings = FALSE)

cognition <- readr::read_csv(cognition_file, show_col_types = FALSE)
covariates <- readr::read_csv(covariate_file, show_col_types = FALSE)

check_columns(cognition, c(id_var, exam_var, age_var, gcc_var, cog_valid_var), "cognition_file")
check_columns(covariates, c(id_var, exam_var, covariates_part1), "covariate_file")

cog <- cognition |>
  dplyr::rename(
    idno = dplyr::all_of(id_var),
    exam = dplyr::all_of(exam_var),
    age = dplyr::all_of(age_var),
    GCC = dplyr::all_of(gcc_var),
    cogvalid = dplyr::all_of(cog_valid_var)
  ) |>
  dplyr::left_join(
    covariates |> dplyr::rename(idno = dplyr::all_of(id_var), exam = dplyr::all_of(exam_var)),
    by = c("idno", "exam")
  ) |>
  dplyr::filter(cogvalid == 1, !is.na(GCC), !is.na(age)) |>
  dplyr::mutate(
    idno = factor(idno),
    GCC_z = as.numeric(scale(GCC))
  )

for (v in covariates_part1) cog[[v]] <- clean_factor(cog[[v]])
cog <- add_age_spline(cog, age_center = age_center, df = spline_df)

linear_formula <- as.formula(paste("GCC_z ~ age_centered_70 +", paste(covariates_part1, collapse = " + "), "+ (1 | idno)"))
spline_formula <- as.formula(paste("GCC_z ~ age_ns1 + age_ns2 + age_ns3 +", paste(covariates_part1, collapse = " + "), "+ (1 | idno)"))

fit_linear <- fit_lmer(linear_formula, cog, reml = FALSE)
fit_spline <- fit_lmer(spline_formula, cog, reml = FALSE)
lrt <- anova(fit_linear, fit_spline)
nonlinearity_p <- lrt$`Pr(>Chisq)`[2]
selected_model <- if (!is.na(nonlinearity_p) && nonlinearity_p < 0.05) fit_spline else fit_linear
selected_age_form <- if (identical(selected_model, fit_spline)) "natural spline" else "linear"

age_only_formula <- if (selected_age_form == "natural spline") {
  GCC_z ~ age_ns1 + age_ns2 + age_ns3 + (1 | idno)
} else {
  GCC_z ~ age_centered_70 + (1 | idno)
}
age_only_model <- fit_lmer(age_only_formula, cog, reml = FALSE)
cog$GCC_expected_for_age <- predict(age_only_model, newdata = cog, re.form = NA)
cog$GCC_age_gap <- cog$GCC_z - cog$GCC_expected_for_age

age_grid <- 55:95
pred <- marginal_age_prediction(selected_model, cog, age_grid, age_center)

p <- ggplot2::ggplot(cog, ggplot2::aes(age, GCC_z, group = idno)) +
  ggplot2::geom_line(color = "gray80", alpha = 0.35, linewidth = 0.25) +
  ggplot2::geom_point(ggplot2::aes(color = factor(exam)), alpha = 0.60, size = 1) +
  ggplot2::geom_line(data = pred, ggplot2::aes(age, estimate, group = 1), inherit.aes = FALSE, linewidth = 1.2, color = "black") +
  ggplot2::labs(
    title = "Longitudinal trajectories of global cognitive composite across age",
    x = "Age, years",
    y = "Standardized GCC",
    color = "Exam"
  ) +
  ggplot2::theme_classic(base_size = 12)
ggplot2::ggsave(file.path(output_dir, "figures", "part1_gcc_by_age.png"), p, width = 7, height = 5, dpi = 300)

model_selection <- data.frame(
  outcome = "GCC",
  selected_age_form = selected_age_form,
  nonlinearity_p = nonlinearity_p,
  n_participants = dplyr::n_distinct(cog$idno),
  n_observations = nrow(cog)
)

write_xlsx(
  list(model_selection = model_selection, fixed_effects = extract_fixed_terms(selected_model)),
  file.path(output_dir, "part1_cognitive_aging_results.xlsx")
)

saveRDS(
  list(data = cog, selected_model = selected_model, age_only_model = age_only_model, model_selection = model_selection),
  file.path(output_dir, "part1_cognitive_aging_results.rds")
)
