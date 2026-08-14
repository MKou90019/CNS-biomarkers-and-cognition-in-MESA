check_columns <- function(data, vars, data_name) {
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop(data_name, " is missing: ", paste(missing_vars, collapse = ", "))
  }
}

clean_factor <- function(x) {
  factor(ifelse(is.na(x), "Missing", as.character(x)))
}

label_biomarker <- function(x) {
  dplyr::case_when(
    x == "NEFL" ~ "NfL",
    x == "abeta42_40_ratio" ~ "Aβ42/40",
    TRUE ~ x
  )
}

fit_lmer <- function(formula, data, weights = NULL, reml = FALSE) {
  tryCatch(
    lmerTest::lmer(formula, data = data, weights = weights, REML = reml),
    error = function(e) NULL
  )
}

extract_fixed_terms <- function(model) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) |>
    dplyr::filter(term != "(Intercept)")
}

add_age_spline <- function(data, age_center = 70, df = 3) {
  data$age_centered_70 <- data$age - age_center
  basis <- splines::ns(data$age_centered_70, df = df)
  colnames(basis) <- paste0("age_ns", seq_len(ncol(basis)))
  dplyr::bind_cols(data, as.data.frame(basis))
}

marginal_age_prediction <- function(model, data, age_grid, age_center = 70) {
  out <- lapply(age_grid, function(a) {
    nd <- data
    nd$age <- a
    nd$age_centered_70 <- a - age_center
    if (all(c("age_ns1", "age_ns2", "age_ns3") %in% names(nd))) {
      basis <- splines::ns(data$age_centered_70, df = 3)
      ns_terms <- predict(basis, newx = nd$age_centered_70)
      nd$age_ns1 <- ns_terms[, 1]
      nd$age_ns2 <- ns_terms[, 2]
      nd$age_ns3 <- ns_terms[, 3]
    }
    data.frame(age = a, estimate = mean(predict(model, newdata = nd, re.form = NA), na.rm = TRUE))
  })
  dplyr::bind_rows(out)
}

make_volcano <- function(data, beta_col, p_col, q_col, title, xlab, file) {
  p <- data |>
    dplyr::mutate(
      significant = ifelse(.data[[q_col]] < 0.05, "Q-value < 0.05", "Q-value >= 0.05"),
      neg_log10_p = -log10(.data[[p_col]])
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[beta_col]], y = neg_log10_p)) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    ggplot2::geom_point(ggplot2::aes(color = significant), alpha = 0.85, size = 2) +
    ggrepel::geom_text_repel(
      data = function(d) dplyr::filter(d, .data[[q_col]] < 0.05),
      ggplot2::aes(label = display_label),
      size = 3,
      max.overlaps = 20
    ) +
    ggplot2::scale_color_manual(values = c("Q-value < 0.05" = "#b2182b", "Q-value >= 0.05" = "gray70"), name = "Significance") +
    ggplot2::labs(title = title, x = xlab, y = "-log10(P-value)") +
    ggplot2::theme_classic(base_size = 12)
  ggplot2::ggsave(file, p, width = 7, height = 6, dpi = 300)
  p
}

write_xlsx <- function(named_tables, file) {
  wb <- openxlsx::createWorkbook()
  for (nm in names(named_tables)) {
    sheet <- substr(nm, 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, named_tables[[nm]])
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}
