## Copy this file to config.R and edit paths locally.
## Do not commit config.R to GitHub.

output_dir <- "results"

cognition_file <- "path/to/analysis_ready_cognition_long.csv"
covariate_file <- "path/to/analysis_ready_covariates_long.csv"
alamar_file <- "path/to/analysis_ready_alamar_long.csv"
domain_file <- "path/to/alamar_domain_mapping.csv"

id_var <- "idno"
exam_var <- "exam"
age_var <- "age"
gcc_var <- "GCC"
cog_valid_var <- "cogvalid"
kit_var <- "kit_sheet"
outlier_var <- "sample_outlier"

abeta40_var <- "Abeta40"
abeta42_var <- "Abeta42"
abeta_ratio_var <- "abeta42_40_ratio"
apoe4_protein_var <- "APOE4"

covariates_part1 <- c(
  "gender", "race", "edu_r", "site", "egfr_binary",
  "obesity", "hypertension", "diabetes", "apoe_e4_carrier"
)

covariates_part2_part3 <- c(
  "gender", "race", "edu_r", "site", "egfr_binary",
  "obesity", "hypertension", "diabetes", "apoe4_carrier_from_npq", "kit_sheet"
)

age_center <- 70
spline_df <- 3
ipw_lower_quantile <- 0.01
ipw_upper_quantile <- 0.99
