required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "ggrepel",
  "lme4",
  "lmerTest",
  "broom.mixed",
  "openxlsx",
  "splines"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Please install required packages: ", paste(missing_packages, collapse = ", "))
}

invisible(lapply(required_packages, library, character.only = TRUE))
