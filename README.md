# Alamar CNS Biomarker Aging and Cognitive Aging

This folder contains simplified R code for the three main manuscript analyses:

1. cognitive aging;
2. Alamar CNS biomarker aging and trajectory clustering;
3. associations between cognition age gaps and biomarker age gaps.

Participant-level MESA data are not included. The scripts require analysis-ready datasets created in an approved secure environment.

## Files

```text
packages.R
functions.R
config_template.R
01_cognitive_aging.R
02_biomarker_aging_and_clustering.R
03_cognition_biomarker_age_gap.R
PUBLIC_RELEASE_CHECKLIST.md
.gitignore
```

## How to Run

All analyses should be performed in R version 4.4.1.

1. Copy `config_template.R` to `config.R`.
2. Edit `config.R` with local paths to approved analysis-ready datasets.
3. Run the scripts in order:

```r
source("01_cognitive_aging.R")
source("02_biomarker_aging_and_clustering.R")
source("03_cognition_biomarker_age_gap.R")
```

## Required Inputs

The code assumes four analysis-ready files:

- cognition long file: one row per participant-visit with age, exam/window, GCC, and cognitive validity indicator;
- covariate long file: one row per participant-visit with demographic and clinical covariates;
- Alamar biomarker file: one row per participant-visit with biomarker columns and plate/kit information;
- domain mapping file: one row per biomarker with the assigned Alamar protein domain.

Expected variables are listed in `config_template.R`. Dataset derivation, cleaning, and restricted-data linkage scripts are not included.

## Main Methods Encoded

Part 1 fits repeated GCC models using age as the time scale, compares linear age and natural spline age functions using likelihood ratio tests, and derives cognition age gaps from age-only mixed-effects models.

Part 2 fits repeated standardized biomarker models using age as the time scale, compares linear and spline age functions, generates marginal age trajectories, and clusters biomarkers using Euclidean distance and complete-linkage hierarchical clustering of re-standardized predicted trajectories from ages 55 to 95.

Part 3 links cognitive aging and biomarker aging by modeling standardized biomarker age gaps in relation to standardized GCC age gaps, using attrition inverse probability weights.

## Data Availability

MESA participant-level data are not provided in this repository. Qualified investigators must obtain access through the MESA Coordinating Center and follow all applicable data-use agreements.
