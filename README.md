# DOAS: Double Observer Adjusted Survey 🔭

<!-- badges: start -->
<!-- badges: end -->

The **DOAS** R package provides a statistical framework to estimate true animal abundance by combining standard, low-effort block counts with a targeted subset of intensive Double Observer (DO) surveys. 

By leveraging mark-recapture principles on the Double Observer sectors, `DOAS` can dynamically estimate and correct for imperfect detection probabilities across the entire study area.

## Installation

You can install the development version of DOAS from GitHub using the `devtools` package:

```r
# install.packages("devtools")
devtools::install_github("achazhardenberg/DOAS", build_vignettes = TRUE)
```

## Overview of Features

The package includes both analysis functions for real data and simulation functions for study design and power analysis:

- **Estimators for Field Data:**
  - `doas()`: The core abundance estimator for applying Double Observer detection adjustments to block count survey data.
  - `chapman_estimate()`: Calculates independent observer mark-recapture estimates.
  - `dependent_zippin_estimate()`: Calculates primary-secondary (dependent) observer removal estimates.
  - `bayesian_conjugate_estimate()`: Bound-safe abundance estimation using Beta-Binomial posteriors.
  
- **Simulation and Power Analysis:**
  - `simulate_pop()`: Generates geographically structured true animal populations.
  - `simulate_survey()`: Simulates field sampling with realistic intra-group and inter-observer detection errors.
  - `doas_power()`: Conducts comprehensive Monte Carlo power analyses to optimize observer deployment efforts.

## Quick Example

### Formatting your data

Data passed to `doas()` must be summarized at the sector level. If a sector received a Double Observer survey, include the `groups_both`, `groups_obs1`, and `groups_obs2` columns. Otherwise, leave them as `NA`.

```r
library(DOAS)

survey_data <- data.frame(
  sector_id = c("A", "B", "C", "D"),
  total_count = c(55, 30, 80, 100),
  detected_groups = c(8, 5, 12, 10),
  groups_both = c(6, 3, NA, NA),
  groups_obs1 = c(1, 1, NA, NA),
  groups_obs2 = c(1, 1, NA, NA)
)

# Run the hybrid estimator
results <- doas(survey_data, method = "hybrid", obs_type = "independent")

summary(results)
```

## Documentation & Tutorials

The package comes with a comprehensive, step-by-step tutorial vignette that covers both analyzing real-world datasets and running programmatic simulations. 

You can view the vignette in R after installation by running:
```r
vignette("DOAS-tutorial", package = "DOAS")
```
