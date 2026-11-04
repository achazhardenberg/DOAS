# DOAS: Double Observer Adjusted Survey

<!-- badges: start -->
<!-- badges: end -->

The **DOAS** R package provides a unified Frequentist and Hierarchical Bayesian statistical framework to estimate wildlife population abundance by combining standard, low-effort single-observer block counts across a study area with targeted Double Observer (DO) surveys in a subset of spatial sectors.

By calibrating detection probability on the Double Observer sectors—and optionally modeling sector-level detection covariates—`DOAS` corrects for imperfect detection ($p < 1$) across the entire landscape while substantially reducing field personnel requirements.

## Installation

You can install the development version of `DOAS` from GitHub using `devtools`:

```r
# install.packages("devtools")
devtools::install_github("achazhardenberg/DOAS", build_vignettes = TRUE)
```

## Overview of Features

### 1. Flexible Data Ingestion & Validation
- **`doas_data()`**: Validates and standardizes survey data from three field recording formats into a `doas_data` object (with a dedicated `print()` summary):
  - `format = "sector"`: Wide sector-level summaries (`total_count`, `detected_groups`, `groups_both`, `groups_obs1`, `groups_obs2`). Automatically infers `detected_groups` or `groups_both` if one is omitted in DO sectors and verifies arithmetic consistency.
  - `format = "group"`: Long group-level capture histories (one row per detected animal group with `obs1` and `obs2` binary detection indicators).
  - `format = "sighting"`: Long observer-sighting records reconciled by `group_id` and `observer`.

### 2. Frequentist & Hierarchical Bayesian Abundance Estimators
- **`doas()`**: Core abundance estimator controlled by three orthogonal design arguments:
  - **`method = c("frequentist", "bayesian")`**: Choose between Frequentist estimation (with analytical multivariate Delta-method or non-parametric spatial block `"bootstrap"` confidence intervals) and multi-chain Hierarchical Bayesian MCMC estimation (reporting posterior means, 95% credible intervals, Gelman–Rubin $\hat{R}$, and effective sample size $n_{\text{eff}}$).
  - **`obs_type = c("independent", "dependent")`**: Supports both simultaneous independent observers (Chapman mark–recapture with symmetric observer pool efficiency $\bar{p} = (p_1 + p_2)/2$) and dependent primary–secondary observer protocols (Zippin removal).
  - **`pooling = c("hybrid", "global")`**: `"hybrid"` (default) retains direct mark–recapture estimates inside surveyed DO sectors and extrapolates only single-observer Block sectors; `"global"` applies the calibrated single-observer detection rate across all sectors.
  - **`formula = ~ covariate`**: Models sector-level variation in detection probability using covariates (e.g., survey effort, terrain ruggedness, visibility).
- **Low-Level Building Blocks**:
  - `chapman_estimate()` & `chapman_variance()`: Bias-corrected independent double-observer estimator and Seber (1982) variance.
  - `dependent_zippin_estimate()` & `dependent_zippin_variance()`: Primary–secondary observer removal estimator and variance.
  - `bayesian_conjugate_estimate()`: Bound-safe conjugate Beta–Binomial Monte Carlo abundance estimator for individual DO sectors.

### 3. Information-Theoretic Model Selection
- **`QAIC()`**: Computes the Quasi-Akaike Information Criterion ($\hat{c}$, $K$, $\Delta\text{QAIC}$, and Akaike weights) to rank competing frequentist `doas` covariate models while accounting for overdispersion.
- **`WAIC()`**: Computes the Widely Applicable Information Criterion ($\text{lppd}$, $p_{\text{WAIC}}$, $\Delta\text{WAIC}$, and model weights) from the pointwise posterior log-likelihood matrix of Bayesian `doas` models.

### 4. Simulation & Prospective Power Analysis
- **`simulate_pop()`**: Generates spatially structured animal populations partitioned into groups across sectors.
- **`simulate_survey()`**: Simulates partial Double Observer and Block count surveys with group-size scaling, counting error, inter-observer skill differences (`obs_effect`), inter-sector heterogeneity (`detvar`), and sector covariate effects (`cov_effect`).
- **`doas_power()`**: Conducts Monte Carlo power and accuracy analyses across varying numbers of Double Observer sectors (`surv_areas_vec`) for any combination of `method`, `obs_type`, `pooling`, and `formula`.

## Quick Start Example

```r
library(DOAS)

# 1. Prepare and validate sector-level survey data
raw_data <- data.frame(
  sector_id       = c("A", "B", "C", "D", "E"),
  total_count     = c(55, 30, 42, 80, 100),
  detected_groups = c(8, 5, 6, 12, 10),
  groups_both     = c(6, 3, 4, NA, NA),
  groups_obs1     = c(1, 1, 1, NA, NA),
  groups_obs2     = c(1, 1, 1, NA, NA),
  effort          = c(1.2, -0.5, 0.3, -0.8, 0.9)
)

survey_dat <- doas_data(raw_data, format = "sector")
print(survey_dat)

# 2. Frequentist Hybrid Estimation (Null vs. Covariate Model)
fit_freq_null <- doas(
  survey_dat,
  method   = "frequentist",
  obs_type = "independent",
  pooling  = "hybrid"
)

fit_freq_cov <- doas(
  survey_dat,
  method   = "frequentist",
  obs_type = "independent",
  pooling  = "hybrid",
  formula  = ~ effort
)

summary(fit_freq_null)
QAIC(fit_freq_null, fit_freq_cov)

# 3. Multi-Chain Hierarchical Bayesian Estimation & WAIC Selection
fit_bayes_null <- doas(
  survey_dat,
  method   = "bayesian",
  obs_type = "independent",
  pooling  = "hybrid",
  chains   = 3,
  iter     = 2000,
  burnin   = 500
)

fit_bayes_cov <- doas(
  survey_dat,
  method   = "bayesian",
  obs_type = "independent",
  pooling  = "hybrid",
  formula  = ~ effort,
  chains   = 3,
  iter     = 2000,
  burnin   = 500
)

summary(fit_bayes_null)
WAIC(fit_bayes_null, fit_bayes_cov)
```

## Documentation & Tutorials

A step-by-step tutorial vignette covering field data preparation, frequentist and Bayesian estimation, covariate selection, and prospective power analysis is included with the package:

```r
vignette("DOAS-tutorial", package = "DOAS")
```

## Citation & Methods

The Double Observer Adjusted Survey framework was originally developed in the publication below. If you use `DOAS` in your research, please cite:

> Panaccio, M., Brambilla, A., Bassano, B., Smith, T., & von Hardenberg, A. (2024). A new double observer-based census framework to improve abundance estimations in mountain ungulates and other gregarious species with a reduced effort. *Ecological Solutions and Evidence*, 5(4), e12405. https://doi.org/10.1002/2688-8319.12405
