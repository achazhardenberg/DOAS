# Double Observer Adjusted Survey (DOAS) Abundance Estimation

Estimates total population abundance using either a full Double Observer
(DO) survey across all sectors or a partial Double Observer Adjusted
Survey (DOAS) combining DO surveys in a subset of sectors with
single-observer block counts (and optionally unsurveyed sectors imputed
via `n_formula`). Supports both Frequentist and multi-chain Bayesian
MCMC inference across independent (Lincoln-Petersen/Chapman) and
dependent (2-pass removal/Zippin) double-observer field protocols.

## Usage

``` r
doas(
  data,
  method = c("frequentist", "bayesian"),
  obs_type = c("independent", "dependent"),
  pooling = c("hybrid", "global"),
  p_formula = NULL,
  n_formula = NULL,
  ci_method = c("lognormal", "bootstrap"),
  R = 1000,
  chains = 4,
  iter = 2000,
  burnin = NULL
)
```

## Arguments

- data:

  A dataframe or
  [`doas_data`](https://achazhardenberg.github.io/DOAS/reference/doas_data.md)
  object containing survey data. Each row represents a sector: -
  `sector_id`: Unique identifier for the sector. - `total_count`: Total
  individuals counted across all groups in the sector (`NA` if
  unsurveyed). - `detected_groups`: Total groups detected in the
  sector. - `groups_both`, `groups_obs1`, `groups_obs2`: Double Observer
  group counts (non-NA for DO sectors, `NA` for Block/Unsurveyed
  sectors).

- method:

  Character string specifying the statistical inference paradigm: -
  `"frequentist"` (default): Uses analytic Chapman (when
  `obs_type = "independent"`) or Zippin removal (when
  `obs_type = "dependent"`) estimators and Delta-method variance. -
  `"bayesian"`: Uses a multi-chain Collapsed Poisson-Gamma-Beta Gibbs
  Sampler (MCMC) to sample the full joint posterior of detection
  probabilities and abundance.

- obs_type:

  Character string specifying the field observer protocol: -
  `"independent"` (default): Both observers scan independently
  (Lincoln-Petersen / Chapman geometry). - `"dependent"`: Primary
  observer calls out detected groups and Secondary observer records only
  missed groups (2-pass removal / Zippin geometry).

- pooling:

  Character string specifying how DO sector estimates are combined with
  Block sectors: - `"hybrid"` (default): Uses sector-specific DO
  estimates for DO sectors and extrapolates only Block sectors using the
  pooled/modeled detection probability. - `"global"`: Extrapolates all
  surveyed sectors using the pooled or covariate-modeled detection
  probability.

- p_formula:

  Optional one-sided formula (e.g., `~ effort * ruggedness`) to model
  detection probability \\\hat{p}\\ across sectors. For
  `method = "frequentist"`, fits a quasi-binomial GLM; for
  `method = "bayesian"`, runs a coupled multi-chain Bayesian logistic
  regression via MCMC.

- n_formula:

  Optional one-sided formula (e.g.,
  `~ elevation + habitat + offset(log(area))`) to model sector-level
  abundance \\\lambda_i\\ as a function of environmental covariates
  (Quasi-Poisson log-link GLM in frequentist mode; Bayesian Poisson
  log-linear regression in Bayesian mode). Enables
  [`predict.doas`](https://achazhardenberg.github.io/DOAS/reference/predict.doas.md)
  for new/unsurveyed sectors and habitat response curves.

- ci_method:

  Character string specifying how to calculate Confidence/Credible
  Intervals: - `"lognormal"` (default): Uses the analytic Chao (1987)
  log-normal approximation for frequentist estimators, or direct
  posterior credible intervals from the MCMC draws when
  `method = "bayesian"`. - `"bootstrap"`: Uses a non-parametric
  stratified sector resampling bootstrap.

- R:

  Integer indicating the number of bootstrap replicates (only used if
  `ci_method = "bootstrap"`). Default is 1000.

- chains:

  Integer number of Markov chains when `method = "bayesian"` (default
  4).

- iter:

  Integer number of total MCMC iterations per chain when
  `method = "bayesian"` (default 2000).

- burnin:

  Integer number of warm-up/burn-in iterations per chain when
  `method = "bayesian"` (default `floor(iter / 4)`).

## Value

A list of class \`doas\` containing the estimate, confidence/credible
intervals, sector-level estimates, fitted submodels, and (when Bayesian)
MCMC convergence diagnostics.
