# Simulate a Population

Simulates a population of individuals split into groups across multiple
sectors. Supports sector-level abundance covariate modeling via
`n_formula` and named coefficient vectors `n_betas`, or user-provided
sector covariates.

## Usage

``` r
simulate_pop(
  pop_size,
  mean_group,
  tot_areas,
  area_sizes = NULL,
  n_formula = NULL,
  n_betas = NULL,
  covariates = NULL,
  n_cov_effect = 0
)
```

## Arguments

- pop_size:

  Average total population size (Poisson distributed).

- mean_group:

  Average group size (Poisson distributed).

- tot_areas:

  Total number of survey sectors.

- area_sizes:

  Vector of relative probabilities/sizes for assigning groups to
  sectors. If NULL, probabilities are computed from `n_formula` and
  `n_betas` (or `n_cov_effect`).

- n_formula:

  Optional one-sided formula for sector abundance covariates (e.g.
  `~ elevation + forest` or `~ elevation * forest`).

- n_betas:

  Optional named numeric vector of log-linear coefficients for terms in
  `n_formula` (e.g., `c(elevation = 0.5, forest = -0.3)`).

- covariates:

  Optional data frame of sector-level covariates containing `tot_areas`
  rows. If omitted, any variables specified in `n_formula` are simulated
  as standard normal variates.

- n_cov_effect:

  Optional numeric log-linear slope for a single simulated `habitat`
  covariate (retained for backward compatibility; equivalent to
  `n_formula = ~ habitat, n_betas = c(habitat = n_cov_effect)`).

## Value

A list containing:

- tot_pop:

  Total simulated individuals (scalar).

- sim_pop:

  List of groups per sector.

- covariates:

  Data frame of sector-level covariates (one row per sector).

- habitat:

  Numeric vector of the `habitat` covariate (if present).
