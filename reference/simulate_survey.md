# Simulate a Survey and format for DOAS

Simulates field observer detections across survey sectors with support
for detection covariates modeled via `p_formula` and named coefficients
`p_betas`.

## Usage

``` r
simulate_survey(
  sim_pop,
  det_p1,
  intra_group_det,
  det_group,
  obs_effect,
  detvar,
  surv_areas,
  p_formula = NULL,
  p_betas = NULL,
  covariates = NULL,
  cov_effect = 0
)
```

## Arguments

- sim_pop:

  A simulated population from
  [`simulate_pop`](https://achazhardenberg.github.io/DOAS/reference/simulate_pop.md)
  (either the returned list or `pop$sim_pop`).

- det_p1:

  Base detection probability for Observer 1.

- intra_group_det:

  Mean factor for count errors (e.g., 0.9 for undercounting).

- det_group:

  Base detection factor for group size scaling.

- obs_effect:

  Observer 2 relative efficiency (e.g., 0.75 means Obs 2 is 75 percent
  as effective).

- detvar:

  Standard deviation for detection probability across sectors.

- surv_areas:

  Number of sectors to receive Double Observer surveys.

- p_formula:

  Optional one-sided formula for detection probability (e.g.
  `~ effort + elevation`).

- p_betas:

  Optional named numeric vector of logit-scale coefficients for
  `p_formula` (e.g., `c(effort = 0.6, elevation = -0.3)`).

- covariates:

  Optional data frame of sector-level covariates. If omitted, extracted
  from `sim_pop`.

- cov_effect:

  Optional numeric log-odds slope for a single simulated `covariate`
  (retained for backward compatibility; equivalent to
  `p_formula = ~ covariate, p_betas = c(covariate = cov_effect)`).

## Value

A data frame with one row per sector, containing observed counts, Double
Observer tallies, and all simulated sector-level covariates ready for
input to
[`doas`](https://achazhardenberg.github.io/DOAS/reference/doas.md).
