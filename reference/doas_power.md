# DOAS Power Analysis

Runs multiple simulations to evaluate the statistical performance (bias,
precision, coverage) of the DOAS estimator under varying sampling
efforts. Supports multi-covariate simulation and modeling on detection
(`p_formula`, `p_betas`) and abundance (`n_formula`, `n_betas`).

## Usage

``` r
doas_power(
  n_sims = 100,
  surv_areas = c(2, 4, 6),
  pop_size = 1000,
  mean_group = 10,
  tot_areas = 15,
  area_sizes = NULL,
  det_p1 = 0.6,
  intra_group_det = 0.9,
  det_group = 0.9,
  obs_effect = 0.75,
  detvar = 0.05,
  p_formula = NULL,
  p_betas = NULL,
  n_formula = NULL,
  n_betas = NULL,
  covariates = NULL,
  cov_effect = 0,
  n_cov_effect = 0,
  method = c("ML", "Bayes", "frequentist", "bayesian"),
  obs_type = c("independent", "dependent"),
  pooling = c("hybrid", "global"),
  ci_method = c("lognormal", "bootstrap")
)
```

## Arguments

- n_sims:

  Number of simulations to run per scenario.

- surv_areas:

  Vector of Double Observer survey efforts to test (number of DO
  sectors).

- pop_size:

  Average total population size for
  [`simulate_pop`](https://achazhardenberg.github.io/DOAS/reference/simulate_pop.md).

- mean_group:

  Average group size for
  [`simulate_pop`](https://achazhardenberg.github.io/DOAS/reference/simulate_pop.md).

- tot_areas:

  Total number of survey sectors.

- area_sizes:

  Vector of probabilities/sizes for assigning groups to sectors. If
  NULL, uniform or covariate-driven probabilities are used.

- det_p1:

  Base detection probability for Observer 1.

- intra_group_det:

  Mean factor for count errors.

- det_group:

  Base detection factor for group size scaling.

- obs_effect:

  Observer 2 relative efficiency.

- detvar:

  Standard deviation for detection probability across sectors.

- p_formula:

  Optional detection formula passed to simulation and
  [`doas`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
  (e.g., `~ effort + elevation`).

- p_betas:

  Optional named numeric vector of logit coefficients for `p_formula` in
  [`simulate_survey`](https://achazhardenberg.github.io/DOAS/reference/simulate_survey.md).

- n_formula:

  Optional abundance formula passed to simulation and
  [`doas`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
  (e.g., `~ elevation + forest`).

- n_betas:

  Optional named numeric vector of log-linear coefficients for
  `n_formula` in
  [`simulate_pop`](https://achazhardenberg.github.io/DOAS/reference/simulate_pop.md).

- covariates:

  Optional data frame of sector-level covariates.

- cov_effect:

  Numeric log-odds effect of a single simulated detection covariate
  (`covariate`) (retained for backward compatibility).

- n_cov_effect:

  Numeric log-linear effect of a single simulated abundance covariate
  (`habitat`) (retained for backward compatibility).

- method:

  Statistical inference paradigm (`"ML"`, `"Bayes"`, `"frequentist"`, or
  `"bayesian"`). Default is `"ML"`.

- obs_type:

  Observer field protocol (`"independent"` or `"dependent"`).

- pooling:

  Sector pooling strategy (`"hybrid"` or `"global"`).

- ci_method:

  Confidence interval method (`"lognormal"` or `"bootstrap"`).

## Value

A data frame summarizing performance metrics for each scenario in
`surv_areas`: - `surv_areas`: Number of DO sectors surveyed -
`successful_sims`: Number of simulations that returned valid estimates -
`mean_true_N`: Average simulated actual population size - `mean_est_N`:
Average estimated population size - `rel_bias`: Relative bias ((mean_est
/ mean_true) - 1) - `rmse`: Root Mean Squared Error - `coverage`:
Proportion of simulations where true pop fell within the 95% CI -
`mean_cv`: Mean Coefficient of Variation (SE / Estimate)
