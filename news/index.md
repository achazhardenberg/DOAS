# Changelog

## DOAS 1.1.0 (2026-11-05)

### Minor Release: Bayesian Diagnostics & Posterior Predictive Checks with `bayesplot`

This release expands the Bayesian capabilities of `DOAS` by integrating
seamlessly with the [`bayesplot`](https://mc-stan.org/bayesplot/)
ecosystem for graphical Posterior Predictive Checks (PPC) and MCMC
diagnostics.

#### Posterior Predictive Simulation & Model Checking

- **[`posterior_predict()`](https://achazhardenberg.github.io/DOAS/reference/posterior_predict.md)
  S3 Method**:
  - Simulates $`S \times M`$ replicated datasets ($`y_{\text{rep}}`$)
    from the posterior predictive distribution of a fitted `doas`
    Bayesian model (`method = "Bayes"`).
  - Double-observer sectors simulate from the multinomial
    capture-recapture process conditional on posterior draws of latent
    abundance and detection probabilities; single-observer sectors
    simulate from a binomial detection process.
  - Attaches observed survey counts as `attr(yrep, "y")`, making the
    output immediately compatible with all `bayesplot::ppc_*` functions.
  - Supports target predictions for total unique counts
    (`type = "count"`), joint detections (`type = "both"`), and
    observer-specific tallies (`type = "obs1"`, `type = "obs2"`).
  - Automatically handles unsurveyed sectors with
    `drop_unsurveyed = TRUE`.
- **[`pp_check()`](https://achazhardenberg.github.io/DOAS/reference/pp_check.md)
  S3 Method**:
  - Implements a direct S3 method for
    [`bayesplot::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
    on `doas` objects.
  - Supported graphical checks include:
    - `"dens_overlay"`: Density overlay comparing empirical counts
      against multiple $`y_{\text{rep}}`$ draws.
    - `"stat"`: Empirical distribution of any test statistic
      (e.g. `stat = "sd"`, `stat = "max"`, `stat = "mean"`) with the
      observed statistic marked.
    - `"stat_2d"`: Bivariate scatterplot of two test statistics
      (e.g. `stat = c("mean", "sd")`).
    - `"bars"`: Bar plot comparing observed vs. simulated count
      frequencies (ideal for inspecting zero-inflation and clumped
      distributions).
    - `"intervals"`, `"ribbon"`, `"rootogram"`: Sector-level predictive
      intervals and rootograms.

#### MCMC Diagnostic Arrays

- **[`as.array()`](https://rdrr.io/r/base/array.html) &
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) S3 Methods**:
  - Converts multi-chain MCMC draws into a 3D array of dimensions
    `[iterations, chains, parameters]` or a 2D matrix.
  - Extracts chain-specific draws for total population abundance
    (`total_abundance`), mean detection probability (`p_hat`), detection
    covariate coefficients (`p_*`), and abundance covariate coefficients
    (`n_*`).
  - Enables all `bayesplot::mcmc_*` diagnostics: `mcmc_trace()`,
    `mcmc_intervals()`, `mcmc_areas()`, and `mcmc_dens_overlay()`.

#### Package Infrastructure & Testing

- Added `bayesplot` and `ggplot2` to `Suggests:` in `DESCRIPTION`.
- Updated `_pkgdown.yml` reference index to include Bayesian diagnostic
  methods under a dedicated “Model Selection & Diagnostics” section.
- Expanded the unit test suite to 164 passing tests covering all new S3
  generics, argument checks, dependent/independent protocols, and error
  handling.
- Verified clean `R CMD check` with zero errors, zero warnings, zero
  notes.

------------------------------------------------------------------------

## DOAS 1.0.0 (2026-11-04)

### Major Release: The Complete DOAS Framework

This major release marks the culmination of 9 months of development
since the initial package release in March 2026. `DOAS` is now a mature,
fully-featured, and mathematically rigorous package for wildlife census
design, combining Double Observer (DO) mark-recapture/removal counts
with single-observer block surveys and environmental covariate modeling.

#### Core Architecture & API Refinements

- **Canonical Inference Paradigms (`method`)**:
  - **`method = "ML"` (default)**: Maximum Likelihood / analytical
    closed-form Chapman (independent observers) or Zippin removal
    (dependent observers) base estimators with multivariate Delta-method
    variance propagation across sectors. Fully backward-compatible with
    `method = "frequentist"`.
  - **`method = "Bayes"`**: Multi-chain Collapsed Poisson-Gamma-Beta
    Gibbs Sampler (Data Augmentation MCMC) sampling the full joint
    posterior of latent undetected groups $`f_0`$ and observer detection
    probabilities. Fully backward-compatible with `method = "bayesian"`.
- **Observer Protocols (`obs_type`)**:
  - `"independent"` (default): Two observers scan independently
    (Lincoln-Petersen / Chapman geometry).
  - `"dependent"`: Primary observer points out detected groups and
    secondary observer records only missed groups (2-pass removal /
    Zippin geometry).
- **Flexible Sector Pooling (`pooling`)**:
  - `"hybrid"` (default): Uses sector-specific double-observer estimates
    for DO sectors and extrapolates single-observer Block sectors using
    the modeled detection probability.
  - `"global"`: Extrapolates all surveyed sectors using the
    pooled/modeled detection probability.
- **Boundary-Safe Confidence & Credible Intervals (`ci_method`)**:
  - `"lognormal"` (default): Implements Chao (1987) log-normal bounds on
    the unobserved component $`f_0 = \hat{N} - N_{\text{obs}}`$,
    strictly guaranteeing
    $`\text{CI}_{\text{lower}} \ge N_{\text{obs}}`$. For Bayesian
    models, evaluates exact 95% posterior credible intervals directly
    from MCMC draws.
  - `"bootstrap"`: Stratified non-parametric sector resampling
    bootstrap.
- **Solitary Species & Individual Counts**: Full native support for
  solitary species (carnivores, territorial birds, solitary ungulates)
  where group size $`\equiv 1`$ and `detected_groups == total_count`.

#### Decoupled Covariate Modeling & Prediction

- **Detection Covariate Modeling (`p_formula`)**: Models primary
  observer detection probability $`\hat{p}_i`$ as a function of
  environmental or survey covariates
  (e.g. `p_formula = ~ effort + ruggedness`). Fits a quasi-binomial
  GLM (ML) or a coupled Metropolis-within-Gibbs logistic regression
  (Bayes).
- **Abundance Covariate Modeling (`n_formula`)**: Models sector-level
  abundance $`\lambda_i`$ as a function of habitat covariates
  (e.g. `n_formula = ~ elevation + forest`). Fits a quasi-Poisson
  GLM (ML) or a Poisson log-linear regression (Bayes).
- **Area & Effort Offsets
  ([`offset()`](https://rdrr.io/r/stats/offset.html))**: Supports
  `offset(log(area))` or `offset(log(transect_length))` in `n_formula`
  to account for variable sector sizes or transect lengths, yielding
  density estimates per unit area or effort.
- **S3 Prediction Method (`predict.doas`)**: Predicts expected
  individual abundance (`type = "abundance"`), detection probability
  (`type = "detection"`), or group count (`type = "groups"`) for
  surveyed sectors or new data (`newdata`).
- **Unsurveyed Sector Imputation**: Automatically imputes abundance and
  standard errors for completely unsurveyed sectors (`total_count = NA`)
  using the fitted `n_formula` habitat response model, propagating both
  process variance and parameter covariance into the total population
  estimate.

#### Multi-Chain Bayesian Gibbs Sampler & Diagnostics

- **Conjugate Beta Updates**: Employs exact conjugate Beta full
  conditionals for detection probabilities. For dependent observers, the
  full conditional
  $`p \mid x_1, x_2, f_0 \sim \text{Beta}(x_1 + x_2 + \alpha, x_2 + 2 f_0 + \beta)`$
  remains strictly defined on $`(0, 1)`$ even in sparse regimes where
  the secondary observer detects as many groups as the primary
  ($`x_1 \le x_2`$).
- **Convergence Diagnostics**: Automatically computes Gelman-Rubin
  potential scale reduction factor ($`\hat{R}`$) and Effective Sample
  Size (ESS) with lag-1 autocorrelation correction across all chains for
  total abundance, sector abundances, and detection probabilities.

#### Data Ingestion & Formatting (`doas_data`)

- **Multi-Format Ingestion**: Ingests and standardizes raw field sheets
  from:
  - `"sector"`: Wide sector-by-sector summary tables.
  - `"group"`: One row per detected group (long format).
  - `"sighting"`: One row per observer sighting record linked by
    `group_id`.
- **Group Size Reconciliation**: Supports `"max"` (default), `"mean"`,
  or `"obs1"` reconciliation rules for sighting records.
- **Validation & Integrity Checks**: Verifies tally consistency
  (`detected_groups == groups_both + groups_obs1 + groups_obs2`) and
  prevents misclassification of DO sectors with 0 secondary detections.

#### Model Selection (`QAIC` & `WAIC`)

- **Quasi-Akaike Information Criterion (`QAIC`)**: Model selection for
  frequentist quasi-binomial detection GLMs and quasi-Poisson abundance
  GLMs, using a common variance inflation factor
  $`\hat{c} = \max(1, \hat{c}_{\text{global}})`$ estimated from the most
  parameterized candidate model (Burnham & Anderson 2002).
- **Watanabe-Akaike Information Criterion (`WAIC`)**: Out-of-sample
  predictive accuracy evaluation for Bayesian MCMC models using
  log-sum-exp stabilized log pointwise predictive density and
  $`p_{\text{WAIC}}`$ variance penalties (Gelman et al. 2014).

#### Simulation Engine & Power Analysis

- **Population Simulation (`simulate_pop`)**: Generates cluster-sampled
  animal populations with environmental covariates and
  Poisson-distributed group sizes.
- **Survey Simulation (`simulate_survey`)**: Simulates field observer
  detections with detection covariates, group-size scaling, and observer
  count error.
- **Power Analysis (`doas_power`)**: Evaluates Relative Bias, Root Mean
  Squared Error (RMSE), and 95% CI Coverage across varying DO survey
  effort levels to optimize field sampling design.

#### Documentation & Infrastructure

- Full `pkgdown` website with automated GitHub Actions continuous
  deployment.
- Comprehensive tutorial vignette (`vignettes/DOAS-tutorial.Rmd`)
  demonstrating ungulate and solitary surveys, offsets, and power
  simulations.
- Official package hex logo in SVG and PNG formats under `man/figures/`.
- 133 passing `testthat` unit tests and clean `R CMD check` status
  (`Status: OK`).

------------------------------------------------------------------------

## DOAS 0.1.0 (2026-03-04)

### Initial Release

- Initial implementation of the Double Observer Adjusted Survey (DOAS)
  framework accompanying Panaccio et al. (2024).
- Core Chapman and Zippin estimators for independent and dependent
  double observers.
- Basic simulation functions for cluster-sampled populations and survey
  detections.
- Bootstrap and log-normal confidence intervals.
