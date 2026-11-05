# Package index

## Main Estimation Workflow

Top-level functions to structure field data and estimate abundance.

- [`doas()`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
  [`summary(`*`<doas>`*`)`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
  : Double Observer Adjusted Survey (DOAS) Abundance Estimation
- [`doas_data()`](https://achazhardenberg.github.io/DOAS/reference/doas_data.md)
  : Format and Validate Survey Data for DOAS
- [`predict(`*`<doas>`*`)`](https://achazhardenberg.github.io/DOAS/reference/predict.doas.md)
  : Predict Method for DOAS Models
- [`print(`*`<doas_data>`*`)`](https://achazhardenberg.github.io/DOAS/reference/print.doas_data.md)
  : Print Method for Validated DOAS Survey Data

## Simulation & Power Analysis

Tools to simulate cluster-sampled populations, survey protocols, and run
power analyses.

- [`simulate_pop()`](https://achazhardenberg.github.io/DOAS/reference/simulate_pop.md)
  : Simulate a Population
- [`simulate_survey()`](https://achazhardenberg.github.io/DOAS/reference/simulate_survey.md)
  : Simulate a Survey and format for DOAS
- [`doas_power()`](https://achazhardenberg.github.io/DOAS/reference/doas_power.md)
  : DOAS Power Analysis

## Model Selection & Diagnostics

Information criteria and Bayesian posterior predictive checks (PPC) via
bayesplot.

- [`QAIC()`](https://achazhardenberg.github.io/DOAS/reference/QAIC.md) :
  Quasi-Akaike Information Criterion (QAIC)
- [`WAIC()`](https://achazhardenberg.github.io/DOAS/reference/WAIC.md) :
  Widely Applicable Information Criterion (WAIC)
- [`pp_check()`](https://achazhardenberg.github.io/DOAS/reference/pp_check.md)
  : Graphical Posterior Predictive Checks for DOAS Bayesian Models
- [`posterior_predict()`](https://achazhardenberg.github.io/DOAS/reference/posterior_predict.md)
  : Draw from Posterior Predictive Distribution for DOAS Bayesian Models
- [`as.array(`*`<doas>`*`)`](https://achazhardenberg.github.io/DOAS/reference/as.array.doas.md)
  : Convert DOAS Bayesian MCMC Draws to Array
- [`as.matrix(`*`<doas>`*`)`](https://achazhardenberg.github.io/DOAS/reference/as.matrix.doas.md)
  : Convert DOAS Bayesian MCMC Draws to Matrix

## Base Estimators

Underlying mathematical estimators for mark-recapture and removal
setups.

- [`chapman_estimate()`](https://achazhardenberg.github.io/DOAS/reference/chapman_estimate.md)
  : Calculate Chapman-adjusted Lincoln-Petersen estimate
- [`chapman_variance()`](https://achazhardenberg.github.io/DOAS/reference/chapman_variance.md)
  : Calculate Variance of the Chapman estimator
- [`dependent_zippin_estimate()`](https://achazhardenberg.github.io/DOAS/reference/dependent_zippin_estimate.md)
  : Calculate Dependent Double Observer (Zippin) Estimate
- [`dependent_zippin_variance()`](https://achazhardenberg.github.io/DOAS/reference/dependent_zippin_variance.md)
  : Calculate Variance of the Dependent Double Observer (Zippin)
  estimator
- [`bayesian_conjugate_estimate()`](https://achazhardenberg.github.io/DOAS/reference/bayesian_conjugate_estimate.md)
  : Calculate Bayesian Multi-Chain MCMC Double Observer Estimate
