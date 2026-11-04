# Package index

## Main Estimation Workflow

Top-level functions to structure field data and estimate abundance.

- [`doas()`](https://achazhardenberg.github.io/DOAS/reference/doas.md) :
  Double Observer Adjusted Survey (DOAS) Abundance Estimation
- [`doas_data()`](https://achazhardenberg.github.io/DOAS/reference/doas_data.md)
  : Format and Validate Survey Data for DOAS
- [`predict(`*`<doas>`*`)`](https://achazhardenberg.github.io/DOAS/reference/predict.doas.md)
  : Predict Method for DOAS Models
- [`summary(`*`<doas>`*`)`](https://achazhardenberg.github.io/DOAS/reference/summary.doas.md)
  : Summary for doas objects
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

## Model Selection

Information criteria for frequentist and Bayesian models.

- [`QAIC()`](https://achazhardenberg.github.io/DOAS/reference/QAIC.md) :
  Quasi-Akaike Information Criterion (QAIC)
- [`QAIC(`*`<doas>`*`)`](https://achazhardenberg.github.io/DOAS/reference/QAIC.doas.md)
  : QAIC Method for Frequentist DOAS Models
- [`WAIC()`](https://achazhardenberg.github.io/DOAS/reference/WAIC.md) :
  Widely Applicable Information Criterion (WAIC)
- [`WAIC(`*`<doas>`*`)`](https://achazhardenberg.github.io/DOAS/reference/WAIC.doas.md)
  : WAIC Method for Bayesian DOAS Models

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
