# Quasi-Akaike Information Criterion (QAIC)

Computes the Quasi-Akaike Information Criterion (QAIC; Burnham &
Anderson 2002) for frequentist statistical models that handle
overdispersion (such as quasi-binomial detection GLM models via
`p_formula` or quasi-poisson abundance GLM models via `n_formula` in
[`doas`](https://achazhardenberg.github.io/DOAS/reference/doas.md)). For
Bayesian MCMC models (`method = "bayesian"`), use
[`WAIC`](https://achazhardenberg.github.io/DOAS/reference/WAIC.md)
instead.

## Usage

``` r
QAIC(object, ...)
```

## Arguments

- object:

  A fitted model object (e.g., a frequentist \`doas\` model).

- ...:

  Optionally more fitted model objects for model comparison.

## Value

A numeric QAIC value, or a comparison data frame if multiple models are
passed.
