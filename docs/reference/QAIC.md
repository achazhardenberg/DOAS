# Quasi-Akaike Information Criterion (QAIC)

Computes the Quasi-Akaike Information Criterion (QAIC; Burnham &
Anderson 2002) for frequentist statistical models that handle
overdispersion (such as quasi-binomial detection GLM models via
`p_formula` or quasi-poisson abundance GLM models via `n_formula` in
[`doas`](https://achazhardenberg.github.io/DOAS/reference/doas.md)). For
Bayesian MCMC models (`method = "bayesian"`), use
[`WAIC`](https://achazhardenberg.github.io/DOAS/reference/WAIC.md)
instead.

Extracts the QAIC for frequentist DOAS models that used a detection
covariate formula (`p_formula`) or an abundance covariate formula
(`n_formula`). Following Burnham & Anderson (2002), when comparing
multiple candidate models, a single common variance inflation factor
\\\hat{c} = \max(1, \hat{c}\_{\text{global}})\\ is estimated from the
most parameterized model in the candidate set.

## Usage

``` r
QAIC(object, ...)

# S3 method for class 'doas'
QAIC(object, ...)
```

## Arguments

- object:

  A fitted \`doas\` object with `method = "frequentist"`.

- ...:

  Optionally more fitted \`doas\` objects.

## Value

A numeric QAIC value, or a comparison data frame if multiple models are
passed.

A numeric QAIC value, or a data frame if multiple models are passed.
