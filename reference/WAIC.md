# Widely Applicable Information Criterion (WAIC)

Computes the Watanabe-Akaike (Widely Applicable) Information Criterion
(WAIC; Watanabe 2010; Gelman, Hwang & Vehtari 2014) from the full
pointwise posterior log-likelihood draws of Bayesian MCMC
[`doas`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
models (`method = "bayesian"`). Can compare both detection covariate
models (`p_formula = ~ ...`), abundance covariate models
(`n_formula = ~ ...`), and the null model.

## Usage

``` r
WAIC(object, ...)
```

## Arguments

- object:

  A fitted Bayesian \`doas\` model object (`method = "bayesian"`).

- ...:

  Optionally more fitted Bayesian \`doas\` model objects for model
  comparison.

## Value

If a single model is provided, a numeric scalar WAIC value (with
attributes `lppd`, `p_WAIC`, and `SE`). If multiple models are provided,
a sorted comparison data frame with columns `lppd`, `p_WAIC`, `WAIC`,
`SE`, and `delta_WAIC`.
