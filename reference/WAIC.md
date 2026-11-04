# Widely Applicable Information Criterion (WAIC)

Computes the Watanabe-Akaike (Widely Applicable) Information Criterion
(WAIC; Watanabe 2010; Gelman, Hwang & Vehtari 2014) from the full
pointwise posterior log-likelihood draws of Bayesian MCMC
[`doas`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
models (`method = "Bayes"`). Can compare both detection covariate models
(`p_formula = ~ ...`), abundance covariate models (`n_formula = ~ ...`),
and the null model.

Evaluates the log pointwise predictive density (`lppd`) and the
effective number of parameters (`p_WAIC`, computed via the posterior
variance of the pointwise log-likelihood across MCMC draws as
recommended by Gelman et al. 2014).

## Usage

``` r
WAIC(object, ...)

# S3 method for class 'doas'
WAIC(object, ...)
```

## Arguments

- object:

  A fitted \`doas\` object with `method = "Bayes"` (or `"bayesian"`).

- ...:

  Optionally more fitted Bayesian \`doas\` objects.

## Value

If a single model is provided, a numeric scalar WAIC value (with
attributes `lppd`, `p_WAIC`, and `SE`). If multiple models are provided,
a sorted comparison data frame with columns `lppd`, `p_WAIC`, `WAIC`,
`SE`, and `delta_WAIC`.

A numeric WAIC value (for a single model) or a comparison data frame
(for multiple models).
