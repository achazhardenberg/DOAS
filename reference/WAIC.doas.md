# WAIC Method for Bayesian DOAS Models

Evaluates the log pointwise predictive density (`lppd`) and the
effective number of parameters (`p_WAIC`, computed via the posterior
variance of the pointwise log-likelihood across MCMC draws as
recommended by Gelman et al. 2014).

## Usage

``` r
# S3 method for class 'doas'
WAIC(object, ...)
```

## Arguments

- object:

  A fitted \`doas\` object with `method = "bayesian"`.

- ...:

  Optionally more fitted Bayesian \`doas\` objects.

## Value

A numeric WAIC value (for a single model) or a comparison data frame
(for multiple models).

## Examples

``` r
if (FALSE) { # \dontrun{
mod_null   <- doas(data, method = "bayesian")
mod_effort <- doas(data, method = "bayesian", p_formula = ~ effort)
WAIC(mod_null, mod_effort)
} # }
```
