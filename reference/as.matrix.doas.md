# Convert DOAS Bayesian MCMC Draws to Matrix

Convert DOAS Bayesian MCMC Draws to Matrix

## Usage

``` r
# S3 method for class 'doas'
as.matrix(x, ...)
```

## Arguments

- x:

  A fitted `doas` model object with `method = "Bayes"`.

- ...:

  Additional arguments passed to
  [`as.array.doas`](https://achazhardenberg.github.io/DOAS/reference/as.array.doas.md).

## Value

A matrix with rows corresponding to pooled draws and columns to
parameters.
