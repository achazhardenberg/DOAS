# Convert DOAS Bayesian MCMC Draws to Array

Formats posterior MCMC draws of model parameters into a 3D array of
dimensions `[iterations, chains, parameters]`, compatible with MCMC
diagnostic functions in `bayesplot` (e.g. `mcmc_trace`,
`mcmc_intervals`).

## Usage

``` r
# S3 method for class 'doas'
as.array(x, ...)
```

## Arguments

- x:

  A fitted `doas` model object with `method = "Bayes"`.

- ...:

  Additional arguments (currently unused).

## Value

A 3-dimensional array with dimensions
`[iterations, chains, parameters]`.
