# Draw from Posterior Predictive Distribution for DOAS Bayesian Models

Simulates replicated datasets from the posterior predictive distribution
of a fitted
[`doas`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
Bayesian model. For double-observer sectors, observations are simulated
from the multinomial capture-recapture process conditional on posterior
draws of latent abundance and observer detection probabilities. For
single-observer sectors, observations are simulated from a binomial
detection process.

## Usage

``` r
posterior_predict(object, ...)

# S3 method for class 'doas'
posterior_predict(
  object,
  n_draws = NULL,
  type = c("count", "detected", "both", "joint", "obs1", "obs2"),
  drop_unsurveyed = TRUE,
  ...
)
```

## Arguments

- object:

  A fitted `doas` model object with `method = "Bayes"`.

- ...:

  Additional arguments (currently unused).

- n_draws:

  Number of posterior predictive draws to simulate. If `NULL` (default),
  all available post-warmup draws across chains are used.

- type:

  Character string indicating what quantity to predict:

  - `"count"` or `"detected"` (default): Total unique individuals/groups
    detected at each sector (comparable to `object$data$total_count`).

  - `"both"` or `"joint"`: Joint detections \\y\_{12}\\ in
    double-observer sectors.

  - `"obs1"`: Observer 1 detections (\\y_1 + y\_{12}\\ in DO sectors).

  - `"obs2"`: Observer 2 detections (\\y_2 + y\_{12}\\ in DO sectors).

- drop_unsurveyed:

  Logical. If `TRUE` (default), unsurveyed sectors (which have no field
  observations) are omitted from `yrep` and `attr(yrep, "y")`.

## Value

A matrix of class `c("doas_yrep", "matrix")` of dimension
`n_draws x n_sectors`, with attributes `"y"` (observed vector),
`"sector_id"` (sector identifiers), and `"target"` (predicted quantity).

## Details

The returned matrix `yrep` has dimensions `n_draws x n_sectors` and
includes the observed data vector as `attr(yrep, "y")`, making it
directly compatible with `bayesplot::ppc_*` functions.
