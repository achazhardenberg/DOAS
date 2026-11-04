# Predict Method for DOAS Models

Predicts sector-level abundance (`type = "abundance"`), group count
(`type = "groups"`), or single-observer detection probability
(`type = "detection"`) from a fitted
[`doas`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
object, either for the original surveyed sectors (`newdata = NULL`) or
for a new data frame of covariates (`newdata`).

## Usage

``` r
# S3 method for class 'doas'
predict(
  object,
  newdata = NULL,
  type = c("abundance", "detection", "groups"),
  ...
)
```

## Arguments

- object:

  A fitted object of class `"doas"`.

- newdata:

  Optional data frame containing covariates required by `p_formula`
  (when `type = "detection"`) or `n_formula` (when `type = "abundance"`
  or `"groups"`). If `NULL`, predictions are returned for the original
  sectors.

- type:

  Character string specifying the prediction target: - `"abundance"`
  (default): Expected individual abundance \\\hat{N}\_i\\. -
  `"detection"`: Single-observer detection probability \\\hat{p}\_i \in
  (0, 1)\\. - `"groups"`: Expected group abundance \\\hat{G}\_i\\.

- ...:

  Additional arguments (currently unused).

## Value

A data frame with columns `fit`, `se`, `ci_lower`, and `ci_upper` (plus
`sector_id` if available).
