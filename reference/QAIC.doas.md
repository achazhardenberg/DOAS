# QAIC Method for Frequentist DOAS Models

Extracts the QAIC for frequentist DOAS models that used a detection
covariate formula (`p_formula`) or an abundance covariate formula
(`n_formula`). Following Burnham & Anderson (2002), when comparing
multiple candidate models, a single common variance inflation factor
\\\hat{c} = \max(1, \hat{c}\_{\text{global}})\\ is estimated from the
most parameterized model in the candidate set.

## Usage

``` r
# S3 method for class 'doas'
QAIC(object, ...)
```

## Arguments

- object:

  A fitted \`doas\` object with `method = "frequentist"`.

- ...:

  Optionally more fitted \`doas\` objects.

## Value

A numeric QAIC value, or a data frame if multiple models are passed.

## Examples

``` r
if (FALSE) { # \dontrun{
mod1 <- doas(data, method = "frequentist", p_formula = ~ effort)
mod2 <- doas(data, method = "frequentist", p_formula = ~ effort + habitat)
QAIC(mod1, mod2)
} # }
```
