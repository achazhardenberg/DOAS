# DOAS Tutorial: Analyzing Real Survey Data and Running Simulations

## Introduction

The **Double Observer Adjusted Survey (DOAS)** framework reduces the
effort required to estimate wild animal abundance by combining standard
block counts with a subset of Double Observer (DO) counts.

This tutorial covers the complete workflow of the `DOAS` package. You
will learn how to: 1. **Analyze real survey data** using the
[`doas()`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
function. 2. **Model covariates on detection probability (`p_formula`)
and abundance (`n_formula`)**, predict onto new data with
[`predict()`](https://rdrr.io/r/stats/predict.html), and compare models
with
[`QAIC()`](https://achazhardenberg.github.io/DOAS/reference/QAIC.md) /
[`WAIC()`](https://achazhardenberg.github.io/DOAS/reference/WAIC.md). 3.
**Simulate animal populations and surveys** using
[`simulate_pop()`](https://achazhardenberg.github.io/DOAS/reference/simulate_pop.md)
and
[`simulate_survey()`](https://achazhardenberg.github.io/DOAS/reference/simulate_survey.md).
4. **Run power analyses** using
[`doas_power()`](https://achazhardenberg.github.io/DOAS/reference/doas_power.md)
to plan your field sampling efforts. 5. **Understand the underlying base
estimators**.

------------------------------------------------------------------------

## 1. Analyzing Survey Data (`doas`)

Your data needs to be summarized at the **sector level**. The final
`data.frame` passed to the `doas` function must contain exactly one row
per survey sector.

### Required Columns:

1.  **`sector_id`**: A unique identifier for the sector.
2.  **`total_count`**: The raw total number of individuals counted in
    that sector.
3.  **`detected_groups`**: The total number of unique groups detected in
    that sector.

### Double Observer (DO) Columns (Optional):

For sectors that received a Double Observer survey, provide the
following columns (leave them as `NA` for standard Block sectors): 4.
**`groups_both`**: How many groups were seen by **both** observers. 5.
**`groups_obs1`**: How many groups were seen by **Observer 1 only**. 6.
**`groups_obs2`**: How many groups were seen by **Observer 2 only**.

> **Note:** The
> [`doas()`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
> function automatically calculates the mean group size internally
> (`total_count / detected_groups`). You do not need to provide it.

### Example Data

Let’s generate a 24-sector dataset using
[`simulate_pop()`](https://achazhardenberg.github.io/DOAS/reference/simulate_pop.md)
and
[`simulate_survey()`](https://achazhardenberg.github.io/DOAS/reference/simulate_survey.md),
where 12 sectors receive Double Observer surveys and 12 receive
single-observer block counts, with a simulated detection covariate
(`effort`) and a simulated abundance covariate (`habitat`):

``` r

library(DOAS)

set.seed(12)
sim_pop_example <- simulate_pop(
  pop_size   = 3500,
  mean_group = 8,
  tot_areas  = 24,
  n_formula  = ~ habitat,
  n_betas    = c(habitat = 0.5)
)

my_survey_data <- simulate_survey(
  sim_pop         = sim_pop_example,
  det_p1          = 0.55,
  intra_group_det = 1.0,
  det_group       = 0.95,
  obs_effect      = 0.90,
  detvar          = 0.02,
  surv_areas      = 12,
  p_formula       = ~ effort,
  p_betas         = c(effort = 1.0)
)

head(my_survey_data)
#>   sector_id survey_type total_count groups_both groups_obs1 groups_obs2
#> 1         1       Block          31          NA          NA          NA
#> 2         2       Block          55          NA          NA          NA
#> 3         3          DO          58           1           5           4
#> 4         4       Block          57          NA          NA          NA
#> 5         5          DO          20           0           0           1
#> 6         6          DO          97           5           3           5
#>   detected_groups    habitat     effort
#> 1               1 -1.4805676 -0.3560288
#> 2              10  1.5771695 -1.3522702
#> 3              10 -0.9567445 -0.5747150
#> 4               7 -0.9200052 -0.3996702
#> 5               1 -1.9976421 -0.5177994
#> 6              13 -0.2722960  0.5118662
```

> **Understanding DO vs. Block Sectors:** In this survey, 12 out of 24
> sectors were randomly selected for Double Observer surveys
> (`surv_areas = 12`): - **Block sectors** (e.g. sectors 1, 2, and 4;
> `survey_type == "Block"`): Only a **single observer** surveyed the
> sector. That observer detected the groups (`detected_groups = 1` for
> Sector 1, counting `total_count = 31` animals). The double-observer
> columns (`groups_both`, `groups_obs1`, `groups_obs2`) are `NA` because
> only one observer surveyed those sectors. - **Double Observer (DO)
> sectors** (e.g. sectors 3, 5, and 6; `survey_type == "DO"`): **Two
> observers** surveyed simultaneously, recording which groups were
> spotted by both observers vs. only Observer 1 or Observer 2. `DOAS`
> uses these DO sectors to estimate detection probability and adjust the
> single-observer block counts.

### Formatting & Inspecting Raw Field Sheets (`doas_data`)

If your field data is recorded in **long format** (for example, one row
per detected group or one row per observer sighting) or uses custom
column names, use
[`doas_data()`](https://achazhardenberg.github.io/DOAS/reference/doas_data.md)
first to aggregate and visually verify the sector-level table before
running
[`doas()`](https://achazhardenberg.github.io/DOAS/reference/doas.md):

``` r

raw_field_sheet <- data.frame(
  valley = c("Valley_A", "Valley_A", "Valley_A", "Valley_B", "Valley_B", "Valley_C"),
  herd_size = c(35, 10, 10, 20, 10, 80),
  seen_by = c("both", "obs1", "obs2", "both", "obs1", "block"),
  effort = c(10, 10, 10, 6, 6, 12)
)

inspected_data <- doas_data(
  raw_field_sheet,
  format = "group",
  sector = "valley",
  count = "herd_size",
  seen_by = "seen_by",
  do_sectors = c("Valley_A", "Valley_B"),
  covariates = "effort"
)

print(inspected_data)
#> DOAS Validated Survey Data: 3 sectors (2 DO, 1 Block) | 165 total individuals in 6 groups
#> --------------------------------------------------------------------------------
#>   sector_id survey_type total_count detected_groups groups_both groups_obs1
#> 1  Valley_A          DO          55               3           1           1
#> 2  Valley_B          DO          30               2           1           1
#> 3  Valley_C       Block          80               1          NA          NA
#>   groups_obs2 mean_group_size effort
#> 1           1        18.33333     10
#> 2           0        15.00000      6
#> 3          NA        80.00000     12
```

### Running the Estimator

[`doas()`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
separates the estimation setup into three intuitive arguments: -
**`method`**: `"ML"` (default, Maximum Likelihood / analytical Chapman
or Zippin removal) or `"Bayes"` (multi-chain Bayesian MCMC collapsed
Gibbs sampler). - **`obs_type`**: `"independent"` (default, observers
scan independently) or `"dependent"` (Observer 1 points out detected
groups and Observer 2 records only missed groups). - **`pooling`**:
`"hybrid"` (default, uses sector-specific DO estimates in DO sectors and
extrapolates Block sectors) or `"global"` (extrapolates all sectors
using the pooled primary detection probability).

``` r

# Default ML estimation with Independent observers and Hybrid pooling
results <- doas(
  my_survey_data,
  method = "ML",
  obs_type = "independent",
  pooling = "hybrid"
)

summary(results)
#> 
#>   Double Observer Adjusted Survey (DOAS) Estimation  
#> ====================================================
#> 
#> Inference Method: Maximum Likelihood / Analytical (Chapman)
#> Sector Pooling:   Hybrid
#> Total sectors: 24
#> Double Observer (DO) sectors: 12 (Independent)
#> 
#> Mean detection probability (p_hat): 0.535 (Var = 0.00207)
#> 
#> Estimated Total Abundance: 3134.2
#> Standard Error (Posterior SD): 190.7
#> 95% Log-normal CI: [2823.1 - 3580.4]
#> ----------------------------------------------------
#> Sector-level Estimates:
#>  sector_id total_count est_groups est_abundance se_abundance ci_lower ci_upper
#>          1          31       1.87         57.99        39.86    34.28   253.32
#>          2          55      18.71        102.88        23.86    74.02   175.52
#>          3          58      20.00        116.00        48.53    71.91   299.82
#>          4          57      13.09        106.62        28.95    74.17   200.35
#>          5          20       1.00         20.00         0.00    20.00    20.00
#>          6          97      15.50        115.65        18.11   100.78   189.15
#>          7         110      13.29        121.79        13.89   111.89   183.40
#>          8          75      11.22        140.29        40.86    96.16   276.45
#>          9          89      29.93        166.48        31.73   124.80   256.65
#>         10         117      18.80        157.11        33.26   126.72   282.59
#>         11          27      11.22         50.50        14.71    34.62    99.52
#>         12          21       7.48         39.28        13.81    25.90    89.23
#>         13          14       3.74         26.19        12.83    16.24    80.21
#>         14          63      16.50        115.50        46.30    74.87   295.18
#>         15         109      28.06        203.89        39.89   152.03   318.21
#>         16          65      16.50        119.17        47.77    77.25   304.55
#>         17         176      31.31        204.08        20.55   183.78   277.37
#>         18         109      33.67        203.89        37.10   154.31   307.72
#>         19          95      23.00        168.08        56.61   114.06   375.12
#>         20          34       3.00         34.00         0.00    34.00    34.00
#>         21          90      18.71        168.35        39.04   121.13   287.21
#>         22         353      34.44        392.14        30.48   363.15   503.90
#>         23         194      35.71        216.47        16.67   200.14   276.26
#>         24          47      14.96         87.91        22.49    61.95   159.00
#> ----------------------------------------------------
```

You can view the specific sector-by-sector estimates directly from the
results object:

``` r

knitr::kable(results$sector_estimates, digits = 2)
```

| sector_id | total_count | est_groups | est_abundance | se_abundance | ci_lower | ci_upper |
|----------:|------------:|-----------:|--------------:|-------------:|---------:|---------:|
|         1 |          31 |       1.87 |         57.99 |        39.86 |    34.28 |   253.32 |
|         2 |          55 |      18.71 |        102.88 |        23.86 |    74.02 |   175.52 |
|         3 |          58 |      20.00 |        116.00 |        48.53 |    71.91 |   299.82 |
|         4 |          57 |      13.09 |        106.62 |        28.95 |    74.17 |   200.35 |
|         5 |          20 |       1.00 |         20.00 |         0.00 |    20.00 |    20.00 |
|         6 |          97 |      15.50 |        115.65 |        18.11 |   100.78 |   189.15 |
|         7 |         110 |      13.29 |        121.79 |        13.89 |   111.89 |   183.40 |
|         8 |          75 |      11.22 |        140.29 |        40.86 |    96.16 |   276.45 |
|         9 |          89 |      29.93 |        166.48 |        31.73 |   124.80 |   256.65 |
|        10 |         117 |      18.80 |        157.11 |        33.26 |   126.72 |   282.59 |
|        11 |          27 |      11.22 |         50.50 |        14.71 |    34.62 |    99.52 |
|        12 |          21 |       7.48 |         39.28 |        13.81 |    25.90 |    89.23 |
|        13 |          14 |       3.74 |         26.19 |        12.83 |    16.24 |    80.21 |
|        14 |          63 |      16.50 |        115.50 |        46.30 |    74.87 |   295.18 |
|        15 |         109 |      28.06 |        203.89 |        39.89 |   152.03 |   318.21 |
|        16 |          65 |      16.50 |        119.17 |        47.77 |    77.25 |   304.55 |
|        17 |         176 |      31.31 |        204.08 |        20.55 |   183.78 |   277.37 |
|        18 |         109 |      33.67 |        203.89 |        37.10 |   154.31 |   307.72 |
|        19 |          95 |      23.00 |        168.08 |        56.61 |   114.06 |   375.12 |
|        20 |          34 |       3.00 |         34.00 |         0.00 |    34.00 |    34.00 |
|        21 |          90 |      18.71 |        168.35 |        39.04 |   121.13 |   287.21 |
|        22 |         353 |      34.44 |        392.14 |        30.48 |   363.15 |   503.90 |
|        23 |         194 |      35.71 |        216.47 |        16.67 |   200.14 |   276.26 |
|        24 |          47 |      14.96 |         87.91 |        22.49 |    61.95 |   159.00 |

### Solitary Species and Individual-Level Counts (Group Size = 1)

While DOAS was originally developed for gregarious mountain ungulates,
it applies directly to **solitary species** (e.g., roe deer, carnivores,
or territorial birds) where observers detect single animals rather than
social groups.

When each detection is an individual animal, group size is identically
$`1`$, so **the number of detected groups equals total count**
(`detected_groups = total_count`).
[`doas_data()`](https://achazhardenberg.github.io/DOAS/reference/doas_data.md)
handles this automatically: if you omit `detected_groups`, it sets
`detected_groups = total_count` and `mean_group_size = 1`:

``` r

# Survey data for a solitary species across 8 sectors (4 DO, 4 Block)
# Omitting detected_groups defaults every observation to an individual animal (group_size = 1)
solitary_data <- data.frame(
  sector_id   = paste0("Sector_", 1:8),
  total_count = c(14, 22, 19, 28, 12, 31, 16, 25),
  groups_both = c(9, 15, 13, 19, NA, NA, NA, NA),
  groups_obs1 = c(3, 4, 3, 5, NA, NA, NA, NA),
  groups_obs2 = c(2, 3, 3, 4, NA, NA, NA, NA),
  ruggedness  = c(0.8, -0.4, 0.2, -0.9, 1.1, -0.7, 0.5, -0.3)
)

validated_solitary <- doas_data(solitary_data, format = "sector")
print(validated_solitary)
#> DOAS Validated Survey Data: 8 sectors (4 DO, 4 Block) | 167 total individuals in 167 groups
#> --------------------------------------------------------------------------------
#>   sector_id survey_type total_count detected_groups groups_both groups_obs1
#> 1  Sector_1          DO          14              14           9           3
#> 2  Sector_2          DO          22              22          15           4
#> 3  Sector_3          DO          19              19          13           3
#> 4  Sector_4          DO          28              28          19           5
#> 5  Sector_5       Block          12              12          NA          NA
#> 6  Sector_6       Block          31              31          NA          NA
#> 7  Sector_7       Block          16              16          NA          NA
#> 8  Sector_8       Block          25              25          NA          NA
#>   groups_obs2 mean_group_size ruggedness
#> 1           2               1        0.8
#> 2           3               1       -0.4
#> 3           3               1        0.2
#> 4           4               1       -0.9
#> 5          NA               1        1.1
#> 6          NA               1       -0.7
#> 7          NA               1        0.5
#> 8          NA               1       -0.3

# Fit DOAS with detection covariate
fit_solitary <- doas(
  validated_solitary,
  method = "frequentist",
  p_formula = ~ruggedness
)

summary(fit_solitary)
#> 
#>   Double Observer Adjusted Survey (DOAS) Estimation  
#> ====================================================
#> 
#> Inference Method: Maximum Likelihood / Analytical (Chapman)
#> Sector Pooling:   Hybrid
#> Detection Formula (p_formula): p_hat ~ruggedness
#> Total sectors: 8
#> Double Observer (DO) sectors: 4 (Independent)
#> 
#> Detection Probability Model (p_formula GLM):
#>                Estimate Std. Error  t value     Pr(>|t|)
#> (Intercept)  1.42475070 0.03471109 41.04598 0.0005930241
#> ruggedness  -0.06188514 0.05290776 -1.16968 0.3626599027
#> 
#> Estimated Total Abundance: 190.1
#> Standard Error (Posterior SD): 5.4
#> 95% Log-normal CI: [181.6 - 203.4]
#> ----------------------------------------------------
#> Sector-level Estimates:
#>  sector_id total_count est_groups est_abundance se_abundance ci_lower ci_upper
#>   Sector_1          14      14.60         14.60         0.92    14.07    19.20
#>   Sector_2          22      22.75         22.75         1.02    22.10    27.60
#>   Sector_3          19      19.64         19.64         0.94    19.08    24.23
#>   Sector_4          28      29.00         29.00         1.20    28.16    34.34
#>   Sector_5          12      15.09         15.09         1.99    12.98    21.78
#>   Sector_6          31      38.14         38.14         2.98    34.26    46.66
#>   Sector_7          16      19.97         19.97         2.23    17.42    27.10
#>   Sector_8          25      30.90         30.90         2.71    27.51    38.90
#> ----------------------------------------------------
```

Notice that for solitary species, `est_groups` and `est_abundance` in
`fit_solitary$sector_estimates` are identical, as expected.

------------------------------------------------------------------------

## 2. Modeling Covariates on Detection (`p_formula`) and Abundance (`n_formula`)

If detection probability varies across sectors due to survey effort,
visibility, habitat, or mean group size (`mean_group_size`, which is
computed automatically inside
[`doas()`](https://achazhardenberg.github.io/DOAS/reference/doas.md)),
you can pass a one-sided R formula to `p_formula`. You can also model
sector-level abundance as a function of habitat covariates via
`n_formula` and generate predictions with
[`predict()`](https://rdrr.io/r/stats/predict.html):

- **ML (`method = "ML"`, default):** Fits a quasi-binomial logistic
  regression (`p_hat ~ covariates`, weighted by `G_est`) on the DO
  sectors (`p_formula`) and/or a quasi-Poisson log-link regression on
  detection-adjusted abundance (`n_formula`). Competing models can be
  compared using
  [`QAIC()`](https://achazhardenberg.github.io/DOAS/reference/QAIC.md).
  (Also accepts `method = "frequentist"`).
- **Bayes (`method = "Bayes"`):** Runs a multi-chain
  Metropolis-within-Gibbs regression coupled to the posterior draws of
  latent group abundance $`G_i`$, reporting posterior means, standard
  deviations, 95% credible intervals, $`\hat{R}`$, and ESS for every
  regression coefficient. (Also accepts `method = "bayesian"`).

``` r

# ML detection covariate model and abundance covariate model
mod_effort <- doas(my_survey_data, method = "ML", p_formula = ~effort, n_formula = ~habitat)
mod_group <- doas(my_survey_data, method = "ML", p_formula = ~mean_group_size, n_formula = ~habitat)

summary(mod_effort)
#> 
#>   Double Observer Adjusted Survey (DOAS) Estimation  
#> ====================================================
#> 
#> Inference Method: Maximum Likelihood / Analytical (Chapman)
#> Sector Pooling:   Hybrid
#> Detection Formula (p_formula): p_hat ~effort
#> Abundance Formula (n_formula): N_hat ~habitat
#> Total sectors: 24
#> Double Observer (DO) sectors: 12 (Independent)
#> 
#> Detection Probability Model (p_formula GLM):
#>               Estimate Std. Error   t value   Pr(>|t|)
#> (Intercept) 0.07503245  0.1610263 0.4659641 0.65122509
#> effort      0.84717674  0.3173853 2.6692377 0.02352331
#> 
#> Sector Abundance Model (n_formula GLM):
#>              Estimate Std. Error   t value     Pr(>|t|)
#> (Intercept) 4.8795479 0.06411012 76.111978 3.829349e-28
#> habitat     0.4586001 0.06164967  7.438808 1.931972e-07
#> 
#> Estimated Total Abundance: 3244.5
#> Standard Error (Posterior SD): 235.6
#> 95% Log-normal CI: [2868.1 - 3806.8]
#> ----------------------------------------------------
#> Sector-level Estimates:
#>  sector_id total_count est_groups est_abundance se_abundance ci_lower ci_upper
#>          1          31       2.25         69.88        52.68    36.26   318.61
#>          2          55      39.17        215.44        94.08   110.28   520.66
#>          3          58      20.00        116.00        48.53    71.91   299.82
#>          4          57      16.11        131.19        40.26    84.41   257.81
#>          5          20       1.00         20.00         0.00    20.00    20.00
#>          6          97      15.50        115.65        18.11   100.78   189.15
#>          7         110      13.29        121.79        13.89   111.89   183.40
#>          8          75      14.78        184.79        63.59   113.26   390.04
#>          9          89      31.35        174.37        33.47   129.68   268.13
#>         10         117      18.80        157.11        33.26   126.72   282.59
#>         11          27      18.74         84.31        34.72    46.16   198.45
#>         12          21      11.21         58.86        26.13    32.14   149.67
#>         13          14       8.89         62.20        46.02    23.95   247.49
#>         14          63      16.50        115.50        46.30    74.87   295.18
#>         15         109      18.79        136.52        21.23   116.21   214.02
#>         16          65      16.50        119.17        47.77    77.25   304.55
#>         17         176      31.31        204.08        20.55   183.78   277.37
#>         18         109      23.74        143.74        22.42   119.93   219.43
#>         19          95      23.00        168.08        56.61   114.06   375.12
#>         20          34       3.00         34.00         0.00    34.00    34.00
#>         21          90      14.22        127.95        25.48   101.48   215.44
#>         22         353      34.44        392.14        30.48   363.15   503.90
#>         23         194      35.71        216.47        16.67   200.14   276.26
#>         24          47      12.81         75.24        17.54    56.22   133.50
#> ----------------------------------------------------

# Predict detection probability or expected abundance for new sectors
predict(mod_effort, newdata = data.frame(effort = c(-1, 0, 1), habitat = c(-1, 0, 1)), type = "detection")
#>         fit         se  ci_lower  ci_upper
#> 1 0.3160154 0.07680297 0.1871589 0.4810813
#> 2 0.5187493 0.04019996 0.4401437 0.5964377
#> 3 0.7154920 0.07256428 0.5556502 0.8349180
predict(mod_effort, newdata = data.frame(effort = c(-1, 0, 1), habitat = c(-1, 0, 1)), type = "abundance")
#>         fit        se  ci_lower ci_upper
#> 1  83.17508  8.337939  68.33802 101.2335
#> 2 131.57116  8.435043 116.03502 149.1875
#> 3 208.12689 15.812700 179.33122 241.5463

# Compare competing detection covariate models using Quasi-AIC (QAIC)
QAIC(mod_effort, mod_group)
#>            df c_hat     QAIC delta_QAIC
#> mod_effort  2 1.413 40.34824   0.000000
#> mod_group   2 1.413 45.82942   5.481188
```

You can also fit covariate models in the multi-chain Bayesian framework
(`method = "Bayes"`) and compare them using the Widely Applicable
Information Criterion
([`WAIC()`](https://achazhardenberg.github.io/DOAS/reference/WAIC.md)):

``` r

set.seed(42)
mod_bayes_null <- doas(
  my_survey_data,
  method = "Bayes",
  obs_type = "independent",
  chains = 4,
  iter = 1000
)

mod_bayes_cov <- doas(
  my_survey_data,
  method = "Bayes",
  obs_type = "independent",
  p_formula = ~effort,
  n_formula = ~habitat,
  chains = 4,
  iter = 1000
)

summary(mod_bayes_cov)
#> 
#>   Double Observer Adjusted Survey (DOAS) Estimation  
#> ====================================================
#> 
#> Inference Method: Bayesian MCMC Independent
#> Sector Pooling:   Hybrid
#> Detection Formula (p_formula): p_hat ~effort
#> Abundance Formula (n_formula): N_hat ~habitat
#> Total sectors: 24
#> Double Observer (DO) sectors: 12 (Independent)
#> MCMC Sampler: 4 chains, 1000 total iterations (250 burn-in), 3000 post-warmup draws
#> Convergence: R-hat = 1.001 | Effective Sample Size (ESS) = 905
#> 
#> Bayesian Logistic Regression Posterior (Detection p_formula):
#>              Mean    SD CI_2.5 CI_97.5  Rhat ESS
#> (Intercept) 0.141 0.173 -0.184   0.479 1.006 297
#> effort      0.648 0.320  0.048   1.259 1.007 300
#> 
#> Bayesian Log-Linear Regression Posterior (Abundance n_formula):
#>              Mean    SD CI_2.5 CI_97.5  Rhat ESS
#> (Intercept) 4.841 0.076  4.702   4.996 1.006 226
#> habitat     0.483 0.068  0.350   0.617 1.009 271
#> 
#> Estimated Total Abundance: 3146.3
#> Standard Error (Posterior SD): 223.0
#> 95% Bayesian Credible Interval: [2783.1 - 3650.4]
#> ----------------------------------------------------
#> Sector-level Estimates:
#>  sector_id total_count est_groups est_abundance se_abundance ci_lower ci_upper
#>          1          31       2.12         65.66        47.38    31.28   186.28
#>          2          55      33.33        183.29        83.94    77.40   396.40
#>          3          58      15.11         87.65        24.87    58.00   150.80
#>          4          57      15.07        122.69        37.43    65.13   211.70
#>          5          20       1.31         26.19        14.15    20.00    60.00
#>          6          97      16.23        121.09        21.27    97.00   171.62
#>          7         110      14.28        130.90        19.60   110.00   183.33
#>          8          75      13.62        170.23        57.92    87.06   312.06
#>          9          89      30.47        169.47        32.87   116.61   244.55
#>         10         117      18.81        157.18        32.65   117.00   234.21
#>         11          27      16.49         74.19        31.13    35.56   152.56
#>         12          21      10.06         52.82        24.21    21.00   115.24
#>         13          14       7.45         52.16        39.47    14.75   161.75
#>         14          63      13.40         93.83        27.70    63.00   161.00
#>         15         109      20.52        149.08        33.29   109.00   231.80
#>         16          65      13.43         97.01        27.88    65.00   166.11
#>         17         176      32.04        208.85        22.90   176.00   267.26
#>         18         109      25.60        155.00        29.45   114.95   230.00
#>         19          95      19.57        143.03        35.64    95.00   233.85
#>         20          34       4.12         46.73        18.92    34.00    90.67
#>         21          90      15.03        135.31        30.38    90.00   206.73
#>         22         353      35.52        404.42        35.33   353.00   489.65
#>         23         194      36.69        222.42        19.25   194.00   266.75
#>         24          47      13.12         77.08        18.34    47.44   117.94
#> ----------------------------------------------------

# Compare Bayesian models using WAIC
WAIC(mod_bayes_null, mod_bayes_cov)
#>                     lppd   p_WAIC     WAIC       SE delta_WAIC
#> mod_bayes_cov  -23.40172 10.24014 67.28373 5.127923   0.000000
#> mod_bayes_null -24.32518 13.27370 75.19777 8.252770   7.914038
```

### Accounting for Sampling Effort, Transect Length, and Area Offsets (`offset()`)

In field surveys with transects or management sectors of unequal size
(e.g. transects of different lengths or sectors with different surface
areas), longer transects or larger sectors naturally intercept more
animals simply because **more area was sampled**.

In `DOAS`, you can model this directly by including an
[`offset()`](https://rdrr.io/r/stats/offset.html) in `n_formula`:
``` math
\log(\lambda_j) = \log(\text{effort}_j) + \beta_0 + \beta_1 \cdot \text{habitat}_j
```

The [`offset()`](https://rdrr.io/r/stats/offset.html) fixes the
coefficient to exactly $`1.0`$ on the log scale, meaning that expected
abundance scales proportionally with sampling effort, and the remaining
regression coefficients represent **density per unit effort or area**.

> **Note on Detection vs. Abundance:** - If transect length or sector
> size alters the **number of animals present in the sampled unit**, put
> `offset(log(transect_length))` in `n_formula`. - If transect length or
> survey duration alters **how well observers spot animals**
> (e.g. observer fatigue on long walks or search speed), include it as
> an unconstrained covariate in `p_formula = ~ fatigue`. -
> [`offset()`](https://rdrr.io/r/stats/offset.html) is fully supported
> in **both frequentist and Bayesian models** (as well as in
> [`predict.doas()`](https://achazhardenberg.github.io/DOAS/reference/predict.doas.md)!).

``` r

# Survey across 8 transects of variable length with categorical habitat
transect_data <- data.frame(
  sector_id       = paste0("Transect_", 1:8),
  total_count     = c(48, 22, 65, 80, 18, 30, 42, 55),
  groups_both     = c(5, 3, 8, 7, NA, NA, NA, NA),
  groups_obs1     = c(2, 1, 2, 2, NA, NA, NA, NA),
  groups_obs2     = c(1, 1, 1, 1, NA, NA, NA, NA),
  detected_groups = c(8, 5, 11, 10, 3, 4, 6, 8),
  transect_length = c(5.0, 2.5, 7.2, 6.8, 2.0, 3.5, 4.5, 5.5), # km
  habitat         = factor(rep(c("Forest", "Open"), 4)),
  fatigue         = c(0.8, -0.4, 1.2, 0.9, -0.8, -0.2, 0.3, 0.7)
)

d_tr <- doas_data(transect_data, format = "sector")

# Frequentist DOAS with detection covariate and abundance offset
fit_offset_freq <- doas(
  d_tr,
  method    = "frequentist",
  p_formula = ~ fatigue,
  n_formula = ~ offset(log(transect_length)) + habitat
)

# Bayesian DOAS with detection covariate and abundance offset
fit_offset_bayes <- doas(
  d_tr,
  method    = "bayesian",
  p_formula = ~ fatigue,
  n_formula = ~ offset(log(transect_length)) + habitat,
  chains    = 2,
  iter      = 500,
  burnin    = 100
)

summary(fit_offset_freq)
#> 
#>   Double Observer Adjusted Survey (DOAS) Estimation  
#> ====================================================
#> 
#> Inference Method: Maximum Likelihood / Analytical (Chapman)
#> Sector Pooling:   Hybrid
#> Detection Formula (p_formula): p_hat ~fatigue
#> Abundance Formula (n_formula): N_hat ~offset(log(transect_length)) + habitat
#> Total sectors: 8
#> Double Observer (DO) sectors: 4 (Independent)
#> 
#> Detection Probability Model (p_formula GLM):
#>              Estimate Std. Error  t value   Pr(>|t|)
#> (Intercept) 1.2318894  0.1524096 8.082757 0.01496397
#> fatigue     0.3166266  0.1668506 1.897666 0.19817114
#> 
#> Sector Abundance Model (n_formula GLM):
#>              Estimate Std. Error   t value     Pr(>|t|)
#> (Intercept) 2.3404735 0.06072699 38.540913 2.037867e-08
#> habitatOpen 0.1105331 0.08403455  1.315329 2.364266e-01
#> 
#> Estimated Total Abundance: 406.5
#> Standard Error (Posterior SD): 20.9
#> 95% Log-normal CI: [380.1 - 467.8]
#> ----------------------------------------------------
#> Sector-level Estimates:
#>   sector_id total_count est_groups est_abundance se_abundance ci_lower ci_upper
#>  Transect_1          48       8.33         50.00         4.00    48.17    72.04
#>  Transect_2          22       5.25         23.10         2.46    22.08    37.16
#>  Transect_3          65      11.22         66.31         3.08    65.09    84.18
#>  Transect_4          80      10.25         82.00         4.47    80.15   107.57
#>  Transect_5          18       4.13         24.77         7.69    19.13    58.35
#>  Transect_6          30       5.24         39.32         9.72    31.74    80.08
#>  Transect_7          42       7.59         53.14        10.02    44.46    92.37
#>  Transect_8          55       9.87         67.86        10.51    58.16   107.25
#> ----------------------------------------------------
summary(fit_offset_bayes)
#> 
#>   Double Observer Adjusted Survey (DOAS) Estimation  
#> ====================================================
#> 
#> Inference Method: Bayesian MCMC Independent
#> Sector Pooling:   Hybrid
#> Detection Formula (p_formula): p_hat ~fatigue
#> Abundance Formula (n_formula): N_hat ~offset(log(transect_length)) + habitat
#> Total sectors: 8
#> Double Observer (DO) sectors: 4 (Independent)
#> MCMC Sampler: 2 chains, 500 total iterations (100 burn-in), 800 post-warmup draws
#> Convergence: R-hat = 1.030 | Effective Sample Size (ESS) = 214
#> 
#> Bayesian Logistic Regression Posterior (Detection p_formula):
#>              Mean    SD CI_2.5 CI_97.5  Rhat ESS
#> (Intercept) 1.550 0.790  0.169   3.277 1.014  13
#> fatigue     0.048 0.847 -1.553   1.536 1.007  13
#> 
#> Bayesian Log-Linear Regression Posterior (Abundance n_formula):
#>              Mean    SD CI_2.5 CI_97.5  Rhat ESS
#> (Intercept) 2.323 0.128  2.078   2.599 1.009  50
#> habitatOpen 0.118 0.148 -0.174   0.395 1.000  81
#> 
#> Estimated Total Abundance: 409.0
#> Standard Error (Posterior SD): 37.3
#> 95% Bayesian Credible Interval: [366.0 - 515.6]
#> ----------------------------------------------------
#> Sector-level Estimates:
#>   sector_id total_count est_groups est_abundance se_abundance ci_lower ci_upper
#>  Transect_1          48       8.32         49.91         4.30       48    60.00
#>  Transect_2          22       5.21         22.92         2.29       22    30.80
#>  Transect_3          65      11.29         66.71         3.67       65    76.82
#>  Transect_4          80      10.32         82.57         5.24       80    96.00
#>  Transect_5          18       4.51         27.08        17.88       18    83.83
#>  Transect_6          30       5.27         39.49        13.40       30    75.00
#>  Transect_7          42       7.49         52.44        11.48       42    76.77
#>  Transect_8          55       9.82         67.54        12.40       55   102.51
#> ----------------------------------------------------

# Predict expected abundance for a new 10 km transect in Forest vs Open habitat
new_transects <- data.frame(
  transect_length = c(10, 10),
  habitat = factor(c("Forest", "Open"), levels = c("Forest", "Open")),
  fatigue = c(0, 0)
)
predict(fit_offset_freq, newdata = new_transects, type = "abundance")
#>        fit       se  ci_lower ci_upper
#> 1 103.8615 6.307198  92.20678 116.9894
#> 2 116.0002 6.738041 103.51766 129.9879
predict(fit_offset_bayes, newdata = new_transects, type = "abundance")
#>        fit       se ci_lower ci_upper
#> 1 102.8859 13.40665 79.89243 134.5608
#> 2 115.6821 14.51252 91.55605 156.0262
```

### Imputing Completely Unsurveyed Sectors

Occasionally, severe weather, inaccessible terrain, or logistical
constraints prevent some sectors or transects from being surveyed at all
(`total_count = NA`).

If you provide an abundance model via `n_formula`,
[`doas()`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
automatically imputes abundance and standard errors for unsurveyed
sectors using the fitted habitat response curve and propagates this
uncertainty into the total population estimate:

``` r

# Transect 8 could not be surveyed due to weather
unsurveyed_sheet <- transect_data
unsurveyed_sheet$total_count[8] <- NA
unsurveyed_sheet$detected_groups[8] <- NA

d_unsurveyed <- doas_data(unsurveyed_sheet, format = "sector")
print(d_unsurveyed)
#> DOAS Validated Survey Data: 8 sectors (4 DO, 4 Block) | 305 total individuals in 47 groups
#> --------------------------------------------------------------------------------
#>    sector_id survey_type total_count detected_groups groups_both groups_obs1
#> 1 Transect_1          DO          48               8           5           2
#> 2 Transect_2          DO          22               5           3           1
#> 3 Transect_3          DO          65              11           8           2
#> 4 Transect_4          DO          80              10           7           2
#> 5 Transect_5       Block          18               3          NA          NA
#> 6 Transect_6       Block          30               4          NA          NA
#> 7 Transect_7       Block          42               6          NA          NA
#> 8 Transect_8  Unsurveyed          NA              NA          NA          NA
#>   groups_obs2 mean_group_size transect_length habitat fatigue
#> 1           1        6.000000             5.0  Forest     0.8
#> 2           1        4.400000             2.5    Open    -0.4
#> 3           1        5.909091             7.2  Forest     1.2
#> 4           1        8.000000             6.8    Open     0.9
#> 5          NA        6.000000             2.0  Forest    -0.8
#> 6          NA        7.500000             3.5    Open    -0.2
#> 7          NA        7.000000             4.5  Forest     0.3
#> 8          NA              NA             5.5    Open     0.7

# Fit DOAS: n_formula imputes Transect 8 based on habitat and transect length
fit_imputed <- doas(
  d_unsurveyed,
  n_formula = ~ offset(log(transect_length)) + habitat
)

# Check sector-level estimates: Transect 8 is imputed with appropriate SE and 95% CI
fit_imputed$sector_estimates[8, ]
#>    sector_id total_count est_groups est_abundance se_abundance ci_lower
#> 8 Transect_8          NA   9.410191      61.06613     7.318747 48.32212
#>   ci_upper
#> 8 77.17113
```

------------------------------------------------------------------------

## 3. Simulating Populations and Surveys

To test the method or plan study designs, the DOAS package includes
functions to simulate realistic cluster-sampled populations.

### Simulating a Population (`simulate_pop`)

[`simulate_pop()`](https://achazhardenberg.github.io/DOAS/reference/simulate_pop.md)
distributes a simulated population into logical sectors and groups. You
can optionally model true sector abundance as a function of
environmental covariates via `n_formula` and a named vector of
log-linear slopes `n_betas`:

``` r

set.seed(42)
# Simulate ~1000 individuals across 15 sectors, with abundance driven by elevation
simulated_population <- simulate_pop(
  pop_size   = 1000,
  mean_group = 10,
  tot_areas  = 15,
  n_formula  = ~ elevation,
  n_betas    = c(elevation = 0.5)
)

# Total simulated individuals and generated sector covariates
simulated_population$tot_pop
#> [1] 1020
head(simulated_population$covariates)
#>   sector_id  elevation
#> 1         1  1.3709584
#> 2         2 -0.5646982
#> 3         3  0.3631284
#> 4         4  0.6328626
#> 5         5  0.4042683
#> 6         6 -0.1061245
```

### Simulating a DOAS Survey (`simulate_survey`)

[`simulate_survey()`](https://achazhardenberg.github.io/DOAS/reference/simulate_survey.md)
takes the simulated population and applies detection probabilities to
mimic field data collection. You can also specify detection covariates
via `p_formula` and `p_betas` (which can reuse covariates from the
population like `elevation` and/or introduce new survey-specific
covariates like `effort`):

``` r

# Simulate field observers surveying the population
survey_data <- simulate_survey(
  sim_pop         = simulated_population,
  det_p1          = 0.6,  # Baseline detection probability (~60%)
  intra_group_det = 1.0,  # Observers count individuals within groups accurately
  det_group       = 0.9,  # Group size scaling factor
  obs_effect      = 0.75, # Observer 2 relative efficiency
  detvar          = 0.05, # Inter-sector detection variance
  surv_areas      = 5,    # 5 DO sectors, 10 Block sectors
  p_formula       = ~ elevation + effort,
  p_betas         = c(elevation = -0.3, effort = 0.6)
)

head(survey_data)
#>   sector_id survey_type total_count groups_both groups_obs1 groups_obs2
#> 1         1          DO          76           1           2           3
#> 2         2          DO          53           1           0           1
#> 3         3       Block          31          NA          NA          NA
#> 4         4       Block          43          NA          NA          NA
#> 5         5          DO          38           2           1           0
#> 6         6       Block          29          NA          NA          NA
#>   detected_groups  elevation     effort
#> 1               6  1.3709584 -0.4650218
#> 2               2 -0.5646982  0.5734078
#> 3               2  0.3631284  0.2549361
#> 4               1  0.6328626  2.6721442
#> 5               3  0.4042683 -0.1344950
#> 6               3 -0.1061245  0.1115741
```

This output is automatically formatted and ready to be plugged into
[`doas()`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
with matching covariate formulas:

``` r

# Fit Frequentist DOAS with detection and abundance formulas
sim_results <- doas(
  survey_data,
  method    = "frequentist",
  p_formula = ~ elevation + effort,
  n_formula = ~ elevation
)
summary(sim_results)
#> 
#>   Double Observer Adjusted Survey (DOAS) Estimation  
#> ====================================================
#> 
#> Inference Method: Maximum Likelihood / Analytical (Chapman)
#> Sector Pooling:   Hybrid
#> Detection Formula (p_formula): p_hat ~elevation + effort
#> Abundance Formula (n_formula): N_hat ~elevation
#> Total sectors: 15
#> Double Observer (DO) sectors: 5 (Independent)
#> 
#> Detection Probability Model (p_formula GLM):
#>                 Estimate Std. Error     t value  Pr(>|t|)
#> (Intercept)  0.558145484  0.7249189  0.76994193 0.5218407
#> elevation   -0.652253943  0.5343064 -1.22074883 0.3465695
#> effort      -0.008493942  0.4775497 -0.01778651 0.9874240
#> 
#> Sector Abundance Model (n_formula GLM):
#>              Estimate Std. Error   t value     Pr(>|t|)
#> (Intercept) 3.7108806 0.14161798 26.203457 1.227930e-12
#> elevation   0.6123455 0.09752624  6.278777 2.838688e-05
#> 
#> Estimated Total Abundance: 993.9
#> Standard Error (Posterior SD): 158.1
#> 95% Log-normal CI: [784.1 - 1438.9]
#> ----------------------------------------------------
#> Sector-level Estimates:
#>  sector_id total_count est_groups est_abundance se_abundance ci_lower ci_upper
#>          1          76       9.00        114.00        40.06    82.98   282.82
#>          2          53       2.00         53.00         0.00    53.00    53.00
#>          3          31       3.45         53.53        27.88    34.41   179.95
#>          4          43       1.88         81.04        73.79    46.29   483.25
#>          5          38       3.00         38.00         0.00    38.00    38.00
#>          6          29       4.60         44.50        19.42    31.31   133.04
#>          7          34       7.56         85.70        59.62    42.53   347.56
#>          8           0       0.00          0.00         0.00     0.00     0.00
#>          9          37      12.52        115.78        76.11    53.04   423.85
#>         10          34       6.20         52.70        21.23    37.14   145.29
#>         11          56      14.00         98.00        38.34    65.13   249.18
#>         12         116      13.00        188.50        70.04   130.76   472.02
#>         13          23       1.23         28.31        14.24    23.31   114.09
#>         14           8       1.48         11.81         7.42     8.33    52.30
#>         15          19       3.05         28.99        14.48    20.24    99.37
#> ----------------------------------------------------
```

### Advanced: The Multi-Chain Bayesian Approach

When observer detections are sparse—or in dependent-observer surveys
where the secondary observer occasionally spots as many groups as the
primary observer—classic frequentist estimators can suffer from
small-sample instability. Passing `method = "bayesian"` runs a
multi-chain Collapsed Poisson-Gamma-Beta Gibbs Sampler (supporting both
`obs_type = "independent"` and `obs_type = "dependent"`) with
Gelman-Rubin $`\hat{R}`$ and Effective Sample Size (ESS) diagnostics.

``` r

bayes_results <- doas(
  survey_data,
  method = "bayesian",
  obs_type = "independent",
  pooling = "global"
)
summary(bayes_results)
#> 
#>   Double Observer Adjusted Survey (DOAS) Estimation  
#> ====================================================
#> 
#> Inference Method: Bayesian MCMC Independent
#> Sector Pooling:   Global
#> Total sectors: 15
#> Double Observer (DO) sectors: 5 (Independent)
#> MCMC Sampler: 4 chains, 2000 total iterations (500 burn-in), 6000 post-warmup draws
#> Convergence: R-hat = 1.003 | Effective Sample Size (ESS) = 2049
#> 
#> Mean detection probability (p_hat): 0.420 (Var = 0.00863)
#> 
#> Estimated Total Abundance: 1203.2
#> Standard Error (Posterior SD): 342.8
#> 95% Bayesian Credible Interval: [762.9 - 2077.8]
#> ----------------------------------------------------
#> Sector-level Estimates:
#>  sector_id total_count est_groups est_abundance se_abundance ci_lower ci_upper
#>          1          76       8.80        111.49        54.56    76.00   261.05
#>          2          53       3.77         99.97        73.91    53.00   306.80
#>          3          31       5.03         77.96        49.59    31.00   216.39
#>          4          43       2.51        108.14        92.36    43.09   344.09
#>          5          38       6.29         79.64        41.70    38.00   188.51
#>          6          29       7.54         72.93        40.30    29.17   174.17
#>          7          34       7.54         85.51        46.23    34.00   203.44
#>          8           0       0.00          0.00         0.00     0.00     0.00
#>          9          37      10.06         93.05        45.14    37.19   203.69
#>         10          34      10.06         85.51        42.22    34.00   186.78
#>         11          56      11.32         79.22        32.95    56.00   165.97
#>         12         116      11.32        164.10        66.82   116.00   342.81
#>         13          23       2.51         57.84        48.56    23.41   184.41
#>         14           8       2.51         20.12        17.30     8.00    63.89
#>         15          19       5.03         47.78        30.88    19.00   123.21
#> ----------------------------------------------------
```

------------------------------------------------------------------------

## 4. Power Analysis (`doas_power`)

If you are planning a field season, you can use
[`doas_power()`](https://achazhardenberg.github.io/DOAS/reference/doas_power.md)
to understand how many DO surveys you need to achieve reliable
confidence intervals and minimize bias.

[`doas_power()`](https://achazhardenberg.github.io/DOAS/reference/doas_power.md)
runs multiple simulated surveys under the hood and reports metrics like
Relative Bias, Root Mean Squared Error (RMSE), and 95% CI Coverage.

``` r

# Run a quick power analysis testing 2 vs 5 DO sectors
# Note: For real power analysis, increase n_sims (e.g., 500)
power_results <- doas_power(
  n_sims = 10,
  surv_areas = c(2, 5), # Test scenarios: 2 DO sectors vs 5 DO sectors
  pop_size = 1000,
  mean_group = 10,
  tot_areas = 15,
  intra_group_det = 1.0,
  det_group = 0.9, # Low group-size heterogeneity
  detvar = 0.05, # Low inter-sector detection variance
  method = "frequentist",
  obs_type = "independent",
  pooling = "global"
)

knitr::kable(power_results, digits = 3)
```

|  | surv_areas | successful_sims | mean_true_N | mean_est_N | rel_bias | rmse | coverage | mean_cv |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 10 | 1003.2 | 1072.463 | 0.069 | 381.001 | 1 | 0.411 |
| 5 | 5 | 10 | 1002.0 | 1007.459 | 0.005 | 216.800 | 1 | 0.247 |

Higher effort (more `surv_areas`) generally yields improved CI coverage
and lower Coefficient of Variation (CV).

------------------------------------------------------------------------

## 5. Understanding Base Estimators

The [`doas()`](https://achazhardenberg.github.io/DOAS/reference/doas.md)
function leverages several group-abundance estimators under the hood.
For advanced applications, these are exported for you to use directly:

### Chapman Estimator (Frequentist Independent Observers)

A bias-adjusted Lincoln-Petersen estimator.

``` r

chapman_estimate(both = 10, obs1_only = 2, obs2_only = 3)
#> [1] 15.54545
chapman_variance(both = 10, obs1_only = 2, obs2_only = 3)
#> [1] 0.7520661
```

### Zippin Estimator (Frequentist Dependent Observers)

A removal estimator for primary/secondary observer setups.

``` r

dependent_zippin_estimate(both = 10, obs1_only = 2, obs2_only = 3)
#> [1] 16
dependent_zippin_variance(both = 10, obs1_only = 2, obs2_only = 3)
#> [1] 2.962963
```

### Bayesian Multi-Chain Gibbs Estimator (Independent & Dependent Observers)

Uses a multi-chain Collapsed Poisson-Gamma-Beta Gibbs Sampler to sample
latent undetected groups and observer detection probabilities.

``` r

set.seed(123)
bayesian_res <- bayesian_conjugate_estimate(
  both = 10,
  obs1_only = 2,
  obs2_only = 3,
  obs_type = "independent",
  iter = 1000
)
c(Mean = bayesian_res$mean, Variance = bayesian_res$var, Rhat = bayesian_res$rhat)
#>      Mean  Variance      Rhat 
#> 15.693667  1.138873  1.000000
```
