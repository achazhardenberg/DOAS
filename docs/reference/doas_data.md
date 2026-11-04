# Format and Validate Survey Data for DOAS

Transforms raw field data from sector-wide, group-level long, or
observer-sighting long formats into a validated, sector-level
\`doas_data\` data frame ready for inspection prior to running
[`doas`](https://achazhardenberg.github.io/DOAS/reference/doas.md).

## Usage

``` r
doas_data(
  data,
  format = c("sector", "group", "sighting"),
  sector = "sector_id",
  count = NULL,
  groups = "detected_groups",
  both = "groups_both",
  obs1 = NULL,
  obs2 = NULL,
  seen_by = NULL,
  observer = "observer",
  group_id = "group_id",
  primary_obs = NULL,
  count_rule = c("max", "mean", "obs1"),
  do_sectors = NULL,
  all_sectors = NULL,
  covariates = NULL
)
```

## Arguments

- data:

  A data frame containing survey records.

- format:

  Character string specifying the input data structure: - `"sector"`
  (default): One row per survey sector (wide format). - `"group"`: One
  row per detected animal group (long format). - `"sighting"`: One row
  per observer sighting record, linked by `group_id` (long format).

- sector:

  Character name of the sector identifier column (default
  `"sector_id"`).

- count:

  Character name of the individual count column (default `"total_count"`
  for `format = "sector"`, or `"group_size"` for `format = "group"` and
  `"sighting"`).

- groups:

  Character name of the total detected groups column for
  `format = "sector"` (default `"detected_groups"`). Optional for DO
  sectors if `both`, `obs1`, and `obs2` are supplied.

- both:

  Character name of the column for groups seen by both observers in
  `format = "sector"` (default `"groups_both"`). Optional if `groups`,
  `obs1`, and `obs2` are supplied.

- obs1:

  Character name of the Observer 1 column (default `"groups_obs1"` for
  `format = "sector"`, or binary indicator `"obs1"` for
  `format = "group"`).

- obs2:

  Character name of the Observer 2 column (default `"groups_obs2"` for
  `format = "sector"`, or binary indicator `"obs2"` for
  `format = "group"`).

- seen_by:

  Optional character name of a single categorical column in
  `format = "group"` indicating which observer(s) detected the group
  (values: `"both"`, `"obs1"`, `"obs2"`, `"block"`). If provided, used
  instead of binary `obs1`/`obs2` columns.

- observer:

  Character name of the observer ID column when `format = "sighting"`
  (default `"observer"`).

- group_id:

  Character name of the unique group ID column when
  `format = "sighting"` (default `"group_id"`).

- primary_obs:

  Value in the `observer` column identifying the Primary Observer
  (Observer 1) when `format = "sighting"` (default `1`, also matches
  `"1"`, `"obs1"`, `"Obs1"` if `NULL`).

- count_rule:

  Character rule for reconciling group size when both observers record a
  sighting of the same `group_id` in `format = "sighting"`: `"max"`
  (default), `"mean"`, or `"obs1"`.

- do_sectors:

  Optional vector of sector IDs that received Double Observer (DO)
  surveys. Strongly recommended for `format = "group"` and
  `format = "sighting"` so DO sectors where Observer 2 happened to
  detect 0 groups are not misclassified as Block sectors.

- all_sectors:

  Optional vector of all surveyed sector IDs (including empty sectors
  with 0 detections).

- covariates:

  Optional character vector of sector-level covariate column names to
  retain.

## Value

A validated sector-level data frame of class
`c("doas_data", "data.frame")` with columns `sector_id`, `survey_type`,
`total_count`, `detected_groups`, `groups_both`, `groups_obs1`,
`groups_obs2`, `mean_group_size`, and any requested `covariates`.

## Examples

``` r
# 1. Group-level long format (1 row per detected group)
raw_groups <- data.frame(
  valley = c("A", "A", "A", "B", "B", "C"),
  size   = c(12, 5, 8, 10, 4, 15),
  who    = c("both", "obs1", "obs2", "both", "obs1", "block"),
  effort = c(8, 8, 8, 6, 6, 5)
)
df_sec <- doas_data(
  raw_groups,
  format = "group",
  sector = "valley",
  count = "size",
  seen_by = "who",
  do_sectors = c("A", "B"),
  covariates = "effort"
)
print(df_sec)
#> DOAS Validated Survey Data: 3 sectors (2 DO, 1 Block) | 54 total individuals in 6 groups
#> --------------------------------------------------------------------------------
#>   sector_id survey_type total_count detected_groups groups_both groups_obs1
#> 1         A          DO          25               3           1           1
#> 2         B          DO          14               2           1           1
#> 3         C       Block          15               1          NA          NA
#>   groups_obs2 mean_group_size effort
#> 1           1        8.333333      8
#> 2           0        7.000000      6
#> 3          NA       15.000000      5
```
