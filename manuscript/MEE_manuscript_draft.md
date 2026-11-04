# DOAS: An R Package for Frequentist and Hierarchical Bayesian Double Observer Adjusted Surveys

**Authors:** Matteo Panaccio, Alice Brambilla, Bruno Bassano, Tessa Smith, Achaz von Hardenberg

**Target Journal:** *Methods in Ecology and Evolution* (Application Article)

---

## Abstract

1. **Background:** Ground-based block counts and aerial surveys are widely used to monitor gregarious wildlife populations (e.g., mountain ungulates), yet they routinely underestimate true abundance because detection probability is imperfect ($p < 1$) and varies across space, observers, and group sizes. The Double Observer Adjusted Survey (DOAS) framework addresses the trade-off between survey accuracy and field logistics by pairing independent or dependent double-observer counts in a subset of spatial sectors with single-observer block counts across the remaining sectors.
2. **Objectives & Implementation:** We present `DOAS`, an open-source R package that provides a complete, validated analytical and simulation workflow for Double Observer Adjusted Surveys. The package structures field data at three levels of resolution (`format = "sector"`, `"group"`, or `"sighting"`) via `doas_data()` and supports four core estimation architectures controlled by orthogonal arguments in `doas()`: paradigm (`method = c("frequentist", "bayesian")`), observer configuration (`obs_type = c("independent", "dependent")`), and spatial extrapolation strategy (`pooling = c("hybrid", "global")`).
3. **Methodological Advancements:** Beyond operationalizing the original DOAS framework (Panaccio et al., 2024), `DOAS` introduces four key methodological advances: (i) pooled capture-history estimation of the global detection rate across Double Observer (DO) sectors, which eliminates small-sample Seber (1982) boundary bias and ad-hoc sector truncation when sector-level counts are sparse; (ii) a multi-chain hierarchical Bayesian Markov chain Monte Carlo (MCMC) engine with Gelman–Rubin ($\hat{R}$) and effective sample size ($n_{\text{eff}}$) diagnostics that shrinks noisy sector-level detection rates toward a global hyperparameter; (iii) sector-level covariate modeling (`formula = ~ covariate`) with full multivariate Delta-method variance propagation, non-parametric spatial block bootstrapping, Quasi-AIC (`QAIC()`) for overdispersed frequentist Generalized Linear Models (GLMs), and Widely Applicable Information Criterion (`WAIC()`) for Bayesian models; and (iv) a corrected simulation engine (`simulate_pop()`, `simulate_survey()`, and `doas_power()`) that resolves an observer-efficiency asymmetry in the original simulation script of Panaccio et al. (2024).
4. **Conclusions:** By combining rigorous frequentist and Bayesian inference with prospective power analysis (`doas_power()`), `DOAS` enables wildlife managers and quantitative ecologists to design cost-effective partial double-observer surveys and obtain unbiased abundance estimates with well-calibrated uncertainty intervals.

**Keywords:** abundance estimation, capture–recapture, Chapman estimator, double-observer survey, hierarchical Bayesian model, imperfect detection, power analysis, wildlife monitoring

---

## 1. Introduction

Reliable estimates of population size are foundational to wildlife conservation, harvest management, and the assessment of demographic trends (Williams et al., 2002; Nichols et al., 2009). For open-habitat and mountain ungulates—such as Alpine ibex (*Capra ibex*), chamois (*Rupicapra* spp.), and bighorn sheep (*Ovis canadensis*)—monitoring programs have historically relied on coordinated ground-based block counts (also termed standardized minimum counts) in which field personnel simultaneously scan predefined valleys or management sectors (Largo et al., 2008; Brambilla et al., 2020). Although block counts are logistically tractable over large, rugged landscapes, they rest on the untenable assumption that detection probability is either perfect ($p = 1$) or constant over time and space. In practice, visibility varies with terrain ruggedness, vegetation cover, group size, and observer experience, causing unadjusted block counts to systematically underestimate true abundance and confound population trends with variation in detectability (Caughley, 1974; Pollock et al., 2002).

Two-sample mark–recapture and removal protocols—specifically the independent Double Observer (DO) method (Nichols et al., 2000; Forsyth & Hickling, 1997) and the dependent primary–secondary observer method (Cook & Jacobson, 1979; Suryawanshi et al., 2012)—explicitly estimate detection probability during visual surveys without requiring physical capture or individual markings. However, deploying two trained observers across every sector of a large protected area doubles personnel costs and is rarely feasible within the narrow phenological windows required for synchronized censuses.

To resolve this logistical bottleneck, Panaccio et al. (2024) introduced the **Double Observer Adjusted Survey (DOAS)** framework. Under a DOAS design, only a representative subset of $n_{\text{DO}}$ spatial sectors is surveyed by a double-observer pair, while the remaining $n_{\text{Block}}$ sectors are surveyed by a single observer following standard block-count protocols. The detection probability calibrated within the DO sectors is then used to adjust the single-observer counts across the unsurveyed Block sectors.

While Panaccio et al. (2024) demonstrated the conceptual efficiency of the DOAS design using custom simulation scripts, applying the methodology to empirical monitoring datasets—especially those characterized by small sector sample sizes, dependent observer protocols, or environmental covariates—requires a standardized, statistically rigorous software implementation. Furthermore, closer mathematical examination of the original simulation scripts revealed two important statistical considerations when applying DOAS in practice: (1) averaging unweighted sector-level Chapman detection estimates across sparse DO sectors can introduce small-sample boundary bias when individual sectors contain few groups, and (2) simulating single-observer block sectors with a single fixed observer efficiency rather than the mean efficiency of the interchangeable observer pool created an artificial positive bias in the original simulation benchmarks whenever observers differed in skill.

Here, we present `DOAS`, a comprehensive R package that formalizes, corrects, and substantially extends the Double Observer Adjusted Survey framework. `DOAS` provides:
1. A defensive data-ingestion pipeline (`doas_data()`) that validates and aggregates sighting-level, group-level, or sector-level field records;
2. Unified Frequentist and multi-chain Hierarchical Bayesian estimators (`doas()`) supporting both independent (Chapman/Lincoln–Petersen) and dependent (Zippin removal) observer configurations under `"hybrid"` or `"global"` spatial pooling;
3. Sector-level detection covariate modeling via standard R formulas, paired with multivariate Delta-method variance propagation, non-parametric spatial bootstrap confidence intervals, `QAIC()` for frequentist quasi-binomial models, and `WAIC()` for Bayesian MCMC models; and
4. A corrected, covariate-aware stochastic simulator (`simulate_pop()`, `simulate_survey()`) and prospective survey design engine (`doas_power()`).

---

## 2. Package Architecture and Data Pipeline

### 2.1 Standardizing Field Data with `doas_data()`

Field teams record double-observer surveys at different levels of aggregation depending on field protocols. Rather than performing silent data transformations inside the model-fitting function, `DOAS` provides an explicit preprocessing constructor, `doas_data()`, which converts raw data frames into a validated S3 object of class `"doas_data"` with a dedicated `print.doas_data()` summary method. `doas_data()` supports three input formats via the `format` argument:

* **`format = "sector"` (Wide sector-level summaries):** One row per spatial sector containing total individuals counted (`total_count`), groups seen by both observers (`groups_both`, $B_m$), groups seen only by Observer 1 (`groups_obs1`, $S_{1,m}$), groups seen only by Observer 2 (`groups_obs2`, $S_{2,m}$), and total distinct groups detected (`detected_groups`, $D_m$). Sectors where $B_m, S_{1,m}, S_{2,m}$ are `NA` are automatically classified as single-observer `"Block"` sectors. In `"DO"` sectors, if either `detected_groups` or `groups_both` is omitted (`NA`), `doas_data()` automatically computes the missing column from the identity $D_m = B_m + S_{1,m} + S_{2,m}$, and raises an informative error if both are supplied but arithmetically inconsistent.
* **`format = "group"` (Long group-level encounter histories):** One row per detected animal group, recording the sector ID, group size, and binary detection indicators (`obs1`, `obs2` $\in \{0, 1\}$) for DO sectors (or `NA` for Block sectors), along with optional sector-level covariates.
* **`format = "sighting"` (Long observer-sighting records):** One row per observer sighting event (referencing a `group_id` and `observer` label), reconciled into group-level capture histories and aggregated to sector-level counts.

---

## 3. Statistical Framework and Methodological Refinements

### 3.1 Independent vs. Dependent Observer Configurations

Let a study area be partitioned into $M = n_{\text{DO}} + n_{\text{Block}}$ non-overlapping spatial sectors. In each DO sector $m \in \{1, \dots, n_{\text{DO}}\}$ containing $G_m$ true groups, two observers survey the sector and detect $B_m$ groups jointly, $S_{1,m}$ groups seen exclusively by Observer 1, and $S_{2,m}$ groups seen exclusively by Observer 2, yielding $D_m = B_m + S_{1,m} + S_{2,m}$ distinct observed groups and a mean observed group size $\bar{s}_m = C_m / D_m$, where $C_m$ is the total number of individuals counted in sector $m$.

`DOAS` supports two field protocols via `obs_type`:

1. **Independent Observers (`obs_type = "independent"`):** Both observers scan the sector independently without communicating during the survey. Within sector $m$, the bias-corrected Chapman (1951) estimator of true group abundance is:
   $$\hat{G}_m = \frac{(B_m + S_{1,m} + 1)(B_m + S_{2,m} + 1)}{B_m + 1} - 1$$
   with finite-sample variance given by Seber (1982):
   $$\widehat{\text{Var}}(\hat{G}_m) = \frac{(B_m + S_{1,m} + 1)(B_m + S_{2,m} + 1) S_{1,m} S_{2,m}}{(B_m + 1)^2 (B_m + 2)}$$
   Because rangers or observers assigned to single-observer Block sectors are drawn interchangeably from the same personnel pool as the two-observer teams, the expected single-observer detection probability is the symmetric mean of Observer 1 and Observer 2:
   $$\bar{p} = \frac{p_1 + p_2}{2}$$

2. **Dependent Observers (`obs_type = "dependent"`):** A primary observer (Observer 1) scans and announces detected groups ($n_{1,m} = B_m + S_{1,m}$), while a secondary observer (Observer 2) records missed groups ($n_{2,m} = S_{2,m}$). Following the two-pass Zippin (1956) removal model with $S_{1,m} = 0$, the primary observer's detection probability is estimated as $\hat{p}_{1,m} = 1 - \frac{n_{2,m}}{n_{1,m} + 1}$, yielding $\hat{G}_m = \frac{n_{1,m} + n_{2,m}}{1 - (1 - \hat{p}_{1,m})^2}$.

### 3.2 Spatial Extrapolation: Hybrid vs. Global Pooling

The `pooling` argument governs how DO sector information is combined with Block sectors:
* **`pooling = "hybrid"` (Default):** Preserves the direct mark–recapture group abundance estimates $\hat{G}_m$ inside the $n_{\text{DO}}$ Double Observer sectors, and extrapolates only the $n_{\text{Block}}$ single-observer sectors using the calibrated single-observer detection probability $\hat{p}_i$:
  $$\hat{N}_{\text{hybrid}} = \sum_{m \in \text{DO}} \hat{G}_m \bar{s}_m + \sum_{b \in \text{Block}} \frac{D_b}{\hat{p}_b} \bar{s}_b$$
* **`pooling = "global"`:** Applies the calibrated single-observer detection probability $\hat{p}_i$ uniformly across all $M$ sectors (using the single-observer equivalent group count $B_m + 0.5(S_{1,m} + S_{2,m})$ in DO sectors for independent observers, or $n_{1,m}$ for dependent observers).

In both modes, analytical variance is propagated across all extrapolated sectors by accounting for both binomial detection variance within each sector and the shared estimation covariance of $\hat{p}$ via the multivariate Delta method:
$$\widehat{\text{Var}}\!\left(\sum_{b} \frac{C_b}{\hat{p}_b}\right) = \sum_{b} C_b \bar{s}_b \frac{1 - \hat{p}_b}{\hat{p}_b^2} + \mathbf{g}_{\boldsymbol{\beta}}^\top \widehat{\boldsymbol{\Sigma}}_{\boldsymbol{\beta}} \mathbf{g}_{\boldsymbol{\beta}}$$
where $\mathbf{g}_{\boldsymbol{\beta}} = \nabla_{\boldsymbol{\beta}} \sum_{b} (C_b / \hat{p}_b)$ and $\widehat{\boldsymbol{\Sigma}}_{\boldsymbol{\beta}}$ is the parameter covariance matrix. Confidence intervals are constructed on the log-normal scale for the unobserved remainder $\hat{N} - C_{\text{total}}$ (Burnham et al., 1987; Chao, 1987), guaranteeing that the lower confidence bound never falls below the minimum number of individuals physically counted in the field ($C_{\text{total}} = \sum_{i=1}^M C_i$).

### 3.3 Methodological Refinements Relative to Panaccio et al. (2024)

Developing and testing the `DOAS` package uncovered two important properties of the original implementation in Panaccio et al. (2024) (`DOAS_script.R`) that we resolve here:

1. **Pooled Capture Histories vs. Unweighted Mean of Sector-Level $\hat{p}_m$:**
   In `DOAS_script.R`, the global single-observer detection rate was computed as the unweighted arithmetic mean of sector-specific Chapman detection rates across the $n_{\text{DO}}$ sectors: $\hat{p}_{\text{legacy}} = \frac{1}{n_{\text{DO}}} \sum_{m=1}^{n_{\text{DO}}} \frac{(B_m + S_{1,m}) + (B_m + S_{2,m})}{2\hat{G}_m}$, combined with an ad-hoc truncation rule setting $\hat{p}_m = 0$ whenever $B_m + S_{1,m} \le 1$ or $B_m + S_{2,m} = 0$. In small or low-density sectors where $S_{2,m} = 0$ (Observer 2 detects no unique groups beyond those seen by Observer 1), the sector-level Chapman formula yields $\hat{G}_m = B_m + S_{1,m}$ and thus forces $\hat{p}_{1,m} = 1.0$, creating a systematic positive boundary bias in $\hat{p}_m$ (or a downward spike if truncated to zero). In `DOAS`, the intercept-only frequentist estimator instead pools the sufficient capture histories across all $n_{\text{DO}}$ sectors ($B_{\text{sum}} = \sum B_m$, $S_{1,\text{sum}} = \sum S_{1,m}$, $S_{2,\text{sum}} = \sum S_{2,m}$) before computing the pooled Chapman abundance $\hat{G}_{\text{pooled}}$ and the symmetric mean observer detection rate:
   $$\hat{p}_{\text{pooled}} = \frac{2 B_{\text{sum}} + S_{1,\text{sum}} + S_{2,\text{sum}}}{2 \hat{G}_{\text{pooled}}}$$
   This eliminates small-sector ratio bias and removes the need for arbitrary truncation rules.

2. **Symmetric Observer Pool Efficiency in Block-Count Simulations:**
   In the simulation study of Panaccio et al. (2024), DOAS abundance estimates exhibited a modest positive bias ($\hat{N} / N_{\text{true}} > 1$) whenever Observer 2 had lower detection efficiency than Observer 1 ($\text{obs.effect} < 1$). Inspection of `DOAS_script.R` revealed that this overestimate was an artifact of the simulation function (`block.count()`) rather than an intrinsic flaw of the DOAS estimator: `block.count()` simulated all single-observer Block sectors using Observer 1's baseline detection rate ($p_{\text{block}} = p_1$), whereas the estimator divided those Block counts by the average detection probability of Observer 1 and Observer 2 ($\bar{p} = \frac{p_1 + p_2}{2} = p_1 \frac{1 + \text{obs.effect}}{2}$). Consequently, simulated Block-sector abundance estimates were inflated by a factor of $\frac{p_1}{\bar{p}} = \frac{2}{1 + \text{obs.effect}}$ (corresponding to $+14.3\%$ when $\text{obs.effect} = 0.75$ and $+33.3\%$ when $\text{obs.effect} = 0.50$). Because field rangers in Gran Paradiso National Park and similar monitoring programs are assigned across DO and Block sectors from the same personnel pool, the expected detection efficiency in a single-observer Block sector is $\bar{p} = p_1 \frac{1 + \text{obs.effect}}{2}$. Correcting `simulate_survey()` in `DOAS` to use this symmetric observer-pool efficiency eliminates the apparent overestimation across all observer-effect levels (Section 5).

---

## 4. Hierarchical Bayesian Estimation and Covariate Selection

### 4.1 Multi-Chain Hierarchical Bayesian MCMC (`method = "bayesian"`)

When the number of surveyed DO sectors or the number of groups per sector is small, frequentist point estimates can remain sensitive to sampling variability. Setting `method = "bayesian"` in `doas()` activates a pure-R multi-chain Markov chain Monte Carlo (MCMC) sampler that combines data augmentation for the unobserved group counts $U_m = G_m - D_m \sim \text{NegBin}(D_m + 1, 1 - q_m)$ with conjugate Beta updates (or Metropolis–Hastings logit-normal random-effect updates when sector covariates are included).

For each MCMC iteration, posterior draws of total abundance $N^{(t)}$ respect the physical floor $C_{\text{total}}$, and `doas()` automatically computes multi-chain convergence diagnostics across all parameters, including the split Gelman–Rubin potential scale reduction factor ($\hat{R}$; Gelman & Rubin, 1992) and effective sample size ($n_{\text{eff}}$).

### 4.2 Covariate Modeling (`formula`) and Model Selection (`QAIC` and `WAIC`)

Both Frequentist and Bayesian estimators accept a one-sided R formula (e.g., `formula = ~ effort + visibility`) to model variation in sector-level detection probability on the logit scale:
$$\text{logit}(p_i) = \mathbf{x}_i^\top \boldsymbol{\beta}$$
To compare candidate covariate models, `DOAS` provides two information criteria tailored to each paradigm:
* **`QAIC()` for Frequentist Models:** Because grouped ecological counts frequently exhibit overdispersion relative to the binomial likelihood, `QAIC()` estimates the variance inflation factor $\hat{c}$ from the Pearson residuals of the global (most parameterized) model in a candidate set and computes the Quasi-Akaike Information Criterion (Burnham & Anderson, 2002):
  $$\text{QAIC} = -\frac{2 \log \mathcal{L}(\hat{\boldsymbol{\beta}})}{\hat{c}} + 2K$$
  returning a ranked comparison table with $\hat{c}$, parameter count $K$, $\Delta\text{QAIC}$, and Akaike model weights.
* **`WAIC()` for Bayesian Models:** Uses the pointwise log-likelihood matrix stored across posterior draws of the DO sectors to compute the Widely Applicable Information Criterion ($\text{WAIC} = -2(\text{lppd} - p_{\text{WAIC}})$; Watanabe, 2010; Vehtari et al., 2017) along with $\Delta\text{WAIC}$ and model weights.

---

## 5. Simulation Benchmark: Original Script vs. `DOAS` Package

To illustrate the impact of the simulation and estimation refinements described in Section 3.3, we conducted a Monte Carlo benchmark (`manuscript/benchmark_MEE_paper.R`) across a 38-sector virtual landscape ($N_{\text{true}} \approx 3500$, mean group size $= 10$, base detection $p_1 = 0.60$, group-size scaling $\text{det\_group} = 0.50$, inter-sector SD $\text{detvar} = 0.10$). We varied the relative efficiency of Observer 2 ($\text{obs\_effect} \in \{0.50, 0.75, 1.00\}$) and the number of Double Observer sectors ($n_{\text{DO}} \in \{5, 10, 19\}$ out of 38 sectors), running 50 replicates per combination.

| Observer Effect (`obs_effect`) | DO Sectors ($n_{\text{DO}}$) | Legacy Script (`DOAS_script.R` Asymmetric Sim) | Legacy Unweighted Estimator (Corrected Sim) | `DOAS` Frequentist (`pooling = "hybrid"`) | `DOAS` Frequentist (`pooling = "global"`) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| **0.50** | 5 | 1.252 | 1.003 | **0.973** | **0.987** |
| **0.50** | 10 | 1.187 | 0.968 | **0.978** | **1.006** |
| **0.50** | 19 | 1.128 | 0.941 | **0.946** | **0.995** |
| **0.75** | 5 | 1.058 | 0.974 | **1.006** | **1.011** |
| **0.75** | 10 | 1.073 | 0.989 | **1.025** | **1.040** |
| **0.75** | 19 | 1.013 | 0.955 | **0.989** | **1.018** |
| **1.00** | 5 | 1.017 | 1.004 | **0.979** | **0.981** |
| **1.00** | 10 | 0.978 | 1.013 | **1.033** | **1.038** |
| **1.00** | 19 | 0.972 | 0.981 | **1.001** | **1.015** |

*Table 1: Mean relative abundance accuracy ($\hat{N} / N_{\text{true}}$; ideal $= 1.000$) across 50 simulated surveys per scenario. Column 3 shows the original `DOAS_script.R` pipeline where Block sectors were simulated with $p_1$ instead of $\bar{p}$, producing the theoretical inflation factor $\approx 2/(1 + \text{obs\_effect})$. Column 4 applies the original unweighted sector-mean estimator to the corrected symmetric simulation, and Columns 5–6 show the `DOAS` package pooled Frequentist Hybrid and Global estimators.*

Furthermore, under sparse, low-density conditions ($N_{\text{true}} \approx 800$ across 38 sectors, $p_1 = 0.35$, $n_{\text{DO}} = 6$), where individual DO sectors frequently contain $0\text{--}2$ detected groups, the unweighted sector-mean estimator in `DOAS_script.R` suffers from severe instability due to zeroed sectors (`obs1 <= 1 || obs2 == 0`), whereas the pooled `DOAS` frequentist and hierarchical Bayesian estimators remain stable and well-calibrated.

---

## 6. Worked Workflow and Prospective Power Analysis

### 6.1 Data Preparation, Estimation, and Model Selection

The following reproducible workflow demonstrates how field biologists can format raw survey data, fit competing frequentist and Bayesian DOAS models, and compare covariate structures:

```r
library(DOAS)

# 1. Simulate a 38-sector population and survey with a sector covariate effect
set.seed(123)
pop <- simulate_pop(pop_size = 3000, mean_group = 10, tot_areas = 38)
raw_survey <- simulate_survey(
  sim_pop = pop$sim_pop, det_p1 = 0.65, intra_group_det = 0.95,
  det_group = 0.5, obs_effect = 0.80, detvar = 0.10,
  surv_areas = 10, cov_effect = 0.6
)

# 2. Validate and inspect the survey table using doas_data()
survey_dat <- doas_data(raw_survey, format = "sector")
print(survey_dat)

# 3. Fit Frequentist null and covariate models and rank via QAIC
fit_null <- doas(survey_dat, method = "frequentist", obs_type = "independent", pooling = "hybrid")
fit_cov  <- doas(survey_dat, method = "frequentist", obs_type = "independent",
                 pooling = "hybrid", formula = ~ covariate)
QAIC(fit_null, fit_cov)
summary(fit_cov)

# 4. Fit Hierarchical Bayesian models and rank via WAIC
bayes_null <- doas(survey_dat, method = "bayesian", obs_type = "independent",
                   pooling = "hybrid", chains = 3, iter = 2000, burnin = 500)
bayes_cov  <- doas(survey_dat, method = "bayesian", obs_type = "independent",
                   pooling = "hybrid", formula = ~ covariate,
                   chains = 3, iter = 2000, burnin = 500)
WAIC(bayes_null, bayes_cov)
summary(bayes_cov)
```

### 6.2 Prospective Survey Design with `doas_power()`

Before deploying field teams, wildlife managers must determine how many sectors require Double Observer coverage ($n_{\text{DO}}$) to achieve target precision (e.g., Coefficient of Variation $\text{CV} < 10\%$) and unbiasedness. The `doas_power()` function automates this prospective power analysis across any combination of `method`, `obs_type`, `pooling`, and `formula`:

```r
# Evaluate survey accuracy across 4, 8, 12, and 16 Double Observer sectors
power_res <- doas_power(
  n_sims = 100,
  pop_size = 3000,
  mean_group = 10,
  tot_areas = 38,
  surv_areas_vec = c(4, 8, 12, 16),
  det_p1 = 0.60,
  intra_group_det = 0.95,
  det_group = 0.50,
  obs_effect = 0.80,
  detvar = 0.15,
  method = "frequentist",
  obs_type = "independent",
  pooling = "hybrid"
)

print(power_res$summary)
plot(power_res)
```

---

## 7. Conclusions

The `DOAS` R package translates the Double Observer Adjusted Survey framework into an accessible, modular, and statistically rigorous toolkit for ecology and wildlife management. By resolving small-sample pooling biases, correcting single-observer block-count simulation assumptions, and introducing hierarchical Bayesian MCMC estimation alongside `QAIC` and `WAIC` model selection, `DOAS` empowers practitioners to account for imperfect detection across large landscapes at a fraction of the cost of full double-observer censuses.

---

## Data and Code Availability

The `DOAS` R package is open-source under the GPL-3 license and includes complete documentation, unit tests, and vignettes. Reproducible benchmark scripts for all analyses presented in this paper are provided in `manuscript/benchmark_MEE_paper.R`.

---

## References

* Brambilla, A., Von Hardenberg, A., Nelli, L., & Bassano, B. (2020). Distribution, status, and recent population dynamics of Alpine ibex *Capra ibex* in Europe. *Mammal Review*, 50(3), 267–277.
* Burnham, K. P., & Anderson, D. R. (2002). *Model Selection and Multimodel Inference: A Practical Information-Theoretical Approach* (2nd ed.). Springer-Verlag.
* Burnham, K. P., Anderson, D. R., White, G. C., Brownie, C., & Pollock, K. H. (1987). *Design and Analysis Methods for Fish Survival Experiments Based on Release-Recapture*. American Fisheries Society Monograph 5.
* Caughley, G. (1974). Bias in aerial survey. *The Journal of Wildlife Management*, 38(4), 921–933.
* Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. *Biometrics*, 43(4), 783–791.
* Chapman, D. G. (1951). Some properties of the hypergeometric distribution with applications to zoological sample censuses. *University of California Publications in Statistics*, 1, 131–160.
* Cook, R. D., & Jacobson, J. O. (1979). A design for estimating visibility bias in aerial surveys. *Biometrics*, 35(4), 735–742.
* Forsyth, D. M., & Hickling, G. J. (1997). An improved technique for indexing abundance of Himalayan thar (*Hemitragus jemlahicus*). *New Zealand Journal of Ecology*, 21(1), 97–101.
* Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. *Statistical Science*, 7(4), 457–472.
* Largo, E., Gaillard, J. M., Festa-Bianchet, M., Toïgo, C., Bassano, B., Cortot, H., ... & Von Hardenberg, A. (2008). Can ground counts reliably monitor ibex *Capra ibex* populations? *Wildlife Biology*, 14(4), 489–499.
* Nichols, J. D., Hines, J. E., Sauer, J. R., Fallon, F. W., Fallon, J. E., & Heglund, P. J. (2000). A double-observer approach for estimating detection probability and abundance from point counts. *The Auk*, 117(2), 393–408.
* Nichols, J. D., Thomas, L., & Conn, P. B. (2009). Inferences about landbird abundance from count data: recent advances and future directions. In *Modeling Demographic Processes in Marked Populations* (pp. 201–235). Springer.
* Panaccio, M., Brambilla, A., Bassano, B., Smith, T., & von Hardenberg, A. (2024). A new double observer-based census framework to improve abundance estimations in mountain ungulates and other gregarious species with a reduced effort. *Ecological Solutions and Evidence*, 5(4), e12405. https://doi.org/10.1002/2688-8319.12405
* Pollock, K. H., Nichols, J. D., Simons, T. R., Farnsworth, G. L., Bailey, L. L., & Sauer, J. R. (2002). Large scale wildlife monitoring studies: statistical methods for design and analysis. *Environmetrics*, 13(2), 105–119.
* Seber, G. A. F. (1982). *The Estimation of Animal Abundance and Related Parameters* (2nd ed.). Blackburn Press.
* Suryawanshi, K. R., Bhatnagar, Y. V., & Mishra, C. (2012). Standardizing the double-observer survey method for estimating mountain ungulateprey of the endangered snow leopard. *Oecologia*, 169(3), 581–590.
* Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. *Statistics and Computing*, 27(5), 1413–1432.
* Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. *Journal of Machine Learning Research*, 11, 3571–3594.
* Williams, B. K., Nichols, J. D., & Conroy, M. J. (2002). *Analysis and Management of Animal Populations*. Academic Press.
* Zippin, C. (1956). An evaluation of the removal method of estimating animal populations. *Biometrics*, 12(2), 163–189.
