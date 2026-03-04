devtools::load_all()

set.seed(42)

# Using the EXACT defaults from doas_power() which match the paper's baseline
cat("Running DOAS Power Analysis (PAPER BASELINE) with Chapman (Independent) Estimator...\n")
res_chapman_paper <- doas_power(
    n_sims = 100,
    surv_areas = c(2, 4, 6, 8),
    # The rest are default = paper's baseline
    method = "hybrid",
    ci_method = "lognormal",
    obs_type = "independent"
)

cat("\nRunning DOAS Power Analysis (PAPER BASELINE) with Bayesian Estimator...\n")
res_bayesian_paper <- doas_power(
    n_sims = 100,
    surv_areas = c(2, 4, 6, 8),
    method = "hybrid",
    ci_method = "lognormal",
    obs_type = "bayesian"
)

cat("\n--- CHAPMAN ESTIMATOR PERFORMANCE (PAPER BASELINE) ---\n")
print(round(res_chapman_paper, 3))

cat("\n--- BAYESIAN ESTIMATOR PERFORMANCE (PAPER BASELINE) ---\n")
print(round(res_bayesian_paper, 3))
