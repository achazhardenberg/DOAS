# Paired Power Benchmark: Original Paper (Gran Paradiso) & Sparse Scenarios
devtools::load_all()

set.seed(42)

cat("Running DOAS Power Analysis with Frequentist Chapman Estimator...\n")
res_chapman <- doas_power(
    n_sims = 100,
    surv_areas = c(2, 4, 6, 8),
    pop_size = 500,
    mean_group = 10,
    tot_areas = 15,
    det_p1 = 0.5,
    intra_group_det = 0.9,
    det_group = 0.5,
    obs_effect = 0.5,
    detvar = 0.3,
    method = "frequentist",
    obs_type = "independent",
    pooling = "hybrid",
    ci_method = "lognormal"
)

cat("\nRunning DOAS Power Analysis with Multi-Chain Bayesian Estimator...\n")
res_bayesian <- doas_power(
    n_sims = 100,
    surv_areas = c(2, 4, 6, 8),
    pop_size = 500,
    mean_group = 10,
    tot_areas = 15,
    det_p1 = 0.5,
    intra_group_det = 0.9,
    det_group = 0.5,
    obs_effect = 0.5,
    detvar = 0.3,
    method = "bayesian",
    obs_type = "independent",
    pooling = "hybrid",
    ci_method = "lognormal"
)

cat("\n--- FREQUENTIST CHAPMAN ESTIMATOR PERFORMANCE ---\n")
print(round(res_chapman, 3))

cat("\n--- MULTI-CHAIN BAYESIAN ESTIMATOR PERFORMANCE ---\n")
print(round(res_bayesian, 3))
