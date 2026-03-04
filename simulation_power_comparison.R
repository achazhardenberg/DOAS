devtools::load_all()

set.seed(42)

cat("Running DOAS Power Analysis with Chapman (Independent) Estimator...\n")
res_chapman <- doas_power(
    n_sims = 100,
    surv_areas = c(2, 4, 6, 8),
    pop_size = 500,
    mean_group = 10,
    tot_areas = 15,
    det_p1 = 0.5,
    intra_group_det = 0.9,
    det_group = 0.5,
    obs_effect = 0.5, # Make observer 2 weak to induce sparsity/low detection
    detvar = 0.3,
    method = "hybrid",
    ci_method = "lognormal",
    obs_type = "independent"
)

cat("\nRunning DOAS Power Analysis with Bayesian Estimator...\n")
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
    method = "hybrid",
    ci_method = "lognormal",
    obs_type = "bayesian"
)

cat("\n--- CHAPMAN ESTIMATOR PERFORMANCE ---\n")
print(round(res_chapman, 3))

cat("\n--- BAYESIAN ESTIMATOR PERFORMANCE ---\n")
print(round(res_bayesian, 3))
