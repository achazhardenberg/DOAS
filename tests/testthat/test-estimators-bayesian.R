test_that("doas computes multi-chain bayesian independent and dependent estimates and diagnostics", {
    data <- data.frame(
        sector_id = 1:5,
        total_count = c(39, 28, 67, 20, 100),
        groups_both = c(10, 5, 12, NA, NA),
        groups_obs1 = c(2, 1, 2, NA, NA),
        groups_obs2 = c(1, 1, 1, NA, NA),
        mean_group_size = c(3, 4, 4.5, NA, NA),
        detected_groups = c(13, 7, 15, 5, 20),
        effort = c(1.5, 0.8, 2.0, 0.5, 2.5)
    )

    set.seed(42)
    res_bayes_ind <- doas(
        data,
        method = "bayesian",
        obs_type = "independent",
        pooling = "hybrid",
        chains = 4,
        iter = 1000,
        burnin = 250
    )

    expect_s3_class(res_bayes_ind, "doas")
    expect_equal(res_bayes_ind$method, "bayesian")
    expect_equal(res_bayes_ind$obs_type, "independent")
    expect_true(res_bayes_ind$total_abundance > 0)
    expect_true(res_bayes_ind$se_abundance > 0)
    expect_true(res_bayes_ind$ci_lower >= sum(data$total_count))
    expect_true(res_bayes_ind$ci_upper >= res_bayes_ind$ci_lower)

    # Check multi-chain MCMC diagnostics
    expect_false(is.null(res_bayes_ind$mcmc_diagnostics))
    expect_equal(res_bayes_ind$mcmc_diagnostics$chains, 4)
    expect_equal(res_bayes_ind$mcmc_diagnostics$post_warmup_draws, 3000)
    expect_true(res_bayes_ind$mcmc_diagnostics$rhat_abundance < 1.1)
    expect_true(res_bayes_ind$mcmc_diagnostics$ess_abundance > 100)

    # Check exact additivity of sector estimates to total abundance
    expect_equal(
        sum(res_bayes_ind$sector_estimates$est_abundance),
        res_bayes_ind$total_abundance,
        tolerance = 1e-6
    )

    # Test Bayesian Dependent (2-pass removal)
    res_bayes_dep <- doas(
        data,
        method = "bayesian",
        obs_type = "dependent",
        pooling = "hybrid",
        chains = 4,
        iter = 1000,
        burnin = 250
    )
    expect_s3_class(res_bayes_dep, "doas")
    expect_equal(res_bayes_dep$method, "bayesian")
    expect_equal(res_bayes_dep$obs_type, "dependent")
    expect_true(res_bayes_dep$total_abundance > sum(data$total_count))
    expect_true(res_bayes_dep$mcmc_diagnostics$rhat_abundance < 1.1)

    # Test edge case where x1 <= x2 (Zippin fails with NA variance, Bayesian Dependent succeeds)
    edge_res <- bayesian_conjugate_estimate(
        both = c(0, 1),
        obs1_only = c(2, 1),
        obs2_only = c(3, 2),
        obs_type = "dependent",
        chains = 2,
        iter = 600,
        burnin = 150
    )
    expect_true(all(edge_res$mean > c(5, 4)))
    expect_true(all(edge_res$var > 0))

    # Test Bayesian covariate model (Metropolis-within-Gibbs)
    res_bayes_cov <- doas(
        data,
        method = "bayesian",
        obs_type = "independent",
        p_formula = ~effort,
        chains = 3,
        iter = 800,
        burnin = 200
    )
    expect_s3_class(res_bayes_cov$bayes_coef, "data.frame")
    expect_equal(nrow(res_bayes_cov$bayes_coef), 2)
    expect_true(all(c("Mean", "SD", "CI_2.5", "CI_97.5", "Rhat", "ESS") %in% names(res_bayes_cov$bayes_coef)))

    # Test WAIC on single Bayesian model and comparing null vs covariate Bayesian models
    waic_single <- WAIC(res_bayes_cov)
    expect_true(is.numeric(waic_single))
    expect_false(is.na(waic_single))
    expect_true(!is.null(attr(waic_single, "p_WAIC")))

    waic_comp <- WAIC(res_bayes_ind, res_bayes_cov)
    expect_s3_class(waic_comp, "data.frame")
    expect_equal(nrow(waic_comp), 2)
    expect_true(all(c("lppd", "p_WAIC", "WAIC", "SE", "delta_WAIC") %in% names(waic_comp)))
})
