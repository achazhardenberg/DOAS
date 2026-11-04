test_that("chapman_estimate calculates correctly", {
    both <- c(10, 5)
    obs1 <- c(2, 1)
    obs2 <- c(3, 2)

    # (10+2+1)*(10+3+1)/(10+1) - 1 = 13 * 14 / 11 - 1 = 182 / 11 - 1 = 16.5454... - 1 = 15.5454...
    res <- chapman_estimate(both, obs1, obs2)
    expect_equal(res[1], (13 * 14 / 11) - 1, tolerance = 1e-5)

    # (5+1+1)*(5+2+1)/(5+1) - 1 = 7 * 8 / 6 - 1 = 56 / 6 - 1 = 9.333... - 1 = 8.333...
    expect_equal(res[2], (7 * 8 / 6) - 1, tolerance = 1e-5)
})

test_that("doas computes hybrid and global frequentist estimates correctly", {
    data <- data.frame(
        sector_id = c(1, 2, 3),
        total_count = c(50, 30, 100),
        groups_both = c(10, 5, NA),
        groups_obs1 = c(2, 1, NA),
        groups_obs2 = c(3, 2, NA),
        mean_group_size = c(3, 4, NA),
        detected_groups = c(15, 8, 20)
    )

    res <- doas(data, method = "frequentist", obs_type = "independent", pooling = "hybrid")

    expect_s3_class(res, "doas")
    expect_equal(res$method, "frequentist")
    expect_equal(res$obs_type, "independent")
    expect_equal(res$pooling, "hybrid")
    expect_equal(res$num_do_sectors, 2)
    expect_equal(res$num_total_sectors, 3)
    expect_true(res$total_abundance > 0)

    res_glob <- doas(data, method = "frequentist", obs_type = "independent", pooling = "global")
    expect_equal(res_glob$pooling, "global")
    expect_true(res_glob$total_abundance > 0)

    # Verify exact symmetry between Observer 1 and Observer 2 for independent observers
    data_swapped <- data
    data_swapped$groups_obs1 <- data$groups_obs2
    data_swapped$groups_obs2 <- data$groups_obs1
    res_swapped <- doas(data_swapped, method = "frequentist", obs_type = "independent", pooling = "hybrid")
    expect_equal(res$total_abundance, res_swapped$total_abundance, tolerance = 1e-10)
    expect_equal(res$se_abundance, res_swapped$se_abundance, tolerance = 1e-10)
})

test_that("simulate_pop and simulate_survey run without error", {
    set.seed(42)
    pop <- simulate_pop(pop_size = 500, mean_group = 10, tot_areas = 5)
    expect_type(pop, "list")
    expect_true("tot_pop" %in% names(pop))

    survey <- simulate_survey(
        pop$sim_pop,
        det_p1 = 0.6,
        intra_group_det = 0.9,
        det_group = 0.5,
        obs_effect = 0.75,
        detvar = 0.3,
        surv_areas = 2
    )
    expect_s3_class(survey, "data.frame")
    expect_equal(nrow(survey), 5)
    do_sector_count <- sum(!is.na(survey$groups_both))
    expect_equal(do_sector_count, 2)
})

test_that("doas_power runs without error", {
    set.seed(42)
    power_res <- doas_power(
        n_sims = 2,
        surv_areas = c(2, 4),
        pop_size = 100,
        tot_areas = 5,
        method = "frequentist",
        obs_type = "independent",
        pooling = "hybrid"
    )

    expect_s3_class(power_res, "data.frame")
    expect_equal(nrow(power_res), 2)
    expect_true(all(
        c(
            "surv_areas",
            "mean_true_N",
            "mean_est_N",
            "rel_bias",
            "rmse",
            "coverage",
            "mean_cv"
        ) %in%
            names(power_res)
    ))

    # Test power analysis with simulated sector covariate and formula
    power_cov <- doas_power(
        n_sims = 2,
        surv_areas = 4,
        pop_size = 200,
        tot_areas = 6,
        cov_effect = 0.5,
        formula = ~covariate
    )
    expect_s3_class(power_cov, "data.frame")
    expect_equal(nrow(power_cov), 1)
})

test_that("doas with formula computes covariate models correctly", {
    data <- data.frame(
        sector_id = c(1, 2, 3, 4, 5),
        total_count = c(50, 30, 80, 20, 100),
        groups_both = c(10, 5, 12, NA, NA),
        groups_obs1 = c(2, 1, 2, NA, NA),
        groups_obs2 = c(3, 2, 3, NA, NA),
        mean_group_size = c(3, 4, 4.5, NA, NA),
        detected_groups = c(15, 8, 17, 5, 20),
        effort = c(1.5, 0.8, 2.0, 0.5, 2.5)
    )

    res <- doas(data, formula = ~effort)

    expect_s3_class(res, "doas")
    expect_equal(res$method, "frequentist")
    expect_equal(res$pooling, "hybrid")
    expect_equal(res$num_do_sectors, 3)
    expect_s3_class(res$model, "glm")
    expect_equal(deparse(res$formula), "~effort")
    expect_true(res$total_abundance > 0)
    expect_true(res$se_abundance > 0)
})

test_that("doas computes bootstrap confidence intervals correctly", {
    data <- data.frame(
        sector_id = c(1, 2, 3, 4, 5),
        total_count = c(50, 30, 80, 20, 100),
        groups_both = c(10, 5, 12, NA, NA),
        groups_obs1 = c(2, 1, 2, NA, NA),
        groups_obs2 = c(3, 2, 3, NA, NA),
        mean_group_size = c(3, 4, 4.5, NA, NA),
        detected_groups = c(15, 8, 17, 5, 20)
    )

    set.seed(42)
    res <- doas(data, ci_method = "bootstrap", R = 10)

    expect_s3_class(res, "doas")
    expect_equal(res$ci_method, "bootstrap")
    expect_true(res$total_abundance > 0)
    expect_true(res$se_abundance > 0)
    expect_true(res$ci_lower >= sum(data$total_count))
    expect_true(res$ci_upper >= res$ci_lower)
})

test_that("QAIC is calculated correctly for covariate models", {
    data <- data.frame(
        sector_id = 1:6,
        total_count = c(55, 30, 80, 100, 20, 0),
        groups_both = c(6, 3, 10, 8, NA, NA),
        groups_obs1 = c(1, 1, 2, 2, NA, NA),
        groups_obs2 = c(1, 1, 1, 1, NA, NA),
        mean_group_size = c(6.875, 6.0, 5.0, 8.0, NA, NA),
        detected_groups = c(8, 5, 13, 11, 2, 0),
        effort = c(12, 10, 8, 4, 15, 2),
        habitat = c("A", "A", "B", "B", "B", "A")
    )

    mod1 <- doas(data, formula = ~effort)
    mod2 <- doas(data, formula = ~ effort + habitat)
    mod3 <- doas(data)

    qaic1 <- QAIC(mod1)
    qaic2 <- QAIC(mod2)

    expect_true(is.numeric(qaic1))
    expect_true(is.numeric(qaic2))

    expect_warning(qaic3 <- QAIC(mod3))
    expect_true(is.na(qaic3))

    comp <- QAIC(mod1, mod2)
    expect_s3_class(comp, "data.frame")
    expect_equal(nrow(comp), 2)
    expect_true("QAIC" %in% names(comp))
    expect_equal(sort(rownames(comp)), c("mod1", "mod2"))
})

test_that("doas computes frequentist dependent (Zippin) estimates correctly", {
    data <- data.frame(
        sector_id = 1:5,
        total_count = c(39, 28, 67, 20, 100),
        groups_both = c(10, 5, 12, NA, NA),
        groups_obs1 = c(2, 1, 2, NA, NA),
        groups_obs2 = c(1, 1, 1, NA, NA),
        mean_group_size = c(3, 4, 4.5, NA, NA),
        detected_groups = c(13, 7, 15, 5, 20)
    )

    res_dep <- doas(data, method = "frequentist", obs_type = "dependent")
    res_ind <- doas(data, method = "frequentist", obs_type = "independent")

    expect_s3_class(res_dep, "doas")
    expect_equal(res_dep$method, "frequentist")
    expect_equal(res_dep$obs_type, "dependent")
    expect_false(res_dep$total_abundance == res_ind$total_abundance)
    expect_true(res_dep$total_abundance > 0)
    expect_true(res_dep$se_abundance > 0)
    expect_true(res_dep$ci_lower >= sum(data$total_count))
    expect_true(res_dep$ci_upper >= res_dep$ci_lower)
})

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
        formula = ~effort,
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

test_that("doas_data formats and validates sector, group-long, and sighting-long datasets", {
    # 1. Group-level long format
    raw_groups <- data.frame(
        valley = c("A", "A", "A", "B", "B", "C"),
        size   = c(12, 5, 8, 10, 4, 15),
        who    = c("both", "obs1", "obs2", "both", "obs1", "block"),
        effort = c(8, 8, 8, 6, 6, 5)
    )
    df_grp <- doas_data(
        raw_groups,
        format = "group",
        sector = "valley",
        count = "size",
        seen_by = "who",
        do_sectors = c("A", "B"),
        covariates = "effort"
    )
    expect_s3_class(df_grp, "doas_data")
    expect_equal(nrow(df_grp), 3)
    expect_equal(df_grp$total_count, c(25, 14, 15))
    expect_equal(df_grp$groups_both, c(1, 1, NA))
    expect_equal(df_grp$groups_obs1, c(1, 1, NA))
    expect_equal(df_grp$groups_obs2, c(1, 0, NA))

    # Fit doas on the inspected doas_data object
    res_grp <- doas(df_grp)
    expect_s3_class(res_grp, "doas")
    expect_true(res_grp$total_abundance > sum(df_grp$total_count))

    # 2. Sighting-level long format (1 row per observer sighting, linked by group_id)
    raw_sightings <- data.frame(
        sector_id = c("S1", "S1", "S1", "S2", "S2", "S3"),
        observer  = c("Obs1", "Obs2", "Obs1", "Obs1", "Obs2", "Obs1"),
        group_id  = c("G1", "G1", "G2", "G3", "G3", "G4"),
        group_size = c(10, 10, 6, 8, 8, 20)
    )
    df_sight <- doas_data(
        raw_sightings,
        format = "sighting",
        do_sectors = c("S1", "S2")
    )
    expect_equal(df_sight$groups_both, c(1, 1, NA))
    expect_equal(df_sight$groups_obs1, c(1, 0, NA))
    expect_equal(df_sight$detected_groups, c(2, 1, 1))

    # 3. Sector-level wide format without groups_both (auto-computed from detected_groups - obs1 - obs2)
    sec_no_both <- data.frame(
        sector_id = c("S1", "S2", "S3"),
        total_count = c(50, 30, 40),
        detected_groups = c(10, 6, 8),
        groups_obs1 = c(2, 1, NA),
        groups_obs2 = c(1, 1, NA)
    )
    df_no_both <- doas_data(sec_no_both, format = "sector")
    expect_equal(df_no_both$groups_both, c(7, 4, NA))
})

