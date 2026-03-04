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

test_that("doas computes hybrid estimate correctly", {
    data <- data.frame(
        sector_id = c(1, 2, 3),
        total_count = c(50, 30, 100),
        groups_both = c(10, 5, NA),
        groups_obs1 = c(2, 1, NA),
        groups_obs2 = c(3, 2, NA),
        mean_group_size = c(3, 4, NA),
        detected_groups = c(15, 8, 20)
    )

    res <- doas(data, method = "hybrid")

    expect_s3_class(res, "doas")
    expect_equal(res$method, "hybrid")
    expect_equal(res$num_do_sectors, 2)
    expect_equal(res$num_total_sectors, 3)

    # Check if it has a sensible total abundance computation
    expect_true(res$total_abundance > 0)
})

test_that("simulate_pop and simulate_survey run without error", {
    set.seed(42)
    pop <- simulate_pop(pop_size = 500, mean_group = 10, tot_areas = 5)
    expect_type(pop, "list")
    expect_true("tot_pop" %in% names(pop))

    survey <- simulate_survey(
        pop$ibex_pop,
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
    # Run a very small power analysis just to check mechanics
    power_res <- doas_power(
        n_sims = 2,
        surv_areas = c(2, 4),
        pop_size = 100,
        tot_areas = 5
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
})

test_that("doas with formula computes covariate models correctly", {
    # Simulate data with a distinct covariate
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
    expect_equal(res$method, "hybrid")
    expect_equal(res$num_do_sectors, 3)

    # Check that a GLM was fitted
    expect_s3_class(res$model, "glm")

    # Check that the formula is saved
    expect_equal(deparse(res$formula), "~effort")

    # Check if calculation completes successfully
    expect_true(res$total_abundance > 0)
    expect_true(res$se_abundance > 0)
})

test_that("doas computes bootstrap confidence intervals correctly", {
    # Simulate basic data
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
    # Run a small bootstrap for testing speed
    res <- doas(data, ci_method = "bootstrap", R = 10)

    expect_s3_class(res, "doas")
    expect_equal(res$ci_method, "bootstrap")

    # Check if calculation completes successfully and bounds make sense
    expect_true(res$total_abundance > 0)
    expect_true(res$se_abundance > 0)
    expect_true(res$ci_lower >= sum(data$total_count))
    expect_true(res$ci_upper >= res$ci_lower)
})

test_that("QAIC is calculated correctly for covariate models", {
    # Simulate data
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

    # Run two models
    mod1 <- doas(data, formula = ~effort)
    mod2 <- doas(data, formula = ~ effort + habitat)
    mod3 <- doas(data) # No formula

    # Test QAIC on single models
    qaic1 <- QAIC(mod1)
    qaic2 <- QAIC(mod2)

    expect_true(is.numeric(qaic1))
    expect_true(is.numeric(qaic2))

    # Test QAIC on model with no formula (should warn and return NA)
    expect_warning(qaic3 <- QAIC(mod3))
    expect_true(is.na(qaic3))

    # Test comparing two models
    comp <- QAIC(mod1, mod2)
    expect_s3_class(comp, "data.frame")
    expect_equal(nrow(comp), 2)
    expect_true("QAIC" %in% names(comp))
    expect_equal(sort(rownames(comp)), c("mod1", "mod2"))
})

test_that("doas computes dependent double observer estimates correctly", {
    # Simulate data where obs2 only saw 1 group missed by obs1 in each DO sector
    data <- data.frame(
        sector_id = 1:5,
        total_count = c(39, 28, 67, 20, 100),
        groups_both = c(10, 5, 12, NA, NA),
        groups_obs1 = c(2, 1, 2, NA, NA),
        groups_obs2 = c(1, 1, 1, NA, NA),
        mean_group_size = c(3, 4, 4.5, NA, NA),
        detected_groups = c(13, 7, 15, 5, 20)
    )

    res_dep <- doas(data, obs_type = "dependent")
    res_ind <- doas(data, obs_type = "independent")

    expect_s3_class(res_dep, "doas")

    # Dependent estimate should generally differ from independent for the same dataset
    expect_false(res_dep$total_abundance == res_ind$total_abundance)

    # The output log should not throw NA errors
    expect_true(res_dep$total_abundance > 0)
    expect_true(res_dep$se_abundance > 0)

    # Check bounds
    expect_true(res_dep$ci_lower >= sum(data$total_count))
    expect_true(res_dep$ci_upper >= res_dep$ci_lower)
})

test_that("doas computes bayesian double observer estimates correctly", {
    data <- data.frame(
        sector_id = 1:5,
        total_count = c(39, 28, 67, 20, 100),
        groups_both = c(10, 5, 12, NA, NA),
        groups_obs1 = c(2, 1, 2, NA, NA),
        groups_obs2 = c(1, 1, 1, NA, NA),
        mean_group_size = c(3, 4, 4.5, NA, NA),
        detected_groups = c(13, 7, 15, 5, 20)
    )

    set.seed(42)
    res_bayes <- doas(data, obs_type = "bayesian")

    expect_s3_class(res_bayes, "doas")

    # Check that estimates are generated
    expect_true(res_bayes$total_abundance > 0)
    expect_true(res_bayes$se_abundance > 0)

    # Compare to independent - they should evaluate differently
    res_ind <- doas(data, obs_type = "independent")
    expect_false(res_bayes$total_abundance == res_ind$total_abundance)

    # Check boundaries (Bayesian mean should remain mathematically sound above counts)
    expect_true(res_bayes$total_abundance >= sum(data$total_count))
})

