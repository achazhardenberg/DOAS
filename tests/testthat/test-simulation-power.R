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

test_that("simulate_pop supports habitat abundance covariate effects", {
    set.seed(12)
    pop_cov <- simulate_pop(pop_size = 1000, mean_group = 8, tot_areas = 10, n_cov_effect = 0.8)
    expect_true("habitat" %in% names(pop_cov))
    expect_equal(length(pop_cov$habitat), 10)

    survey_cov <- simulate_survey(
        pop_cov$sim_pop,
        det_p1 = 0.6,
        intra_group_det = 1.0,
        det_group = 0.9,
        obs_effect = 0.8,
        detvar = 0.05,
        surv_areas = 5,
        cov_effect = 0.5
    )
    expect_true("habitat" %in% names(survey_cov))
    expect_true("covariate" %in% names(survey_cov))
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

    # Test power analysis with simulated sector covariate and p_formula
    power_cov <- doas_power(
        n_sims = 2,
        surv_areas = 4,
        pop_size = 200,
        tot_areas = 6,
        cov_effect = 0.5,
        p_formula = ~covariate
    )
    expect_s3_class(power_cov, "data.frame")
    expect_equal(nrow(power_cov), 1)
})
