test_that("doas with p_formula computes covariate models correctly", {
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

    res <- doas(data, p_formula = ~effort)

    expect_s3_class(res, "doas")
    expect_equal(res$method, "frequentist")
    expect_equal(res$pooling, "hybrid")
    expect_equal(res$num_do_sectors, 3)
    expect_s3_class(res$model, "glm")
    expect_equal(deparse(res$p_formula), "~effort")
    expect_true(res$total_abundance > 0)
    expect_true(res$se_abundance > 0)
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

    mod1 <- doas(data, p_formula = ~effort)
    mod2 <- doas(data, p_formula = ~ effort + habitat)
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

test_that("multiple covariates, interactions, n_formula, and predict.doas() work", {
    set.seed(42)
    cov_df <- data.frame(
        sector_id = paste0("S", 1:10),
        total_count = c(60, 45, 80, 50, 70, 40, 90, 55, 65, 75),
        detected_groups = c(10, 8, 12, 9, 11, 7, 14, 9, 10, 11),
        groups_both = c(6, 4, 8, 5, 7, 3, NA, NA, NA, NA),
        groups_obs1 = c(2, 2, 2, 2, 2, 2, NA, NA, NA, NA),
        groups_obs2 = c(2, 2, 2, 2, 2, 2, NA, NA, NA, NA),
        effort = rnorm(10),
        ruggedness = rnorm(10),
        habitat = factor(rep(c("Forest", "Open"), 5))
    )
    d_cov <- doas_data(cov_df)

    # Categorical and interaction terms
    fit_cat <- doas(d_cov, method = "frequentist", p_formula = ~ effort + habitat, n_formula = ~ habitat)
    expect_s3_class(fit_cat$model, "glm")
    expect_true("habitatOpen" %in% names(stats::coef(fit_cat$model)))
    expect_true("habitatOpen" %in% names(stats::coef(fit_cat$n_model)))

    fit_f_int <- doas(d_cov, method = "frequentist", p_formula = ~ effort * ruggedness)
    expect_equal(length(stats::coef(fit_f_int$model)), 4)

    fit_b_int <- doas(
        d_cov,
        method = "bayesian",
        p_formula = ~ effort * ruggedness,
        chains = 2,
        iter = 400,
        burnin = 100
    )
    expect_equal(nrow(fit_b_int$bayes_coef), 4)
    expect_true("effort:ruggedness" %in% rownames(fit_b_int$bayes_coef))

    # Abundance covariate modeling (n_formula), unsurveyed sectors, and predict.doas()
    fit_fn <- doas(
        d_cov,
        method = "frequentist",
        p_formula = ~effort,
        n_formula = ~ruggedness
    )
    expect_s3_class(fit_fn$n_model, "glm")

    new_grid <- data.frame(
        effort = c(-1, 0, 1),
        ruggedness = c(-0.5, 0, 0.5)
    )
    pred_n <- predict(fit_fn, newdata = new_grid, type = "abundance")
    pred_p <- predict(fit_fn, newdata = new_grid, type = "detection")
    pred_g <- predict(fit_fn, newdata = new_grid, type = "groups")
    expect_equal(nrow(pred_n), 3)
    expect_true(all(c("fit", "se", "ci_lower", "ci_upper") %in% names(pred_n)))
    expect_true(all(pred_p$fit > 0 & pred_p$fit < 1))
    expect_true(all(pred_g$fit > 0))

    # Bayesian n_formula + predict()
    fit_bn <- doas(
        d_cov,
        method = "bayesian",
        p_formula = ~effort,
        n_formula = ~ruggedness,
        chains = 2,
        iter = 400,
        burnin = 100
    )
    expect_equal(nrow(fit_bn$bayes_n_coef), 2)
    pred_bn <- predict(fit_bn, newdata = new_grid, type = "abundance")
    expect_equal(nrow(pred_bn), 3)
    expect_true(all(pred_bn$fit > 0))

    # Imputing an unsurveyed sector (total_count = NA) via n_formula
    cov_unsurv <- cov_df
    cov_unsurv$total_count[10] <- NA
    cov_unsurv$detected_groups[10] <- NA
    d_unsurv <- doas_data(cov_unsurv)
    expect_equal(d_unsurv$survey_type[10], "Unsurveyed")
    expect_error(doas(d_unsurv), "n_formula")

    fit_unsurv <- doas(d_unsurv, n_formula = ~ruggedness)
    expect_true(fit_unsurv$sector_estimates$est_abundance[10] > 0)
})
