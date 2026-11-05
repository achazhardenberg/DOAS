test_that("posterior_predict generates valid yrep matrices", {
  set.seed(123)
  toy_data <- data.frame(
    sector_id = 1:8,
    survey_type = c(rep("DO", 5), rep("SO", 3)),
    total_count = c(6, 2, 0, 7, 3, 2, 1, 0),
    detected_groups = c(6, 2, 0, 7, 3, 2, 1, 0),
    groups_both = c(3, 1, 0, 3, 1, NA, NA, NA),
    groups_obs1 = c(2, 1, 0, 2, 1, NA, NA, NA),
    groups_obs2 = c(1, 0, 0, 2, 1, NA, NA, NA),
    mean_group_size = 1
  )

  fit_bayes <- doas(toy_data, method = "Bayes", obs_type = "independent", iter = 1000, burnin = 200)

  yrep <- posterior_predict(fit_bayes, n_draws = 50)
  expect_s3_class(yrep, "doas_yrep")
  expect_true(is.matrix(yrep))
  expect_equal(dim(yrep), c(50, 8))
  expect_equal(length(attr(yrep, "y")), 8)
  expect_equal(as.numeric(attr(yrep, "y")), toy_data$total_count)

  # Test target types
  yrep_both <- posterior_predict(fit_bayes, n_draws = 20, type = "both")
  expect_equal(dim(yrep_both), c(20, 8))
  expect_equal(attr(yrep_both, "target"), "both")

  yrep_obs1 <- posterior_predict(fit_bayes, n_draws = 20, type = "obs1")
  expect_equal(dim(yrep_obs1), c(20, 8))

  yrep_obs2 <- posterior_predict(fit_bayes, n_draws = 20, type = "obs2")
  expect_equal(dim(yrep_obs2), c(20, 8))

  # Test dependent observers
  fit_dep <- doas(toy_data, method = "Bayes", obs_type = "dependent", iter = 1000, burnin = 200)
  yrep_dep <- posterior_predict(fit_dep, n_draws = 30)
  expect_equal(dim(yrep_dep), c(30, 8))
})

test_that("posterior_predict drops unsurveyed sectors correctly", {
  set.seed(42)
  toy_data <- data.frame(
    sector_id = 1:5,
    survey_type = c("DO", "DO", "SO", "Unsurveyed", "Unsurveyed"),
    total_count = c(4, 2, 1, NA, NA),
    detected_groups = c(4, 2, 1, NA, NA),
    groups_both = c(2, 1, NA, NA, NA),
    groups_obs1 = c(1, 1, NA, NA, NA),
    groups_obs2 = c(1, 0, NA, NA, NA),
    habitat = c(1.2, 0.8, 1.5, 0.9, 1.1),
    mean_group_size = 1
  )

  fit_unsurv <- doas(
    toy_data,
    method = "Bayes",
    obs_type = "independent",
    n_formula = ~ habitat,
    iter = 1000,
    burnin = 200
  )

  # Default drop_unsurveyed = TRUE
  yrep <- posterior_predict(fit_unsurv, n_draws = 20)
  expect_equal(ncol(yrep), 3)
  expect_equal(length(attr(yrep, "y")), 3)
  expect_false(any(is.na(yrep)))

  # drop_unsurveyed = FALSE
  yrep_all <- posterior_predict(fit_unsurv, n_draws = 20, drop_unsurveyed = FALSE)
  expect_equal(ncol(yrep_all), 5)
  expect_true(all(is.na(yrep_all[, 4:5])))
})

test_that("pp_check produces ggplot objects when bayesplot is available", {
  skip_if_not_installed("bayesplot")

  set.seed(99)
  toy_data <- data.frame(
    sector_id = 1:6,
    survey_type = rep("DO", 6),
    total_count = c(5, 3, 0, 8, 2, 0),
    detected_groups = c(5, 3, 0, 8, 2, 0),
    groups_both = c(2, 1, 0, 4, 1, 0),
    groups_obs1 = c(2, 1, 0, 2, 1, 0),
    groups_obs2 = c(1, 1, 0, 2, 0, 0),
    mean_group_size = 1
  )

  fit_bayes <- doas(toy_data, method = "Bayes", obs_type = "independent", iter = 1000, burnin = 200)

  p_dens <- pp_check(fit_bayes, type = "dens_overlay", n_draws = 15)
  expect_s3_class(p_dens, "ggplot")

  p_stat <- pp_check(fit_bayes, type = "stat", stat = "sd")
  expect_s3_class(p_stat, "ggplot")

  p_bars <- pp_check(fit_bayes, type = "bars")
  expect_s3_class(p_bars, "ggplot")
})

test_that("as.array and as.matrix format MCMC draws correctly", {
  set.seed(101)
  toy_data <- data.frame(
    sector_id = 1:6,
    survey_type = rep("DO", 6),
    total_count = c(5, 3, 0, 8, 2, 0),
    detected_groups = c(5, 3, 0, 8, 2, 0),
    groups_both = c(2, 1, 0, 4, 1, 0),
    groups_obs1 = c(2, 1, 0, 2, 1, 0),
    groups_obs2 = c(1, 1, 0, 2, 0, 0),
    cov1 = rnorm(6),
    mean_group_size = 1
  )

  fit_bayes <- doas(
    toy_data,
    method = "Bayes",
    obs_type = "independent",
    p_formula = ~ cov1,
    iter = 1000,
    burnin = 200,
    chains = 2
  )

  arr <- as.array(fit_bayes)
  expect_equal(length(dim(arr)), 3)
  expect_equal(dim(arr)[1], 800)  # iter - burnin
  expect_equal(dim(arr)[2], 2)    # chains
  expect_true("total_abundance" %in% dimnames(arr)$parameters)
  expect_true("p_hat" %in% dimnames(arr)$parameters)

  mat <- as.matrix(fit_bayes)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 1600)  # 800 * 2
  expect_equal(ncol(mat), dim(arr)[3])

  # Test compatibility with bayesplot mcmc_trace
  if (requireNamespace("bayesplot", quietly = TRUE)) {
    p_trace <- bayesplot::mcmc_trace(arr)
    expect_s3_class(p_trace, "ggplot")
  }
})

test_that("PPC and array methods reject non-Bayesian models with informative errors", {
  set.seed(42)
  toy_data <- data.frame(
    sector_id = 1:4,
    survey_type = rep("DO", 4),
    total_count = c(5, 3, 2, 4),
    detected_groups = c(5, 3, 2, 4),
    groups_both = c(2, 1, 1, 2),
    groups_obs1 = c(2, 1, 0, 1),
    groups_obs2 = c(1, 1, 1, 1),
    mean_group_size = 1
  )

  fit_ml <- doas(toy_data, method = "ML", obs_type = "independent")

  expect_error(posterior_predict(fit_ml), "only supported for Bayesian models")
  expect_error(pp_check(fit_ml), "only supported for Bayesian models")
  expect_error(as.array(fit_ml), "only supported for Bayesian models")
  expect_error(as.matrix(fit_ml), "only supported for Bayesian models")
})
