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

test_that("dependent_zippin_estimate and variance calculate correctly", {
    both <- 10
    obs1 <- 2
    obs2 <- 3

    # Zippin estimate
    est <- dependent_zippin_estimate(both, obs1, obs2)
    var_est <- dependent_zippin_variance(both, obs1, obs2)

    expect_true(is.numeric(est))
    expect_true(is.numeric(var_est))
    expect_true(est >= both + obs1 + obs2)
    expect_true(var_est > 0)
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
