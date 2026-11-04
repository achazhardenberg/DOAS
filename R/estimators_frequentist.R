#' Calculate Chapman-adjusted Lincoln-Petersen estimate
#'
#' @param both Integer vector of numbers of groups seen by both observers
#' @param obs1_only Integer vector of numbers of groups seen by Observer 1 only
#' @param obs2_only Integer vector of numbers of groups seen by Observer 2 only
#'
#' @return A numeric vector of the estimated total number of groups (\eqn{G_{est}})
#' @export
chapman_estimate <- function(both, obs1_only, obs2_only) {
    ((both + obs1_only + 1) * (both + obs2_only + 1) / (both + 1)) - 1
}

#' Calculate Variance of the Chapman estimator
#'
#' @param both Integer vector of numbers of groups seen by both observers
#' @param obs1_only Integer vector of numbers of groups seen by Observer 1 only
#' @param obs2_only Integer vector of numbers of groups seen by Observer 2 only
#'
#' @return A numeric vector describing the variance of the estimated number of groups
#' @export
chapman_variance <- function(both, obs1_only, obs2_only) {
    num <- (both + obs1_only + 1) *
        (both + obs2_only + 1) *
        obs1_only *
        obs2_only
    den <- ((both + 1)^2) * (both + 2)
    num / den
}

#' Calculate Dependent Double Observer (Zippin) Estimate
#'
#' @param both Integer vector of numbers of groups seen by both observers
#' @param obs1_only Integer vector of numbers of groups seen by Observer 1 (Primary) only
#' @param obs2_only Integer vector of numbers of groups seen by Observer 2 (Secondary) only
#'
#' @return A numeric vector of the estimated total number of groups (\eqn{G_{est}})
#' @export
dependent_zippin_estimate <- function(both, obs1_only, obs2_only) {
    x1 <- both + obs1_only
    x2 <- obs2_only
    ifelse(x1 > x2, (x1^2) / (x1 - x2), x1 + x2)
}

#' Calculate Variance of the Dependent Double Observer (Zippin) estimator
#'
#' @param both Integer vector of numbers of groups seen by both observers
#' @param obs1_only Integer vector of numbers of groups seen by Observer 1 (Primary) only
#' @param obs2_only Integer vector of numbers of groups seen by Observer 2 (Secondary) only
#'
#' @return A numeric vector describing the variance of the estimated number of groups
#' @export
dependent_zippin_variance <- function(both, obs1_only, obs2_only) {
    x1 <- both + obs1_only
    x2 <- obs2_only
    ifelse(x1 > x2, (x1^2 * x2^2 * (x1 + x2)) / ((x1 - x2)^4), NA_real_)
}

#' Compute Chao (1987) Log-Normal Confidence Interval Bounds
#'
#' @param est Numeric point estimate of abundance.
#' @param var_est Numeric variance of the abundance estimate.
#' @param min_count Numeric minimum observed count (lower bound floor).
#' @return Numeric vector of length 2: `c(ci_lower, ci_upper)`.
#' @keywords internal
#' @noRd
compute_lognormal_ci <- function(est, var_est, min_count) {
    se_est <- sqrt(var_est)
    f0 <- est - min_count
    if (is.na(se_est) || is.nan(se_est) || f0 <= 0 || se_est == 0) {
        return(c(min_count, est))
    }
    c_term <- exp(1.96 * sqrt(log(1 + (var_est / (f0^2)))))
    c(min_count + (f0 / c_term), min_count + (f0 * c_term))
}

#' Compute Frequentist DO Sector Group Estimates and Variances
#'
#' @param do_data Data frame of DO sectors.
#' @param obs_type Character: `"independent"` or `"dependent"`.
#' @return Updated `do_data` with `G_est`, `G_var`, and `p_hat`.
#' @keywords internal
#' @noRd
compute_frequentist_do <- function(do_data, obs_type) {
    if (obs_type == "independent") {
        do_data$G_est <- chapman_estimate(
            both = do_data$groups_both,
            obs1_only = do_data$groups_obs1,
            obs2_only = do_data$groups_obs2
        )
        do_data$G_var <- chapman_variance(
            both = do_data$groups_both,
            obs1_only = do_data$groups_obs1,
            obs2_only = do_data$groups_obs2
        )
    } else {
        do_data$G_est <- dependent_zippin_estimate(
            both = do_data$groups_both,
            obs1_only = do_data$groups_obs1,
            obs2_only = do_data$groups_obs2
        )
        do_data$G_var <- dependent_zippin_variance(
            both = do_data$groups_both,
            obs1_only = do_data$groups_obs1,
            obs2_only = do_data$groups_obs2
        )
    }

    if (obs_type == "independent") {
        do_data$p_hat <- ifelse(
            do_data$G_est > 0,
            pmin(
                1,
                (2 * do_data$groups_both + do_data$groups_obs1 + do_data$groups_obs2) /
                    (2 * do_data$G_est)
            ),
            0
        )
    } else {
        do_data$p_hat <- ifelse(
            do_data$G_est > 0,
            pmin(1, (do_data$groups_both + do_data$groups_obs1) / do_data$G_est),
            0
        )
    }
    do_data
}

#' Compute Frequentist Pooled Detection Probability and Variance
#'
#' @param do_data Data frame of DO sectors with `p_hat` computed.
#' @param obs_type Character: `"independent"` or `"dependent"`.
#' @return List with `mean_p_hat` and `var_mean_p_hat`.
#' @importFrom stats var
#' @keywords internal
#' @noRd
compute_frequentist_global_p <- function(do_data, obs_type) {
    b_sum <- sum(do_data$groups_both, na.rm = TRUE)
    s1_sum <- sum(do_data$groups_obs1, na.rm = TRUE)
    s2_sum <- sum(do_data$groups_obs2, na.rm = TRUE)

    if (obs_type == "independent") {
        total_G <- chapman_estimate(b_sum, s1_sum, s2_sum)
        var_total_G <- chapman_variance(b_sum, s1_sum, s2_sum)
        mean_p_hat <- (2 * b_sum + s1_sum + s2_sum) / (2 * total_G)
    } else {
        total_G <- dependent_zippin_estimate(b_sum, s1_sum, s2_sum)
        var_total_G <- dependent_zippin_variance(b_sum, s1_sum, s2_sum)
        mean_p_hat <- (b_sum + s1_sum) / total_G
    }

    emp_var_p <- if (nrow(do_data) > 1) {
        stats::var(do_data$p_hat, na.rm = TRUE) / nrow(do_data)
    } else {
        0
    }
    if (is.na(emp_var_p)) {
        emp_var_p <- 0
    }

    pool_var_p <- if (!is.na(var_total_G) && total_G > 0) {
        (mean_p_hat^2) *
            (var_total_G / (total_G^2)) +
            (mean_p_hat * (1 - mean_p_hat)) / max(total_G, 1)
    } else {
        (mean_p_hat * (1 - mean_p_hat)) / max(total_G, 1)
    }

    var_mean_p_hat <- max(emp_var_p, pool_var_p, na.rm = TRUE)
    if (is.na(var_mean_p_hat)) {
        var_mean_p_hat <- 0
    }

    list(mean_p_hat = mean_p_hat, var_mean_p_hat = var_mean_p_hat)
}

#' Predict Sector-Specific Detection Probability from GLM
#'
#' @param model Fitted `glm` object.
#' @param row Single-row data frame for the sector.
#' @param fallback_var Numeric fallback variance if GLM standard error is `NaN`.
#' @return List with `p_hat` and `var_p_hat`.
#' @importFrom stats predict
#' @keywords internal
#' @noRd
predict_sector_glm_p <- function(model, row, fallback_var) {
    pred <- stats::predict(model, newdata = row, type = "link", se.fit = TRUE)
    link_val <- pred$fit
    link_se <- pred$se.fit

    sector_p_hat <- exp(link_val) / (1 + exp(link_val))
    sector_var_p_hat <- (sector_p_hat * (1 - sector_p_hat))^2 * (link_se^2)

    if (is.na(sector_var_p_hat) || is.nan(sector_var_p_hat)) {
        sector_var_p_hat <- fallback_var
        warning(paste(
            "GLM produced NaN standard errors for sector",
            row$sector_id,
            "(likely due to too few DO sectors). Falling back to basic mean variance."
        ))
    }

    list(p_hat = sector_p_hat, var_p_hat = sector_var_p_hat)
}

#' Compute Frequentist Sector-Level and Total Abundance Estimates
#'
#' @param data Full survey data frame.
#' @param do_data DO sector subset data frame.
#' @param pooling Character: `"hybrid"` or `"global"`.
#' @param model Optional fitted `glm` object.
#' @param mean_p_hat Numeric global detection probability.
#' @param var_mean_p_hat Numeric variance of `mean_p_hat`.
#' @return List with `sector_estimates`, `total_abundance`, `total_variance`, `se_abundance`, `ci_lower`, `ci_upper`.
#' @importFrom stats delete.response terms model.matrix vcov
#' @keywords internal
#' @noRd
compute_frequentist_sectors <- function(
    data,
    do_data,
    obs_type,
    pooling,
    model,
    mean_p_hat,
    var_mean_p_hat
) {
    n_sec <- nrow(data)
    sector_estimates <- data.frame(
        sector_id = data$sector_id,
        total_count = data$total_count,
        est_groups = numeric(n_sec),
        est_abundance = numeric(n_sec),
        se_abundance = numeric(n_sec),
        ci_lower = numeric(n_sec),
        ci_upper = numeric(n_sec)
    )

    total_indep_variance <- 0
    shared_extrap_C_sum <- 0

    X_all <- NULL
    vcov_beta <- NULL
    grad_beta <- NULL
    has_valid_vcov <- FALSE

    if (!is.null(model)) {
        rhs_terms <- stats::delete.response(stats::terms(model))
        X_all <- tryCatch(
            stats::model.matrix(rhs_terms, data = data),
            error = function(e) NULL
        )
        vcov_beta <- tryCatch(
            stats::vcov(model),
            error = function(e) NULL
        )
        if (
            !is.null(X_all) &&
                !is.null(vcov_beta) &&
                !any(is.na(vcov_beta)) &&
                !any(is.nan(vcov_beta))
        ) {
            has_valid_vcov <- TRUE
            grad_beta <- numeric(ncol(X_all))
        }
    }

    for (i in seq_len(n_sec)) {
        row <- data[i, ]

        if (pooling == "hybrid" && row$survey_type == "DO") {
            match_idx <- which(do_data$sector_id == row$sector_id)[1]
            est_groups <- do_data$G_est[match_idx]
            var_groups <- do_data$G_var[match_idx]
            mean_size <- do_data$mean_group_size[match_idx]
            est_abundance <- est_groups * mean_size
            var_abundance <- var_groups * (mean_size^2)
            total_indep_variance <- total_indep_variance +
                ifelse(is.na(var_abundance), 0, var_abundance)
        } else {
            mean_size_block <- row$mean_group_size
            if (!is.null(model)) {
                pred_p <- predict_sector_glm_p(model, row, var_mean_p_hat)
                sector_p_hat <- pred_p$p_hat
                sector_var_p_hat <- pred_p$var_p_hat
            } else {
                sector_p_hat <- mean_p_hat
                sector_var_p_hat <- var_mean_p_hat
            }

            primary_groups <- if (row$survey_type == "DO") {
                if (obs_type == "independent") {
                    row$groups_both + 0.5 * (row$groups_obs1 + row$groups_obs2)
                } else {
                    row$groups_both + row$groups_obs1
                }
            } else {
                row$detected_groups
            }

            est_groups <- primary_groups / sector_p_hat
            est_abundance <- est_groups * mean_size_block

            var_groups_binom <- primary_groups *
                (1 - sector_p_hat) /
                (sector_p_hat^2)
            var_groups_p <- (primary_groups^2) *
                sector_var_p_hat /
                (sector_p_hat^4)
            var_groups <- var_groups_binom + var_groups_p
            var_abundance <- var_groups * (mean_size_block^2)

            if (!is.null(model)) {
                if (has_valid_vcov) {
                    total_indep_variance <- total_indep_variance +
                        (var_groups_binom * (mean_size_block^2))
                    c_i <- primary_groups * mean_size_block
                    grad_beta <- grad_beta -
                        (c_i * ((1 - sector_p_hat) / sector_p_hat)) * X_all[i, ]
                } else {
                    total_indep_variance <- total_indep_variance + var_abundance
                }
            } else {
                total_indep_variance <- total_indep_variance +
                    (var_groups_binom * (mean_size_block^2))
                shared_extrap_C_sum <- shared_extrap_C_sum +
                    (primary_groups * mean_size_block)
            }
        }

        bounds <- compute_lognormal_ci(
            est_abundance,
            var_abundance,
            row$total_count
        )
        sector_estimates$est_groups[i] <- est_groups
        sector_estimates$est_abundance[i] <- est_abundance
        sector_estimates$se_abundance[i] <- sqrt(var_abundance)
        sector_estimates$ci_lower[i] <- bounds[1]
        sector_estimates$ci_upper[i] <- bounds[2]
    }

    total_abundance <- sum(sector_estimates$est_abundance)
    shared_p_variance <- if (!is.null(model) && has_valid_vcov) {
        as.numeric(t(grad_beta) %*% vcov_beta %*% grad_beta)
    } else if (is.null(model) && shared_extrap_C_sum > 0) {
        (shared_extrap_C_sum^2) * var_mean_p_hat / (mean_p_hat^4)
    } else {
        0
    }
    total_variance <- total_indep_variance + shared_p_variance
    tot_bounds <- compute_lognormal_ci(
        total_abundance,
        total_variance,
        sum(data$total_count)
    )

    list(
        sector_estimates = sector_estimates,
        total_abundance = total_abundance,
        total_variance = total_variance,
        se_abundance = sqrt(total_variance),
        ci_lower = tot_bounds[1],
        ci_upper = tot_bounds[2]
    )
}
