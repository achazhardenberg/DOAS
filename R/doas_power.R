#' Run a Single Power Simulation Replicate
#'
#' @keywords internal
#' @noRd
run_single_power_sim <- function(
    sa,
    pop_size,
    mean_group,
    tot_areas,
    area_sizes,
    det_p1,
    intra_group_det,
    det_group,
    obs_effect,
    detvar,
    cov_effect,
    method,
    obs_type,
    pooling,
    formula,
    ci_method
) {
    pop <- simulate_pop(
        pop_size = pop_size,
        mean_group = mean_group,
        tot_areas = tot_areas,
        area_sizes = area_sizes
    )

    survey_data <- simulate_survey(
        sim_pop = pop$sim_pop,
        det_p1 = det_p1,
        intra_group_det = intra_group_det,
        det_group = det_group,
        obs_effect = obs_effect,
        detvar = detvar,
        surv_areas = sa,
        cov_effect = cov_effect
    )

    est <- tryCatch(
        {
            doas(
                survey_data,
                method = method,
                obs_type = obs_type,
                pooling = pooling,
                formula = formula,
                ci_method = ci_method
            )
        },
        error = function(e) NULL
    )

    if (!is.null(est) && !is.na(est$total_abundance)) {
        c(
            true_N = pop$tot_pop,
            est_N = est$total_abundance,
            se_N = est$se_abundance,
            ci_lower = est$ci_lower,
            ci_upper = est$ci_upper
        )
    } else {
        c(
            true_N = NA_real_,
            est_N = NA_real_,
            se_N = NA_real_,
            ci_lower = NA_real_,
            ci_upper = NA_real_
        )
    }
}

#' DOAS Power Analysis
#'
#' Runs multiple simulations to evaluate the statistical performance (bias,
#' precision, coverage) of the DOAS estimator under varying sampling efforts.
#'
#' @param n_sims Number of simulations to run per scenario.
#' @param surv_areas Vector of Double Observer survey efforts to test (number of DO sectors).
#' @param pop_size Average total population size for \code{\link{simulate_pop}}.
#' @param mean_group Average group size for \code{\link{simulate_pop}}.
#' @param tot_areas Total number of survey sectors.
#' @param area_sizes Vector of probabilities/sizes for assigning groups to sectors. If NULL, uniform probabilities are used.
#' @param det_p1 Base detection probability for Observer 1.
#' @param intra_group_det Mean factor for count errors.
#' @param det_group Base detection factor for group size scaling.
#' @param obs_effect Observer 2 relative efficiency.
#' @param detvar Standard deviation for detection probability across sectors.
#' @param cov_effect Numeric log-odds effect of a simulated sector-level covariate (\code{covariate}) on detection probability (default \code{0}).
#' @param method Statistical inference paradigm (\code{"frequentist"} or \code{"bayesian"}).
#' @param obs_type Observer field protocol (\code{"independent"} or \code{"dependent"}).
#' @param pooling Sector pooling strategy (\code{"hybrid"} or \code{"global"}).
#' @param formula Optional formula (e.g., \code{~ covariate} or \code{~ mean_group_size}) passed to \code{\link{doas}}.
#' @param ci_method Confidence interval method (\code{"lognormal"} or \code{"bootstrap"}).
#'
#' @return A data frame summarizing performance metrics for each scenario in \code{surv_areas}:
#'   - \code{surv_areas}: Number of DO sectors surveyed
#'   - \code{successful_sims}: Number of simulations that returned valid estimates
#'   - \code{mean_true_N}: Average simulated actual population size
#'   - \code{mean_est_N}: Average estimated population size
#'   - \code{rel_bias}: Relative bias ((mean_est / mean_true) - 1)
#'   - \code{rmse}: Root Mean Squared Error
#'   - \code{coverage}: Proportion of simulations where true pop fell within the 95\% CI
#'   - \code{mean_cv}: Mean Coefficient of Variation (SE / Estimate)
#' @export
doas_power <- function(
    n_sims = 100,
    surv_areas = c(2, 4, 6),
    pop_size = 1000,
    mean_group = 10,
    tot_areas = 15,
    area_sizes = NULL,
    det_p1 = 0.6,
    intra_group_det = 0.9,
    det_group = 0.9,
    obs_effect = 0.75,
    detvar = 0.05,
    cov_effect = 0,
    method = c("frequentist", "bayesian"),
    obs_type = c("independent", "dependent"),
    pooling = c("hybrid", "global"),
    formula = NULL,
    ci_method = c("lognormal", "bootstrap")
) {
    method <- match.arg(method)
    obs_type <- match.arg(obs_type)
    pooling <- match.arg(pooling)
    ci_method <- match.arg(ci_method)

    results_list <- list()

    for (sa in surv_areas) {
        sim_mat <- vapply(
            seq_len(n_sims),
            function(i) {
                run_single_power_sim(
                    sa = sa,
                    pop_size = pop_size,
                    mean_group = mean_group,
                    tot_areas = tot_areas,
                    area_sizes = area_sizes,
                    det_p1 = det_p1,
                    intra_group_det = intra_group_det,
                    det_group = det_group,
                    obs_effect = obs_effect,
                    detvar = detvar,
                    cov_effect = cov_effect,
                    method = method,
                    obs_type = obs_type,
                    pooling = pooling,
                    formula = formula,
                    ci_method = ci_method
                )
            },
            numeric(5)
        )

        sim_results <- as.data.frame(t(sim_mat))
        valid_sims <- sim_results[!is.na(sim_results$est_N), , drop = FALSE]
        n_valid <- nrow(valid_sims)

        if (n_valid > 0) {
            mean_true <- mean(valid_sims$true_N)
            mean_est <- mean(valid_sims$est_N)
            rel_bias <- (mean_est - mean_true) / mean_true
            rmse <- sqrt(mean((valid_sims$est_N - valid_sims$true_N)^2))
            coverage <- mean(
                valid_sims$true_N >= valid_sims$ci_lower &
                    valid_sims$true_N <= valid_sims$ci_upper
            )
            mean_cv <- mean(valid_sims$se_N / valid_sims$est_N, na.rm = TRUE)

            results_list[[as.character(sa)]] <- data.frame(
                surv_areas = sa,
                successful_sims = n_valid,
                mean_true_N = mean_true,
                mean_est_N = mean_est,
                rel_bias = rel_bias,
                rmse = rmse,
                coverage = coverage,
                mean_cv = mean_cv
            )
        } else {
            warning(sprintf("All simulations failed for surv_areas = %d", sa))
        }
    }

    if (length(results_list) > 0) {
        final_results <- do.call(rbind, results_list)
        rownames(final_results) <- NULL
        return(final_results)
    } else {
        stop("All simulations failed.")
    }
}
