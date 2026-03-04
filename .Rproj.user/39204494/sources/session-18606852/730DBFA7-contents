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
#' @param method Estimation method ("hybrid" or "global").
#' @param ci_method Confidence interval method ("lognormal" or "bootstrap").
#' @param obs_type Observer geometry ("independent", "dependent", or "bayesian").
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
    det_group = 0.5,
    obs_effect = 0.75,
    detvar = 0.3,
    method = c("hybrid", "global"),
    ci_method = c("lognormal", "bootstrap"),
    obs_type = c("independent", "dependent", "bayesian")
) {
    method <- match.arg(method)
    ci_method <- match.arg(ci_method)
    obs_type <- match.arg(obs_type)

    results_list <- list()

    for (sa in surv_areas) {
        sim_results <- data.frame(
            sim_id = 1:n_sims,
            true_N = NA,
            est_N = NA,
            se_N = NA,
            ci_lower = NA,
            ci_upper = NA
        )

        for (i in seq_len(n_sims)) {
            # Simulate population
            pop <- simulate_pop(
                pop_size = pop_size,
                mean_group = mean_group,
                tot_areas = tot_areas,
                area_sizes = area_sizes
            )

            # Simulate survey for this number of DO areas
            survey_data <- simulate_survey(
                sim_pop = pop$sim_pop,
                det_p1 = det_p1,
                intra_group_det = intra_group_det,
                det_group = det_group,
                obs_effect = obs_effect,
                detvar = detvar,
                surv_areas = sa
            )

            # Estimate abundance
            est <- tryCatch(
                {
                    doas(
                        survey_data,
                        method = method,
                        ci_method = ci_method,
                        obs_type = obs_type
                    )
                },
                error = function(e) NULL
            )

            if (!is.null(est) && !is.na(est$total_abundance)) {
                sim_results$true_N[i] <- pop$tot_pop
                sim_results$est_N[i] <- est$total_abundance
                sim_results$se_N[i] <- est$se_abundance
                sim_results$ci_lower[i] <- est$ci_lower
                sim_results$ci_upper[i] <- est$ci_upper
            }
        }

        # Remove failed simulations (e.g. 0 DO counts)
        valid_sims <- sim_results[!is.na(sim_results$est_N), ]
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
