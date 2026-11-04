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
    p_formula,
    p_betas,
    n_formula,
    n_betas,
    covariates,
    cov_effect,
    n_cov_effect,
    method,
    obs_type,
    pooling,
    ci_method
) {
    pop <- simulate_pop(
        pop_size = pop_size,
        mean_group = mean_group,
        tot_areas = tot_areas,
        area_sizes = area_sizes,
        n_formula = n_formula,
        n_betas = n_betas,
        covariates = covariates,
        n_cov_effect = n_cov_effect
    )

    survey_data <- simulate_survey(
        sim_pop = pop,
        det_p1 = det_p1,
        intra_group_det = intra_group_det,
        det_group = det_group,
        obs_effect = obs_effect,
        detvar = detvar,
        surv_areas = sa,
        p_formula = p_formula,
        p_betas = p_betas,
        cov_effect = cov_effect
    )

    est <- tryCatch(
        {
            doas(
                survey_data,
                method = method,
                obs_type = obs_type,
                pooling = pooling,
                p_formula = p_formula,
                n_formula = n_formula,
                ci_method = ci_method
            )
        },
        error = function(e) NULL
    )

    if (!is.null(est) && !is.na(est$total_abundance)) {
        c(
            true_N = as.numeric(pop$tot_pop),
            est_N = as.numeric(est$total_abundance),
            se_N = as.numeric(est$se_abundance),
            ci_lower = as.numeric(est$ci_lower),
            ci_upper = as.numeric(est$ci_upper)
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
#' Supports multi-covariate simulation and modeling on detection (\code{p_formula}, \code{p_betas})
#' and abundance (\code{n_formula}, \code{n_betas}).
#'
#' @param n_sims Number of simulations to run per scenario.
#' @param surv_areas Vector of Double Observer survey efforts to test (number of DO sectors).
#' @param pop_size Average total population size for \code{\link{simulate_pop}}.
#' @param mean_group Average group size for \code{\link{simulate_pop}}.
#' @param tot_areas Total number of survey sectors.
#' @param area_sizes Vector of probabilities/sizes for assigning groups to sectors. If NULL, uniform or covariate-driven probabilities are used.
#' @param det_p1 Base detection probability for Observer 1.
#' @param intra_group_det Mean factor for count errors.
#' @param det_group Base detection factor for group size scaling.
#' @param obs_effect Observer 2 relative efficiency.
#' @param detvar Standard deviation for detection probability across sectors.
#' @param p_formula Optional detection formula passed to simulation and \code{\link{doas}} (e.g., \code{~ effort + elevation}).
#' @param p_betas Optional named numeric vector of logit coefficients for \code{p_formula} in \code{\link{simulate_survey}}.
#' @param n_formula Optional abundance formula passed to simulation and \code{\link{doas}} (e.g., \code{~ elevation + forest}).
#' @param n_betas Optional named numeric vector of log-linear coefficients for \code{n_formula} in \code{\link{simulate_pop}}.
#' @param covariates Optional data frame of sector-level covariates.
#' @param cov_effect Numeric log-odds effect of a single simulated detection covariate (\code{covariate}) (retained for backward compatibility).
#' @param n_cov_effect Numeric log-linear effect of a single simulated abundance covariate (\code{habitat}) (retained for backward compatibility).
#' @param method Statistical inference paradigm (\code{"ML"}, \code{"Bayes"}, \code{"frequentist"}, or \code{"bayesian"}). Default is \code{"ML"}.
#' @param obs_type Observer field protocol (\code{"independent"} or \code{"dependent"}).
#' @param pooling Sector pooling strategy (\code{"hybrid"} or \code{"global"}).
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
    p_formula = NULL,
    p_betas = NULL,
    n_formula = NULL,
    n_betas = NULL,
    covariates = NULL,
    cov_effect = 0,
    n_cov_effect = 0,
    method = c("ML", "Bayes", "frequentist", "bayesian"),
    obs_type = c("independent", "dependent"),
    pooling = c("hybrid", "global"),
    ci_method = c("lognormal", "bootstrap")
) {
    method <- match.arg(method)
    if (method %in% c("frequentist", "ML")) {
        method <- "ML"
    } else {
        method <- "Bayes"
    }
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
                    p_formula = p_formula,
                    p_betas = p_betas,
                    n_formula = n_formula,
                    n_betas = n_betas,
                    covariates = covariates,
                    cov_effect = cov_effect,
                    n_cov_effect = n_cov_effect,
                    method = method,
                    obs_type = obs_type,
                    pooling = pooling,
                    ci_method = ci_method
                )
            },
            numeric(5)
        )

        valid_idx <- which(!is.na(sim_mat["est_N", ]))
        n_valid <- length(valid_idx)

        if (n_valid > 0) {
            true_N <- sim_mat["true_N", valid_idx]
            est_N <- sim_mat["est_N", valid_idx]
            se_N <- sim_mat["se_N", valid_idx]
            ci_l <- sim_mat["ci_lower", valid_idx]
            ci_u <- sim_mat["ci_upper", valid_idx]

            mean_true <- mean(true_N)
            mean_est <- mean(est_N)
            rel_bias <- (mean_est - mean_true) / mean_true
            rmse <- sqrt(mean((est_N - true_N)^2))
            coverage <- mean(true_N >= ci_l & true_N <= ci_u)
            mean_cv <- mean(se_N / est_N, na.rm = TRUE)
        } else {
            mean_true <- NA_real_
            mean_est <- NA_real_
            rel_bias <- NA_real_
            rmse <- NA_real_
            coverage <- NA_real_
            mean_cv <- NA_real_
        }

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
    }

    do.call(rbind, results_list)
}
