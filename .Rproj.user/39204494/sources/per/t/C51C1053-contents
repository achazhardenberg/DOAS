#' Simulate a Population
#'
#' Simulates a population of individuals split into groups across multiple sectors.
#'
#' @param pop_size Average total population size (Poisson distributed).
#' @param mean_group Average group size (Poisson distributed).
#' @param tot_areas Total number of survey sectors.
#' @param area_sizes Vector of probabilities/sizes for assigning groups to sectors. If NULL, uniform probabilities are used.
#'
#' @return A list containing the total population size and a list of groups per sector.
#' @importFrom stats rpois rnorm runif
#' @export
simulate_pop <- function(pop_size, mean_group, tot_areas, area_sizes = NULL) {
    tot_pop <- rpois(1, pop_size)
    if (tot_pop == 0) {
        return(list(tot_pop = 0, sim_pop = vector("list", tot_areas)))
    }

    if (is.null(area_sizes)) {
        area_sizes <- rep(1, tot_areas)
    }

    sim_list <- as.character(1:tot_pop)
    sim_areas <- split(
        sample(sim_list),
        sample(
            1:tot_areas,
            length(sim_list),
            replace = TRUE,
            prob = area_sizes
        )
    )

    sim_pop <- vector("list", tot_areas)
    for (i in 1:tot_areas) {
        area_pop <- unlist(sim_areas[[as.character(i)]] %||% character(0))
        if (length(area_pop) > 0) {
            k <- rpois(1, length(area_pop) / mean_group)
            if (k == 0) {
                k <- 1
            }
            groups <- split(
                sample(area_pop),
                sample(1:k, length(area_pop), replace = TRUE)
            )
            names(groups) <- paste0("g", 1:length(groups))
            sim_pop[[i]] <- groups
        } else {
            sim_pop[[i]] <- list()
        }
    }

    return(list(sim_pop = sim_pop, tot_pop = tot_pop))
}

# Helper to provide %||%
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Simulate a Survey and format for DOAS
#'
#' @param sim_pop A simulated population from \code{\link{simulate_pop}}.
#' @param det_p1 Base detection probability for Observer 1.
#' @param intra_group_det Mean factor for count errors (e.g., 0.9 for undercounting).
#' @param det_group Base detection factor for group size scaling.
#' @param obs_effect Observer 2 relative efficiency (e.g., 0.75 means Obs 2 is 75% as effective).
#' @param detvar Standard deviation for detection probability across sectors.
#' @param surv_areas Number of sectors to receive Double Observer surveys.
#'
#' @return A list containing the total population size and a list of groups per sector.
#' @importFrom truncnorm rtruncnorm
#' @importFrom scales rescale
#' @export
simulate_survey <- function(
    sim_pop,
    det_p1,
    intra_group_det,
    det_group,
    obs_effect,
    detvar,
    surv_areas
) {
    tot_areas <- length(sim_pop)

    if (surv_areas > tot_areas) {
        surv_areas <- tot_areas
    }
    do_sectors <- sample(1:tot_areas, surv_areas, replace = FALSE)

    survey_data <- data.frame(
        sector_id = 1:tot_areas,
        total_count = 0,
        groups_both = NA_integer_,
        groups_obs1 = NA_integer_,
        groups_obs2 = NA_integer_,
        detected_groups = 0,
        stringsAsFactors = FALSE
    )

    for (m in 1:tot_areas) {
        sim_groups <- sim_pop[[m]]
        if (length(sim_groups) > 0) {
            group_sizes <- sapply(sim_groups, length)
            grp_sqrt <- sqrt(group_sizes)
            grp_sqrt <- ifelse(grp_sqrt >= sqrt(15), sqrt(15), grp_sqrt)

            # Group detectabilities based on size
            if (length(grp_sqrt) > 1) {
                group_det <- rescale(grp_sqrt, to = c(det_group, 1))
            } else {
                group_det <- det_group
            }

            detp_factor <- rtruncnorm(
                1,
                a = 0,
                b = 1,
                mean = det_p1,
                sd = detvar
            )

            # Observer 1 detection
            p1 <- detp_factor * group_det
            obs1_det <- runif(length(p1)) < p1

            # For DO sectors, calculate Observer 2
            if (m %in% do_sectors) {
                p2 <- p1 * obs_effect
                obs2_det <- runif(length(p2)) < p2

                # Determine intersections
                idx_both <- which(obs1_det & obs2_det)
                idx_obs1 <- which(obs1_det & !obs2_det)
                idx_obs2 <- which(!obs1_det & obs2_det)

                survey_data$groups_both[m] <- length(idx_both)
                survey_data$groups_obs1[m] <- length(idx_obs1)
                survey_data$groups_obs2[m] <- length(idx_obs2)

                # Find detected groups
                idx_any <- which(obs1_det | obs2_det)

                # Total groups detected
                survey_data$detected_groups[m] <- length(idx_any)

                if (length(idx_any) > 0) {
                    err_factor <- rnorm(1, intra_group_det, 0.05)
                    survey_data$total_count[m] <- round(
                        (sum(group_sizes[idx_any]) * err_factor),
                        0
                    )
                }
            } else {
                # Regular block survey (Observer 1 only)
                idx_obs1 <- which(obs1_det)
                survey_data$detected_groups[m] <- length(idx_obs1)

                if (length(idx_obs1) > 0) {
                    err_factor <- rnorm(1, intra_group_det, 0.05)
                    survey_data$total_count[m] <- round(
                        (sum(group_sizes[idx_obs1]) * err_factor),
                        0
                    )
                }
            }
        }
    }

    return(survey_data)
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
#' @param method Estimation method ("hybrid" or "global").
#' @param ci_method Confidence interval method ("lognormal" or "bootstrap").
#' @param obs_type Observer geometry ("independent", "dependent", or "bayesian").
#'
#' @return A data frame summarizing performance metrics for each scenario in \code{surv_areas}:
#' \describe{
#'   \item{\code{surv_areas}}{Number of DO sectors surveyed}
#'   \item{\code{mean_true_N}}{Average simulated actual population size}
#'   \item{\code{mean_est_N}}{Average estimated population size}
#'   \item{\code{rel_bias}}{Relative bias ((mean_est / mean_true) - 1)}
#'   \item{\code{rmse}}{Root Mean Squared Error}
#'   \item{\code{coverage}}{Proportion of simulations where true pop fell within the 95\% CI}
#'   \item{\code{mean_cv}}{Mean Coefficient of Variation (SE / Estimate)}
#' }
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
