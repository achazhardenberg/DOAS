#' Calculate Chapman-adjusted Lincoln-Petersen estimate
#'
#' @param both Integer vector of numbers of groups seen by both observers
#' @param obs1_only Integer vector of numbers of groups seen by Observer 1 only
#' @param obs2_only Integer vector of numbers of groups seen by Observer 2 only
#'
#' @return A numeric vector of the estimated total number of groups ($G_{est}$)
#' @export
chapman_estimate <- function(both, obs1_only, obs2_only) {
    # Formula: G_est = (((B + S1 + 1) * (B + S2 + 1)) / (B + 1)) - 1
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
#' @return A numeric vector of the estimated total number of groups ($G_{est}$)
#' @export
dependent_zippin_estimate <- function(both, obs1_only, obs2_only) {
    x1 <- both + obs1_only # Total seen by Primary Observer
    x2 <- obs2_only # Total seen by Secondary Observer (missed by Primary)

    # Zippin 2-pass estimate: N = x1^2 / (x1 - x2)
    # If x1 <= x2, the model fails (secondary saw more than primary), fallback to total unique groups (x1 + x2)
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

    # Var(N) = (x1^2 * x2^2 * (x1 + x2)) / (x1 - x2)^4
    # If x1 <= x2, variance is undefined via this method, fallback to NA (will be handled downstream)
    ifelse(x1 > x2, (x1^2 * x2^2 * (x1 + x2)) / ((x1 - x2)^4), NA)
}

#' Calculate Bayesian Conjugate Double Observer Estimate
#'
#' Uses Monte Carlo simulation from Beta posterior distributions (assuming Beta(1,1) priors)
#' to estimate the mean group abundance and its variance.
#'
#' @param both Integer vector of numbers of groups seen by both observers
#' @param obs1_only Integer vector of numbers of groups seen by Observer 1 only
#' @param obs2_only Integer vector of numbers of groups seen by Observer 2 only
#' @param iter Integer number of Monte Carlo iterations (default 10000)
#'
#' @return A list containing the posterior 'mean' estimate vector and 'var' variance vector.
#' @export
bayesian_conjugate_estimate <- function(
    both,
    obs1_only,
    obs2_only,
    iter = 10000
) {
    n_sectors <- length(both)
    est_mean <- numeric(n_sectors)
    est_var <- numeric(n_sectors)

    for (i in 1:n_sectors) {
        b <- both[i]
        s1 <- obs1_only[i]
        s2 <- obs2_only[i]

        # If no animals were seen at all, abundance is 0 with 0 variance
        if (is.na(b) || (b == 0 && s1 == 0 && s2 == 0)) {
            est_mean[i] <- 0
            est_var[i] <- 0
            next
        }

        # Posterior Beta parameters (assuming Beta(1,1) prior)
        # p1: successes = b + s1, failures = s2
        alpha_1 <- b + s1 + 1
        beta_1 <- s2 + 1

        # p2: successes = b + s2, failures = s1
        alpha_2 <- b + s2 + 1
        beta_2 <- s1 + 1

        # Simulate posteriors
        p1_sim <- stats::rbeta(iter, alpha_1, beta_1)
        p2_sim <- stats::rbeta(iter, alpha_2, beta_2)

        # Overall detection probability (at least one observer)
        p_overall_sim <- p1_sim + p2_sim - (p1_sim * p2_sim)

        # Total distinct groups detected
        detected <- b + s1 + s2

        # Posterior abundance estimates
        N_sim <- detected / p_overall_sim

        est_mean[i] <- mean(N_sim)
        est_var[i] <- stats::var(N_sim)
    }

    list(mean = est_mean, var = est_var)
}

#' Double Observer Adjusted Survey (DOAS) Abundance Estimation
#'
#' Estimates total population abundance using a combination of block counts
#' in all sectors and Double Observer (DO) surveys in a subset of sectors.
#'
#' @param data A dataframe containing survey data. Each row should represent a sector.
#'   Required columns:
#'   - \code{sector_id}: Unique identifier for the sector.
#'   - \code{total_count}: Total individuals counted across all groups in the sector.
#'   - \code{detected_groups}: Total groups detected in the sector (for Block sectors this is the final group count; for DO sectors it can be B + S1 + S2 or supplied).
#'
#'   Optional/Conditional columns for Double Observer (DO) sectors:
#'   - \code{groups_both}: Number of groups seen by both observers.
#'   - \code{groups_obs1}: Number of groups seen by Observer 1 only.
#'   - \code{groups_obs2}: Number of groups seen by Observer 2 only.
#'   If the groups_both, groups_obs1, and groups_obs2 columns are present and non-NA for a row, it is treated as a DO survey.
#'   Otherwise it defaults to a standard Block count.
#'
#' @param method Character string specifying the estimation method.
#'   "hybrid" (default) uses the specific DO estimate for DO sectors and adjusts only Block sectors.
#'   "global" adjusts all sectors using the average detection probability.
#' @param formula Optional formula (e.g., \code{~ covariate1 + covariate2}) to model detection probability
#'   \code{\\hat{p}} across sectors. If provided, fits a quasi-binomial GLM on DO sectors and predicts
#'   specific detection probabilities for each Block sector.
#' @param ci_method Character string specifying how to calculate Confidence Intervals.
#'   "lognormal" (default) uses the analytic Delta Method variance and Chao (1987) log-normal approximation.
#'   "bootstrap" uses a non-parametric block resampling bootstrap.
#' @param R Integer indicating the number of bootstrap replicates (only used if \code{ci_method = "bootstrap"}). Default is 1000.
#' @param obs_type Character string specifying the observer relationship.
#'   "independent" (default) uses the Chapman estimator assuming both observers scanned independently.
#'   "dependent" assumes Observer 1 points out all animals, and Observer 2 only records what Observer 1 missed, reducing to a Zippin Removal estimator.
#'   "bayesian" uses Monte Carlo simulation from Conjugate Beta-Binomial posteriors to naturally bound credible estimates avoiding normality assumptions.
#'
#' @return A list of class `doas` containing the estimate and intermediate calculations.
#' @export
doas <- function(
    data,
    method = c("hybrid", "global"),
    formula = NULL,
    ci_method = c("lognormal", "bootstrap"),
    R = 1000,
    obs_type = c("independent", "dependent", "bayesian")
) {
    method <- match.arg(method)
    ci_method <- match.arg(ci_method)
    obs_type <- match.arg(obs_type)

    # Validate required columns
    req_cols <- c("sector_id", "total_count", "detected_groups")
    do_cols <- c("groups_both", "groups_obs1", "groups_obs2")
    if (!all(req_cols %in% names(data))) {
        stop("Missing required basic columns in data.")
    }

    # Calculate mean_group_size internally
    data$mean_group_size <- ifelse(
        data$detected_groups > 0,
        data$total_count / data$detected_groups,
        0
    )

    # Auto-detect survey type
    if (all(do_cols %in% names(data))) {
        data$survey_type <- ifelse(
            !is.na(data$groups_both) &
                !is.na(data$groups_obs1) &
                !is.na(data$groups_obs2),
            "DO",
            "Block"
        )
    } else {
        data$survey_type <- "Block"
        # Pad missing DO columns with NA so the rest of the code works
        for (col in do_cols) {
            data[[col]] <- NA
        }
    }

    do_data <- data[data$survey_type == "DO", ]
    if (nrow(do_data) == 0) {
        stop(
            "No 'DO' sectors found in data. Ensure you have provided the groups_both, groups_obs1, and groups_obs2 columns with non-NA values for at least one sector."
        )
    }

    # Calculate Estimates and Detection Probabilities for DO sectors
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
    } else if (obs_type == "dependent") {
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
    } else if (obs_type == "bayesian") {
        bayes_res <- bayesian_conjugate_estimate(
            both = do_data$groups_both,
            obs1_only = do_data$groups_obs1,
            obs2_only = do_data$groups_obs2,
            iter = 10000
        )
        do_data$G_est <- bayes_res$mean
        do_data$G_var <- bayes_res$var
    }

    # Group detection probability for the PRIMARY observer in each DO sector
    do_data$p_hat <- (do_data$groups_both + do_data$groups_obs1) / do_data$G_est

    # Optional Modeling of p_hat
    model <- NULL
    if (!is.null(formula)) {
        if (!inherits(formula, "formula")) {
            stop("formula must be a valid R formula (e.g. ~ covariate)")
        }

        # Ensure 'p_hat' is the response variable
        model_formula <- update(formula, p_hat ~ .)

        # Fit a quasi-binomial model (bounded [0, 1] but allows variance)
        # Weight by the number of estimated groups to give more weight to larger sectors
        model <- tryCatch(
            {
                eval(bquote(stats::glm(
                    .(model_formula),
                    family = stats::quasibinomial(link = "logit"),
                    data = do_data,
                    weights = G_est
                )))
            },
            error = function(e) {
                stop(paste(
                    "Failed to fit GLM for detection probabilities:",
                    e$message
                ))
            }
        )
    }
    # Calculate mean_p_hat globally to avoid small-sample bias in individual sector estimates
    if (obs_type == "independent") {
        total_G <- chapman_estimate(
            both = sum(do_data$groups_both, na.rm = TRUE),
            obs1_only = sum(do_data$groups_obs1, na.rm = TRUE),
            obs2_only = sum(do_data$groups_obs2, na.rm = TRUE)
        )
    } else if (obs_type == "dependent") {
        total_G <- dependent_zippin_estimate(
            both = sum(do_data$groups_both, na.rm = TRUE),
            obs1_only = sum(do_data$groups_obs1, na.rm = TRUE),
            obs2_only = sum(do_data$groups_obs2, na.rm = TRUE)
        )
    } else if (obs_type == "bayesian") {
        # For bayesian, we could pool the posteriors, but running one aggregate model is preferred
        bayes_res <- bayesian_conjugate_estimate(
            both = sum(do_data$groups_both, na.rm = TRUE),
            obs1_only = sum(do_data$groups_obs1, na.rm = TRUE),
            obs2_only = sum(do_data$groups_obs2, na.rm = TRUE),
            iter = 10000
        )
        total_G <- bayes_res$mean
    }

    mean_p_hat <- sum(do_data$groups_both + do_data$groups_obs1, na.rm = TRUE) /
        total_G

    if (mean_p_hat <= 0 || is.na(mean_p_hat) || mean_p_hat > 1) {
        stop("Calculated mean detection probability is invalid.")
    }

    # Variance of mean detection probability
    var_mean_p_hat <- stats::var(do_data$p_hat, na.rm = TRUE) / nrow(do_data)
    if (is.na(var_mean_p_hat)) {
        var_mean_p_hat <- 0
    }

    results <- list(
        method = method,
        formula = formula,
        model = model,
        mean_p_hat = mean_p_hat,
        var_mean_p_hat = var_mean_p_hat,
        num_do_sectors = nrow(do_data),
        num_total_sectors = nrow(data),
        obs_type = obs_type,
        sector_estimates = data.frame(
            sector_id = data$sector_id,
            total_count = data$total_count,
            est_groups = NA,
            est_abundance = NA,
            se_abundance = NA,
            ci_lower = NA,
            ci_upper = NA
        )
    )

    total_abundance <- 0
    total_variance <- 0
    total_seen <- sum(data$total_count)

    for (i in seq_len(nrow(data))) {
        row <- data[i, ]

        if (method == "hybrid" && row$survey_type == "DO") {
            # Hybrid uses specific DO estimate for this sector
            # Force match to the first occurrence to handle duplicated rows generated during Bootstrap resampling
            match_idx <- which(do_data$sector_id == row$sector_id)[1]
            est_groups <- do_data$G_est[match_idx]
            var_groups <- do_data$G_var[match_idx]
            mean_size <- do_data$mean_group_size[match_idx]
            est_abundance <- est_groups * mean_size
            # For simplicity, assuming mean_size is a constant without variance
            var_abundance <- var_groups * (mean_size^2)
        } else {
            # Adjust block counts using mean detection probability
            # Using average group size for the sector as calculated by total_count / detected_groups
            mean_size_block <- ifelse(
                row$detected_groups > 0,
                row$total_count / row$detected_groups,
                0
            )

            if (!is.null(model)) {
                # Predict p_hat and its standard error on the logit scale
                pred <- stats::predict(
                    model,
                    newdata = row,
                    type = "link",
                    se.fit = TRUE
                )
                link_val <- pred$fit
                link_se <- pred$se.fit

                # Inverse link function (logit to probability)
                sector_p_hat <- exp(link_val) / (1 + exp(link_val))

                # Delta method for variance of p_hat: Var(g(X)) ~ [g'(X)]^2 * Var(X)
                # For exp(x)/(1+exp(x)), derivative is exp(x)/(1+exp(x))^2 = p * (1 - p)
                sector_var_p_hat <- (sector_p_hat * (1 - sector_p_hat))^2 *
                    (link_se^2)

                # Fallback if the GLM lacked degrees of freedom and returned NaN for standard error
                if (is.na(sector_var_p_hat) || is.nan(sector_var_p_hat)) {
                    sector_var_p_hat <- var_mean_p_hat
                    warning(paste(
                        "GLM produced NaN standard errors for sector",
                        row$sector_id,
                        "(likely due to too few DO sectors). Falling back to basic mean variance."
                    ))
                }
            } else {
                sector_p_hat <- mean_p_hat
                sector_var_p_hat <- var_mean_p_hat
            }

            if (row$survey_type == "DO") {
                primary_groups <- row$groups_both + row$groups_obs1
            } else {
                primary_groups <- row$detected_groups
            }

            est_groups <- primary_groups / sector_p_hat
            est_abundance <- est_groups * mean_size_block

            # Variance using Delta method: Var(C / p) ~ C^2 * Var(p) / p^4
            var_groups <- (primary_groups^2) *
                sector_var_p_hat /
                (sector_p_hat^4)
            var_abundance <- var_groups * (mean_size_block^2)
        }

        results$sector_estimates$est_groups[i] <- est_groups
        results$sector_estimates$est_abundance[i] <- est_abundance

        # Sector-level SE and Log-normal CI
        se_abund <- sqrt(var_abundance)
        results$sector_estimates$se_abundance[i] <- se_abund

        f0_sec <- est_abundance - row$total_count
        if (is.na(se_abund) || f0_sec <= 0 || se_abund == 0) {
            results$sector_estimates$ci_lower[i] <- row$total_count
            results$sector_estimates$ci_upper[i] <- est_abundance
        } else {
            C_term_sec <- exp(
                1.96 * sqrt(log(1 + (var_abundance / (f0_sec^2))))
            )
            results$sector_estimates$ci_lower[i] <- row$total_count +
                (f0_sec / C_term_sec)
            results$sector_estimates$ci_upper[i] <- row$total_count +
                (f0_sec * C_term_sec)
        }

        total_abundance <- total_abundance + est_abundance
        total_variance <- total_variance + var_abundance
    }

    results$total_abundance <- total_abundance
    results$total_variance <- total_variance
    results$se_abundance <- sqrt(total_variance)

    total_seen <- sum(data$total_count)
    results$ci_method <- ci_method

    if (ci_method == "lognormal") {
        # Log-normal Confidence Interval (Chao 1987) using analytic Delta Method variance
        f0 <- total_abundance - total_seen
        if (
            is.na(results$se_abundance) || f0 <= 0 || results$se_abundance == 0
        ) {
            results$ci_lower <- total_seen
            results$ci_upper <- total_abundance
        } else {
            C_term <- exp(
                1.96 * sqrt(log(1 + (results$total_variance / (f0^2))))
            )
            results$ci_lower <- total_seen + (f0 / C_term)
            results$ci_upper <- total_seen + (f0 * C_term)
        }
    } else if (ci_method == "bootstrap") {
        # Non-parametric Bootstrap Confidence Intervals
        boot_estimates <- numeric(R)

        do_indices <- which(
            !is.na(data$groups_both) &
                !is.na(data$groups_obs1) &
                !is.na(data$groups_obs2)
        )
        block_indices <- setdiff(seq_len(nrow(data)), do_indices)

        for (b in seq_len(R)) {
            # Resample sectors WITH replacement, keeping design structure
            boot_data_do <- data[
                sample(do_indices, length(do_indices), replace = TRUE),
                ,
                drop = FALSE
            ]
            boot_data_block <- data[
                sample(block_indices, length(block_indices), replace = TRUE),
                ,
                drop = FALSE
            ]
            boot_data <- rbind(boot_data_do, boot_data_block)

            # Re-run doas purely for the total point estimate
            boot_res <- suppressWarnings(tryCatch(
                {
                    doas(
                        boot_data,
                        method = method,
                        formula = formula,
                        ci_method = "lognormal",
                        obs_type = obs_type
                    )
                },
                error = function(e) NULL
            ))

            if (!is.null(boot_res)) {
                boot_estimates[b] <- boot_res$total_abundance
            } else {
                boot_estimates[b] <- NA
            }
        }

        valid_boots <- boot_estimates[!is.na(boot_estimates)]
        if (length(valid_boots) < R * 0.5) {
            warning(
                "More than 50% of bootstrap iterations failed. CI may be unreliable."
            )
        }

        # Override Standard Error and CI using empirical Bootstrap distribution
        results$se_abundance <- stats::sd(valid_boots)
        results$total_variance <- results$se_abundance^2

        # Use simple percentiles for CI. Bound lower by what was actually seen in the original single survey.
        raw_bounds <- stats::quantile(
            valid_boots,
            probs = c(0.025, 0.975),
            na.rm = TRUE
        )
        results$ci_lower <- max(total_seen, raw_bounds[1])
        results$ci_upper <- max(total_seen, raw_bounds[2])
    }

    class(results) <- "doas"
    return(results)
}

#' Summary for doas objects
#'
#' @param object An object of class `doas`
#' @param ... Additional arguments
#'
#' @export
summary.doas <- function(object, ...) {
    cat("\n  Double Observer Adjusted Survey (DOAS) Estimation  \n")
    cat("====================================================\n\n")
    cat(sprintf("Method used: %s\n", tools::toTitleCase(object$method)))
    if (!is.null(object$formula)) {
        cat(sprintf(
            "Covariate Formula: p_hat %s\n",
            paste(deparse(object$formula), collapse = " ")
        ))
    }
    cat(sprintf("Total sectors surveyed: %d\n", object$num_total_sectors))

    obs_label <- if (!is.null(object$obs_type)) {
        object$obs_type
    } else {
        "independent"
    }
    obs_label <- tools::toTitleCase(obs_label)
    cat(sprintf(
        "Double Observer (DO) sectors: %d (%s)\n",
        object$num_do_sectors,
        obs_label
    ))

    if (is.null(object$formula)) {
        cat(sprintf(
            "\nMean detection probability (p_hat): %.3f (Var = %.5f)\n",
            object$mean_p_hat,
            object$var_mean_p_hat
        ))
    } else {
        cat("\nDetection Probability Model (GLM):\n")
        print(stats::coef(summary(object$model)))
    }

    ci_label <- if (object$ci_method == "bootstrap") {
        "Bootstrap Percentile CI"
    } else {
        "Log-normal CI"
    }

    cat(sprintf("\nEstimated Total Abundance: %.1f\n", object$total_abundance))
    cat(sprintf("Standard Error: %.1f\n", object$se_abundance))
    cat(sprintf(
        "95%% %s: [%.1f - %.1f]\n",
        ci_label,
        object$ci_lower,
        object$ci_upper
    ))
    cat("----------------------------------------------------\n")
    cat("Sector-level Estimates:\n")
    # Only round numeric columns to prevent error with character sector_id
    num_cols <- sapply(object$sector_estimates, is.numeric)
    print_df <- object$sector_estimates
    print_df[, num_cols] <- round(print_df[, num_cols], 2)
    print(print_df, row.names = FALSE)
    cat("----------------------------------------------------\n")
}
