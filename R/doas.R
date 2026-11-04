#' Validate and Prepare Survey Data for DOAS
#'
#' Validates required columns, checks consistency between `detected_groups` and
#' `(groups_both, groups_obs1, groups_obs2)` in DO sectors, calculates `mean_group_size`,
#' and identifies `"DO"` versus `"Block"` survey sectors.
#'
#' @param data A sector-level data frame or `doas_data` object.
#' @return A list with `data` (prepared full data frame) and `do_data` (subset of DO sectors).
#' @keywords internal
#' @noRd
prepare_doas_data <- function(data) {
    if (!inherits(data, "doas_data")) {
        if (!all(c("sector_id", "total_count") %in% names(data))) {
            stop(
                "Missing required columns ('sector_id', 'total_count'). Use doas_data() to format and inspect tables with custom column names or long-format group/sighting records before calling doas()."
            )
        }
        cov_cols <- setdiff(
            names(data),
            c(
                "sector_id",
                "survey_type",
                "total_count",
                "detected_groups",
                "groups_both",
                "groups_obs1",
                "groups_obs2",
                "mean_group_size"
            )
        )
        data <- format_sector_wide(
            data = data,
            sector = "sector_id",
            count = "total_count",
            groups = "detected_groups",
            both = "groups_both",
            obs1 = "groups_obs1",
            obs2 = "groups_obs2",
            do_sectors = if ("survey_type" %in% names(data)) {
                data$sector_id[data$survey_type == "DO"]
            } else {
                NULL
            },
            covariates = if (length(cov_cols) > 0) cov_cols else NULL
        )
    }

    do_data <- data[data$survey_type == "DO", , drop = FALSE]
    if (nrow(do_data) == 0) {
        stop(
            "No 'DO' sectors found in data. Ensure you have provided Double Observer group counts for at least one sector."
        )
    }
    if (sum(do_data$groups_both + do_data$groups_obs1, na.rm = TRUE) <= 0) {
        stop(
            "Primary observer (Observer 1) detected 0 groups across all DO sectors; cannot estimate primary detection probability."
        )
    }

    list(data = data, do_data = do_data)
}

#' Fit Quasi-Binomial GLM for Primary Observer Detection Probability
#'
#' @param formula R formula for covariates.
#' @param do_data Data frame of DO sectors with `p_hat` and `G_est`.
#' @return A fitted `glm` object.
#' @importFrom stats glm quasibinomial update
#' @keywords internal
#' @noRd
fit_detection_glm <- function(formula, do_data) {
    if (!inherits(formula, "formula")) {
        stop("formula must be a valid R formula (e.g. ~ covariate)")
    }

    do_data_pos <- do_data[do_data$G_est > 0, , drop = FALSE]
    if (nrow(do_data_pos) < 2) {
        stop("At least 2 DO sectors with positive group detections are required to fit a covariate GLM.")
    }

    model_formula <- stats::update(formula, p_hat ~ .)
    tryCatch(
        {
            eval(bquote(stats::glm(
                .(model_formula),
                family = stats::quasibinomial(link = "logit"),
                data = do_data_pos,
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

#' Run Non-Parametric Stratified Bootstrap for DOAS Confidence Intervals
#'
#' @param data Prepared survey data frame.
#' @param method Character: `"frequentist"` or `"bayesian"`.
#' @param obs_type Character: `"independent"` or `"dependent"`.
#' @param pooling Character: `"hybrid"` or `"global"`.
#' @param formula Optional covariate formula.
#' @param R Integer number of bootstrap replicates.
#' @param chains Integer MCMC chains (reduced for bootstrap speed if Bayesian).
#' @param iter Integer MCMC iterations.
#' @param burnin Integer MCMC burn-in.
#' @return List with `se_abundance`, `total_variance`, `ci_lower`, `ci_upper`.
#' @importFrom stats sd quantile
#' @keywords internal
#' @noRd
run_doas_bootstrap <- function(
    data,
    method,
    obs_type,
    pooling,
    formula,
    R,
    chains,
    iter,
    burnin
) {
    do_indices <- which(data$survey_type == "DO")
    block_indices <- setdiff(seq_len(nrow(data)), do_indices)
    total_seen <- sum(data$total_count)

    boot_single <- function(b) {
        boot_do <- data[
            sample(do_indices, length(do_indices), replace = TRUE),
            ,
            drop = FALSE
        ]
        boot_block <- data[
            sample(block_indices, length(block_indices), replace = TRUE),
            ,
            drop = FALSE
        ]
        boot_data <- rbind(boot_do, boot_block)
        boot_data$sector_id <- seq_len(nrow(boot_data))
        class(boot_data) <- c("doas_data", "data.frame")

        boot_res <- suppressWarnings(tryCatch(
            {
                doas(
                    boot_data,
                    method = method,
                    obs_type = obs_type,
                    pooling = pooling,
                    formula = formula,
                    ci_method = "lognormal",
                    chains = if (method == "bayesian") 2 else chains,
                    iter = if (method == "bayesian") 600 else iter,
                    burnin = if (method == "bayesian") 150 else burnin
                )
            },
            error = function(e) NULL
        ))

        if (!is.null(boot_res)) boot_res$total_abundance else NA_real_
    }

    boot_estimates <- vapply(seq_len(R), boot_single, numeric(1))
    valid_boots <- boot_estimates[!is.na(boot_estimates)]

    if (length(valid_boots) < R * 0.5) {
        warning(
            "More than 50% of bootstrap iterations failed. CI may be unreliable."
        )
    }

    se_abund <- stats::sd(valid_boots)
    raw_bounds <- stats::quantile(
        valid_boots,
        probs = c(0.025, 0.975),
        na.rm = TRUE,
        names = FALSE
    )

    list(
        se_abundance = se_abund,
        total_variance = se_abund^2,
        ci_lower = max(total_seen, raw_bounds[1]),
        ci_upper = max(total_seen, raw_bounds[2])
    )
}

#' Double Observer Adjusted Survey (DOAS) Abundance Estimation
#'
#' Estimates total population abundance using a combination of block counts
#' in all sectors and Double Observer (DO) surveys in a subset of sectors.
#' Supports both Frequentist and multi-chain Bayesian MCMC inference across
#' independent (Lincoln-Petersen/Chapman) and dependent (2-pass removal/Zippin)
#' double-observer field protocols.
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
#'   If \code{groups_both}, \code{groups_obs1}, and \code{groups_obs2} are present and non-NA for a row, it is treated as a DO survey.
#'   Otherwise it defaults to a standard Block count.
#'
#' @param method Character string specifying the statistical inference paradigm:
#'   - \code{"frequentist"} (default): Uses analytic Chapman (when \code{obs_type = "independent"}) or Zippin removal (when \code{obs_type = "dependent"}) estimators and Delta-method variance.
#'   - \code{"bayesian"}: Uses a multi-chain Collapsed Poisson-Gamma-Beta Gibbs Sampler (MCMC) to sample the full joint posterior of detection probabilities and abundance.
#' @param obs_type Character string specifying the field observer protocol:
#'   - \code{"independent"} (default): Both observers scan independently (Lincoln-Petersen / Chapman geometry).
#'   - \code{"dependent"}: Primary observer calls out detected groups and Secondary observer records only missed groups (2-pass removal / Zippin geometry).
#' @param pooling Character string specifying how DO sector estimates are combined with Block sectors:
#'   - \code{"hybrid"} (default): Uses sector-specific DO estimates for DO sectors and extrapolates only Block sectors using the pooled/modeled primary detection probability.
#'   - \code{"global"}: Extrapolates all sectors using the pooled or covariate-modeled primary detection probability.
#' @param formula Optional formula (e.g., \code{~ covariate1 + covariate2}) to model detection probability
#'   \eqn{\hat{p}} across sectors. For \code{method = "frequentist"}, fits a quasi-binomial GLM; for
#'   \code{method = "bayesian"}, runs a coupled multi-chain Bayesian logistic regression via MCMC.
#' @param ci_method Character string specifying how to calculate Confidence/Credible Intervals:
#'   - \code{"lognormal"} (default): Uses the analytic Chao (1987) log-normal approximation for frequentist estimators, or direct posterior credible intervals from the MCMC draws when \code{method = "bayesian"}.
#'   - \code{"bootstrap"}: Uses a non-parametric stratified sector resampling bootstrap.
#' @param R Integer indicating the number of bootstrap replicates (only used if \code{ci_method = "bootstrap"}). Default is 1000.
#' @param chains Integer number of Markov chains when \code{method = "bayesian"} (default 4).
#' @param iter Integer number of total MCMC iterations per chain when \code{method = "bayesian"} (default 2000).
#' @param burnin Integer number of warm-up/burn-in iterations per chain when \code{method = "bayesian"} (default \code{floor(iter / 4)}).
#'
#' @return A list of class `doas` containing the estimate, confidence/credible intervals, sector-level estimates, and (when Bayesian) MCMC convergence diagnostics.
#' @importFrom stats var coef
#' @importFrom tools toTitleCase
#' @export
doas <- function(
    data,
    method = c("frequentist", "bayesian"),
    obs_type = c("independent", "dependent"),
    pooling = c("hybrid", "global"),
    formula = NULL,
    ci_method = c("lognormal", "bootstrap"),
    R = 1000,
    chains = 4,
    iter = 2000,
    burnin = NULL
) {
    method <- match.arg(method)
    obs_type <- match.arg(obs_type)
    pooling <- match.arg(pooling)
    ci_method <- match.arg(ci_method)

    if (is.null(burnin)) {
        burnin <- floor(iter / 4)
    }

    prep <- prepare_doas_data(data)
    data <- prep$data
    do_data <- prep$do_data

    bayes_res <- NULL
    bayes_cov_res <- NULL
    model <- NULL

    if (method == "frequentist") {
        do_data <- compute_frequentist_do(do_data, obs_type)
        glob_p <- compute_frequentist_global_p(do_data, obs_type)
        mean_p_hat <- glob_p$mean_p_hat
        var_mean_p_hat <- glob_p$var_mean_p_hat

        if (!is.null(formula)) {
            model <- fit_detection_glm(formula, do_data)
        }
    } else {
        bayes_res <- bayesian_conjugate_estimate(
            both = do_data$groups_both,
            obs1_only = do_data$groups_obs1,
            obs2_only = do_data$groups_obs2,
            obs_type = obs_type,
            iter = iter,
            chains = chains,
            burnin = burnin
        )
        do_data$G_est <- bayes_res$mean
        do_data$G_var <- bayes_res$var
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

        mu_p1_vec <- as.vector(bayes_res$mu_p1_draws)
        mean_p_hat <- mean(mu_p1_vec)
        var_mean_p_hat <- stats::var(mu_p1_vec)

        if (!is.null(formula)) {
            model <- fit_detection_glm(formula, do_data)
            bayes_cov_res <- fit_bayesian_covariate_mcmc(
                formula = formula,
                do_data = do_data,
                all_data = data,
                obs_type = obs_type,
                G_draws_list = bayes_res$G_draws,
                glm_mod = model
            )
        }
    }

    if (mean_p_hat <= 0 || is.na(mean_p_hat) || mean_p_hat > 1) {
        stop("Calculated mean detection probability is invalid.")
    }

    sec_out <- if (method == "bayesian") {
        compute_bayesian_sectors(
            data = data,
            do_data = do_data,
            obs_type = obs_type,
            pooling = pooling,
            bayes_res = bayes_res,
            bayes_cov_res = bayes_cov_res,
            iter = iter,
            burnin = burnin
        )
    } else {
        compute_frequentist_sectors(
            data = data,
            do_data = do_data,
            obs_type = obs_type,
            pooling = pooling,
            model = model,
            mean_p_hat = mean_p_hat,
            var_mean_p_hat = var_mean_p_hat
        )
    }

    results <- list(
        method = method,
        obs_type = obs_type,
        pooling = pooling,
        formula = formula,
        model = model,
        bayes_coef = if (!is.null(bayes_cov_res)) {
            bayes_cov_res$beta_summary
        } else {
            NULL
        },
        mean_p_hat = mean_p_hat,
        var_mean_p_hat = var_mean_p_hat,
        num_do_sectors = nrow(do_data),
        num_total_sectors = nrow(data),
        sector_estimates = sec_out$sector_estimates,
        total_abundance = sec_out$total_abundance,
        total_variance = sec_out$total_variance,
        se_abundance = sec_out$se_abundance,
        ci_lower = sec_out$ci_lower,
        ci_upper = sec_out$ci_upper,
        ci_method = ci_method,
        mcmc_diagnostics = sec_out$mcmc_diagnostics,
        log_lik = sec_out$log_lik
    )

    if (ci_method == "bootstrap") {
        boot_out <- run_doas_bootstrap(
            data = data,
            method = method,
            obs_type = obs_type,
            pooling = pooling,
            formula = formula,
            R = R,
            chains = chains,
            iter = iter,
            burnin = burnin
        )
        results$se_abundance <- boot_out$se_abundance
        results$total_variance <- boot_out$total_variance
        results$ci_lower <- boot_out$ci_lower
        results$ci_upper <- boot_out$ci_upper
    }

    class(results) <- "doas"
    results
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

    estimator_label <- if (object$method == "frequentist") {
        if (object$obs_type == "independent") {
            "Frequentist Independent (Chapman)"
        } else {
            "Frequentist Dependent (Zippin Removal)"
        }
    } else {
        if (object$obs_type == "independent") {
            "Bayesian MCMC Independent"
        } else {
            "Bayesian MCMC Dependent (2-Pass Removal)"
        }
    }

    cat(sprintf("Inference Method: %s\n", estimator_label))
    cat(sprintf("Sector Pooling:   %s\n", tools::toTitleCase(object$pooling)))
    if (!is.null(object$formula)) {
        cat(sprintf(
            "Covariate Formula: p_hat %s\n",
            paste(deparse(object$formula), collapse = " ")
        ))
    }
    cat(sprintf("Total sectors surveyed: %d\n", object$num_total_sectors))
    cat(sprintf(
        "Double Observer (DO) sectors: %d (%s)\n",
        object$num_do_sectors,
        tools::toTitleCase(object$obs_type)
    ))

    if (!is.null(object$mcmc_diagnostics)) {
        md <- object$mcmc_diagnostics
        cat(sprintf(
            "MCMC Sampler: %d chains, %d total iterations (%d burn-in), %d post-warmup draws\n",
            md$chains,
            md$iter,
            md$burnin,
            md$post_warmup_draws
        ))
        cat(sprintf(
            "Convergence: R-hat = %.3f | Effective Sample Size (ESS) = %d\n",
            md$rhat_abundance,
            md$ess_abundance
        ))
    }

    if (is.null(object$formula)) {
        cat(sprintf(
            "\nMean detection probability (p_hat): %.3f (Var = %.5f)\n",
            object$mean_p_hat,
            object$var_mean_p_hat
        ))
    } else if (!is.null(object$bayes_coef)) {
        cat("\nBayesian Logistic Regression Posterior (Detection Probability):\n")
        print(round(object$bayes_coef, 3))
        waic_val <- suppressWarnings(WAIC(object))
        cat(sprintf(
            "Model Fit: WAIC = %.2f (p_WAIC = %.2f, lppd = %.2f)\n",
            as.numeric(waic_val),
            attr(waic_val, "p_WAIC"),
            attr(waic_val, "lppd")
        ))
    } else {
        cat("\nDetection Probability Model (GLM):\n")
        print(stats::coef(summary(object$model)))
        qaic_val <- suppressWarnings(QAIC(object))
        if (!is.na(qaic_val)) {
            cat(sprintf("Model Fit: QAIC = %.2f\n", qaic_val))
        }
    }

    ci_label <- if (object$ci_method == "bootstrap") {
        "Bootstrap Percentile CI"
    } else if (identical(object$method, "bayesian")) {
        "Bayesian Credible Interval"
    } else {
        "Log-normal CI"
    }

    cat(sprintf("\nEstimated Total Abundance: %.1f\n", object$total_abundance))
    cat(sprintf("Standard Error (Posterior SD): %.1f\n", object$se_abundance))
    cat(sprintf(
        "95%% %s: [%.1f - %.1f]\n",
        ci_label,
        object$ci_lower,
        object$ci_upper
    ))
    cat("----------------------------------------------------\n")
    cat("Sector-level Estimates:\n")
    num_cols <- vapply(object$sector_estimates, is.numeric, logical(1))
    print_df <- object$sector_estimates
    print_df[, num_cols] <- round(print_df[, num_cols], 2)
    print(print_df, row.names = FALSE)
    cat("----------------------------------------------------\n")
}
