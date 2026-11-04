#' Validate and Prepare Survey Data for DOAS
#'
#' Validates required columns, checks consistency between `detected_groups` and
#' `(groups_both, groups_obs1, groups_obs2)` in DO sectors, calculates `mean_group_size`,
#' and identifies `"DO"`, `"Block"`, and `"Unsurveyed"` survey sectors.
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
#' @param p_formula R formula for detection covariates.
#' @param do_data Data frame of DO sectors with `p_hat` and `G_est`.
#' @return A fitted `glm` object.
#' @importFrom stats glm quasibinomial update
#' @keywords internal
#' @noRd
fit_detection_glm <- function(p_formula, do_data) {
    if (!inherits(p_formula, "formula")) {
        stop("p_formula must be a valid R formula (e.g. ~ covariate)")
    }

    do_data_pos <- do_data[do_data$G_est > 0, , drop = FALSE]
    if (nrow(do_data_pos) < 2) {
        stop("At least 2 DO sectors with positive group detections are required to fit a detection covariate GLM.")
    }

    model_formula <- stats::update(p_formula, p_hat ~ .)
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
#' @param p_formula Optional detection covariate formula.
#' @param n_formula Optional abundance covariate formula.
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
    p_formula,
    n_formula,
    R,
    chains,
    iter,
    burnin
) {
    do_indices <- which(data$survey_type == "DO")
    block_indices <- setdiff(seq_len(nrow(data)), do_indices)
    total_seen <- sum(data$total_count, na.rm = TRUE)

    boot_single <- function(b) {
        boot_do <- data[
            sample(do_indices, length(do_indices), replace = TRUE),
            ,
            drop = FALSE
        ]
        boot_block <- if (length(block_indices) > 0) {
            data[
                sample(block_indices, length(block_indices), replace = TRUE),
                ,
                drop = FALSE
            ]
        } else {
            data[integer(0), , drop = FALSE]
        }
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
                    p_formula = p_formula,
                    n_formula = n_formula,
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
#' Estimates total population abundance using either a full Double Observer (DO)
#' survey across all sectors or a partial Double Observer Adjusted Survey (DOAS)
#' combining DO surveys in a subset of sectors with single-observer block counts
#' (and optionally unsurveyed sectors imputed via \code{n_formula}).
#' Supports both Frequentist and multi-chain Bayesian MCMC inference across
#' independent (Lincoln-Petersen/Chapman) and dependent (2-pass removal/Zippin)
#' double-observer field protocols.
#'
#' @param data A dataframe or \code{\link{doas_data}} object containing survey data.
#'   Each row represents a sector:
#'   - \code{sector_id}: Unique identifier for the sector.
#'   - \code{total_count}: Total individuals counted across all groups in the sector (\code{NA} if unsurveyed).
#'   - \code{detected_groups}: Total groups detected in the sector.
#'   - \code{groups_both}, \code{groups_obs1}, \code{groups_obs2}: Double Observer group counts (non-NA for DO sectors, \code{NA} for Block/Unsurveyed sectors).
#' @param method Character string specifying the statistical inference paradigm:
#'   - \code{"frequentist"} (default): Uses analytic Chapman (when \code{obs_type = "independent"}) or Zippin removal (when \code{obs_type = "dependent"}) estimators and Delta-method variance.
#'   - \code{"bayesian"}: Uses a multi-chain Collapsed Poisson-Gamma-Beta Gibbs Sampler (MCMC) to sample the full joint posterior of detection probabilities and abundance.
#' @param obs_type Character string specifying the field observer protocol:
#'   - \code{"independent"} (default): Both observers scan independently (Lincoln-Petersen / Chapman geometry).
#'   - \code{"dependent"}: Primary observer calls out detected groups and Secondary observer records only missed groups (2-pass removal / Zippin geometry).
#' @param pooling Character string specifying how DO sector estimates are combined with Block sectors:
#'   - \code{"hybrid"} (default): Uses sector-specific DO estimates for DO sectors and extrapolates only Block sectors using the pooled/modeled detection probability.
#'   - \code{"global"}: Extrapolates all surveyed sectors using the pooled or covariate-modeled detection probability.
#' @param p_formula Optional one-sided formula (e.g., \code{~ effort * ruggedness}) to model detection probability
#'   \eqn{\hat{p}} across sectors. For \code{method = "frequentist"}, fits a quasi-binomial GLM; for
#'   \code{method = "bayesian"}, runs a coupled multi-chain Bayesian logistic regression via MCMC.
#' @param n_formula Optional one-sided formula (e.g., \code{~ elevation + habitat + offset(log(area))}) to model
#'   sector-level abundance \eqn{\lambda_i} as a function of environmental covariates (Quasi-Poisson log-link GLM in
#'   frequentist mode; Bayesian Poisson log-linear regression in Bayesian mode). Enables \code{\link{predict.doas}}
#'   for new/unsurveyed sectors and habitat response curves.
#' @param ci_method Character string specifying how to calculate Confidence/Credible Intervals:
#'   - \code{"lognormal"} (default): Uses the analytic Chao (1987) log-normal approximation for frequentist estimators, or direct posterior credible intervals from the MCMC draws when \code{method = "bayesian"}.
#'   - \code{"bootstrap"}: Uses a non-parametric stratified sector resampling bootstrap.
#' @param R Integer indicating the number of bootstrap replicates (only used if \code{ci_method = "bootstrap"}). Default is 1000.
#' @param chains Integer number of Markov chains when \code{method = "bayesian"} (default 4).
#' @param iter Integer number of total MCMC iterations per chain when \code{method = "bayesian"} (default 2000).
#' @param burnin Integer number of warm-up/burn-in iterations per chain when \code{method = "bayesian"} (default \code{floor(iter / 4)}).
#'
#' @return A list of class `doas` containing the estimate, confidence/credible intervals, sector-level estimates, fitted submodels, and (when Bayesian) MCMC convergence diagnostics.
#' @importFrom stats var coef
#' @importFrom tools toTitleCase
#' @export
doas <- function(
    data,
    method = c("frequentist", "bayesian"),
    obs_type = c("independent", "dependent"),
    pooling = c("hybrid", "global"),
    p_formula = NULL,
    n_formula = NULL,
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

    if (any(data$survey_type == "Unsurveyed") && is.null(n_formula)) {
        stop(
            "Unsurveyed sectors (total_count = NA) require an abundance formula (n_formula = ~ ...) to predict their abundance."
        )
    }

    bayes_res <- NULL
    bayes_cov_res <- NULL
    model <- NULL

    if (method == "frequentist") {
        do_data <- compute_frequentist_do(do_data, obs_type)
        glob_p <- compute_frequentist_global_p(do_data, obs_type)
        mean_p_hat <- glob_p$mean_p_hat
        var_mean_p_hat <- glob_p$var_mean_p_hat

        if (!is.null(p_formula)) {
            model <- fit_detection_glm(p_formula, do_data)
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

        if (!is.null(p_formula)) {
            model <- fit_detection_glm(p_formula, do_data)
            bayes_cov_res <- fit_bayesian_covariate_mcmc(
                formula = p_formula,
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
            n_formula = n_formula,
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
            n_formula = n_formula,
            mean_p_hat = mean_p_hat,
            var_mean_p_hat = var_mean_p_hat
        )
    }

    results <- list(
        method = method,
        obs_type = obs_type,
        pooling = pooling,
        p_formula = p_formula,
        formula = p_formula,
        n_formula = n_formula,
        model = model,
        p_model = model,
        n_model = sec_out$n_model,
        bayes_coef = if (!is.null(bayes_cov_res)) {
            bayes_cov_res$beta_summary
        } else {
            NULL
        },
        bayes_n_coef = if (!is.null(sec_out$bayes_n_res)) {
            sec_out$bayes_n_res$gamma_summary
        } else {
            NULL
        },
        p_beta_draws = if (!is.null(bayes_cov_res)) {
            bayes_cov_res$beta_draws
        } else {
            NULL
        },
        n_gamma_draws = if (!is.null(sec_out$bayes_n_res)) {
            sec_out$bayes_n_res$gamma_draws
        } else {
            NULL
        },
        mu_p1_draws = if (!is.null(bayes_res)) {
            bayes_res$mu_p1_draws
        } else {
            NULL
        },
        mean_gs_pool = sec_out$mean_gs_pool,
        data = data,
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
        log_lik = sec_out$log_lik,
        n_log_lik = if (!is.null(sec_out$bayes_n_res)) {
            sec_out$bayes_n_res$n_log_lik
        } else {
            NULL
        }
    )

    if (ci_method == "bootstrap") {
        boot_out <- run_doas_bootstrap(
            data = data,
            method = method,
            obs_type = obs_type,
            pooling = pooling,
            p_formula = p_formula,
            n_formula = n_formula,
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

#' Predict Method for DOAS Models
#'
#' Predicts sector-level abundance (\code{type = "abundance"}), group count
#' (\code{type = "groups"}), or single-observer detection probability
#' (\code{type = "detection"}) from a fitted \code{\link{doas}} object, either for
#' the original surveyed sectors (\code{newdata = NULL}) or for a new data frame
#' of covariates (\code{newdata}).
#'
#' @param object A fitted object of class \code{"doas"}.
#' @param newdata Optional data frame containing covariates required by \code{p_formula}
#'   (when \code{type = "detection"}) or \code{n_formula} (when \code{type = "abundance"}
#'   or \code{"groups"}). If \code{NULL}, predictions are returned for the original sectors.
#' @param type Character string specifying the prediction target:
#'   - \code{"abundance"} (default): Expected individual abundance \eqn{\hat{N}_i}.
#'   - \code{"detection"}: Single-observer detection probability \eqn{\hat{p}_i \in (0, 1)}.
#'   - \code{"groups"}: Expected group abundance \eqn{\hat{G}_i}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with columns \code{fit}, \code{se}, \code{ci_lower}, and \code{ci_upper}
#'   (plus \code{sector_id} if available).
#' @importFrom stats predict plogis qnorm quantile sd model.matrix model.frame model.offset delete.response terms
#' @export
predict.doas <- function(
    object,
    newdata = NULL,
    type = c("abundance", "detection", "groups"),
    ...
) {
    type <- match.arg(type)
    eval_df <- if (is.null(newdata)) object$data else newdata
    n_rows <- nrow(eval_df)
    has_sec_id <- "sector_id" %in% names(eval_df)

    if (type == "detection") {
        if (is.null(object$p_formula)) {
            p_val <- object$mean_p_hat
            se_val <- sqrt(object$var_mean_p_hat)
            if (identical(object$method, "bayesian") && !is.null(object$mu_p1_draws)) {
                qs <- stats::quantile(
                    as.vector(object$mu_p1_draws),
                    probs = c(0.025, 0.975),
                    names = FALSE
                )
            } else {
                qs <- c(
                    max(0, p_val - 1.96 * se_val),
                    min(1, p_val + 1.96 * se_val)
                )
            }
            out <- data.frame(
                fit = rep(p_val, n_rows),
                se = rep(se_val, n_rows),
                ci_lower = rep(qs[1], n_rows),
                ci_upper = rep(qs[2], n_rows)
            )
        } else if (identical(object$method, "frequentist")) {
            pred_link <- stats::predict(
                object$p_model,
                newdata = eval_df,
                type = "link",
                se.fit = TRUE
            )
            eta <- as.numeric(pred_link$fit)
            se_eta <- as.numeric(pred_link$se.fit)
            p_fit <- stats::plogis(eta)
            se_p <- p_fit * (1 - p_fit) * se_eta
            out <- data.frame(
                fit = p_fit,
                se = se_p,
                ci_lower = stats::plogis(eta - 1.96 * se_eta),
                ci_upper = stats::plogis(eta + 1.96 * se_eta)
            )
        } else {
            rhs_terms <- stats::delete.response(stats::terms(object$p_model))
            X_new <- stats::model.matrix(
                rhs_terms,
                data = eval_df,
                xlev = object$p_model$xlevels
            )
            b_draws <- object$p_beta_draws
            n_keep <- dim(b_draws)[1]
            n_ch <- dim(b_draws)[2]
            n_par <- dim(b_draws)[3]
            b_mat <- matrix(b_draws, nrow = n_keep * n_ch, ncol = n_par)
            eta_mat <- X_new %*% t(b_mat)
            p_mat <- 1 / (1 + exp(-pmax(pmin(eta_mat, 20), -20)))

            out <- data.frame(
                fit = rowMeans(p_mat),
                se = apply(p_mat, 1, stats::sd),
                ci_lower = apply(p_mat, 1, stats::quantile, probs = 0.025, names = FALSE),
                ci_upper = apply(p_mat, 1, stats::quantile, probs = 0.975, names = FALSE)
            )
        }
    } else {
        # type == "abundance" or "groups"
        scale_div <- if (type == "groups") max(object$mean_gs_pool, 1e-6) else 1.0

        if (is.null(object$n_formula)) {
            if (!is.null(newdata)) {
                stop(
                    "Predicting abundance or groups on 'newdata' requires a fitted abundance model (n_formula = ~ ...)."
                )
            }
            sec <- object$sector_estimates
            if (type == "abundance") {
                out <- data.frame(
                    fit = sec$est_abundance,
                    se = sec$se_abundance,
                    ci_lower = sec$ci_lower,
                    ci_upper = sec$ci_upper
                )
            } else {
                out <- data.frame(
                    fit = sec$est_groups,
                    se = sec$se_abundance / scale_div,
                    ci_lower = sec$ci_lower / scale_div,
                    ci_upper = sec$ci_upper / scale_div
                )
            }
        } else if (identical(object$method, "frequentist")) {
            pred_link <- stats::predict(
                object$n_model,
                newdata = eval_df,
                type = "link",
                se.fit = TRUE
            )
            eta <- as.numeric(pred_link$fit)
            se_eta <- as.numeric(pred_link$se.fit)
            fit_val <- exp(eta) / scale_div
            se_val <- fit_val * se_eta
            out <- data.frame(
                fit = fit_val,
                se = se_val,
                ci_lower = exp(eta - 1.96 * se_eta) / scale_div,
                ci_upper = exp(eta + 1.96 * se_eta) / scale_div
            )
        } else {
            rhs_terms <- stats::delete.response(stats::terms(object$n_model))
            mf_new <- stats::model.frame(
                rhs_terms,
                data = eval_df,
                xlev = object$n_model$xlevels
            )
            Z_new <- stats::model.matrix(rhs_terms, data = mf_new)
            off_new <- stats::model.offset(mf_new)
            if (is.null(off_new)) off_new <- rep(0, nrow(Z_new))

            g_draws <- object$n_gamma_draws
            n_keep <- dim(g_draws)[1]
            n_ch <- dim(g_draws)[2]
            n_par <- dim(g_draws)[3]
            g_mat <- matrix(g_draws, nrow = n_keep * n_ch, ncol = n_par)
            eta_mat <- sweep(Z_new %*% t(g_mat), 1, off_new, "+")
            lam_mat <- exp(pmax(pmin(eta_mat, 15), -15)) / scale_div

            out <- data.frame(
                fit = rowMeans(lam_mat),
                se = apply(lam_mat, 1, stats::sd),
                ci_lower = apply(lam_mat, 1, stats::quantile, probs = 0.025, names = FALSE),
                ci_upper = apply(lam_mat, 1, stats::quantile, probs = 0.975, names = FALSE)
            )
        }
    }

    if (has_sec_id) {
        out <- cbind(
            data.frame(sector_id = eval_df$sector_id, stringsAsFactors = FALSE),
            out
        )
    }
    rownames(out) <- NULL
    out
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
    if (!is.null(object$p_formula)) {
        cat(sprintf(
            "Detection Formula (p_formula): p_hat %s\n",
            paste(deparse(object$p_formula), collapse = " ")
        ))
    }
    if (!is.null(object$n_formula)) {
        cat(sprintf(
            "Abundance Formula (n_formula): N_hat %s\n",
            paste(deparse(object$n_formula), collapse = " ")
        ))
    }
    cat(sprintf("Total sectors: %d\n", object$num_total_sectors))
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

    if (is.null(object$p_formula)) {
        cat(sprintf(
            "\nMean detection probability (p_hat): %.3f (Var = %.5f)\n",
            object$mean_p_hat,
            object$var_mean_p_hat
        ))
    } else if (!is.null(object$bayes_coef)) {
        cat("\nBayesian Logistic Regression Posterior (Detection p_formula):\n")
        print(round(object$bayes_coef, 3))
    } else if (!is.null(object$p_model)) {
        cat("\nDetection Probability Model (p_formula GLM):\n")
        print(stats::coef(summary(object$p_model)))
    }

    if (!is.null(object$n_formula)) {
        if (!is.null(object$bayes_n_coef)) {
            cat("\nBayesian Log-Linear Regression Posterior (Abundance n_formula):\n")
            print(round(object$bayes_n_coef, 3))
        } else if (!is.null(object$n_model)) {
            cat("\nSector Abundance Model (n_formula GLM):\n")
            print(stats::coef(summary(object$n_model)))
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
