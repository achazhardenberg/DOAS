#' Quasi-Akaike Information Criterion (QAIC)
#'
#' Computes the Quasi-Akaike Information Criterion (QAIC; Burnham & Anderson 2002)
#' for frequentist statistical models that handle overdispersion (such as quasi-binomial
#' detection GLM models via \code{p_formula} or quasi-poisson abundance GLM models
#' via \code{n_formula} in \code{\link{doas}}).
#' For Bayesian MCMC models (\code{method = "bayesian"}), use \code{\link{WAIC}} instead.
#'
#' @param object A fitted model object (e.g., a frequentist `doas` model).
#' @param ... Optionally more fitted model objects for model comparison.
#'
#' @return A numeric QAIC value, or a comparison data frame if multiple models are passed.
#' @export
QAIC <- function(object, ...) {
    UseMethod("QAIC")
}

#' QAIC Method for Frequentist DOAS Models
#'
#' Extracts the QAIC for frequentist DOAS models that used a detection covariate formula
#' (\code{p_formula}) or an abundance covariate formula (\code{n_formula}).
#' Following Burnham & Anderson (2002), when comparing multiple candidate models, a single
#' common variance inflation factor \eqn{\hat{c} = \max(1, \hat{c}_{\text{global}})} is
#' estimated from the most parameterized model in the candidate set.
#'
#' @param object A fitted `doas` object with \code{method = "frequentist"}.
#' @param ... Optionally more fitted `doas` objects.
#'
#' @return A numeric QAIC value, or a data frame if multiple models are passed.
#' @importFrom stats dbinom dpois coef
#' @export
#'
#' @examples
#' \dontrun{
#' mod1 <- doas(data, method = "frequentist", p_formula = ~ effort)
#' mod2 <- doas(data, method = "frequentist", p_formula = ~ effort + habitat)
#' QAIC(mod1, mod2)
#' }
QAIC.doas <- function(object, ...) {
    mods <- list(object, ...)

    has_n <- vapply(
        mods,
        function(m) inherits(m, "doas") && !is.null(m$n_formula),
        logical(1)
    )
    has_p <- vapply(
        mods,
        function(m) inherits(m, "doas") && !is.null(m$p_formula),
        logical(1)
    )
    use_n_submodel <- all(has_n) && (
        !any(has_p) ||
            length(unique(vapply(mods, function(m) paste(deparse(m$p_formula), collapse = " "), character(1)))) == 1
    )

    extract_glm <- function(mod) {
        if (!inherits(mod, "doas") || identical(mod$method, "bayesian")) {
            return(NULL)
        }
        if (use_n_submodel) mod$n_model else mod$p_model
    }

    glm_list <- lapply(mods, extract_glm)
    k_vec <- vapply(
        glm_list,
        function(g) if (is.null(g)) NA_real_ else length(stats::coef(g)),
        numeric(1)
    )
    disp_vec <- vapply(
        glm_list,
        function(g) {
            if (is.null(g)) {
                NA_real_
            } else {
                d <- tryCatch(summary(g)$dispersion, error = function(e) NA_real_)
                if (is.nan(d)) NA_real_ else d
            }
        },
        numeric(1)
    )

    # Common c_hat from the most parameterized valid model in the set (Burnham & Anderson 2002), bounded >= 1
    valid_disp <- which(!is.na(disp_vec) & is.finite(disp_vec) & disp_vec > 0)
    common_c_hat <- if (length(valid_disp) > 0) {
        best_idx <- valid_disp[which.max(k_vec[valid_disp])]
        max(1.0, disp_vec[best_idx])
    } else {
        NA_real_
    }

    get_qaic <- function(idx) {
        mod <- mods[[idx]]
        if (!inherits(mod, "doas")) {
            return(NA_real_)
        }
        if (identical(mod$method, "bayesian")) {
            warning(
                "QAIC is intended for frequentist GLM models. For Bayesian MCMC models (method = 'bayesian'), use WAIC() instead."
            )
            return(NA_real_)
        }
        glm_mod <- glm_list[[idx]]
        if (is.null(glm_mod)) {
            warning(
                "DOAS model did not use a covariate formula (p_formula or n_formula). QAIC is undefined."
            )
            return(NA_real_)
        }
        if (is.na(disp_vec[idx]) || is.na(common_c_hat)) {
            warning(
                "QAIC cannot be computed because the model ran out of residual degrees of freedom (c-hat is NaN)."
            )
            return(NA_real_)
        }

        if (use_n_submodel) {
            y <- round(glm_mod$y)
            mu <- pmax(1e-6, glm_mod$fitted.values)
            log_L <- sum(stats::dpois(y, lambda = mu, log = TRUE))
        } else {
            y <- glm_mod$y
            mu <- pmin(1 - 1e-6, pmax(1e-6, glm_mod$fitted.values))
            wt <- glm_mod$prior.weights
            log_L <- sum(ifelse(
                wt > 0,
                stats::dbinom(round(y * wt), round(wt), mu, log = TRUE),
                0
            ))
        }

        K <- k_vec[idx]
        -2 * (log_L / common_c_hat) + 2 * K
    }

    qaic_vals <- vapply(seq_along(mods), get_qaic, FUN.VALUE = numeric(1))

    if (length(mods) == 1) {
        return(qaic_vals)
    }

    call_names <- as.character(match.call()[-1])
    res_df <- data.frame(
        df = k_vec,
        c_hat = rep(round(common_c_hat, 3), length(mods)),
        QAIC = qaic_vals
    )
    rownames(res_df) <- call_names
    res_df$delta_QAIC <- res_df$QAIC - min(res_df$QAIC, na.rm = TRUE)
    res_df[order(res_df$QAIC), ]
}

#' Widely Applicable Information Criterion (WAIC)
#'
#' Computes the Watanabe-Akaike (Widely Applicable) Information Criterion (WAIC;
#' Watanabe 2010; Gelman, Hwang & Vehtari 2014) from the full pointwise posterior
#' log-likelihood draws of Bayesian MCMC \code{\link{doas}} models (\code{method = "bayesian"}).
#' Can compare both detection covariate models (\code{p_formula = ~ ...}), abundance
#' covariate models (\code{n_formula = ~ ...}), and the null model.
#'
#' @param object A fitted Bayesian `doas` model object (\code{method = "bayesian"}).
#' @param ... Optionally more fitted Bayesian `doas` model objects for model comparison.
#'
#' @return If a single model is provided, a numeric scalar WAIC value (with attributes
#'   \code{lppd}, \code{p_WAIC}, and \code{SE}). If multiple models are provided, a sorted
#'   comparison data frame with columns \code{lppd}, \code{p_WAIC}, \code{WAIC}, \code{SE},
#'   and \code{delta_WAIC}.
#' @export
WAIC <- function(object, ...) {
    UseMethod("WAIC")
}

#' WAIC Method for Bayesian DOAS Models
#'
#' Evaluates the log pointwise predictive density (\code{lppd}) and the effective
#' number of parameters (\code{p_WAIC}, computed via the posterior variance of the
#' pointwise log-likelihood across MCMC draws as recommended by Gelman et al. 2014).
#'
#' @param object A fitted `doas` object with \code{method = "bayesian"}.
#' @param ... Optionally more fitted Bayesian `doas` objects.
#'
#' @return A numeric WAIC value (for a single model) or a comparison data frame (for multiple models).
#' @importFrom stats var dbinom
#' @export
#'
#' @examples
#' \dontrun{
#' mod_null   <- doas(data, method = "bayesian")
#' mod_effort <- doas(data, method = "bayesian", p_formula = ~ effort)
#' WAIC(mod_null, mod_effort)
#' }
WAIC.doas <- function(object, ...) {
    mods <- list(object, ...)
    has_n <- vapply(
        mods,
        function(m) inherits(m, "doas") && !is.null(m$n_log_lik),
        logical(1)
    )
    has_p <- vapply(
        mods,
        function(m) inherits(m, "doas") && !is.null(m$p_formula),
        logical(1)
    )
    use_n_submodel <- all(has_n) && (
        !any(has_p) ||
            length(unique(vapply(mods, function(m) paste(deparse(m$p_formula), collapse = " "), character(1)))) == 1
    )

    compute_single_waic <- function(mod) {
        if (!inherits(mod, "doas")) {
            return(c(lppd = NA_real_, p_WAIC = NA_real_, WAIC = NA_real_, SE = NA_real_))
        }
        if (!identical(mod$method, "bayesian") || is.null(mod$log_lik)) {
            warning(
                "WAIC requires posterior MCMC log-likelihood draws from a Bayesian DOAS model (method = 'bayesian'). For frequentist models, use QAIC() instead."
            )
            return(c(lppd = NA_real_, p_WAIC = NA_real_, WAIC = NA_real_, SE = NA_real_))
        }

        ll_mat <- if (use_n_submodel && !is.null(mod$n_log_lik)) {
            mod$n_log_lik
        } else {
            mod$log_lik
        }
        n_obs <- nrow(ll_mat)

        lppd_i <- numeric(n_obs)
        p_waic_i <- numeric(n_obs)

        for (i in seq_len(n_obs)) {
            ll_i <- ll_mat[i, ]
            m_i <- max(ll_i)
            lppd_i[i] <- m_i + log(mean(exp(ll_i - m_i)))
            p_waic_i[i] <- stats::var(ll_i)
        }

        waic_i <- -2 * (lppd_i - p_waic_i)
        lppd_tot <- sum(lppd_i)
        p_waic_tot <- sum(p_waic_i)
        waic_tot <- sum(waic_i)
        se_waic <- if (n_obs > 1) sqrt(n_obs * stats::var(waic_i)) else NA_real_

        c(lppd = lppd_tot, p_WAIC = p_waic_tot, WAIC = waic_tot, SE = se_waic)
    }

    mat <- vapply(mods, compute_single_waic, FUN.VALUE = numeric(4))

    if (length(mods) == 1) {
        val <- unname(mat["WAIC", 1])
        attr(val, "lppd") <- unname(mat["lppd", 1])
        attr(val, "p_WAIC") <- unname(mat["p_WAIC", 1])
        attr(val, "SE") <- unname(mat["SE", 1])
        return(val)
    }

    call_names <- as.character(match.call()[-1])
    res_df <- as.data.frame(t(mat))
    rownames(res_df) <- call_names
    res_df$delta_WAIC <- res_df$WAIC - min(res_df$WAIC, na.rm = TRUE)
    res_df[order(res_df$WAIC), ]
}
