#' Quasi-Akaike Information Criterion (QAIC)
#'
#' Computes the Quasi-Akaike Information Criterion (QAIC) for frequentist statistical
#' models that handle overdispersion (such as quasi-binomial GLM covariate models in \code{\link{doas}}).
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
#' Extracts the QAIC for frequentist DOAS models that used a covariate formula.
#' Standard AIC is undefined for the quasi-binomial GLM used under the hood,
#' so QAIC penalizes the binomial likelihood by the estimated overdispersion
#' parameter (\eqn{\hat{c}}).
#'
#' @param object A fitted `doas` object with \code{method = "frequentist"}.
#' @param ... Optionally more fitted `doas` objects.
#'
#' @return A numeric QAIC value, or a data frame if multiple models are passed.
#' @importFrom stats dbinom coef
#' @export
#'
#' @examples
#' \dontrun{
#' mod1 <- doas(data, method = "frequentist", formula = ~ effort)
#' mod2 <- doas(data, method = "frequentist", formula = ~ effort + habitat)
#' QAIC(mod1, mod2)
#' }
QAIC.doas <- function(object, ...) {
    get_qaic <- function(mod) {
        if (!inherits(mod, "doas")) {
            return(NA_real_)
        }
        if (identical(mod$method, "bayesian")) {
            warning(
                "QAIC is intended for frequentist GLM models. For Bayesian MCMC models (method = 'bayesian'), use WAIC() instead."
            )
        }
        if (is.null(mod$formula)) {
            warning(
                "DOAS model did not use a covariate formula. QAIC is undefined."
            )
            return(NA_real_)
        }

        glm_mod <- mod$model
        y <- glm_mod$y
        mu <- glm_mod$fitted.values
        wt <- glm_mod$prior.weights
        log_L <- sum(ifelse(
            wt > 0,
            wt * log(stats::dbinom(round(y * wt), round(wt), mu)),
            0
        ))

        c_hat <- summary(glm_mod)$dispersion
        K <- length(stats::coef(glm_mod))
        qaic_val <- -2 * (log_L / c_hat) + 2 * K

        if (is.nan(qaic_val) || is.na(qaic_val)) {
            warning(
                "QAIC cannot be computed because the model ran out of residual degrees of freedom (c-hat is NaN)."
            )
            return(NA_real_)
        }

        qaic_val
    }

    mods <- list(object, ...)
    qaic_vals <- vapply(mods, get_qaic, FUN.VALUE = numeric(1))

    if (length(mods) == 1) {
        return(qaic_vals)
    }

    call_names <- as.character(match.call()[-1])
    res_df <- data.frame(
        df = vapply(
            mods,
            function(m) {
                if (inherits(m, "doas") && !is.null(m$formula)) {
                    length(stats::coef(m$model))
                } else {
                    NA_real_
                }
            },
            FUN.VALUE = numeric(1)
        ),
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
#' Can compare both covariate models (\code{formula = ~ ...}) and the null model (\code{formula = NULL}).
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
#' pointwise log-likelihood across MCMC draws as recommended by Gelman et al. 2014)
#' across Double Observer sectors.
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
#' mod_effort <- doas(data, method = "bayesian", formula = ~ effort)
#' WAIC(mod_null, mod_effort)
#' }
WAIC.doas <- function(object, ...) {
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

        ll_mat <- mod$log_lik
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

    mods <- list(object, ...)
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
