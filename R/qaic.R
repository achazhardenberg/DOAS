#' Quasi-Akaike Information Criterion
#'
#' Computes the Quasi-Akaike Information Criterion (QAIC) for statistical models
#' that handle overdispersion (such as quasi-binomial models).
#'
#' @param object A fitted model object (e.g., a `doas` model).
#' @param ... Optionally more fitted model objects.
#'
#' @return A numeric QAIC value, or a data frame if multiple models are passed.
#' @export
#'
QAIC <- function(object, ...) {
    UseMethod("QAIC")
}

#' QAIC method for DOAS models
#'
#' Extracts the QAIC for DOAS models that used a covariate formula.
#' Standard AIC is undefined for the quasi-binomial GLM used under the hood,
#' so QAIC penalizes the binomial likelihood by the estimated overdispersion
#' parameter (c-hat).
#'
#' @param object A fitted `doas` object.
#' @param ... Optionally more fitted `doas` objects.
#'
#' @return A numeric QAIC value, or a data frame if multiple models are passed.
#' @export
#'
#' @examples
#' \dontrun{
#' mod1 <- doas(data, formula = ~ effort)
#' mod2 <- doas(data, formula = ~ effort + habitat)
#' QAIC(mod1, mod2)
#' }
QAIC.doas <- function(object, ...) {
    # helper for a single model
    get_qaic <- function(mod) {
        if (!inherits(mod, "doas")) {
            return(NA)
        }

        # If no formula was used (i.e. just flat mean), QAIC isn't really comparable
        if (is.null(mod$formula)) {
            warning(
                "DOAS model did not use a covariate formula. QAIC is undefined."
            )
            return(NA)
        }

        glm_mod <- mod$model

        # 1. Get Log-Likelihood of the equivalent Binomial model
        # To avoid R `update()` scoping bugs when doas() is run inside other functions or tests,
        # we manually compute the exact Binomial log-likelihood using the cached data.
        y <- glm_mod$y
        mu <- glm_mod$fitted.values
        wt <- glm_mod$prior.weights
        log_L <- sum(ifelse(
            wt > 0,
            wt * log(stats::dbinom(round(y * wt), round(wt), mu)),
            0
        ))

        # 2. Get the overdispersion parameter (c-hat) from the quasibinomial summary
        c_hat <- summary(glm_mod)$dispersion

        # 3. Number of parameters (K)
        K <- length(stats::coef(glm_mod))

        # 4. Compute QAIC
        qaic_val <- -2 * (log_L / c_hat) + 2 * K

        # If c_hat is NaN (e.g. 0 residual degrees of freedom due to too few DO sectors)
        if (is.nan(qaic_val) || is.na(qaic_val)) {
            warning(
                "QAIC cannot be computed because the model ran out of residual degrees of freedom (c-hat is NaN)."
            )
            return(NA)
        }

        return(qaic_val)
    }

    mods <- list(object, ...)
    qaic_vals <- vapply(mods, get_qaic, FUN.VALUE = numeric(1))

    # If only one model was passed, return raw scalar
    if (length(mods) == 1) {
        return(qaic_vals)
    }

    # If multiple models passed, return a dataframe mimicking typical AIC outputs
    call_names <- as.character(match.call()[-1])
    res_df <- data.frame(
        df = vapply(
            mods,
            function(m) {
                if (inherits(m, "doas") && !is.null(m$formula)) {
                    length(stats::coef(m$model))
                } else {
                    NA
                }
            },
            FUN.VALUE = numeric(1)
        ),
        QAIC = qaic_vals
    )
    rownames(res_df) <- call_names

    # Sort by QAIC
    res_df <- res_df[order(res_df$QAIC), ]
    return(res_df)
}
