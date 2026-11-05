#' Draw from Posterior Predictive Distribution for DOAS Bayesian Models
#'
#' Simulates replicated datasets from the posterior predictive distribution of a fitted
#' \code{\link{doas}} Bayesian model. For double-observer sectors, observations are simulated
#' from the multinomial capture-recapture process conditional on posterior draws of latent
#' abundance and observer detection probabilities. For single-observer sectors, observations
#' are simulated from a binomial detection process.
#'
#' The returned matrix \code{yrep} has dimensions \code{n_draws x n_sectors} and includes
#' the observed data vector as \code{attr(yrep, "y")}, making it directly compatible with
#' \code{bayesplot::ppc_*} functions.
#'
#' @param object A fitted \code{doas} model object with \code{method = "Bayes"}.
#' @param n_draws Number of posterior predictive draws to simulate. If \code{NULL} (default),
#'   all available post-warmup draws across chains are used.
#' @param type Character string indicating what quantity to predict:
#'   \itemize{
#'     \item \code{"count"} or \code{"detected"} (default): Total unique individuals/groups detected
#'       at each sector (comparable to \code{object$data$total_count}).
#'     \item \code{"both"} or \code{"joint"}: Joint detections \eqn{y_{12}} in double-observer sectors.
#'     \item \code{"obs1"}: Observer 1 detections (\eqn{y_1 + y_{12}} in DO sectors).
#'     \item \code{"obs2"}: Observer 2 detections (\eqn{y_2 + y_{12}} in DO sectors).
#'   }
#' @param drop_unsurveyed Logical. If \code{TRUE} (default), unsurveyed sectors (which have
#'   no field observations) are omitted from \code{yrep} and \code{attr(yrep, "y")}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A matrix of class \code{c("doas_yrep", "matrix")} of dimension \code{n_draws x n_sectors},
#'   with attributes \code{"y"} (observed vector), \code{"sector_id"} (sector identifiers),
#'   and \code{"target"} (predicted quantity).
#' @importFrom stats rmultinom rbinom rbeta quantile median sd var
#' @rdname posterior_predict
#' @export
posterior_predict <- function(object, ...) {
    UseMethod("posterior_predict")
}

#' @rdname posterior_predict
#' @export
posterior_predict.doas <- function(
    object,
    n_draws = NULL,
    type = c("count", "detected", "both", "joint", "obs1", "obs2"),
    drop_unsurveyed = TRUE,
    ...
) {
    if (!inherits(object, "doas")) {
        stop("object must be of class 'doas'.", call. = FALSE)
    }
    if (object$method != "Bayes") {
        stop("posterior_predict() is only supported for Bayesian models (method = 'Bayes').", call. = FALSE)
    }
    if (is.null(object$mcmc_samples) || is.null(object$mcmc_samples$sector_N)) {
        stop("No posterior draws found in object$mcmc_samples.", call. = FALSE)
    }

    type <- match.arg(type)
    if (type == "detected") type <- "count"
    if (type == "joint") type <- "both"

    data <- object$data
    n_sec <- nrow(data)

    # Determine available draws
    sec_N_list <- object$mcmc_samples$sector_N
    n_keep <- nrow(sec_N_list[[1]])
    n_chains <- ncol(sec_N_list[[1]])
    s_avail <- n_keep * n_chains

    if (is.null(n_draws)) {
        n_draws <- s_avail
        draw_idx <- seq_len(s_avail)
    } else {
        n_draws <- as.integer(n_draws)
        if (n_draws <= 0) {
            stop("n_draws must be a positive integer.", call. = FALSE)
        }
        if (n_draws <= s_avail) {
            draw_idx <- round(seq(1, s_avail, length.out = n_draws))
        } else {
            draw_idx <- sample(s_avail, n_draws, replace = TRUE)
        }
    }

    yrep <- matrix(0, nrow = n_draws, ncol = n_sec)
    sec_ids <- if ("sector_id" %in% names(data)) as.character(data$sector_id) else paste0("sector_", seq_len(n_sec))
    colnames(yrep) <- sec_ids

    has_p_all <- !is.null(object$mcmc_samples$p1_all_draws)
    has_mu_p1 <- !is.null(object$mcmc_samples$mu_p1)

    for (s in seq_len(n_draws)) {
        idx <- draw_idx[s]

        for (i in seq_len(n_sec)) {
            st <- as.character(data$survey_type[i])
            if (st == "Unsurveyed") {
                yrep[s, i] <- NA_real_
                next
            }

            N_i <- max(0, round(as.vector(sec_N_list[[i]])[idx]))

            # Extraction of p1
            p1 <- if (has_p_all) {
                as.vector(object$mcmc_samples$p1_all_draws[i, , ])[idx]
            } else if (has_mu_p1) {
                as.vector(object$mcmc_samples$mu_p1)[idx]
            } else {
                object$mean_p_hat
            }
            p1 <- min(0.9999, max(0.0001, p1))

            if (st == "SO") {
                # Single Observer
                det <- stats::rbinom(1, size = N_i, prob = p1)
                if (type == "count") {
                    yrep[s, i] <- det
                } else if (type == "obs1") {
                    yrep[s, i] <- det
                } else {
                    yrep[s, i] <- 0
                }
            } else {
                # Double Observer
                if (object$obs_type == "independent") {
                    p2 <- p1
                    pi_both <- p1 * p2
                    pi_s1 <- p1 * (1 - p2)
                    pi_s2 <- (1 - p1) * p2
                    pi_miss <- (1 - p1) * (1 - p2)

                    if (N_i > 0) {
                        sim <- stats::rmultinom(1, size = N_i, prob = c(pi_s1, pi_s2, pi_both, pi_miss))
                        s1_c <- sim[1]
                        s2_c <- sim[2]
                        both_c <- sim[3]
                    } else {
                        s1_c <- 0
                        s2_c <- 0
                        both_c <- 0
                    }

                    if (type == "count") {
                        yrep[s, i] <- s1_c + s2_c + both_c
                    } else if (type == "both") {
                        yrep[s, i] <- both_c
                    } else if (type == "obs1") {
                        yrep[s, i] <- s1_c + both_c
                    } else if (type == "obs2") {
                        yrep[s, i] <- s2_c + both_c
                    }
                } else {
                    # Dependent Observer (Removal)
                    x1 <- stats::rbinom(1, size = N_i, prob = p1)
                    x2 <- stats::rbinom(1, size = max(0, N_i - x1), prob = p1)

                    if (type == "count") {
                        yrep[s, i] <- x1 + x2
                    } else if (type == "both" || type == "obs1") {
                        yrep[s, i] <- x1
                    } else if (type == "obs2") {
                        yrep[s, i] <- x2
                    }
                }
            }
        }
    }

    # Extract observed vector
    y_obs <- if (type == "count") {
        data$total_count
    } else if (type == "both") {
        ifelse(is.na(data$groups_both), 0, data$groups_both)
    } else if (type == "obs1") {
        b <- ifelse(is.na(data$groups_both), 0, data$groups_both)
        s1 <- ifelse(is.na(data$groups_obs1), 0, data$groups_obs1)
        ifelse(data$survey_type == "SO", data$total_count, b + s1)
    } else if (type == "obs2") {
        b <- ifelse(is.na(data$groups_both), 0, data$groups_both)
        s2 <- ifelse(is.na(data$groups_obs2), 0, data$groups_obs2)
        ifelse(data$survey_type == "SO", 0, b + s2)
    }
    names(y_obs) <- sec_ids

    if (drop_unsurveyed && any(data$survey_type == "Unsurveyed")) {
        surv_idx <- which(data$survey_type != "Unsurveyed")
        yrep <- yrep[, surv_idx, drop = FALSE]
        y_obs <- y_obs[surv_idx]
    }

    attr(yrep, "y") <- y_obs
    attr(yrep, "target") <- type
    class(yrep) <- c("doas_yrep", "matrix")
    yrep
}

#' Graphical Posterior Predictive Checks for DOAS Bayesian Models
#'
#' Leverages the \code{bayesplot} package to produce graphical posterior predictive
#' checks (PPC) for \code{\link{doas}} Bayesian models. Compares observed survey
#' counts \eqn{y} against replicated datasets \eqn{y_{\text{rep}}} generated by
#' \code{\link{posterior_predict}}.
#'
#' @param object A \code{doas} model object fitted with \code{method = "Bayes"}.
#' @param type Character string specifying the plot type:
#'   \itemize{
#'     \item \code{"dens_overlay"}: Density overlay of observed data and multiple \eqn{y_{\text{rep}}} draws.
#'     \item \code{"stat"}: Distribution of a test statistic across \eqn{y_{\text{rep}}} with the observed statistic marked.
#'     \item \code{"stat_2d"}: Joint scatterplot of two test statistics across \eqn{y_{\text{rep}}}.
#'     \item \code{"bars"}: Bar plot comparing observed vs simulated count frequencies (ideal for zero-inflation).
#'     \item \code{"intervals"}: Medians and central posterior predictive intervals by sector.
#'     \item \code{"ribbon"}: Ribbon plot of predictive intervals across sectors.
#'     \item \code{"rootogram"}: Hanging rootogram comparing observed vs predicted counts.
#'   }
#' @param n_draws Number of posterior predictive draws to display (for types like \code{"dens_overlay"}). Default 50.
#' @param stat Character string or function for the test statistic when \code{type = "stat"}. Default \code{"mean"}. Common choices: \code{"mean"}, \code{"sd"}, \code{"max"}, \code{"min"}, or custom function.
#' @param target Character string specifying the observation target: \code{"count"} (default), \code{"both"}, \code{"obs1"}, \code{"obs2"}.
#' @param ... Additional arguments passed to the underlying \code{bayesplot::ppc_*} function.
#'
#' @return A \code{ggplot} object produced by \code{bayesplot}.
#' @rdname pp_check
#' @export
pp_check <- function(object, ...) {
    UseMethod("pp_check")
}

#' @rdname pp_check
#' @export
pp_check.doas <- function(
    object,
    type = c("dens_overlay", "stat", "stat_2d", "bars", "intervals", "ribbon", "rootogram"),
    n_draws = 50,
    stat = "mean",
    target = c("count", "both", "obs1", "obs2"),
    ...
) {
    if (!inherits(object, "doas")) {
        stop("object must be of class 'doas'.", call. = FALSE)
    }
    if (object$method != "Bayes") {
        stop("pp_check() is only supported for Bayesian models (method = 'Bayes').", call. = FALSE)
    }
    if (!requireNamespace("bayesplot", quietly = TRUE)) {
        stop("Package 'bayesplot' is required for pp_check(). Please install it with: install.packages('bayesplot')", call. = FALSE)
    }

    type <- match.arg(type)
    target <- match.arg(target)

    # Generate yrep
    n_sims <- max(n_draws, 100)
    yrep <- posterior_predict(object, n_draws = n_sims, type = target)
    y <- attr(yrep, "y")

    # Subset draws for plotting if n_draws < nrow(yrep)
    draws_subset <- if (n_draws < nrow(yrep)) yrep[seq_len(n_draws), , drop = FALSE] else yrep

    switch(
        type,
        dens_overlay = bayesplot::ppc_dens_overlay(y = y, yrep = draws_subset, ...),
        stat = bayesplot::ppc_stat(y = y, yrep = yrep, stat = stat, ...),
        stat_2d = bayesplot::ppc_stat_2d(y = y, yrep = yrep, stat = stat, ...),
        bars = bayesplot::ppc_bars(y = y, yrep = yrep, ...),
        intervals = bayesplot::ppc_intervals(y = y, yrep = yrep, ...),
        ribbon = bayesplot::ppc_ribbon(y = y, yrep = yrep, ...),
        rootogram = bayesplot::ppc_rootogram(y = y, yrep = yrep, ...)
    )
}

#' Convert DOAS Bayesian MCMC Draws to Array
#'
#' Formats posterior MCMC draws of model parameters into a 3D array of dimensions
#' \code{[iterations, chains, parameters]}, compatible with MCMC diagnostic functions
#' in \code{bayesplot} (e.g. \code{mcmc_trace}, \code{mcmc_intervals}).
#'
#' @param x A fitted \code{doas} model object with \code{method = "Bayes"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A 3-dimensional array with dimensions \code{[iterations, chains, parameters]}.
#' @export
as.array.doas <- function(x, ...) {
    if (!inherits(x, "doas")) {
        stop("x must be of class 'doas'.", call. = FALSE)
    }
    if (x$method != "Bayes" || is.null(x$mcmc_samples)) {
        stop("as.array() is only supported for Bayesian models (method = 'Bayes').", call. = FALSE)
    }

    samples <- x$mcmc_samples
    tot_mat <- samples$total_N
    if (is.null(tot_mat)) {
        stop("No MCMC samples available in the doas object.", call. = FALSE)
    }
    n_keep <- nrow(tot_mat)
    n_chains <- ncol(tot_mat)

    param_list <- list(total_abundance = tot_mat)

    if (!is.null(samples$mu_p1)) {
        param_list$p_hat <- samples$mu_p1
    }

    if (!is.null(samples$p_beta)) {
        beta_arr <- samples$p_beta
        p_names <- if (!is.null(x$bayes_coef) && "parameter" %in% names(x$bayes_coef)) {
            paste0("p_", x$bayes_coef$parameter)
        } else {
            paste0("p_beta[", seq_len(dim(beta_arr)[3]), "]")
        }
        for (j in seq_len(dim(beta_arr)[3])) {
            param_list[[p_names[j]]] <- beta_arr[, , j]
        }
    }

    if (!is.null(samples$n_gamma)) {
        gamma_arr <- samples$n_gamma
        g_names <- if (!is.null(x$bayes_n_coef) && "parameter" %in% names(x$bayes_n_coef)) {
            paste0("n_", x$bayes_n_coef$parameter)
        } else {
            paste0("n_gamma[", seq_len(dim(gamma_arr)[3]), "]")
        }
        for (j in seq_len(dim(gamma_arr)[3])) {
            param_list[[g_names[j]]] <- gamma_arr[, , j]
        }
    }

    n_params <- length(param_list)
    out_arr <- array(
        0,
        dim = c(n_keep, n_chains, n_params),
        dimnames = list(
            iterations = NULL,
            chains = paste0("chain:", seq_len(n_chains)),
            parameters = names(param_list)
        )
    )

    for (p_idx in seq_len(n_params)) {
        out_arr[, , p_idx] <- param_list[[p_idx]]
    }

    out_arr
}

#' Convert DOAS Bayesian MCMC Draws to Matrix
#'
#' @param x A fitted \code{doas} model object with \code{method = "Bayes"}.
#' @param ... Additional arguments passed to \code{\link{as.array.doas}}.
#'
#' @return A matrix with rows corresponding to pooled draws and columns to parameters.
#' @export
as.matrix.doas <- function(x, ...) {
    arr <- as.array(x, ...)
    d <- dim(arr)
    n_draws <- d[1] * d[2]
    n_params <- d[3]
    mat <- matrix(0, nrow = n_draws, ncol = n_params)
    colnames(mat) <- dimnames(arr)$parameters
    for (p in seq_len(n_params)) {
        mat[, p] <- as.vector(arr[, , p])
    }
    mat
}
