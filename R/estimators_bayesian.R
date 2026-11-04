#' Compute Gelman-Rubin R-hat and Effective Sample Size (ESS) from MCMC chains
#'
#' @param mat Numeric matrix of posterior draws with dimensions `n_keep x chains`.
#' @return A list with `rhat` and `ess`.
#' @importFrom stats var cor
#' @keywords internal
#' @noRd
compute_mcmc_diagnostics <- function(mat) {
    if (!is.matrix(mat)) {
        mat <- as.matrix(mat)
    }
    n_keep <- nrow(mat)
    chains <- ncol(mat)

    if (n_keep < 4 || chains < 1) {
        return(list(rhat = 1.0, ess = n_keep * chains))
    }

    chain_vars <- apply(mat, 2, stats::var)
    W <- mean(chain_vars, na.rm = TRUE)

    if (is.na(W) || W <= 1e-12) {
        return(list(rhat = 1.0, ess = n_keep * chains))
    }

    if (chains > 1) {
        chain_means <- colMeans(mat, na.rm = TRUE)
        B <- n_keep * stats::var(chain_means, na.rm = TRUE)
        V_hat <- ((n_keep - 1) / n_keep) * W + (B / n_keep)
        rhat <- sqrt(max(V_hat / W, 1.0))
    } else {
        rhat <- 1.0
    }

    rho1_vec <- vapply(
        seq_len(chains),
        function(c_idx) {
            x <- mat[, c_idx]
            if (stats::var(x) > 1e-12) {
                r1 <- stats::cor(x[-n_keep], x[-1])
                if (is.na(r1)) 0 else max(0, min(0.999, r1))
            } else {
                0
            }
        },
        numeric(1)
    )
    rho1 <- mean(rho1_vec)
    ess <- round((n_keep * chains) * (1 - rho1) / (1 + rho1))
    ess <- max(1, min(n_keep * chains, ess))

    list(rhat = rhat, ess = ess)
}

#' Run Multi-Chain Gibbs Sampler for Independent Double Observers
#'
#' @keywords internal
#' @noRd
run_gibbs_independent <- function(
    b_act,
    s1_act,
    s2_act,
    chains,
    iter,
    burnin,
    prior_alpha,
    prior_beta
) {
    n_act <- length(b_act)
    n_keep <- iter - burnin
    D_act <- b_act + s1_act + s2_act

    b_mat <- matrix(b_act, nrow = n_act, ncol = chains)
    s1_mat <- matrix(s1_act, nrow = n_act, ncol = chains)
    s2_mat <- matrix(s2_act, nrow = n_act, ncol = chains)
    D_mat <- b_mat + s1_mat + s2_mat

    b_sum <- sum(b_act)
    s1_sum <- sum(s1_act)
    s2_sum <- sum(s2_act)
    D_sum <- sum(D_act)

    init_multipliers <- seq(0.1, 1.0, length.out = chains)
    f0_mat <- outer(D_act, init_multipliers, function(d, m) round(d * m))
    f0_pool <- round(D_sum * init_multipliers)

    b0_sec <- 1 / max(mean(D_act) * 1.5, 2)
    b0_pool <- 1 / max(D_sum * 1.5, 4)
    kappa <- if (n_act > 1) 6 else 0

    G_draws <- lapply(seq_len(n_act), function(i) {
        matrix(0, nrow = n_keep, ncol = chains)
    })
    p1_draws <- lapply(seq_len(n_act), function(i) {
        matrix(0, nrow = n_keep, ncol = chains)
    })
    mu_p1_mat <- matrix(0.5, nrow = n_keep, ncol = chains)

    for (t in seq_len(iter)) {
        mu_p1 <- stats::rbeta(
            chains,
            b_sum + s1_sum + prior_alpha,
            s2_sum + f0_pool + prior_beta
        )
        mu_p2 <- stats::rbeta(
            chains,
            b_sum + s2_sum + prior_alpha,
            s1_sum + f0_pool + prior_beta
        )
        p_ov_pool <- mu_p1 + mu_p2 - (mu_p1 * mu_p2)
        f0_pool <- stats::rnbinom(
            chains,
            size = max(D_sum, 1),
            prob = (p_ov_pool + b0_pool) / (1 + b0_pool)
        )

        if (n_act > 1) {
            a1_p <- rep(mu_p1 * kappa, each = n_act)
            b1_p <- rep((1 - mu_p1) * kappa, each = n_act)
            a2_p <- rep(mu_p2 * kappa, each = n_act)
            b2_p <- rep((1 - mu_p2) * kappa, each = n_act)
        } else {
            a1_p <- prior_alpha
            b1_p <- prior_beta
            a2_p <- prior_alpha
            b2_p <- prior_beta
        }

        p1_mat <- matrix(
            stats::rbeta(
                n_act * chains,
                b_mat + s1_mat + a1_p,
                s2_mat + f0_mat + b1_p
            ),
            nrow = n_act,
            ncol = chains
        )
        p2_mat <- matrix(
            stats::rbeta(
                n_act * chains,
                b_mat + s2_mat + a2_p,
                s1_mat + f0_mat + b2_p
            ),
            nrow = n_act,
            ncol = chains
        )
        p_ov_mat <- p1_mat + p2_mat - (p1_mat * p2_mat)

        f0_mat <- matrix(
            stats::rnbinom(
                n_act * chains,
                size = pmax(D_mat, 1),
                prob = (p_ov_mat + b0_sec) / (1 + b0_sec)
            ),
            nrow = n_act,
            ncol = chains
        )

        if (t > burnin) {
            keep_idx <- t - burnin
            G_cur <- D_mat + f0_mat
            p_mean_mat <- 0.5 * (p1_mat + p2_mat)
            for (k in seq_len(n_act)) {
                G_draws[[k]][keep_idx, ] <- G_cur[k, ]
                p1_draws[[k]][keep_idx, ] <- p_mean_mat[k, ]
            }
            mu_p1_mat[keep_idx, ] <- 0.5 * (mu_p1 + mu_p2)
        }
    }

    list(G_draws = G_draws, p1_draws = p1_draws, mu_p1_mat = mu_p1_mat)
}

#' Run Multi-Chain Gibbs Sampler for Dependent Double Observers (2-Pass Removal)
#'
#' Models 2-pass removal where Primary Observer detects \eqn{x_1 = b + s_1} out of \eqn{G}
#' groups with probability \eqn{p}, and Secondary Observer detects \eqn{x_2 = s_2} out of
#' the remaining \eqn{G - x_1} groups with probability \eqn{p}.
#' The exact conjugate Beta full conditional for \eqn{p} given latent missed groups \eqn{f_0} is
#' \eqn{\text{Beta}(x_1 + x_2 + \alpha, x_2 + 2 f_0 + \beta)}, which remains strictly valid in (0, 1)
#' even when \eqn{x_1 \le x_2}.
#'
#' @keywords internal
#' @noRd
run_gibbs_dependent <- function(
    b_act,
    s1_act,
    s2_act,
    chains,
    iter,
    burnin,
    prior_alpha,
    prior_beta
) {
    n_act <- length(b_act)
    n_keep <- iter - burnin
    x1_act <- b_act + s1_act
    x2_act <- s2_act
    D_act <- x1_act + x2_act

    x1_mat <- matrix(x1_act, nrow = n_act, ncol = chains)
    x2_mat <- matrix(x2_act, nrow = n_act, ncol = chains)
    D_mat <- x1_mat + x2_mat

    x1_sum <- sum(x1_act)
    x2_sum <- sum(x2_act)
    D_sum <- sum(D_act)

    init_multipliers <- seq(0.1, 1.0, length.out = chains)
    f0_mat <- outer(D_act, init_multipliers, function(d, m) round(d * m))
    f0_pool <- round(D_sum * init_multipliers)

    b0_sec <- 1 / max(mean(D_act) * 1.5, 2)
    b0_pool <- 1 / max(D_sum * 1.5, 4)
    kappa <- if (n_act > 1) 6 else 0

    G_draws <- lapply(seq_len(n_act), function(i) {
        matrix(0, nrow = n_keep, ncol = chains)
    })
    p1_draws <- lapply(seq_len(n_act), function(i) {
        matrix(0, nrow = n_keep, ncol = chains)
    })
    mu_p1_mat <- matrix(0.5, nrow = n_keep, ncol = chains)

    for (t in seq_len(iter)) {
        # Step 1: Pooled 2-pass removal Gibbs update
        # Likelihood: p^(x1 + x2) * (1 - p)^(x2 + 2 * f0)
        mu_p1 <- stats::rbeta(
            chains,
            D_sum + prior_alpha,
            x2_sum + 2 * f0_pool + prior_beta
        )
        p_ov_pool <- 1 - ((1 - mu_p1)^2)
        f0_pool <- stats::rnbinom(
            chains,
            size = max(D_sum, 1),
            prob = (p_ov_pool + b0_pool) / (1 + b0_pool)
        )

        if (n_act > 1) {
            a1_p <- rep(mu_p1 * kappa, each = n_act)
            b1_p <- rep((1 - mu_p1) * kappa, each = n_act)
        } else {
            a1_p <- prior_alpha
            b1_p <- prior_beta
        }

        # Step 2: Sector-specific detection probability p_i given f0_i
        p1_mat <- matrix(
            stats::rbeta(
                n_act * chains,
                D_mat + a1_p,
                x2_mat + 2 * f0_mat + b1_p
            ),
            nrow = n_act,
            ncol = chains
        )
        p_ov_mat <- 1 - ((1 - p1_mat)^2)

        # Step 3: Sector-specific latent missed groups f0_i given p_i
        f0_mat <- matrix(
            stats::rnbinom(
                n_act * chains,
                size = pmax(D_mat, 1),
                prob = (p_ov_mat + b0_sec) / (1 + b0_sec)
            ),
            nrow = n_act,
            ncol = chains
        )

        if (t > burnin) {
            keep_idx <- t - burnin
            G_cur <- D_mat + f0_mat
            for (k in seq_len(n_act)) {
                G_draws[[k]][keep_idx, ] <- G_cur[k, ]
                p1_draws[[k]][keep_idx, ] <- p1_mat[k, ]
            }
            mu_p1_mat[keep_idx, ] <- mu_p1
        }
    }

    list(G_draws = G_draws, p1_draws = p1_draws, mu_p1_mat = mu_p1_mat)
}

#' Calculate Bayesian Multi-Chain MCMC Double Observer Estimate
#'
#' Uses a multi-chain Collapsed Poisson-Gamma-Beta Gibbs Sampler (Data Augmentation MCMC)
#' to sample the joint posterior distribution of latent undetected groups \eqn{f_{0,i}}
#' and observer detection probabilities across sectors. Supports both `"independent"`
#' (two independent observers) and `"dependent"` (2-pass primary/secondary removal) protocols.
#'
#' @param both Integer vector of numbers of groups seen by both observers
#' @param obs1_only Integer vector of numbers of groups seen by Observer 1 only
#' @param obs2_only Integer vector of numbers of groups seen by Observer 2 only
#' @param obs_type Character string specifying observer protocol: `"independent"` (default) or `"dependent"`.
#' @param iter Integer number of total MCMC iterations per chain (default 2000)
#' @param chains Integer number of independent Markov chains (default 4)
#' @param burnin Integer number of warm-up/burn-in iterations per chain (default `floor(iter / 4)`)
#' @param prior_alpha Numeric Beta hyperprior shape1 parameter for detection probability (default 1)
#' @param prior_beta Numeric Beta hyperprior shape2 parameter for detection probability (default 0.5)
#'
#' @return A list containing the posterior `mean`, `median`, `var`, 95% credible interval
#'   vectors (`ci_lower`, `ci_upper`), multi-chain convergence diagnostics (`rhat`, `ess`),
#'   and full posterior sample arrays (`G_draws`, `p1_draws`, `mu_p1_draws`).
#' @importFrom stats rbeta rnbinom var sd cor median quantile
#' @export
bayesian_conjugate_estimate <- function(
    both,
    obs1_only,
    obs2_only,
    obs_type = c("independent", "dependent"),
    iter = 2000,
    chains = 4,
    burnin = NULL,
    prior_alpha = 1,
    prior_beta = 0.5
) {
    obs_type <- match.arg(obs_type)
    if (is.null(burnin)) {
        burnin <- floor(iter / 4)
    }
    if (burnin >= iter) {
        stop("burnin must be strictly less than iter.")
    }
    n_keep <- iter - burnin
    n_sectors <- length(both)
    total_draws <- n_keep * chains

    est_mean <- numeric(n_sectors)
    est_median <- numeric(n_sectors)
    est_var <- numeric(n_sectors)
    ci_lower <- numeric(n_sectors)
    ci_upper <- numeric(n_sectors)
    rhat_vec <- rep(1.0, n_sectors)
    ess_vec <- rep(total_draws, n_sectors)

    b_clean <- ifelse(is.na(both), 0, both)
    s1_clean <- ifelse(is.na(obs1_only), 0, obs1_only)
    s2_clean <- ifelse(is.na(obs2_only), 0, obs2_only)
    D_vec <- b_clean + s1_clean + s2_clean
    valid_mask <- !is.na(both) & (D_vec > 0)

    G_draws_list <- lapply(seq_len(n_sectors), function(i) {
        matrix(0, nrow = n_keep, ncol = chains)
    })
    p1_draws_list <- lapply(seq_len(n_sectors), function(i) {
        matrix(0, nrow = n_keep, ncol = chains)
    })
    mu_p1_mat <- matrix(0.5, nrow = n_keep, ncol = chains)

    if (!any(valid_mask)) {
        return(list(
            mean = est_mean,
            median = est_median,
            var = est_var,
            ci_lower = ci_lower,
            ci_upper = ci_upper,
            rhat = rhat_vec,
            ess = ess_vec,
            G_draws = G_draws_list,
            p1_draws = p1_draws_list,
            mu_p1_draws = mu_p1_mat,
            mu_p1_rhat = 1.0,
            mu_p1_ess = total_draws
        ))
    }

    idx_valid <- which(valid_mask)
    gibbs_sampler <- if (obs_type == "independent") {
        run_gibbs_independent
    } else {
        run_gibbs_dependent
    }

    gibbs_out <- gibbs_sampler(
        b_act = b_clean[idx_valid],
        s1_act = s1_clean[idx_valid],
        s2_act = s2_clean[idx_valid],
        chains = chains,
        iter = iter,
        burnin = burnin,
        prior_alpha = prior_alpha,
        prior_beta = prior_beta
    )

    mu_p1_mat <- gibbs_out$mu_p1_mat

    for (k in seq_along(idx_valid)) {
        orig_i <- idx_valid[k]
        g_mat <- gibbs_out$G_draws[[k]]
        G_draws_list[[orig_i]] <- g_mat
        p1_draws_list[[orig_i]] <- gibbs_out$p1_draws[[k]]

        g_vec <- as.vector(g_mat)
        est_mean[orig_i] <- mean(g_vec)
        est_median[orig_i] <- stats::median(g_vec)
        est_var[orig_i] <- stats::var(g_vec)
        qs <- stats::quantile(g_vec, probs = c(0.025, 0.975), names = FALSE)
        ci_lower[orig_i] <- max(D_vec[orig_i], qs[1])
        ci_upper[orig_i] <- max(D_vec[orig_i], qs[2])

        diag_i <- compute_mcmc_diagnostics(g_mat)
        rhat_vec[orig_i] <- diag_i$rhat
        ess_vec[orig_i] <- diag_i$ess
    }

    mu_diag <- compute_mcmc_diagnostics(mu_p1_mat)

    list(
        mean = est_mean,
        median = est_median,
        var = est_var,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        rhat = rhat_vec,
        ess = ess_vec,
        G_draws = G_draws_list,
        p1_draws = p1_draws_list,
        mu_p1_draws = mu_p1_mat,
        mu_p1_rhat = mu_diag$rhat,
        mu_p1_ess = mu_diag$ess
    )
}

#' Fit Multi-Chain Bayesian Logistic Regression for Detection Covariates
#'
#' @importFrom stats model.matrix coef vcov rnorm dnorm runif quantile sd update
#' @keywords internal
#' @noRd
fit_bayesian_covariate_mcmc <- function(
    formula,
    do_data,
    all_data,
    obs_type,
    G_draws_list,
    glm_mod,
    prior_sd = 2.5
) {
    rhs_formula <- stats::update(formula, NULL ~ .)
    X_do <- stats::model.matrix(rhs_formula, data = do_data)
    X_all <- stats::model.matrix(rhs_formula, data = all_data)

    n_params <- ncol(X_do)
    n_do <- nrow(X_do)
    n_all <- nrow(X_all)
    n_keep <- nrow(G_draws_list[[1]])
    chains <- ncol(G_draws_list[[1]])

    y_obs <- if (obs_type == "independent") {
        do_data$groups_both + 0.5 * (do_data$groups_obs1 + do_data$groups_obs2)
    } else {
        do_data$groups_both + do_data$groups_obs1
    }

    beta_init <- stats::coef(glm_mod)
    beta_init[is.na(beta_init)] <- 0

    cov_init <- tryCatch(
        stats::vcov(glm_mod),
        error = function(e) diag(0.25, n_params)
    )
    if (any(is.na(cov_init)) || any(is.nan(cov_init))) {
        cov_init <- diag(0.25, n_params)
    } else {
        cov_init <- cov_init + diag(1e-4, n_params)
    }

    scale_factor <- (2.38^2) / max(n_params, 1)
    chol_prop <- tryCatch(
        chol(cov_init * scale_factor),
        error = function(e) diag(0.3, n_params)
    )

    log_post_fn <- function(beta_vec, G_vec) {
        eta <- pmax(-20, pmin(20, as.vector(X_do %*% beta_vec)))
        ll <- sum(y_obs * eta - G_vec * log1p(exp(eta)))
        lp <- sum(stats::dnorm(beta_vec, mean = 0, sd = prior_sd, log = TRUE))
        ll + lp
    }

    beta_draws <- array(0, dim = c(n_keep, chains, n_params))
    p1_all_draws <- array(0.5, dim = c(n_all, n_keep, chains))
    warmup_mh <- 200

    for (c_idx in seq_len(chains)) {
        beta_cur <- beta_init +
            as.vector(t(chol_prop) %*% stats::rnorm(n_params, 0, 0.5))
        G_init_vec <- vapply(
            seq_len(n_do),
            function(i) max(y_obs[i], G_draws_list[[i]][1, c_idx]),
            numeric(1)
        )
        lp_cur <- log_post_fn(beta_cur, G_init_vec)

        for (w in seq_len(warmup_mh)) {
            beta_prop <- beta_cur +
                as.vector(t(chol_prop) %*% stats::rnorm(n_params))
            lp_prop <- log_post_fn(beta_prop, G_init_vec)
            if (
                is.finite(lp_prop) && log(stats::runif(1)) < (lp_prop - lp_cur)
            ) {
                beta_cur <- beta_prop
                lp_cur <- lp_prop
            }
        }

        for (s in seq_len(n_keep)) {
            G_s_vec <- vapply(
                seq_len(n_do),
                function(i) max(y_obs[i], G_draws_list[[i]][s, c_idx]),
                numeric(1)
            )
            lp_cur <- log_post_fn(beta_cur, G_s_vec)
            beta_prop <- beta_cur +
                as.vector(t(chol_prop) %*% stats::rnorm(n_params))
            lp_prop <- log_post_fn(beta_prop, G_s_vec)
            if (
                is.finite(lp_prop) && log(stats::runif(1)) < (lp_prop - lp_cur)
            ) {
                beta_cur <- beta_prop
                lp_cur <- lp_prop
            }
            beta_draws[s, c_idx, ] <- beta_cur
            eta_all <- pmax(-20, pmin(20, as.vector(X_all %*% beta_cur)))
            p1_all_draws[, s, c_idx] <- 1 / (1 + exp(-eta_all))
        }
    }

    param_names <- colnames(X_do)
    beta_summary <- data.frame(
        Mean = numeric(n_params),
        SD = numeric(n_params),
        CI_2.5 = numeric(n_params),
        CI_97.5 = numeric(n_params),
        Rhat = numeric(n_params),
        ESS = numeric(n_params),
        row.names = param_names
    )

    for (p_idx in seq_len(n_params)) {
        p_mat <- matrix(beta_draws[, , p_idx], nrow = n_keep, ncol = chains)
        p_vec <- as.vector(p_mat)
        diag_p <- compute_mcmc_diagnostics(p_mat)
        qs <- stats::quantile(p_vec, probs = c(0.025, 0.975), names = FALSE)
        beta_summary$Mean[p_idx] <- mean(p_vec)
        beta_summary$SD[p_idx] <- stats::sd(p_vec)
        beta_summary$CI_2.5[p_idx] <- qs[1]
        beta_summary$CI_97.5[p_idx] <- qs[2]
        beta_summary$Rhat[p_idx] <- diag_p$rhat
        beta_summary$ESS[p_idx] <- diag_p$ess
    }

    list(
        p1_all_draws = p1_all_draws,
        beta_summary = beta_summary,
        beta_draws = beta_draws
    )
}

#' Fit Multi-Chain Bayesian Log-Linear Regression for Sector Abundance Covariates (n_formula)
#'
#' @importFrom stats model.matrix model.frame model.offset delete.response terms coef vcov rnorm dnorm dpois runif quantile sd update quasipoisson glm
#' @keywords internal
#' @noRd
fit_bayesian_abundance_mcmc <- function(
    n_formula,
    surv_data,
    all_data,
    N_surv_draws_list,
    n_glm_mod,
    prior_sd = 2.5
) {
    rhs_terms <- stats::delete.response(stats::terms(n_glm_mod))
    mf_surv <- stats::model.frame(rhs_terms, data = surv_data, xlev = n_glm_mod$xlevels)
    Z_surv <- stats::model.matrix(rhs_terms, data = mf_surv)
    off_surv <- stats::model.offset(mf_surv)
    if (is.null(off_surv)) off_surv <- rep(0, nrow(Z_surv))

    mf_all <- stats::model.frame(rhs_terms, data = all_data, xlev = n_glm_mod$xlevels)
    Z_all <- stats::model.matrix(rhs_terms, data = mf_all)
    off_all <- stats::model.offset(mf_all)
    if (is.null(off_all)) off_all <- rep(0, nrow(Z_all))

    n_params <- ncol(Z_surv)
    n_surv <- nrow(Z_surv)
    n_all <- nrow(Z_all)
    n_keep <- nrow(N_surv_draws_list[[1]])
    chains <- ncol(N_surv_draws_list[[1]])

    gamma_init <- stats::coef(n_glm_mod)
    gamma_init[is.na(gamma_init)] <- 0

    cov_init <- tryCatch(
        stats::vcov(n_glm_mod),
        error = function(e) diag(0.1, n_params)
    )
    if (
        any(is.na(cov_init)) ||
            any(is.nan(cov_init)) ||
            nrow(cov_init) != n_params
    ) {
        cov_init <- diag(0.1, n_params)
    } else {
        cov_init <- cov_init + diag(1e-4, n_params)
    }

    scale_factor <- (2.38^2) / max(n_params, 1)
    chol_prop <- tryCatch(
        chol(cov_init * scale_factor),
        error = function(e) diag(0.15, n_params)
    )

    log_post_fn <- function(gamma_vec, N_vec) {
        eta <- pmax(-15, pmin(15, as.vector(Z_surv %*% gamma_vec) + off_surv))
        ll <- sum(N_vec * eta - exp(eta))
        lp <- sum(stats::dnorm(gamma_vec, mean = 0, sd = prior_sd, log = TRUE))
        ll + lp
    }

    gamma_draws <- array(0, dim = c(n_keep, chains, n_params))
    lambda_all_draws <- array(1, dim = c(n_all, n_keep, chains))
    n_log_lik <- matrix(0, nrow = n_surv, ncol = n_keep * chains)
    warmup_mh <- 200

    for (c_idx in seq_len(chains)) {
        gamma_cur <- gamma_init +
            as.vector(t(chol_prop) %*% stats::rnorm(n_params, 0, 0.25))
        N_init_vec <- vapply(
            seq_len(n_surv),
            function(i) max(0, N_surv_draws_list[[i]][1, c_idx]),
            numeric(1)
        )
        lp_cur <- log_post_fn(gamma_cur, N_init_vec)

        for (w in seq_len(warmup_mh)) {
            gamma_prop <- gamma_cur +
                as.vector(t(chol_prop) %*% stats::rnorm(n_params))
            lp_prop <- log_post_fn(gamma_prop, N_init_vec)
            if (is.finite(lp_prop) && log(stats::runif(1)) < (lp_prop - lp_cur)) {
                gamma_cur <- gamma_prop
                lp_cur <- lp_prop
            }
        }

        for (s in seq_len(n_keep)) {
            N_s_vec <- vapply(
                seq_len(n_surv),
                function(i) max(0, N_surv_draws_list[[i]][s, c_idx]),
                numeric(1)
            )
            lp_cur <- log_post_fn(gamma_cur, N_s_vec)
            gamma_prop <- gamma_cur +
                as.vector(t(chol_prop) %*% stats::rnorm(n_params))
            lp_prop <- log_post_fn(gamma_prop, N_s_vec)
            if (is.finite(lp_prop) && log(stats::runif(1)) < (lp_prop - lp_cur)) {
                gamma_cur <- gamma_prop
                lp_cur <- lp_prop
            }
            gamma_draws[s, c_idx, ] <- gamma_cur
            eta_all <- pmax(-15, pmin(15, as.vector(Z_all %*% gamma_cur) + off_all))
            lambda_all_draws[, s, c_idx] <- exp(eta_all)

            eta_surv <- pmax(-15, pmin(15, as.vector(Z_surv %*% gamma_cur) + off_surv))
            draw_col <- (c_idx - 1) * n_keep + s
            n_log_lik[, draw_col] <- stats::dpois(
                round(N_s_vec),
                lambda = pmax(1e-6, exp(eta_surv)),
                log = TRUE
            )
        }
    }

    param_names <- colnames(Z_surv)
    gamma_summary <- data.frame(
        Mean = numeric(n_params),
        SD = numeric(n_params),
        CI_2.5 = numeric(n_params),
        CI_97.5 = numeric(n_params),
        Rhat = numeric(n_params),
        ESS = numeric(n_params),
        row.names = param_names
    )

    for (p_idx in seq_len(n_params)) {
        p_mat <- matrix(gamma_draws[, , p_idx], nrow = n_keep, ncol = chains)
        p_vec <- as.vector(p_mat)
        diag_p <- compute_mcmc_diagnostics(p_mat)
        qs <- stats::quantile(p_vec, probs = c(0.025, 0.975), names = FALSE)
        gamma_summary$Mean[p_idx] <- mean(p_vec)
        gamma_summary$SD[p_idx] <- stats::sd(p_vec)
        gamma_summary$CI_2.5[p_idx] <- qs[1]
        gamma_summary$CI_97.5[p_idx] <- qs[2]
        gamma_summary$Rhat[p_idx] <- diag_p$rhat
        gamma_summary$ESS[p_idx] <- diag_p$ess
    }

    list(
        lambda_all_draws = lambda_all_draws,
        gamma_summary = gamma_summary,
        gamma_draws = gamma_draws,
        n_log_lik = n_log_lik
    )
}

#' Compute Bayesian Sector-Level and Total Posterior Abundance Estimates
#'
#' @keywords internal
#' @noRd
compute_bayesian_sectors <- function(
    data,
    do_data,
    obs_type,
    pooling,
    bayes_res,
    bayes_cov_res,
    n_formula = NULL,
    iter,
    burnin
) {
    n_sec <- nrow(data)
    n_keep <- nrow(bayes_res$mu_p1_draws)
    n_ch <- ncol(bayes_res$mu_p1_draws)
    total_N_draws_mat <- matrix(0, nrow = n_keep, ncol = n_ch)

    sector_estimates <- data.frame(
        sector_id = data$sector_id,
        total_count = data$total_count,
        est_groups = numeric(n_sec),
        est_abundance = numeric(n_sec),
        se_abundance = numeric(n_sec),
        ci_lower = numeric(n_sec),
        ci_upper = numeric(n_sec)
    )

    N_sec_draws_all <- vector("list", n_sec)

    for (i in seq_len(n_sec)) {
        row <- data[i, ]
        if (identical(row$survey_type, "Unsurveyed")) {
            next
        }
        mean_size_i <- row$mean_group_size

        if (pooling == "hybrid" && row$survey_type == "DO") {
            match_idx <- which(do_data$sector_id == row$sector_id)[1]
            G_sec_mat <- bayes_res$G_draws[[match_idx]]
            G_point <- mean(G_sec_mat)
        } else {
            primary_groups <- if (row$survey_type == "DO") {
                if (obs_type == "independent") {
                    row$groups_both + 0.5 * (row$groups_obs1 + row$groups_obs2)
                } else {
                    row$groups_both + row$groups_obs1
                }
            } else {
                row$detected_groups
            }

            p1_sec_mat <- if (!is.null(bayes_cov_res)) {
                matrix(
                    bayes_cov_res$p1_all_draws[i, , ],
                    nrow = n_keep,
                    ncol = n_ch
                )
            } else {
                bayes_res$mu_p1_draws
            }

            p_eff_mat <- pmax(p1_sec_mat, 0.05)
            G_expected_mat <- primary_groups / p_eff_mat
            G_point <- mean(G_expected_mat)

            if (primary_groups > 0) {
                f0_draw_mat <- matrix(
                    stats::rnbinom(
                        n_keep * n_ch,
                        size = max(1, round(primary_groups)),
                        prob = pmin(0.9999, pmax(0.0001, p_eff_mat))
                    ),
                    nrow = n_keep,
                    ncol = n_ch
                )
                f0_mean_mat <- round(primary_groups) * (1 - p_eff_mat) / p_eff_mat
                G_sec_mat <- G_expected_mat + (f0_draw_mat - f0_mean_mat)
            } else {
                G_sec_mat <- G_expected_mat
            }
        }

        N_sec_mat <- G_sec_mat * mean_size_i
        N_point <- G_point * mean_size_i
        N_sec_mat <- matrix(
            pmax(N_sec_mat - mean(N_sec_mat) + N_point, row$total_count),
            nrow = n_keep,
            ncol = n_ch
        )
        N_sec_draws_all[[i]] <- N_sec_mat
        total_N_draws_mat <- total_N_draws_mat + N_sec_mat

        N_sec_vec <- as.vector(N_sec_mat)
        sec_qs <- stats::quantile(
            N_sec_vec,
            probs = c(0.025, 0.975),
            names = FALSE
        )

        sector_estimates$est_groups[i] <- G_point
        sector_estimates$est_abundance[i] <- N_point
        sector_estimates$se_abundance[i] <- stats::sd(N_sec_vec)
        sector_estimates$ci_lower[i] <- max(row$total_count, sec_qs[1])
        sector_estimates$ci_upper[i] <- max(row$total_count, sec_qs[2])
    }

    n_model <- NULL
    bayes_n_res <- NULL
    surv_mask <- data$survey_type != "Unsurveyed"
    mean_gs_pool <- sum(data$total_count[surv_mask], na.rm = TRUE) /
        max(1, sum(data$detected_groups[surv_mask], na.rm = TRUE))

    if (!is.null(n_formula)) {
        if (!inherits(n_formula, "formula")) {
            stop("n_formula must be a valid R formula (e.g. ~ elevation + habitat)")
        }
        surv_df <- data[surv_mask, , drop = FALSE]
        surv_df$est_abundance <- sector_estimates$est_abundance[surv_mask]
        if (nrow(surv_df) < 2) {
            stop("At least 2 surveyed sectors are required to fit an abundance model (n_formula).")
        }
        n_mod_formula <- stats::update(n_formula, est_abundance ~ .)
        n_model <- tryCatch(
            eval(bquote(stats::glm(
                .(n_mod_formula),
                family = stats::quasipoisson(link = "log"),
                data = surv_df
            ))),
            error = function(e) {
                stop(paste("Failed to fit GLM for sector abundance (n_formula):", e$message))
            }
        )
        N_surv_draws_list <- N_sec_draws_all[which(surv_mask)]
        bayes_n_res <- fit_bayesian_abundance_mcmc(
            n_formula = n_formula,
            surv_data = surv_df,
            all_data = data,
            N_surv_draws_list = N_surv_draws_list,
            n_glm_mod = n_model
        )

        unsurv_idx <- which(!surv_mask)
        if (length(unsurv_idx) > 0) {
            for (u_i in unsurv_idx) {
                lam_u_mat <- matrix(
                    bayes_n_res$lambda_all_draws[u_i, , ],
                    nrow = n_keep,
                    ncol = n_ch
                )
                N_u_mat <- matrix(
                    stats::rpois(n_keep * n_ch, lambda = pmax(0, lam_u_mat)),
                    nrow = n_keep,
                    ncol = n_ch
                )
                total_N_draws_mat <- total_N_draws_mat + N_u_mat
                N_u_vec <- as.vector(N_u_mat)
                N_u_point <- mean(lam_u_mat)
                u_qs <- stats::quantile(N_u_vec, probs = c(0.025, 0.975), names = FALSE)

                sector_estimates$est_groups[u_i] <- N_u_point / max(mean_gs_pool, 1e-6)
                sector_estimates$est_abundance[u_i] <- N_u_point
                sector_estimates$se_abundance[u_i] <- stats::sd(N_u_vec)
                sector_estimates$ci_lower[u_i] <- u_qs[1]
                sector_estimates$ci_upper[u_i] <- u_qs[2]
            }
        }
    }

    tot_N_vec <- as.vector(total_N_draws_mat)
    total_abundance <- mean(tot_N_vec)
    total_variance <- stats::var(tot_N_vec)
    total_seen <- sum(data$total_count, na.rm = TRUE)
    tot_qs <- stats::quantile(tot_N_vec, probs = c(0.025, 0.975), names = FALSE)

    tot_diag <- compute_mcmc_diagnostics(total_N_draws_mat)
    mcmc_diagnostics <- list(
        chains = n_ch,
        iter = iter,
        burnin = burnin,
        post_warmup_draws = n_keep * n_ch,
        rhat_abundance = tot_diag$rhat,
        ess_abundance = tot_diag$ess,
        rhat_p_hat = bayes_res$mu_p1_rhat,
        ess_p_hat = bayes_res$mu_p1_ess,
        sector_rhat = bayes_res$rhat,
        sector_ess = bayes_res$ess
    )

    # Compute pointwise log-likelihood matrix (n_do x total_draws) for WAIC
    n_do <- nrow(do_data)
    total_draws <- n_keep * n_ch
    log_lik_mat <- matrix(0, nrow = n_do, ncol = total_draws)
    for (j in seq_len(n_do)) {
        sec_id <- do_data$sector_id[j]
        all_idx <- which(data$sector_id == sec_id)[1]
        y_j <- if (obs_type == "independent") {
            round(do_data$groups_both[j] + 0.5 * (do_data$groups_obs1[j] + do_data$groups_obs2[j]))
        } else {
            do_data$groups_both[j] + do_data$groups_obs1[j]
        }
        G_j_vec <- pmax(y_j, round(as.vector(bayes_res$G_draws[[j]])))

        p1_j_vec <- if (!is.null(bayes_cov_res)) {
            as.vector(bayes_cov_res$p1_all_draws[all_idx, , ])
        } else {
            as.vector(bayes_res$mu_p1_draws)
        }
        p1_j_vec <- pmin(1 - 1e-6, pmax(1e-6, p1_j_vec))
        log_lik_mat[j, ] <- stats::dbinom(
            x = y_j,
            size = G_j_vec,
            prob = p1_j_vec,
            log = TRUE
        )
    }

    list(
        sector_estimates = sector_estimates,
        total_abundance = total_abundance,
        total_variance = total_variance,
        se_abundance = sqrt(total_variance),
        ci_lower = max(total_seen, tot_qs[1]),
        ci_upper = max(total_seen, tot_qs[2]),
        mcmc_diagnostics = mcmc_diagnostics,
        log_lik = log_lik_mat,
        n_model = n_model,
        bayes_n_res = bayes_n_res,
        mean_gs_pool = mean_gs_pool
    )
}
