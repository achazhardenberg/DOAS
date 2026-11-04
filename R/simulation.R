# Helper to provide %||%
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Simulate Group Detections Within a Single Sector
#'
#' @param sim_groups List of simulated groups in the sector.
#' @param is_do Logical indicating whether this sector is a Double Observer (DO) sector.
#' @param det_p1 Base detection probability for Observer 1.
#' @param intra_group_det Mean factor for count errors.
#' @param det_group Base detection factor for group size scaling.
#' @param obs_effect Observer 2 relative efficiency.
#' @param detvar Standard deviation for detection probability across sectors.
#' @param det_p1_sec Resolved baseline detection probability for Observer 1 in this sector (accounting for p_formula / p_betas).
#' @return A list with `total_count`, `groups_both`, `groups_obs1`, `groups_obs2`, and `detected_groups`.
#' @importFrom stats rnorm runif plogis qlogis pnorm qnorm
#' @keywords internal
#' @noRd
simulate_sector_detection <- function(
    sim_groups,
    is_do,
    det_p1,
    intra_group_det,
    det_group,
    obs_effect,
    detvar,
    det_p1_sec = det_p1
) {
    if (length(sim_groups) == 0) {
        return(list(
            total_count = 0,
            groups_both = if (is_do) 0L else NA_integer_,
            groups_obs1 = if (is_do) 0L else NA_integer_,
            groups_obs2 = if (is_do) 0L else NA_integer_,
            detected_groups = 0L
        ))
    }

    group_sizes <- vapply(sim_groups, length, integer(1))
    grp_sqrt <- pmin(sqrt(group_sizes), sqrt(15))

    # Group detectabilities based on size (anchored at 1 = sqrt(1) as in original DOAS_script.R)
    grp_sqrt_anchored <- c(grp_sqrt, 1)
    rng <- range(grp_sqrt_anchored)
    if (rng[2] > rng[1]) {
        group_det <- (grp_sqrt - rng[1]) / (rng[2] - rng[1]) * (1 - det_group) + det_group
    } else {
        group_det <- rep(det_group, length(grp_sqrt))
    }

    detp_factor <- if (detvar > 0) {
        u_low <- stats::pnorm(0, mean = det_p1_sec, sd = detvar)
        u_high <- stats::pnorm(1, mean = det_p1_sec, sd = detvar)
        stats::qnorm(
            stats::runif(1, min = u_low, max = u_high),
            mean = det_p1_sec,
            sd = detvar
        )
    } else {
        det_p1_sec
    }

    p1 <- detp_factor * group_det
    obs1_det <- stats::runif(length(p1)) < p1

    if (is_do) {
        p2 <- p1 * obs_effect
        obs2_det <- stats::runif(length(p2)) < p2

        idx_both <- which(obs1_det & obs2_det)
        idx_obs1 <- which(obs1_det & !obs2_det)
        idx_obs2 <- which(!obs1_det & obs2_det)
        idx_any <- which(obs1_det | obs2_det)

        total_count <- if (length(idx_any) > 0) {
            err_factor <- stats::rnorm(1, intra_group_det, 0.05)
            round(sum(group_sizes[idx_any]) * err_factor, 0)
        } else {
            0
        }

        list(
            total_count = total_count,
            groups_both = length(idx_both),
            groups_obs1 = length(idx_obs1),
            groups_obs2 = length(idx_obs2),
            detected_groups = length(idx_any)
        )
    } else {
        # Block sectors are surveyed by rangers from the same observer pool (mean efficiency (1 + obs_effect) / 2)
        p_block <- p1 * ((1 + obs_effect) / 2)
        obs_blk_det <- stats::runif(length(p_block)) < p_block
        idx_blk <- which(obs_blk_det)
        total_count <- if (length(idx_blk) > 0) {
            err_factor <- stats::rnorm(1, intra_group_det, 0.05)
            round(sum(group_sizes[idx_blk]) * err_factor, 0)
        } else {
            0
        }

        list(
            total_count = total_count,
            groups_both = NA_integer_,
            groups_obs1 = NA_integer_,
            groups_obs2 = NA_integer_,
            detected_groups = length(idx_blk)
        )
    }
}

#' Simulate a Population
#'
#' Simulates a population of individuals split into groups across multiple sectors.
#' Supports sector-level abundance covariate modeling via \code{n_formula} and named
#' coefficient vectors \code{n_betas}, or user-provided sector covariates.
#'
#' @param pop_size Average total population size (Poisson distributed).
#' @param mean_group Average group size (Poisson distributed).
#' @param tot_areas Total number of survey sectors.
#' @param area_sizes Vector of relative probabilities/sizes for assigning groups to sectors.
#'   If NULL, probabilities are computed from \code{n_formula} and \code{n_betas} (or \code{n_cov_effect}).
#' @param n_formula Optional one-sided formula for sector abundance covariates (e.g. \code{~ elevation + forest}
#'   or \code{~ elevation * forest}).
#' @param n_betas Optional named numeric vector of log-linear coefficients for terms in \code{n_formula}
#'   (e.g., \code{c(elevation = 0.5, forest = -0.3)}).
#' @param covariates Optional data frame of sector-level covariates containing \code{tot_areas} rows.
#'   If omitted, any variables specified in \code{n_formula} are simulated as standard normal variates.
#' @param n_cov_effect Optional numeric log-linear slope for a single simulated \code{habitat} covariate
#'   (retained for backward compatibility; equivalent to \code{n_formula = ~ habitat, n_betas = c(habitat = n_cov_effect)}).
#'
#' @return A list containing:
#'   \item{tot_pop}{Total simulated individuals (scalar).}
#'   \item{sim_pop}{List of groups per sector.}
#'   \item{covariates}{Data frame of sector-level covariates (one row per sector).}
#'   \item{habitat}{Numeric vector of the \code{habitat} covariate (if present).}
#' @importFrom stats rpois rnorm runif model.matrix
#' @export
simulate_pop <- function(
    pop_size,
    mean_group,
    tot_areas,
    area_sizes = NULL,
    n_formula = NULL,
    n_betas = NULL,
    covariates = NULL,
    n_cov_effect = 0
) {
    # 1. Backward compatibility: if n_cov_effect != 0 and n_formula is NULL, set n_formula and n_betas
    if (n_cov_effect != 0 && is.null(n_formula)) {
        n_formula <- ~ habitat
        n_betas <- c(habitat = n_cov_effect)
    }

    # 2. Extract needed variable names from n_formula
    vars_needed <- unique(all.vars(n_formula))
    if (length(vars_needed) == 0 && n_cov_effect == 0 && is.null(covariates)) {
        # Default single habitat covariate for convenience
        vars_needed <- "habitat"
    }

    # 3. Build or validate covariates data frame
    if (is.null(covariates)) {
        cov_list <- list(sector_id = seq_len(tot_areas))
        for (v in vars_needed) {
            cov_list[[v]] <- stats::rnorm(tot_areas, mean = 0, sd = 1)
        }
        cov_df <- as.data.frame(cov_list, stringsAsFactors = FALSE)
    } else {
        if (!is.data.frame(covariates) || nrow(covariates) != tot_areas) {
            stop(sprintf("'covariates' must be a data frame with exactly %d rows.", tot_areas))
        }
        cov_df <- covariates
        if (!"sector_id" %in% names(cov_df)) {
            cov_df$sector_id <- seq_len(tot_areas)
        }
        for (v in vars_needed) {
            if (!v %in% names(cov_df)) {
                cov_df[[v]] <- stats::rnorm(tot_areas, mean = 0, sd = 1)
            }
        }
    }

    # 4. Compute area_sizes from n_formula and n_betas if not supplied
    if (is.null(area_sizes)) {
        if (!is.null(n_formula) && !is.null(n_betas)) {
            mm_n <- stats::model.matrix(n_formula, data = cov_df)
            active_betas <- n_betas[intersect(names(n_betas), colnames(mm_n))]
            if (length(active_betas) > 0) {
                eta_n <- as.vector(mm_n[, names(active_betas), drop = FALSE] %*% active_betas)
            } else {
                eta_n <- rep(0, tot_areas)
            }
            area_sizes <- exp(eta_n)
        } else if ("habitat" %in% names(cov_df) && n_cov_effect != 0) {
            area_sizes <- exp(n_cov_effect * cov_df[["habitat"]])
        } else {
            area_sizes <- rep(1, tot_areas)
        }
    }

    tot_pop <- stats::rpois(1, pop_size)
    if (tot_pop == 0) {
        empty_pop <- vector("list", tot_areas)
        attr(empty_pop, "covariates") <- cov_df
        if ("habitat" %in% names(cov_df)) attr(empty_pop, "habitat") <- cov_df[["habitat"]]
        return(list(
            tot_pop = 0,
            sim_pop = empty_pop,
            covariates = cov_df,
            habitat = cov_df[["habitat"]]
        ))
    }

    sim_list <- as.character(seq_len(tot_pop))
    sim_areas <- split(
        sample(sim_list),
        sample(
            seq_len(tot_areas),
            length(sim_list),
            replace = TRUE,
            prob = area_sizes
        )
    )

    sim_pop <- lapply(seq_len(tot_areas), function(i) {
        area_pop <- unlist(sim_areas[[as.character(i)]] %||% character(0))
        if (length(area_pop) > 0) {
            k <- max(1L, stats::rpois(1, length(area_pop) / mean_group))
            groups <- split(
                sample(area_pop),
                sample(seq_len(k), length(area_pop), replace = TRUE)
            )
            names(groups) <- paste0("g", seq_along(groups))
            groups
        } else {
            list()
        }
    })

    attr(sim_pop, "covariates") <- cov_df
    if ("habitat" %in% names(cov_df)) attr(sim_pop, "habitat") <- cov_df[["habitat"]]

    list(
        sim_pop = sim_pop,
        tot_pop = tot_pop,
        covariates = cov_df,
        habitat = cov_df[["habitat"]]
    )
}

#' Simulate a Survey and format for DOAS
#'
#' Simulates field observer detections across survey sectors with support for
#' detection covariates modeled via \code{p_formula} and named coefficients \code{p_betas}.
#'
#' @param sim_pop A simulated population from \code{\link{simulate_pop}} (either the returned list or \code{pop$sim_pop}).
#' @param det_p1 Base detection probability for Observer 1.
#' @param intra_group_det Mean factor for count errors (e.g., 0.9 for undercounting).
#' @param det_group Base detection factor for group size scaling.
#' @param obs_effect Observer 2 relative efficiency (e.g., 0.75 means Obs 2 is 75 percent as effective).
#' @param detvar Standard deviation for detection probability across sectors.
#' @param surv_areas Number of sectors to receive Double Observer surveys.
#' @param p_formula Optional one-sided formula for detection probability (e.g. \code{~ effort + elevation}).
#' @param p_betas Optional named numeric vector of logit-scale coefficients for \code{p_formula}
#'   (e.g., \code{c(effort = 0.6, elevation = -0.3)}).
#' @param covariates Optional data frame of sector-level covariates. If omitted, extracted from \code{sim_pop}.
#' @param cov_effect Optional numeric log-odds slope for a single simulated \code{covariate}
#'   (retained for backward compatibility; equivalent to \code{p_formula = ~ covariate, p_betas = c(covariate = cov_effect)}).
#'
#' @return A data frame with one row per sector, containing observed counts,
#'   Double Observer tallies, and all simulated sector-level covariates ready for input to \code{\link{doas}}.
#' @importFrom stats qlogis plogis model.matrix rnorm
#' @export
simulate_survey <- function(
    sim_pop,
    det_p1,
    intra_group_det,
    det_group,
    obs_effect,
    detvar,
    surv_areas,
    p_formula = NULL,
    p_betas = NULL,
    covariates = NULL,
    cov_effect = 0
) {
    # 1. Extract sim_pop list and existing covariates
    if (is.list(sim_pop) && "sim_pop" %in% names(sim_pop)) {
        pop_covs <- sim_pop$covariates
        hab_vec <- sim_pop$habitat
        sim_pop <- sim_pop$sim_pop
    } else {
        pop_covs <- attr(sim_pop, "covariates")
        hab_vec <- attr(sim_pop, "habitat")
    }

    tot_areas <- length(sim_pop)

    # 2. Backward compatibility: cov_effect != 0 maps to p_formula = ~ covariate
    if (cov_effect != 0 && is.null(p_formula)) {
        p_formula <- ~ covariate
        p_betas <- c(covariate = cov_effect)
    }

    # 3. Assemble full sector covariates table
    if (!is.null(covariates)) {
        cov_df <- covariates
    } else if (!is.null(pop_covs)) {
        cov_df <- pop_covs
    } else {
        cov_df <- data.frame(sector_id = seq_len(tot_areas))
        if (!is.null(hab_vec)) cov_df$habitat <- hab_vec
    }

    # Ensure all variables required by p_formula are present
    p_vars <- all.vars(p_formula)
    for (v in p_vars) {
        if (!v %in% names(cov_df)) {
            cov_df[[v]] <- stats::rnorm(tot_areas, mean = 0, sd = 1)
        }
    }
    # Backward compatibility: if cov_effect was used or 'covariate' is not present but expected
    if ("covariate" %in% p_vars && !"covariate" %in% names(cov_df)) {
        cov_df$covariate <- stats::rnorm(tot_areas, mean = 0, sd = 1)
    }
    if (cov_effect != 0 && !"covariate" %in% names(cov_df)) {
        cov_df$covariate <- stats::rnorm(tot_areas, mean = 0, sd = 1)
    }

    # 4. Calculate sector-specific det_p1_sec using p_formula and p_betas
    base_logit <- stats::qlogis(pmin(0.99, pmax(0.01, det_p1)))
    if (!is.null(p_formula) && !is.null(p_betas)) {
        mm_p <- stats::model.matrix(p_formula, data = cov_df)
        active_betas <- p_betas[intersect(names(p_betas), colnames(mm_p))]
        if (length(active_betas) > 0) {
            eta_p <- as.vector(mm_p[, names(active_betas), drop = FALSE] %*% active_betas)
        } else {
            eta_p <- rep(0, tot_areas)
        }
        det_p1_vec <- stats::plogis(base_logit + eta_p)
    } else {
        det_p1_vec <- rep(det_p1, tot_areas)
    }

    surv_areas <- min(surv_areas, tot_areas)
    do_sectors <- sample(seq_len(tot_areas), surv_areas, replace = FALSE)

    sec_list <- lapply(seq_len(tot_areas), function(m) {
        res <- simulate_sector_detection(
            sim_groups = sim_pop[[m]],
            is_do = m %in% do_sectors,
            det_p1 = det_p1,
            intra_group_det = intra_group_det,
            det_group = det_group,
            obs_effect = obs_effect,
            detvar = detvar,
            det_p1_sec = det_p1_vec[m]
        )
        data.frame(
            sector_id = m,
            survey_type = if (m %in% do_sectors) "DO" else "Block",
            total_count = res$total_count,
            groups_both = res$groups_both,
            groups_obs1 = res$groups_obs1,
            groups_obs2 = res$groups_obs2,
            detected_groups = res$detected_groups,
            stringsAsFactors = FALSE
        )
    })

    res_df <- do.call(rbind, sec_list)

    # Attach all covariate columns (except duplicate sector_id if present)
    cov_cols_to_add <- setdiff(names(cov_df), "sector_id")
    for (col_name in cov_cols_to_add) {
        res_df[[col_name]] <- cov_df[[col_name]]
    }

    res_df
}
