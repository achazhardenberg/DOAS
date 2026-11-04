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
#' @param cov_effect Numeric log-odds slope of the sector-level covariate on detection probability (default 0).
#' @return A list with `total_count`, `groups_both`, `groups_obs1`, `groups_obs2`, `detected_groups`, and `covariate`.
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
    cov_effect = 0
) {
    cov_val <- stats::rnorm(1, mean = 0, sd = 1)

    if (length(sim_groups) == 0) {
        return(list(
            total_count = 0,
            groups_both = if (is_do) 0L else NA_integer_,
            groups_obs1 = if (is_do) 0L else NA_integer_,
            groups_obs2 = if (is_do) 0L else NA_integer_,
            detected_groups = 0L,
            covariate = cov_val
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

    det_p1_sec <- if (cov_effect != 0) {
        base_logit <- stats::qlogis(pmin(0.99, pmax(0.01, det_p1)))
        stats::plogis(base_logit + cov_effect * cov_val)
    } else {
        det_p1
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
            detected_groups = length(idx_any),
            covariate = cov_val
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
            detected_groups = length(idx_blk),
            covariate = cov_val
        )
    }
}

#' Simulate a Population
#'
#' Simulates a population of individuals split into groups across multiple sectors.
#'
#' @param pop_size Average total population size (Poisson distributed).
#' @param mean_group Average group size (Poisson distributed).
#' @param tot_areas Total number of survey sectors.
#' @param area_sizes Vector of probabilities/sizes for assigning groups to sectors. If NULL,
#'   probabilities are proportional to \code{exp(n_cov_effect * habitat)} (uniform when \code{n_cov_effect = 0}).
#' @param n_cov_effect Optional numeric log-linear effect of a simulated standard-normal sector
#'   habitat covariate (\code{habitat}) on true sector abundance (default \code{0}).
#'
#' @return A list containing \code{tot_pop} (total simulated individuals), \code{sim_pop}
#'   (list of groups per sector), and \code{habitat} (simulated sector-level habitat covariate vector).
#' @importFrom stats rpois rnorm runif
#' @export
simulate_pop <- function(
    pop_size,
    mean_group,
    tot_areas,
    area_sizes = NULL,
    n_cov_effect = 0
) {
    habitat_vec <- stats::rnorm(tot_areas, mean = 0, sd = 1)
    tot_pop <- stats::rpois(1, pop_size)
    if (tot_pop == 0) {
        empty_pop <- vector("list", tot_areas)
        attr(empty_pop, "habitat") <- habitat_vec
        return(list(tot_pop = 0, sim_pop = empty_pop, habitat = habitat_vec))
    }

    if (is.null(area_sizes)) {
        area_sizes <- exp(n_cov_effect * habitat_vec)
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

    attr(sim_pop, "habitat") <- habitat_vec
    list(sim_pop = sim_pop, tot_pop = tot_pop, habitat = habitat_vec)
}

#' Simulate a Survey and format for DOAS
#'
#' @param sim_pop A simulated population from \code{\link{simulate_pop}} (either the full list or \code{pop$sim_pop}).
#' @param det_p1 Base detection probability for Observer 1.
#' @param intra_group_det Mean factor for count errors (e.g., 0.9 for undercounting).
#' @param det_group Base detection factor for group size scaling.
#' @param obs_effect Observer 2 relative efficiency (e.g., 0.75 means Obs 2 is 75 percent as effective).
#' @param detvar Standard deviation for detection probability across sectors.
#' @param surv_areas Number of sectors to receive Double Observer surveys.
#' @param cov_effect Optional numeric log-odds effect of a simulated standard-normal sector covariate
#'   (\code{covariate}) on detection probability (default \code{0}, meaning no covariate effect).
#'
#' @return A data frame with one row per sector, containing observed counts,
#'   Double Observer tallies, a sector detection covariate (\code{covariate}), and a sector
#'   abundance covariate (\code{habitat}) ready for input to \code{\link{doas}}.
#' @export
simulate_survey <- function(
    sim_pop,
    det_p1,
    intra_group_det,
    det_group,
    obs_effect,
    detvar,
    surv_areas,
    cov_effect = 0
) {
    if (is.list(sim_pop) && "sim_pop" %in% names(sim_pop)) {
        hab_vec <- sim_pop$habitat
        sim_pop <- sim_pop$sim_pop
    } else {
        hab_vec <- attr(sim_pop, "habitat")
    }

    tot_areas <- length(sim_pop)
    if (is.null(hab_vec) || length(hab_vec) != tot_areas) {
        hab_vec <- stats::rnorm(tot_areas, mean = 0, sd = 1)
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
            cov_effect = cov_effect
        )
        data.frame(
            sector_id = m,
            total_count = res$total_count,
            groups_both = res$groups_both,
            groups_obs1 = res$groups_obs1,
            groups_obs2 = res$groups_obs2,
            detected_groups = res$detected_groups,
            covariate = res$covariate,
            habitat = hab_vec[m],
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, sec_list)
}
