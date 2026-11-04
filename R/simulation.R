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
#' @importFrom stats rnorm runif plogis qlogis
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
    if (length(unique(grp_sqrt_anchored)) > 1) {
        group_det <- scales::rescale(
            grp_sqrt_anchored,
            to = c(det_group, 1)
        )[seq_along(grp_sqrt)]
    } else {
        group_det <- rep(det_group, length(grp_sqrt))
    }

    det_p1_sec <- if (cov_effect != 0) {
        base_logit <- stats::qlogis(pmin(0.99, pmax(0.01, det_p1)))
        stats::plogis(base_logit + cov_effect * cov_val)
    } else {
        det_p1
    }

    detp_factor <- truncnorm::rtruncnorm(
        1,
        a = 0,
        b = 1,
        mean = det_p1_sec,
        sd = detvar
    )

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
#' @param area_sizes Vector of probabilities/sizes for assigning groups to sectors. If NULL, uniform probabilities are used.
#'
#' @return A list containing the total population size and a list of groups per sector.
#' @importFrom stats rpois rnorm runif
#' @export
simulate_pop <- function(pop_size, mean_group, tot_areas, area_sizes = NULL) {
    tot_pop <- stats::rpois(1, pop_size)
    if (tot_pop == 0) {
        return(list(tot_pop = 0, sim_pop = vector("list", tot_areas)))
    }

    if (is.null(area_sizes)) {
        area_sizes <- rep(1, tot_areas)
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

    list(sim_pop = sim_pop, tot_pop = tot_pop)
}

#' Simulate a Survey and format for DOAS
#'
#' @param sim_pop A simulated population from \code{\link{simulate_pop}}.
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
#'   Double Observer tallies, and a sector-level \code{covariate} column ready for input to \code{\link{doas}}.
#' @importFrom truncnorm rtruncnorm
#' @importFrom scales rescale
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
    tot_areas <- length(sim_pop)
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
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, sec_list)
}
