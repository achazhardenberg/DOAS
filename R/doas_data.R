#' Format and Validate Survey Data for DOAS
#'
#' Transforms raw field data from sector-wide, group-level long, or observer-sighting
#' long formats into a validated, sector-level `doas_data` data frame ready for
#' inspection prior to running \code{\link{doas}}.
#'
#' @param data A data frame containing survey records.
#' @param format Character string specifying the input data structure:
#'   - \code{"sector"} (default): One row per survey sector (wide format).
#'   - \code{"group"}: One row per detected animal group (long format).
#'   - \code{"sighting"}: One row per observer sighting record, linked by \code{group_id} (long format).
#' @param sector Character name of the sector identifier column (default \code{"sector_id"}).
#' @param count Character name of the individual count column (default \code{"total_count"} for
#'   \code{format = "sector"}, or \code{"group_size"} for \code{format = "group"} and \code{"sighting"}).
#' @param groups Character name of the total detected groups column for \code{format = "sector"}
#'   (default \code{"detected_groups"}). Optional for DO sectors if \code{both}, \code{obs1}, and \code{obs2} are supplied.
#' @param both Character name of the column for groups seen by both observers in \code{format = "sector"}
#'   (default \code{"groups_both"}). Optional if \code{groups}, \code{obs1}, and \code{obs2} are supplied.
#' @param obs1 Character name of the Observer 1 column (default \code{"groups_obs1"} for \code{format = "sector"},
#'   or binary indicator \code{"obs1"} for \code{format = "group"}).
#' @param obs2 Character name of the Observer 2 column (default \code{"groups_obs2"} for \code{format = "sector"},
#'   or binary indicator \code{"obs2"} for \code{format = "group"}).
#' @param seen_by Optional character name of a single categorical column in \code{format = "group"}
#'   indicating which observer(s) detected the group (values: \code{"both"}, \code{"obs1"}, \code{"obs2"}, \code{"block"}).
#'   If provided, used instead of binary \code{obs1}/\code{obs2} columns.
#' @param observer Character name of the observer ID column when \code{format = "sighting"} (default \code{"observer"}).
#' @param group_id Character name of the unique group ID column when \code{format = "sighting"} (default \code{"group_id"}).
#' @param primary_obs Value in the \code{observer} column identifying the Primary Observer (Observer 1)
#'   when \code{format = "sighting"} (default \code{1}, also matches \code{"1"}, \code{"obs1"}, \code{"Obs1"} if \code{NULL}).
#' @param count_rule Character rule for reconciling group size when both observers record a sighting of the
#'   same \code{group_id} in \code{format = "sighting"}: \code{"max"} (default), \code{"mean"}, or \code{"obs1"}.
#' @param do_sectors Optional vector of sector IDs that received Double Observer (DO) surveys. Strongly
#'   recommended for \code{format = "group"} and \code{format = "sighting"} so DO sectors where Observer 2
#'   happened to detect 0 groups are not misclassified as Block sectors.
#' @param all_sectors Optional vector of all surveyed sector IDs (including empty sectors with 0 detections).
#' @param covariates Optional character vector of sector-level covariate column names to retain.
#'
#' @return A validated sector-level data frame of class \code{c("doas_data", "data.frame")} with columns
#'   \code{sector_id}, \code{survey_type}, \code{total_count}, \code{detected_groups}, \code{groups_both},
#'   \code{groups_obs1}, \code{groups_obs2}, \code{mean_group_size}, and any requested \code{covariates}.
#' @importFrom stats na.omit
#' @export
#'
#' @examples
#' # 1. Group-level long format (1 row per detected group)
#' raw_groups <- data.frame(
#'   valley = c("A", "A", "A", "B", "B", "C"),
#'   size   = c(12, 5, 8, 10, 4, 15),
#'   who    = c("both", "obs1", "obs2", "both", "obs1", "block"),
#'   effort = c(8, 8, 8, 6, 6, 5)
#' )
#' df_sec <- doas_data(
#'   raw_groups,
#'   format = "group",
#'   sector = "valley",
#'   count = "size",
#'   seen_by = "who",
#'   do_sectors = c("A", "B"),
#'   covariates = "effort"
#' )
#' print(df_sec)
doas_data <- function(
    data,
    format = c("sector", "group", "sighting"),
    sector = "sector_id",
    count = NULL,
    groups = "detected_groups",
    both = "groups_both",
    obs1 = NULL,
    obs2 = NULL,
    seen_by = NULL,
    observer = "observer",
    group_id = "group_id",
    primary_obs = NULL,
    count_rule = c("max", "mean", "obs1"),
    do_sectors = NULL,
    all_sectors = NULL,
    covariates = NULL
) {
    format <- match.arg(format)
    count_rule <- match.arg(count_rule)

    if (!is.data.frame(data)) {
        stop("data must be a data.frame.")
    }
    if (!sector %in% names(data)) {
        stop(sprintf("Sector column '%s' not found in data.", sector))
    }
    if (!is.null(covariates) && !all(covariates %in% names(data))) {
        missing_cov <- setdiff(covariates, names(data))
        stop(sprintf(
            "Covariate column(s) not found in data: %s",
            paste(missing_cov, collapse = ", ")
        ))
    }

    out <- if (format == "sector") {
        format_sector_wide(
            data = data,
            sector = sector,
            count = if (is.null(count)) "total_count" else count,
            groups = groups,
            both = both,
            obs1 = if (is.null(obs1)) "groups_obs1" else obs1,
            obs2 = if (is.null(obs2)) "groups_obs2" else obs2,
            do_sectors = do_sectors,
            covariates = covariates
        )
    } else if (format == "group") {
        format_group_long(
            data = data,
            sector = sector,
            count = if (is.null(count)) "group_size" else count,
            obs1 = if (is.null(obs1)) "obs1" else obs1,
            obs2 = if (is.null(obs2)) "obs2" else obs2,
            seen_by = seen_by,
            do_sectors = do_sectors,
            all_sectors = all_sectors,
            covariates = covariates
        )
    } else {
        format_sighting_long(
            data = data,
            sector = sector,
            count = if (is.null(count)) "group_size" else count,
            observer = observer,
            group_id = group_id,
            primary_obs = primary_obs,
            count_rule = count_rule,
            do_sectors = do_sectors,
            all_sectors = all_sectors,
            covariates = covariates
        )
    }

    class(out) <- c("doas_data", "data.frame")
    out
}

#' Format Sector-Level Wide Data
#'
#' @keywords internal
#' @noRd
format_sector_wide <- function(
    data,
    sector,
    count,
    groups,
    both,
    obs1,
    obs2,
    do_sectors,
    covariates
) {
    if (!count %in% names(data)) {
        stop(sprintf("Count column '%s' not found in data.", count))
    }
    if (any(duplicated(data[[sector]]))) {
        stop(
            "Duplicate sector IDs found in 'sector' format. If your data has multiple rows per sector (e.g., one row per group or sighting), use doas_data(..., format = 'group') or format = 'sighting'."
        )
    }

    n_sec <- nrow(data)
    has_groups <- !is.null(groups) && (groups %in% names(data))
    has_both <- !is.null(both) && (both %in% names(data))
    has_obs1 <- !is.null(obs1) && (obs1 %in% names(data))
    has_obs2 <- !is.null(obs2) && (obs2 %in% names(data))

    det_vec <- if (has_groups) as.numeric(data[[groups]]) else rep(NA_real_, n_sec)
    both_vec <- if (has_both) as.numeric(data[[both]]) else rep(NA_real_, n_sec)
    obs1_vec <- if (has_obs1) as.numeric(data[[obs1]]) else rep(NA_real_, n_sec)
    obs2_vec <- if (has_obs2) as.numeric(data[[obs2]]) else rep(NA_real_, n_sec)

    is_do <- if (!is.null(do_sectors)) {
        data[[sector]] %in% do_sectors
    } else {
        !is.na(obs1_vec) & !is.na(obs2_vec) & (!is.na(both_vec) | !is.na(det_vec))
    }

    tot_cnt <- as.numeric(data[[count]])
    is_unsurv <- !is_do & is.na(tot_cnt)

    for (i in seq_len(n_sec)) {
        if (is_do[i]) {
            if (is.na(obs1_vec[i]) || is.na(obs2_vec[i])) {
                stop(sprintf(
                    "DO sector '%s' is missing 'obs1' or 'obs2' group counts.",
                    as.character(data[[sector]][i])
                ))
            }
            if (!is.na(both_vec[i]) && is.na(det_vec[i])) {
                det_vec[i] <- both_vec[i] + obs1_vec[i] + obs2_vec[i]
            } else if (is.na(both_vec[i]) && !is.na(det_vec[i])) {
                both_vec[i] <- det_vec[i] - obs1_vec[i] - obs2_vec[i]
                if (both_vec[i] < 0) {
                    stop(sprintf(
                        "In DO sector '%s', groups_obs1 + groups_obs2 (%d) exceeds detected_groups (%d). Ensure obs1 and obs2 represent groups seen ONLY by that observer.",
                        as.character(data[[sector]][i]),
                        as.integer(obs1_vec[i] + obs2_vec[i]),
                        as.integer(det_vec[i])
                    ))
                }
            } else if (!is.na(both_vec[i]) && !is.na(det_vec[i])) {
                expected_det <- both_vec[i] + obs1_vec[i] + obs2_vec[i]
                if (abs(det_vec[i] - expected_det) > 1e-6) {
                    stop(sprintf(
                        "In DO sector '%s', detected_groups (%g) does not equal groups_both + groups_obs1 + groups_obs2 (%g).",
                        as.character(data[[sector]][i]),
                        det_vec[i],
                        expected_det
                    ))
                }
            } else {
                stop(sprintf(
                    "DO sector '%s' must provide either 'groups_both' or 'detected_groups'.",
                    as.character(data[[sector]][i])
                ))
            }
        } else if (is_unsurv[i]) {
            det_vec[i] <- NA_real_
            both_vec[i] <- NA_real_
            obs1_vec[i] <- NA_real_
            obs2_vec[i] <- NA_real_
        } else {
            if (is.na(det_vec[i])) {
                if (!has_groups) {
                    # Individual-level survey (group size = 1): each individual is 1 group
                    det_vec[i] <- tot_cnt[i]
                } else {
                    stop(sprintf(
                        "Block sector '%s' is missing 'detected_groups'.",
                        as.character(data[[sector]][i])
                    ))
                }
            }
            both_vec[i] <- NA_real_
            obs1_vec[i] <- NA_real_
            obs2_vec[i] <- NA_real_
        }
    }

    surv_type <- ifelse(is_do, "DO", ifelse(is_unsurv, "Unsurveyed", "Block"))
    out <- data.frame(
        sector_id = data[[sector]],
        survey_type = surv_type,
        total_count = tot_cnt,
        detected_groups = det_vec,
        groups_both = both_vec,
        groups_obs1 = obs1_vec,
        groups_obs2 = obs2_vec,
        mean_group_size = ifelse(
            is_unsurv,
            NA_real_,
            ifelse(!is.na(det_vec) & det_vec > 0, tot_cnt / det_vec, 0)
        ),
        stringsAsFactors = FALSE
    )

    if (is.null(covariates)) {
        core_cols <- c(
            sector,
            count,
            groups,
            both,
            obs1,
            obs2,
            "sector_id",
            "survey_type",
            "total_count",
            "detected_groups",
            "groups_both",
            "groups_obs1",
            "groups_obs2",
            "mean_group_size"
        )
        extra_cols <- setdiff(names(data), core_cols)
        if (length(extra_cols) > 0) {
            covariates <- extra_cols
        }
    }

    if (!is.null(covariates)) {
        for (cov_col in covariates) {
            out[[cov_col]] <- data[[cov_col]]
        }
    }

    out
}

#' Format Group-Level Long Data (1 Row per Detected Group)
#'
#' @keywords internal
#' @noRd
format_group_long <- function(
    data,
    sector,
    count,
    obs1,
    obs2,
    seen_by,
    do_sectors,
    all_sectors,
    covariates
) {
    if (!count %in% names(data)) {
        if (identical(count, "group_size")) {
            # Individual-level long format: default each row to group_size = 1
            data[[count]] <- 1
        } else {
            stop(sprintf("Group size column '%s' not found in data.", count))
        }
    }

    use_seen_by <- !is.null(seen_by)
    if (use_seen_by && !seen_by %in% names(data)) {
        stop(sprintf("Column '%s' not found in data.", seen_by))
    }
    if (!use_seen_by && !obs1 %in% names(data)) {
        stop(sprintf(
            "Either 'seen_by' or binary observer columns ('%s', '%s') must be present in data.",
            obs1,
            obs2
        ))
    }

    sec_ids <- if (!is.null(all_sectors)) {
        unique(c(as.character(all_sectors), as.character(data[[sector]])))
    } else {
        unique(data[[sector]])
    }

    rows_list <- lapply(sec_ids, function(s_id) {
        sub <- data[data[[sector]] == s_id, , drop = FALSE]
        sizes <- as.numeric(sub[[count]])
        active_mask <- !is.na(sizes) & (sizes > 0)
        sub_act <- sub[active_mask, , drop = FALSE]
        sizes_act <- sizes[active_mask]

        if (use_seen_by) {
            tags <- tolower(trimws(as.character(sub_act[[seen_by]])))
            has_do_tags <- any(tags %in% c("both", "11", "obs2", "01", "2"))
            is_do <- if (!is.null(do_sectors)) {
                s_id %in% do_sectors
            } else {
                has_do_tags || (nrow(sub_act) > 0 && all(tags %in% c("both", "obs1", "obs2", "11", "10", "01")))
            }

            if (is_do) {
                b_cnt <- sum(tags %in% c("both", "11"))
                s1_cnt <- sum(tags %in% c("obs1", "10", "1", "primary"))
                s2_cnt <- sum(tags %in% c("obs2", "01", "2", "secondary"))
                det_cnt <- b_cnt + s1_cnt + s2_cnt
            } else {
                b_cnt <- NA_real_
                s1_cnt <- NA_real_
                s2_cnt <- NA_real_
                det_cnt <- length(sizes_act)
            }
        } else {
            o1 <- as.logical(sub_act[[obs1]])
            o2 <- if (obs2 %in% names(sub_act)) {
                as.logical(sub_act[[obs2]])
            } else {
                rep(NA, nrow(sub_act))
            }

            is_do <- if (!is.null(do_sectors)) {
                s_id %in% do_sectors
            } else {
                nrow(sub_act) > 0 && !all(is.na(o2))
            }

            if (is_do) {
                o2[is.na(o2)] <- FALSE
                o1[is.na(o1)] <- FALSE
                b_cnt <- sum(o1 & o2)
                s1_cnt <- sum(o1 & !o2)
                s2_cnt <- sum(!o1 & o2)
                det_cnt <- b_cnt + s1_cnt + s2_cnt
            } else {
                b_cnt <- NA_real_
                s1_cnt <- NA_real_
                s2_cnt <- NA_real_
                det_cnt <- sum(o1, na.rm = TRUE)
            }
        }

        tot_cnt <- sum(sizes_act, na.rm = TRUE)
        df_row <- data.frame(
            sector_id = s_id,
            survey_type = if (is_do) "DO" else "Block",
            total_count = tot_cnt,
            detected_groups = det_cnt,
            groups_both = b_cnt,
            groups_obs1 = s1_cnt,
            groups_obs2 = s2_cnt,
            mean_group_size = if (det_cnt > 0) tot_cnt / det_cnt else 0,
            stringsAsFactors = FALSE
        )

        if (!is.null(covariates)) {
            for (cov_col in covariates) {
                df_row[[cov_col]] <- if (nrow(sub) > 0) sub[[cov_col]][1] else NA
            }
        }
        df_row
    })

    do.call(rbind, rows_list)
}

#' Format Observer-Sighting Long Data (1 Row per Sighting, Linked by group_id)
#'
#' @keywords internal
#' @noRd
format_sighting_long <- function(
    data,
    sector,
    count,
    observer,
    group_id,
    primary_obs,
    count_rule,
    do_sectors,
    all_sectors,
    covariates
) {
    if (!count %in% names(data) && identical(count, "group_size")) {
        data[[count]] <- 1
    }
    for (col in c(count, observer, group_id)) {
        if (!col %in% names(data)) {
            stop(sprintf("Required column '%s' not found in data.", col))
        }
    }

    obs_vals <- unique(stats::na.omit(data[[observer]]))
    if (is.null(primary_obs)) {
        candidates <- c("1", "obs1", "observer1", "primary", "a", "obs_1")
        match_idx <- which(tolower(as.character(obs_vals)) %in% candidates)
        primary_obs <- if (length(match_idx) > 0) obs_vals[match_idx[1]] else obs_vals[1]
    }

    sec_ids <- if (!is.null(all_sectors)) {
        unique(c(as.character(all_sectors), as.character(data[[sector]])))
    } else {
        unique(data[[sector]])
    }

    rows_list <- lapply(sec_ids, function(s_id) {
        sub <- data[data[[sector]] == s_id, , drop = FALSE]
        sizes <- as.numeric(sub[[count]])
        active_mask <- !is.na(sizes) & (sizes > 0) & !is.na(sub[[group_id]])
        sub_act <- sub[active_mask, , drop = FALSE]

        is_do <- if (!is.null(do_sectors)) {
            s_id %in% do_sectors
        } else {
            any(as.character(sub[[observer]]) != as.character(primary_obs), na.rm = TRUE)
        }

        if (nrow(sub_act) == 0) {
            df_row <- data.frame(
                sector_id = s_id,
                survey_type = if (is_do) "DO" else "Block",
                total_count = 0,
                detected_groups = 0,
                groups_both = if (is_do) 0 else NA_real_,
                groups_obs1 = if (is_do) 0 else NA_real_,
                groups_obs2 = if (is_do) 0 else NA_real_,
                mean_group_size = 0,
                stringsAsFactors = FALSE
            )
        } else {
            grp_splits <- split(sub_act, sub_act[[group_id]])
            b_cnt <- 0
            s1_cnt <- 0
            s2_cnt <- 0
            tot_cnt <- 0

            for (g_df in grp_splits) {
                obs_in_g <- as.character(g_df[[observer]])
                saw_1 <- any(obs_in_g == as.character(primary_obs))
                saw_2 <- any(obs_in_g != as.character(primary_obs))

                if (saw_1 && saw_2) {
                    b_cnt <- b_cnt + 1
                } else if (saw_1) {
                    s1_cnt <- s1_cnt + 1
                } else if (saw_2) {
                    s2_cnt <- s2_cnt + 1
                }

                g_counts <- as.numeric(g_df[[count]])
                g_size <- if (count_rule == "max") {
                    max(g_counts, na.rm = TRUE)
                } else if (count_rule == "mean") {
                    mean(g_counts, na.rm = TRUE)
                } else {
                    if (saw_1) {
                        g_counts[obs_in_g == as.character(primary_obs)][1]
                    } else {
                        g_counts[1]
                    }
                }
                tot_cnt <- tot_cnt + g_size
            }

            det_cnt <- b_cnt + s1_cnt + s2_cnt
            df_row <- data.frame(
                sector_id = s_id,
                survey_type = if (is_do) "DO" else "Block",
                total_count = tot_cnt,
                detected_groups = det_cnt,
                groups_both = if (is_do) b_cnt else NA_real_,
                groups_obs1 = if (is_do) s1_cnt else NA_real_,
                groups_obs2 = if (is_do) s2_cnt else NA_real_,
                mean_group_size = if (det_cnt > 0) tot_cnt / det_cnt else 0,
                stringsAsFactors = FALSE
            )
        }

        if (!is.null(covariates)) {
            for (cov_col in covariates) {
                df_row[[cov_col]] <- if (nrow(sub) > 0) sub[[cov_col]][1] else NA
            }
        }
        df_row
    })

    do.call(rbind, rows_list)
}

#' Print Method for Validated DOAS Survey Data
#'
#' @param x An object of class `doas_data`
#' @param ... Additional arguments passed to `print.data.frame`
#' @export
print.doas_data <- function(x, ...) {
    n_tot <- nrow(x)
    n_do <- sum(x$survey_type == "DO", na.rm = TRUE)
    n_blk <- n_tot - n_do
    cat(sprintf(
        "DOAS Validated Survey Data: %d sectors (%d DO, %d Block) | %g total individuals in %g groups\n",
        n_tot,
        n_do,
        n_blk,
        sum(x$total_count, na.rm = TRUE),
        sum(x$detected_groups, na.rm = TRUE)
    ))
    cat("--------------------------------------------------------------------------------\n")
    NextMethod("print")
    invisible(x)
}
