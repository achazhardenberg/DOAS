test_that("doas_data formats group-level long format correctly", {
    raw_groups <- data.frame(
        valley = c("A", "A", "A", "B", "B", "C"),
        size   = c(12, 5, 8, 10, 4, 15),
        who    = c("both", "obs1", "obs2", "both", "obs1", "block"),
        effort = c(8, 8, 8, 6, 6, 5)
    )
    df_grp <- doas_data(
        raw_groups,
        format = "group",
        sector = "valley",
        count = "size",
        seen_by = "who",
        do_sectors = c("A", "B"),
        covariates = "effort"
    )
    expect_s3_class(df_grp, "doas_data")
    expect_equal(nrow(df_grp), 3)
    expect_equal(df_grp$total_count, c(25, 14, 15))
    expect_equal(df_grp$groups_both, c(1, 1, NA))
    expect_equal(df_grp$groups_obs1, c(1, 1, NA))
    expect_equal(df_grp$groups_obs2, c(1, 0, NA))

    res_grp <- doas(df_grp)
    expect_s3_class(res_grp, "doas")
    expect_true(res_grp$total_abundance > sum(df_grp$total_count))
})

test_that("doas_data formats sighting-level long format correctly", {
    raw_sightings <- data.frame(
        sector_id = c("S1", "S1", "S1", "S2", "S2", "S3"),
        observer  = c("Obs1", "Obs2", "Obs1", "Obs1", "Obs2", "Obs1"),
        group_id  = c("G1", "G1", "G2", "G3", "G3", "G4"),
        group_size = c(10, 10, 6, 8, 8, 20)
    )
    df_sight <- doas_data(
        raw_sightings,
        format = "sighting",
        do_sectors = c("S1", "S2")
    )
    expect_equal(df_sight$groups_both, c(1, 1, NA))
    expect_equal(df_sight$groups_obs1, c(1, 0, NA))
    expect_equal(df_sight$detected_groups, c(2, 1, 1))
})

test_that("doas_data formats sector-level wide format and computes missing groups_both", {
    sec_no_both <- data.frame(
        sector_id = c("S1", "S2", "S3"),
        total_count = c(50, 30, 40),
        detected_groups = c(10, 6, 8),
        groups_obs1 = c(2, 1, NA),
        groups_obs2 = c(1, 1, NA)
    )
    df_no_both <- doas_data(sec_no_both, format = "sector")
    expect_equal(df_no_both$groups_both, c(7, 4, NA))
})

test_that("doas_data supports individual-level counts (solitary species)", {
    ind_sec <- data.frame(
        sector_id   = c("S1", "S2", "S3", "S4"),
        total_count = c(10, 8, 15, 12),
        groups_both = c(6, 4, NA, NA),
        groups_obs1 = c(2, 2, NA, NA),
        groups_obs2 = c(2, 2, NA, NA)
    )
    d_ind <- doas_data(ind_sec, format = "sector")
    expect_equal(d_ind$detected_groups, c(10, 8, 15, 12))
    expect_equal(d_ind$mean_group_size, c(1, 1, 1, 1))
    res_ind_only <- doas(d_ind)
    expect_equal(
        res_ind_only$sector_estimates$est_abundance,
        res_ind_only$sector_estimates$est_groups
    )

    ind_long <- data.frame(
        sector_id = c("S1", "S1", "S1", "S2"),
        seen_by   = c("both", "obs1", "obs2", "block")
    )
    d_ind_long <- doas_data(ind_long, format = "group", seen_by = "seen_by", do_sectors = "S1")
    expect_equal(d_ind_long$total_count, c(3, 1))
    expect_equal(d_ind_long$mean_group_size, c(1, 1))
})
