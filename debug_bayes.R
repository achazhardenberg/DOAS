devtools::load_all(".")

set.seed(42)
simulated_population <- simulate_pop(
    pop_size = 1000,
    mean_group = 10,
    tot_areas = 30
)

survey_data <- simulate_survey(
    sim_pop = simulated_population$sim_pop,
    det_p1 = 0.6,
    intra_group_det = 1.0,
    det_group = 0.5,
    obs_effect = 0.75,
    detvar = 0.3,
    surv_areas = 10
)

cat("\n=== INDEPENDENT ===\n")
res_ind <- doas(survey_data, method = "global", obs_type = "independent")
cat(sprintf(
    "True N: %d | Estimated N: %.2f | p_hat: %.3f\n",
    simulated_population$tot_pop,
    res_ind$total_abundance,
    res_ind$mean_p_hat
))

cat("\n=== BAYESIAN ===\n")
res_bay <- doas(survey_data, method = "global", obs_type = "bayesian")
cat(sprintf(
    "True N: %d | Estimated N: %.2f | p_hat: %.3f\n",
    simulated_population$tot_pop,
    res_bay$total_abundance,
    res_bay$mean_p_hat
))

# Let's manually compute what we expect
do_data <- survey_data[!is.na(survey_data$groups_both), ]
cat(sprintf("\nDO sectors: %d\n", nrow(do_data)))
cat(sprintf("Total groups_both: %d\n", sum(do_data$groups_both)))
cat(sprintf("Total groups_obs1: %d\n", sum(do_data$groups_obs1)))
cat(sprintf("Total groups_obs2: %d\n", sum(do_data$groups_obs2)))
cat(sprintf(
    "Total detected (both+obs1+obs2): %d\n",
    sum(do_data$groups_both + do_data$groups_obs1 + do_data$groups_obs2)
))
cat(sprintf("Total detected_groups: %d\n", sum(do_data$detected_groups)))

# Check what Chapman gives pooled
chap_pool <- chapman_estimate(
    both = sum(do_data$groups_both),
    obs1_only = sum(do_data$groups_obs1),
    obs2_only = sum(do_data$groups_obs2)
)
cat(sprintf("\nPooled Chapman Total G estimate: %.2f\n", chap_pool))
cat(sprintf(
    "Primary-obs count: %d\n",
    sum(do_data$groups_both + do_data$groups_obs1)
))
cat(sprintf(
    "Expected p_hat (primary): %.3f\n",
    sum(do_data$groups_both + do_data$groups_obs1) / chap_pool
))

# Check what Bayesian gives pooled
set.seed(42)
bay_pool <- bayesian_conjugate_estimate(
    both = sum(do_data$groups_both),
    obs1_only = sum(do_data$groups_obs1),
    obs2_only = sum(do_data$groups_obs2),
    iter = 10000
)
cat(sprintf("\nPooled Bayesian Total G estimate: %.2f\n", bay_pool$mean))
cat(sprintf(
    "Expected p_hat (primary, Bayesian): %.3f\n",
    sum(do_data$groups_both + do_data$groups_obs1) / bay_pool$mean
))
