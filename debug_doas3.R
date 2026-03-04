library(DOAS)

set.seed(42)
simulated_population <- simulate_pop(
    pop_size = 1000,
    mean_group = 10,
    tot_areas = 15
)

# NO heterogeneity
survey_data <- simulate_survey(
    sim_pop = simulated_population$sim_pop,
    det_p1 = 0.6,
    intra_group_det = 1.0,
    det_group = 1.0,
    obs_effect = 0.75,
    detvar = 0.0,
    surv_areas = 15
)

sim_results <- doas(survey_data, method = "hybrid", obs_type = "independent")

do_data <- sim_results$sector_estimates
do_data$p_hat <- survey_data$detected_groups / do_data$est_groups

cat("\nSector details:\n")
for (i in 1:15) {
    true_groups = length(simulated_population$sim_pop[[i]])
    cat(sprintf(
        "Sector %2d: True_G=%2d Det_G=%2d Both=%2d Obs1=%2d Obs2=%2d | Est_G=%.2f p_hat=%.2f\n",
        i,
        true_groups,
        survey_data$detected_groups[i],
        survey_data$groups_both[i],
        survey_data$groups_obs1[i],
        survey_data$groups_obs2[i],
        do_data$est_groups[i],
        do_data$p_hat[i]
    ))
}

cat(sprintf("\nMean p_hat: %.3f\n", sim_results$mean_p_hat))
