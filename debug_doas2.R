devtools::load_all(".")

set.seed(42)
simulated_population <- simulate_pop(
    pop_size = 1000,
    mean_group = 10,
    tot_areas = 15
)

# NO heterogeneity: det_group = 1.0, detvar = 0.0
survey_data <- simulate_survey(
    sim_pop = simulated_population$sim_pop,
    det_p1 = 0.6,
    intra_group_det = 1.0, # No internal group undercounting either
    det_group = 1.0,
    obs_effect = 0.75,
    detvar = 0.0,
    surv_areas = 15 # all sectors are DO to remove extrapolation error
)

sim_results <- doas(survey_data, method = "global", obs_type = "independent")

cat(sprintf(
    "\nTotal True Groups: %d\n",
    sum(sapply(simulated_population$sim_pop, length))
))
cat(sprintf("Total Detected Groups: %d\n", sum(survey_data$detected_groups)))
cat(sprintf(
    "Total Estimated Groups: %.2f\n",
    sum(sim_results$sector_estimates$est_groups)
))

cat(sprintf("\nTotal True N: %d\n", simulated_population$tot_pop))
cat(sprintf("Total Estimated N: %.2f\n", sim_results$total_abundance))

cat(sprintf("\np_hat: %.3f\n", mean(sim_results$results$p_hat, na.rm = T)))
