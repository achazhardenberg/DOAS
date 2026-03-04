library(DOAS)

set.seed(42)
simulated_population <- simulate_pop(
  pop_size = 1000,
  mean_group = 10,
  tot_areas = 15
)

survey_data <- simulate_survey(
  sim_pop = simulated_population$sim_pop,
  det_p1 = 0.6,
  intra_group_det = 0.9,
  det_group = 0.5,
  obs_effect = 0.75,
  detvar = 0.2,
  surv_areas = 5
)

sim_results <- doas(survey_data, method = "hybrid", obs_type = "bayesian")

cat("\nTrue groups vs Detected groups:\n")
for (i in 1:15) {
  true_groups = length(simulated_population$sim_pop[[i]])
  if (true_groups > 0) {
    true_inds = sum(sapply(simulated_population$sim_pop[[i]], length))
  } else {
    true_inds = 0
  }
  cat(sprintf(
    "Sector %2d: True_G=%2d True_N=%4d | Det_G=%2d Det_N=%4d | Est_N=%.2f\n",
    i,
    true_groups,
    true_inds,
    survey_data$detected_groups[i],
    survey_data$total_count[i],
    sim_results$sector_estimates$est_abundance[i]
  ))
}

# The problem
cat(sprintf("\nTotal True N: %d\n", simulated_population$tot_pop))
cat(sprintf(
  "Total Estimated N: %.2f\n",
  sim_results$total_abundance
))

cat(sprintf(
  "Total True Groups: %d\n",
  sum(sapply(simulated_population$sim_pop, length))
))
cat(sprintf(
  "Total Estimated Groups: %.2f\n",
  sum(sim_results$sector_estimates$est_groups)
))

true_mean_group <- sum(sapply(simulated_population$sim_pop, function(x) {
  sum(sapply(x, length))
})) /
  sum(sapply(simulated_population$sim_pop, length))
obs_mean_group <- sum(survey_data$total_count) /
  sum(survey_data$detected_groups)
cat(sprintf("True Mean Group Size: %.2f\n", true_mean_group))
cat(sprintf("Observed Mean Group Size: %.2f\n", obs_mean_group))
