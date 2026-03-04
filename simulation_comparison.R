# DOAS Estimator Comparison: Chapman vs Bayesian
#
# This script simulates a dataset where observers have low detection rates
# and occasionally miss animals entirely (sparse data).
# It then runs both the classic Chapman (Independent) estimator and the
# new Cartesian Bayesian estimator to compare how they handle the uncertainty.

# Load the development version of the package
devtools::load_all()

# Simulate a dataset specifically designed to be "sparse" (low counts)
set.seed(123)
n_sectors <- 10

# True Abundance per sector (True N) = 50
true_N <- rep(50, n_sectors)

# Assume detection probabilities are low (e.g., tough terrain like dense forest)
p1 <- 0.3 # Observer 1 detection probability
p2 <- 0.2 # Observer 2 detection probability

# Simulate observer counts based on True N and detection probs
detected_by_1 <- rbinom(n_sectors, true_N, p1)
detected_by_2 <- rbinom(n_sectors, true_N, p2)

# The overlap (seen by both) is out of the animals detected by obs 1
groups_both <- rbinom(n_sectors, detected_by_1, p2)

# Calculate the exclusive groups seen by each observer
groups_obs1 <- detected_by_1 - groups_both
groups_obs2 <- detected_by_2 - groups_both

# Prevent negative numbers due to random variation in the simulation bounds
groups_obs2[groups_obs2 < 0] <- 0

# Total unique groups detected
total_detected <- groups_both + groups_obs1 + groups_obs2

# Construct the survey dataframe
sim_data <- data.frame(
    sector_id = 1:n_sectors,
    total_count = total_detected * 2, # Assume average group size of 2
    groups_both = groups_both,
    groups_obs1 = groups_obs1,
    groups_obs2 = groups_obs2,
    mean_group_size = rep(2, n_sectors),
    detected_groups = total_detected
)

cat("=================================================================\n")
cat("                  RAW SIMULATED DO DATA                          \n")
cat("=================================================================\n")
print(sim_data[, c(
    "sector_id",
    "groups_both",
    "groups_obs1",
    "groups_obs2",
    "detected_groups"
)])
cat("\nNote that in Sector 6, Observer 1 saw ZERO exclusive groups.\n\n")

# Run the 2 estimators
res_chapman <- doas(sim_data, obs_type = "independent")
res_bayesian <- doas(sim_data, obs_type = "bayesian")

# Extract and Compare Results
comparison <- data.frame(
    Sector = 1:n_sectors,
    True_N = true_N,
    Detected_N = sim_data$total_count,
    Chapman_Est = round(res_chapman$sector_estimates$est_abundance, 1),
    Bayesian_Est = round(res_bayesian$sector_estimates$est_abundance, 1)
)

cat("=================================================================\n")
cat("               ESTIMATOR COMPARISON TABLE                        \n")
cat("=================================================================\n")
print(comparison)

cat("\n=================================================================\n")
cat(sprintf(
    "TOTAL CHAPMAN ESTIMATE:  %.1f (SE: %.1f)\n",
    res_chapman$total_abundance,
    res_chapman$se_abundance
))
cat(sprintf(
    "TOTAL BAYESIAN ESTIMATE: %.1f (SE: %.1f)\n",
    res_bayesian$total_abundance,
    res_bayesian$se_abundance
))
cat("=================================================================\n\n")
