# DOAS Estimator Comparison: Frequentist (Chapman) vs Multi-Chain Bayesian MCMC
#
# Simulates a sparse 10-sector Double Observer survey and compares the
# Frequentist (Chapman) and Multi-Chain Bayesian MCMC estimators on the
# exact same simulated dataset.

devtools::load_all()

set.seed(123)
n_sectors <- 10

# True Abundance per sector (True N) = 50 groups * 2 individuals = 100
true_G <- rep(50, n_sectors)
p1 <- 0.3 # Observer 1 detection probability
p2 <- 0.2 # Observer 2 detection probability

detected_by_1 <- rbinom(n_sectors, true_G, p1)
detected_by_2 <- rbinom(n_sectors, true_G, p2)
groups_both <- rbinom(n_sectors, detected_by_1, p2)
groups_obs1 <- detected_by_1 - groups_both
groups_obs2 <- pmax(0, detected_by_2 - groups_both)
total_detected <- groups_both + groups_obs1 + groups_obs2

sim_data <- data.frame(
    sector_id = 1:n_sectors,
    total_count = total_detected * 2,
    groups_both = groups_both,
    groups_obs1 = groups_obs1,
    groups_obs2 = groups_obs2,
    detected_groups = total_detected
)

res_chapman <- doas(
    sim_data,
    method = "frequentist",
    obs_type = "independent",
    pooling = "hybrid"
)
res_bayesian <- doas(
    sim_data,
    method = "bayesian",
    obs_type = "independent",
    pooling = "hybrid"
)

comparison <- data.frame(
    Sector = 1:n_sectors,
    True_N = true_G * 2,
    Detected_N = sim_data$total_count,
    Chapman_Est = round(res_chapman$sector_estimates$est_abundance, 1),
    Bayesian_Est = round(res_bayesian$sector_estimates$est_abundance, 1)
)

print(comparison)
summary(res_chapman)
summary(res_bayesian)
