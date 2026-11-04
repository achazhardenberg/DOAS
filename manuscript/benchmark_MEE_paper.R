# ==============================================================================
# Benchmark Script for Methods in Ecology and Evolution (MEE) Application Paper
# Comparing the original published script (DOAS_script.R) against the DOAS package
# ==============================================================================
# Note: This file resides in manuscript/ (excluded from CRAN build via .Rbuildignore)
# and provides reproducible simulation benchmarks for the MEE manuscript.

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(DOAS)
}

# ------------------------------------------------------------------------------
# 1. Wrapper for the Original Paper's Estimator (DOAS_script.R lines 156-209)
#    Applied directly to a standard DOAS sector-level survey data frame so that
#    legacy and new estimators can be compared on the EXACT same simulated data.
# ------------------------------------------------------------------------------
legacy_doas_estimate <- function(df) {
  is_do <- !is.na(df$groups_both)
  do_rows <- df[is_do, , drop = FALSE]
  n_do <- nrow(do_rows)
  
  # Sector-by-sector detection probability (DOAS_script.R lines 156-167)
  group_detpSR <- numeric(n_do)
  for (m in seq_len(n_do)) {
    B  <- do_rows$groups_both[m]
    S1 <- do_rows$groups_obs1[m]
    S2 <- do_rows$groups_obs2[m]
    obs1_count <- B + S1
    obs2_count <- B + S2
    
    # Original paper condition: if (length(obs1.det) > 1 && length(obs2.det) > 0)
    if (obs1_count > 1 && obs2_count > 0) {
      est_groups <- (((obs1_count + 1) * (obs2_count + 1)) / (B + 1)) - 1
      group_detpSR[m] <- ((obs1_count / est_groups) + (obs2_count / est_groups)) / 2
    } else {
      group_detpSR[m] <- 0
    }
  }
  
  # Unweighted mean across DO sectors (DOAS_script.R line 187)
  estSR_groupdetp <- mean(group_detpSR, na.rm = TRUE)
  if (is.na(estSR_groupdetp) || estSR_groupdetp <= 0) {
    return(list(N_hat = NA_real_, p_hat = NA_real_, sector_p = group_detpSR))
  }
  
  # Extrapolate all sectors using estSR_groupdetp (DOAS_script.R lines 198-204)
  # In DO sectors, total_count is the union of Obs1 & Obs2, so single-observer
  # equivalent groups = B + 0.5 * (S1 + S2) (or Chapman in hybrid mode)
  mean_size <- ifelse(df$detected_groups > 0, df$total_count / df$detected_groups, 0)
  single_obs_groups <- ifelse(
    is_do,
    df$groups_both + 0.5 * (df$groups_obs1 + df$groups_obs2),
    df$detected_groups
  )
  N_legacy <- sum((single_obs_groups / estSR_groupdetp) * mean_size, na.rm = TRUE)
  
  list(
    N_hat = N_legacy,
    p_hat = estSR_groupdetp,
    sector_p = group_detpSR
  )
}

# ------------------------------------------------------------------------------
# 2. Simulator Comparing Original Asymmetric Block Count vs. Corrected Simulator
# ------------------------------------------------------------------------------
simulate_legacy_asymmetric_survey <- function(sim_pop, surv_areas, det_p1, det_group,
                                              intra_group_det, obs_effect, detvar) {
  # Simulates survey exactly as DOAS::simulate_survey(), except Block sectors
  # use p1 directly (as in DOAS_script.R's block.count()) instead of the average
  # observer efficiency p1 * (1 + obs_effect) / 2.
  df_sym <- simulate_survey(
    sim_pop = sim_pop,
    surv_areas = surv_areas,
    det_p1 = det_p1,
    det_group = det_group,
    intra_group_det = intra_group_det,
    obs_effect = obs_effect,
    detvar = detvar
  )
  
  # Simulate Block sectors with obs_effect = 1.0 (so p_block = p1, matching DOAS_script.R)
  df_p1_block <- simulate_survey(
    sim_pop = sim_pop,
    surv_areas = 0,
    det_p1 = det_p1,
    det_group = det_group,
    intra_group_det = intra_group_det,
    obs_effect = 1.0,
    detvar = detvar
  )
  
  block_idx <- which(is.na(df_sym$groups_both))
  df_sym[block_idx, ] <- df_p1_block[block_idx, ]
  df_sym
}

# ------------------------------------------------------------------------------
# 3. Benchmark Runner: Compare Legacy vs. DOAS Package Estimators
# ------------------------------------------------------------------------------
run_mee_benchmark <- function(n_sims = 50,
                              pop_size = 3500,
                              mean_group = 10,
                              tot_areas = 38,
                              surv_areas_vec = c(5, 10, 19),
                              det_p1 = 0.6,
                              obs_effect_vec = c(0.5, 0.75, 1.0),
                              det_group = 0.5,
                              intra_group_det = 1.0,
                              detvar = 0.1,
                              include_bayesian = FALSE,
                              seed = 42) {
  set.seed(seed)
  results <- list()
  idx <- 1
  
  for (obs_eff in obs_effect_vec) {
    for (s_areas in surv_areas_vec) {
      for (sim in seq_len(n_sims)) {
        pop <- simulate_pop(
          pop_size = pop_size,
          mean_group = mean_group,
          tot_areas = tot_areas
        )
        true_N <- pop$tot_pop
        
        # A. Original Asymmetric Simulation + Legacy Estimator (DOAS_script.R)
        df_asym <- simulate_legacy_asymmetric_survey(
          sim_pop = pop$sim_pop, surv_areas = s_areas, det_p1 = det_p1,
          det_group = det_group, intra_group_det = intra_group_det,
          obs_effect = obs_eff, detvar = detvar
        )
        leg_asym <- legacy_doas_estimate(df_asym)$N_hat / true_N
        
        # B. Corrected Symmetric Simulation (DOAS::simulate_survey)
        df_corr <- simulate_survey(
          sim_pop = pop$sim_pop, surv_areas = s_areas, det_p1 = det_p1,
          det_group = det_group, intra_group_det = intra_group_det,
          obs_effect = obs_eff, detvar = detvar
        )
        
        # 1. Legacy Estimator on Corrected Simulation
        leg_corr <- legacy_doas_estimate(df_corr)$N_hat / true_N
        
        # 2. DOAS Frequentist Hybrid (Pooled p_hat, Chapman in DO sectors)
        fit_freq_hyb <- doas(
          df_corr, method = "frequentist", obs_type = "independent",
          pooling = "hybrid"
        )
        freq_hyb <- fit_freq_hyb$total_abundance / true_N
        
        # 3. DOAS Frequentist Global (Pooled p_hat applied to all sectors)
        fit_freq_glob <- doas(
          df_corr, method = "frequentist", obs_type = "independent",
          pooling = "global"
        )
        freq_glob <- fit_freq_glob$total_abundance / true_N
        
        # 4. Optional: DOAS Bayesian Hybrid
        bayes_hyb <- NA_real_
        if (include_bayesian) {
          fit_bayes <- doas(
            df_corr, method = "bayesian", obs_type = "independent",
            pooling = "hybrid", iter = 1000, burnin = 250, chains = 2
          )
          bayes_hyb <- fit_bayes$total_abundance / true_N
        }
        
        results[[idx]] <- data.frame(
          obs_effect = obs_eff,
          surv_areas = s_areas,
          sim = sim,
          legacy_asymmetric_sim = leg_asym,
          legacy_corrected_sim  = leg_corr,
          doas_freq_hybrid      = freq_hyb,
          doas_freq_global      = freq_glob,
          doas_bayes_hybrid     = bayes_hyb
        )
        idx <- idx + 1
      }
    }
  }
  
  raw_df <- do.call(rbind, results)
  
  # Summarize Mean N_hat / N_true across simulations
  summary_df <- aggregate(
    cbind(legacy_asymmetric_sim, legacy_corrected_sim, doas_freq_hybrid, doas_freq_global) ~
      obs_effect + surv_areas,
    data = raw_df,
    FUN = function(x) round(mean(x, na.rm = TRUE), 3)
  )
  
  list(raw = raw_df, summary = summary_df)
}

# Run a verification benchmark if executed directly
if (sys.nframe() == 0) {
  cat("Running MEE Benchmark Comparison (n_sims = 50 per scenario)...\n\n")
  bench <- run_mee_benchmark(n_sims = 50)
  print(bench$summary, row.names = FALSE)
}
