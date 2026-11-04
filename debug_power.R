devtools::load_all(".")

cat("=== detvar=0.0 (no heterogeneity) ===\n")
p0 <- doas_power(
  n_sims = 200,
  pop_size = 1000,
  mean_group = 10,
  tot_areas = 30,
  det_p1 = 0.6,
  intra_group_det = 1.0,
  det_group = 0.5,
  obs_effect = 0.75,
  detvar = 0.0,
  surv_areas = c(5, 10, 20),
  method = "global",
  obs_type = "independent"
)
print(p0)

cat("\n=== detvar=0.3 (high heterogeneity), independent ===\n")
p1 <- doas_power(
  n_sims = 200,
  pop_size = 1000,
  mean_group = 10,
  tot_areas = 30,
  det_p1 = 0.6,
  intra_group_det = 1.0,
  det_group = 0.5,
  obs_effect = 0.75,
  detvar = 0.3,
  surv_areas = c(5, 10, 20),
  method = "global",
  obs_type = "independent"
)
print(p1)
