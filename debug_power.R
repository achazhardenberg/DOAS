devtools::load_all(".")

p_res <- doas_power(
  n_sims = 100,
  pop_size = 1000,
  mean_group = 10,
  tot_areas = 15,
  det_p1 = 0.6,
  intra_group_det = 1.0,
  det_group = 1.0,
  obs_effect = 0.75,
  detvar = 0.0,
  surv_areas = c(5, 10, 15),
  method = "hybrid",
  obs_type = "bayesian"
)
print(p_res)
