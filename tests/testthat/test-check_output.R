library(testthat)        # load testthat package
library(cities)       # load our package

# Test whether the output is a data frame
test_that("data_generator() returns a list data frame", {
  total_data = 3
  reference_id = 1
  threshold = NA
  timepoints = c(0,24,48,72,96,120,144)
  IR_display = TRUE
  delta_adjustment_in = c(0,1)

  n_patient_ctrl = 120
  n_patient_expt = 150
  n_patient_vector = c(n_patient_ctrl, n_patient_expt)
  n_total = sum(n_patient_vector)

  mean_control = c(0,0,0,0,0,0,0)
  mean_treatment = c(0,0.1,0.2,0.4,0.6,0.8,1)
  mean_list = list(mean_control, mean_treatment)

  sigma_ar_vec = c(1, 1)
  pacf_list = list(c(-0.2, 0.4),
                   c(-0.2, 0.4))

  beta_list = list(c(1.25, 1.25),
                   c(1.25, 1.25))
  covariate_df = NA

  # LoE & EE
  up_good = "Up"
  p_loe_max = 0.75
  z_l_loe = -7
  z_u_loe = -1
  p_ee_max = 0.1
  z_l_ee = 4
  z_u_ee = 10

  # Admin & AE

  p_admin_ctrl = 0.02
  p_admin_expt = 0.02
  p_admin = c(p_admin_ctrl, p_admin_expt)

  prob_ae_ctrl = 0.7
  prob_ae_expt = 0.9
  prob_ae = c(prob_ae_ctrl, prob_ae_expt)

  rate_dc_ae_ctrl = 0.1
  rate_dc_ae_expt = 0.1
  rate_dc_ae = c(rate_dc_ae_ctrl, rate_dc_ae_expt)

  starting_seed_val = 1
  static_output = TRUE


  data_out = data_generator(n_patient_vector = n_patient_vector,
                            p_loe_max = p_loe_max, z_l_loe = z_l_loe, z_u_loe = z_u_loe,
                            p_ee_max = p_ee_max, z_l_ee = z_l_ee, z_u_ee = z_u_ee, timepoints = timepoints,
                            pacf_list = pacf_list, sigma_ar_vec = sigma_ar_vec, mean_list = mean_list,
                            beta_list = beta_list, p_admin = p_admin, rate_dc_ae = rate_dc_ae,
                            prob_ae = prob_ae, seed_val = starting_seed_val, reference_id = reference_id,
                            plot_po = FALSE, up_good = up_good, threshold = threshold,
                            delta_adjustment_in = delta_adjustment_in,
                            covariate_df = covariate_df)
  expect_type(data_out, "list")
})

# Test whether the output is a data frame with NA for beta_list
test_that("data_generator() returns a list data frame", {
  total_data = 3
  reference_id = 1
  threshold = NA
  timepoints = c(0,24,48,72,96,120,144)
  IR_display = TRUE
  delta_adjustment_in = c(0,1)

  n_patient_ctrl = 120
  n_patient_expt = 150
  n_patient_vector = c(n_patient_ctrl, n_patient_expt)
  n_total = sum(n_patient_vector)

  mean_control = c(0,0,0,0,0,0,0)
  mean_treatment = c(0,0.1,0.2,0.4,0.6,0.8,1)
  mean_list = list(mean_control, mean_treatment)

  sigma_ar_vec = c(1, 1)
  pacf_list = list(c(-0.2, 0.4),
                   c(-0.2, 0.4))

  beta_list = NA
  covariate_df = NA

  # LoE & EE
  up_good = "Up"
  p_loe_max = 0.75
  z_l_loe = -7
  z_u_loe = -1
  p_ee_max = 0.1
  z_l_ee = 4
  z_u_ee = 10

  # Admin & AE

  p_admin_ctrl = 0.02
  p_admin_expt = 0.02
  p_admin = c(p_admin_ctrl, p_admin_expt)

  prob_ae_ctrl = 0.7
  prob_ae_expt = 0.9
  prob_ae = c(prob_ae_ctrl, prob_ae_expt)

  rate_dc_ae_ctrl = 0.1
  rate_dc_ae_expt = 0.1
  rate_dc_ae = c(rate_dc_ae_ctrl, rate_dc_ae_expt)

  starting_seed_val = 1
  static_output = TRUE


  data_out = data_generator(n_patient_vector = n_patient_vector,
                            p_loe_max = p_loe_max, z_l_loe = z_l_loe, z_u_loe = z_u_loe,
                            p_ee_max = p_ee_max, z_l_ee = z_l_ee, z_u_ee = z_u_ee, timepoints = timepoints,
                            pacf_list = pacf_list, sigma_ar_vec = sigma_ar_vec, mean_list = mean_list,
                            beta_list = beta_list, p_admin = p_admin, rate_dc_ae = rate_dc_ae,
                            prob_ae = prob_ae, seed_val = starting_seed_val, reference_id = reference_id,
                            plot_po = FALSE, up_good = up_good, threshold = threshold,
                            delta_adjustment_in = delta_adjustment_in,
                            covariate_df = covariate_df)
  expect_type(data_out, "list")
})

# Test whether the output is a data frame with custom covariances and covariates
test_that("data_generator() returns a list data frame", {
  p_admin_expt = 0.02
  p_admin_ctrl = 0.03

  prob_ae_expt = 0.9
  rate_dc_ae_expt = 0.1
  prob_ae_ctrl = 0.7
  rate_dc_ae_ctrl = 0.1

  n_patient_expt = 120
  n_patient_expt2 = 100
  n_patient_ctrl = 150

  mean_treatment = c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.1, 1.2, 1.3, 1.4, 1.5)
  mean_treatment2 = 2*c(0, 0.1, 0.2, 0.4, 0.6, 0.8, 1,  1.1, 1.2, 1.3, 1.4, 1.5)
  mean_control = rep(0, length(mean_treatment))

  beta_control = c(1.25, 1.25, 1)
  beta_expt = c(1.25, 1.25, 1)
  beta_expt2 = c(1.25, 1.25, 1)

  p_loe_max = 0.6
  z_l_loe = -5
  z_u_loe = -1.2
  p_ee_max = 0.15
  z_l_ee = 2.5
  z_u_ee = 4.5
  timepoints = c(0:(length(mean_treatment)-1))*24

  pacf_vec = c(0.5, -0.2)
  sigma_ar = 1

  up_good = "Up"
  delta_adjustment_in = NA
  threshold = NA

  mean_list = list(mean_control, mean_treatment, mean_treatment2)

  p_admin = c(p_admin_ctrl, p_admin_expt, p_admin_expt)
  rate_dc_ae = c(rate_dc_ae_ctrl, rate_dc_ae_expt, rate_dc_ae_expt)
  prob_ae = c(prob_ae_ctrl, prob_ae_expt,prob_ae_expt)
  n_patient_vector = c(n_patient_ctrl, n_patient_expt, n_patient_expt2)

  sigma_ar_vec = NA

  total_data = 3
  reference_id = 1

  starting_seed_val = 1

  A = matrix(runif(length(timepoints)^2)*2-1, ncol=length(timepoints))
  Sigma = t(A) %*% A
  pacf_list = list(Sigma,
                   Sigma,
                   Sigma)

  n_total = sum(n_patient_vector)
  beta_list = list(beta_control, beta_expt, beta_expt2)

  covariate_df = data.frame(continuous_1 = (rnorm(n = n_total, mean = 0, sd = 1)),
                            continuous_2 = (rnorm(n = n_total, mean = 0, sd = 1)),
                            binary_1 = rbinom(n = n_total, size = 1, prob = 0.5))

  data_out = data_generator(n_patient_vector = n_patient_vector,
                            p_loe_max = p_loe_max, z_l_loe = z_l_loe, z_u_loe = z_u_loe,
                            p_ee_max = p_ee_max, z_l_ee = z_l_ee, z_u_ee = z_u_ee, timepoints = timepoints,
                            pacf_list = pacf_list, sigma_ar_vec = sigma_ar_vec, mean_list = mean_list,
                            beta_list = beta_list, p_admin = p_admin, rate_dc_ae = rate_dc_ae,
                            prob_ae = prob_ae, seed_val = starting_seed_val, reference_id = reference_id,
                            plot_po = FALSE, up_good = up_good, threshold = threshold,
                            delta_adjustment_in = delta_adjustment_in,
                            covariate_df = covariate_df)
  expect_type(data_out, "list")
})
