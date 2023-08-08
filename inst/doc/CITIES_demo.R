## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 

## -----------------------------------------------------------------------------
library(cities)
library(dplyr)
library(tidyr)
library(plotly)
library(ggplot2)
library(ggthemes)


## ----echo = T, results = 'hide'-----------------------------------------------
total_data = 50
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

## ----echo=TRUE, results='hide', out.width="50%", fig.keep='all'---------------
mean_out = plot_means(n_patient_vector = n_patient_vector,
                      timepoints = timepoints,
                      pacf_list = pacf_list,
                      sigma_ar_vec = sigma_ar_vec,
                      mean_list = mean_list,
                      beta_list = beta_list,
                      reference_id = reference_id,
                      seed_val = starting_seed_val,
                      total_data = total_data,
                      threshold = threshold,
                      covariate_df = covariate_df,
                      static_output = static_output)


## ----echo=TRUE, results='hide', out.width="50%", fig.keep='all'---------------
plot_loe_ee (mean_list = mean_list,
             ref_grp = reference_id,
             stdev_vec = sigma_ar_vec,
             p_loe_max = p_loe_max,
             z_l_loe = z_l_loe,
             z_u_loe = z_u_loe,
             p_ee_max = p_ee_max,
             z_l_ee = z_l_ee,
             z_u_ee = z_u_ee,
             up_good = up_good,
             greyscale = FALSE,
             static_output = TRUE)

## ----echo=TRUE, results='hide', fig.keep='all',out.width="50%", message = FALSE----
data_out = data_generator_loop(n_patient_vector,
                               p_loe_max, 
                               z_l_loe,
                               z_u_loe,
                               p_ee_max,
                               z_l_ee,
                               z_u_ee,
                               timepoints,
                               pacf_list,
                               sigma_ar_vec,
                               mean_list,
                               beta_list,
                               p_admin,
                               rate_dc_ae,
                               prob_ae,
                               starting_seed_val,
                               reference_id, 
                               plot_po = FALSE,
                               up_good,
                               threshold,
                               total_data,
                               delta_adjustment_in,
                               covariate_df)

estimates_out = plot_estimates(data_out = data_out,  
                               total_data = total_data,
                               timepoints = timepoints,
                               reference_id = reference_id,
                               IR_display = IR_display,
                               normal_output = TRUE,
                               static_output = static_output)


## -----------------------------------------------------------------------------
display_table = estimates_out %>%
  mutate(mean_se = paste0(mean, " (", round(se, 2) , ")")) %>%
  dplyr::select(-se, -Arm, -mean) %>%
  pivot_wider(names_from = Estimand, values_from = mean_se)
display_table[,c(1,2,5,3,4,6,7)]

## -----------------------------------------------------------------------------
dc_out = plot_dc(data_out = data_out, 
                 total_data = total_data, 
                 timepoints = timepoints,
                 static_output = static_output)

dc_out1 = dc_out %>% 
  ungroup() %>%
  filter(Timepoints == max(timepoints),
         Reason != "OVERALL") %>%
  select(Arm, Reason, Value) %>%
  pivot_wider(names_from = Arm,
              values_from = Value) %>%
  arrange(factor(Reason, levels = c("AE", "LOE", "EE", "ADMIN")))
rbind(dc_out1,c("OVERALL",colSums(dc_out1[,-1])))

## ----echo = T, results = 'hide'-----------------------------------------------
total_data = 100
reference_id = 1
threshold = NA
timepoints = c(0,6,12,18,26)
IR_display = FALSE
delta_adjustment_in = NA

n_patient_ctrl = 195
n_patient_expt = 192
n_patient_vector = c(n_patient_ctrl, n_patient_expt)
n_total = sum(n_patient_vector)

mean_control = c(8, 8, 7.98, 7.97, 7.94)
mean_treatment = c(8, 7.45, 7.26, 7.21, 7.16)
mean_list = list(mean_control, mean_treatment)

sigma_ar_vec = c(0.8, 0.8)
pacf_list = c(0.5, 0.5)

beta_list = NA
covariate_df = NA

# LoE & EE
up_good = "Down" 
p_loe_max = 0.25 
z_l_loe = 1
z_u_loe = 4
p_ee_max = 0
z_l_ee = 0
z_u_ee = 0

# Admin & AE

p_admin_ctrl = 0.03 #(2+3+1+15+1+4)/192 = 0.1354167
p_admin_expt = 0.02 #(5+2+9+3)/197 = 0.0964467
p_admin = c(p_admin_ctrl, p_admin_expt)

prob_ae_ctrl = 0.53 #101/192
prob_ae_expt = 0.6 #118/197
prob_ae = c(prob_ae_ctrl, prob_ae_expt)

rate_dc_ae_ctrl = 0.01 # 2/192
rate_dc_ae_expt = 0.02 # 4/197
rate_dc_ae = c(rate_dc_ae_ctrl, rate_dc_ae_expt)

starting_seed_val = 1
static_output = TRUE

## ----echo=TRUE, results='hide', fig.keep='all', out.width="50%", message = FALSE----
mean_out = plot_means(n_patient_vector = n_patient_vector,
                      timepoints = timepoints,
                      pacf_list = pacf_list,
                      sigma_ar_vec = sigma_ar_vec,
                      mean_list = mean_list,
                      beta_list = beta_list,
                      reference_id = reference_id,
                      seed_val = starting_seed_val,
                      total_data = total_data,
                      threshold = threshold,
                      covariate_df = covariate_df,
                      static_output = static_output)


## ----echo=TRUE, results='hide', fig.keep='all',out.width="50%", message = FALSE----
plot_loe_ee (mean_list = mean_list,
             ref_grp = reference_id,
             stdev_vec = sigma_ar_vec,
             p_loe_max = p_loe_max,
             z_l_loe = z_l_loe,
             z_u_loe = z_u_loe,
             p_ee_max = p_ee_max,
             z_l_ee = z_l_ee,
             z_u_ee = z_u_ee,
             up_good = up_good,
             greyscale = FALSE,
             static_output = TRUE)

## ----echo=TRUE, results='hide', fig.keep='all',out.width="50%", message = FALSE----
data_out = data_generator_loop(n_patient_vector,
                               p_loe_max, 
                               z_l_loe,
                               z_u_loe,
                               p_ee_max,
                               z_l_ee,
                               z_u_ee,
                               timepoints,
                               pacf_list,
                               sigma_ar_vec,
                               mean_list,
                               beta_list,
                               p_admin,
                               rate_dc_ae,
                               prob_ae,
                               starting_seed_val,
                               reference_id, 
                               plot_po = FALSE,
                               up_good,
                               threshold,
                               total_data,
                               delta_adjustment_in,
                               covariate_df)

estimates_out = plot_estimates(data_out = data_out,  
                               total_data = total_data,
                               timepoints = timepoints,
                               reference_id = reference_id,
                               IR_display = IR_display,
                               normal_output = TRUE,
                               static_output = static_output)


## -----------------------------------------------------------------------------
display_table = estimates_out %>%
  mutate(mean_se = paste0(mean, " (", round(se, 2) , ")")) %>%
  dplyr::select(-se, -Arm, -mean) %>%
  pivot_wider(names_from = Estimand, values_from = mean_se)
display_table[,c(1,2,5,3,4)]

## -----------------------------------------------------------------------------
dc_out = plot_dc(data_out = data_out, 
                 total_data = total_data, 
                 timepoints = timepoints,
                 static_output = static_output)

## -----------------------------------------------------------------------------
dc_out1 = dc_out %>% 
  ungroup() %>%
  filter(Timepoints == max(timepoints),
         Reason != "OVERALL") %>%
  select(Arm, Reason, Value) %>%
  pivot_wider(names_from = Arm,
              values_from = Value) %>%
  arrange(factor(Reason, levels = c("AE", "LOE", "EE", "ADMIN")))
rbind(dc_out1,c("OVERALL",colSums(dc_out1[,-1])))
  

## ----echo = T, results = 'hide'-----------------------------------------------
total_data = 100
reference_id = 1
threshold = NA
timepoints = c(0, 12, 24, 36, 52, 64, 76)
IR_display = FALSE
delta_adjustment_in = NA

n_patient_ctrl = 126
n_patient_expt = 131
n_patient_vector = c(n_patient_ctrl, n_patient_expt)
n_total = sum(n_patient_vector)

sigma_ar_vec = c(10, 10)
pacf_list = c(0.5, 0.5)

mean_control = c(106, 105.87, 104.81, 102.85, 99.31, 97.67, 95.97)
mean_treatment = c(106, 106.35, 106.06, 105.29, 102.98, 101.09, 99.17)
mean_list = list(mean_control, mean_treatment)


beta_list = NA
covariate_df = NA

# LoE & EE
up_good = "Up" 
p_loe_max = 0 
z_l_loe = 0
z_u_loe = 0
p_ee_max = 0
z_l_ee = 0
z_u_ee = 0

# Admin & AE
p_admin_ctrl = 0.03#
p_admin_expt = 0.01 #
p_admin = c(p_admin_ctrl, p_admin_expt)

prob_ae_ctrl = 0.9 # 
prob_ae_expt = 0.9 # 
prob_ae = c(prob_ae_ctrl, prob_ae_expt)

rate_dc_ae_ctrl = 0.07 # 
rate_dc_ae_expt = 0.3 # 
rate_dc_ae = c(rate_dc_ae_ctrl, rate_dc_ae_expt)

starting_seed_val = 1
static_output = TRUE

## ----echo=TRUE, results='hide', fig.keep='all', out.width="50%",message = FALSE----
mean_out = plot_means(n_patient_vector = n_patient_vector,
                      timepoints = timepoints,
                      pacf_list = pacf_list,
                      sigma_ar_vec = sigma_ar_vec,
                      mean_list = mean_list,
                      beta_list = beta_list,
                      reference_id = reference_id,
                      seed_val = starting_seed_val,
                      total_data = total_data,
                      threshold = threshold,
                      covariate_df = covariate_df,
                      static_output = static_output)


## ----echo=TRUE, results='hide', fig.keep='all',out.width="50%", message = FALSE----
plot_loe_ee (mean_list = mean_list,
             ref_grp = reference_id,
             stdev_vec = sigma_ar_vec,
             p_loe_max = p_loe_max,
             z_l_loe = z_l_loe,
             z_u_loe = z_u_loe,
             p_ee_max = p_ee_max,
             z_l_ee = z_l_ee,
             z_u_ee = z_u_ee,
             up_good = up_good,
             greyscale = FALSE,
             static_output = TRUE)

## ----echo=TRUE, results='hide', fig.keep='all', out.width="50%",message = FALSE----
data_out = data_generator_loop(n_patient_vector,
                               p_loe_max, 
                               z_l_loe,
                               z_u_loe,
                               p_ee_max,
                               z_l_ee,
                               z_u_ee,
                               timepoints,
                               pacf_list,
                               sigma_ar_vec,
                               mean_list,
                               beta_list,
                               p_admin,
                               rate_dc_ae,
                               prob_ae,
                               starting_seed_val,
                               reference_id, 
                               plot_po = FALSE,
                               up_good,
                               threshold,
                               total_data,
                               delta_adjustment_in,
                               covariate_df)

estimates_out = plot_estimates(data_out = data_out,  
                               total_data = total_data,
                               timepoints = timepoints,
                               IR_display = IR_display,
                               reference_id = reference_id,
                               static_output = static_output)

## -----------------------------------------------------------------------------
display_table = estimates_out %>%
  mutate(mean_se = paste0(mean, " (", round(se, 2) , ")")) %>%
  dplyr::select(-se, -Arm, -mean) %>%
  pivot_wider(names_from = Estimand, values_from = mean_se)
display_table[,c(1,2,5,3,4)]

## ---- out.width="50%"---------------------------------------------------------
dc_out = plot_dc(data_out = data_out, 
                 total_data = total_data, 
                 timepoints = timepoints,
                 static_output = static_output)

## -----------------------------------------------------------------------------
dc_out1 = dc_out %>% 
  ungroup() %>%
  filter(Timepoints == max(timepoints),
         Reason != "OVERALL") %>%
  select(Arm, Reason, Value) %>%
  pivot_wider(names_from = Arm,
              values_from = Value) %>%
  arrange(factor(Reason, levels = c("AE", "LOE", "EE", "ADMIN")))
rbind(dc_out1,c("OVERALL",colSums(dc_out1[,-1])))
  

## ----echo = T, results = 'hide'-----------------------------------------------

p_admin_expt = 0.005
p_admin_ctrl = 0.005

prob_ae_expt = 0.6
rate_dc_ae_expt = 0.1
prob_ae_ctrl = 0.3
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

p_loe_max = 0.3 
z_l_loe = -5
z_u_loe = -3
p_ee_max = 0.1
z_l_ee = 6
z_u_ee = 7.5
timepoints = c(0:(length(mean_treatment)-1))*24

up_good = "Up" 
delta_adjustment_in = NA
threshold = NA

mean_list = list(mean_control, mean_treatment, mean_treatment2)

p_admin = c(p_admin_ctrl, p_admin_expt, p_admin_expt)
rate_dc_ae = c(rate_dc_ae_ctrl, rate_dc_ae_expt, rate_dc_ae_expt)
prob_ae = c(prob_ae_ctrl, prob_ae_expt,prob_ae_expt)
n_patient_vector = c(n_patient_ctrl, n_patient_expt, n_patient_expt2)

sigma_ar_vec = NA

total_data = 20
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

starting_seed_val = 1
static_output = TRUE
IR_display = TRUE

## ----echo=TRUE, results='hide', fig.keep='all', out.width="50%", message = FALSE----
mean_out = plot_means(n_patient_vector = n_patient_vector,
                      timepoints = timepoints,
                      pacf_list = pacf_list,
                      sigma_ar_vec = sigma_ar_vec,
                      mean_list = mean_list,
                      beta_list = beta_list,
                      reference_id = reference_id,
                      seed_val = starting_seed_val,
                      total_data = total_data,
                      threshold = threshold,
                      covariate_df = covariate_df,
                      static_output = static_output)


## ----echo=TRUE, results='hide', fig.keep='all', out.width="50%", message = FALSE----
plot_loe_ee (mean_list = mean_list,
             ref_grp = reference_id,
             stdev_vec = sigma_ar_vec,
             p_loe_max = p_loe_max,
             z_l_loe = z_l_loe,
             z_u_loe = z_u_loe,
             p_ee_max = p_ee_max,
             z_l_ee = z_l_ee,
             z_u_ee = z_u_ee,
             up_good = up_good,
             greyscale = FALSE,
             static_output = TRUE)

## ----echo=TRUE, results='hide', fig.keep='all', out.width="50%", message = FALSE----
data_out = data_generator_loop(n_patient_vector,
                               p_loe_max, 
                               z_l_loe,
                               z_u_loe,
                               p_ee_max,
                               z_l_ee,
                               z_u_ee,
                               timepoints,
                               pacf_list,
                               sigma_ar_vec,
                               mean_list,
                               beta_list,
                               p_admin,
                               rate_dc_ae,
                               prob_ae,
                               starting_seed_val,
                               reference_id, 
                               plot_po = FALSE,
                               up_good,
                               threshold,
                               total_data,
                               delta_adjustment_in,
                               covariate_df)

estimates_out = plot_estimates(data_out = data_out,  
                               total_data = total_data,
                               timepoints = timepoints,
                               IR_display = IR_display,
                               reference_id = reference_id,
                               static_output = static_output)

## -----------------------------------------------------------------------------
estimates_out %>%
  mutate(mean_se = paste0(mean, " (", round(se, 2) , ")")) %>%
  dplyr::select(-se, -mean) %>%
  pivot_wider(names_from = Estimand, values_from = mean_se)

## ---- out.width="50%"---------------------------------------------------------
dc_out = plot_dc(data_out = data_out, 
                 total_data = total_data, 
                 timepoints = timepoints,
                 static_output = static_output)

## -----------------------------------------------------------------------------
dc_out1 = dc_out %>% 
  ungroup() %>%
  filter(Timepoints == max(timepoints),
         Reason != "OVERALL") %>%
  select(Arm, Reason, Value) %>%
  pivot_wider(names_from = Arm,
              values_from = Value) %>%
  arrange(factor(Reason, levels = c("AE", "LOE", "EE", "ADMIN")))
rbind(dc_out1,c("OVERALL",colSums(dc_out1[,-1])))
  

