#' @title data_generator
#'
#' @description Helper function to simulate single clinical trial
#'
#' @param n_patient_vector Vector of number of patients
#' @param p_loe_max The maximum probability of discontinuing due to LoE
#' @param z_l_loe The lower (or left) threshold of the LoE curve
#' @param z_u_loe The upper (or right) threshold of the LoE curve
#' @param p_ee_max The maximum probability of discontinuing due to EE
#' @param z_l_ee The lower (or left) threshold of the EE curve
#' @param z_u_ee The upper (or right) threshold of the EE curve
#' @param timepoints Vector of timepoints (e.g. weeks, days, time indices)
#' @param pacf_list List of pacf vectors
#' @param sigma_ar_vec Vector of variances per arm associated with list of pacf vectors
#' @param mean_list List of vectors of means per arm
#' @param beta_list List of vectors of beta coefficients per arm.
#' All vectors must have the same length and must be the same
#' as the number of columns for the covariate_df.
#' @param p_admin Vector of probabilities of discontinuing due to admin reasons
#' @param rate_dc_ae Vector of probabilities of observing at least one adverse
#' event
#' @param prob_ae Vector of proportions of discontinuing due to adverse event
#' @param seed_val Starting seed value
#' @param reference_id ID for pairwise comparisons, e.g. for three arms,
#' if reference_id=1, then arms 2 and 3 will be compared only to arm 1
#' @param plot_po TRUE, if plotting data only. Otherwise, set to FALSE
#' @param up_good "Up" if higher outcome values indicate better responses
#' @param threshold Value to dichotomize continuous outcomes on
#' @param delta_adjustment_in Vector of delta adjustment values or NA if none.
#' E.g. (2,3,1) when reference_id = 1 means no delta adjustment on arm 1
#' (even though 2 was supplied, but since arm 1 is the reference arm, this will
#' be defaulted to 0 regardless), 3 on arm 2 and 1 on arm 3.
#' @param covariate_df Matrix or dataframe of covariates. Set NA if using default
#' covariates, which comprises one continuous (standard normal) and binary
#' (bernoulli with prob 0.5) covariates. Rows correspond to the total number of
#' subjects. Order matters. For instance, if you want to simulate a trial
#' with 3 arms, each of size 30,50 and 80, then covariate_df would have 30+50+80
#' rows such that the first 30 rows are covariates for arm 1, the next 50 rows
#' are covariates for arm 2 and the last 80 rows are covariates for arm 3.
#'
#' @return List of dataframes of estimands and simulated data, including delta
#' adjusted ones if requested:
#'   \item{estimand_mean}{List of means of the FULL, S_{++}, S_{*+} and
#'   PP estimands}
#'   \item{estimand_sd}{List of standard deviations of the FULL, S_{++}, S_{*+}
#'   and PP estimands}
#'   \item{dc_mean_list}{List of proportions of discontinuations}
#'   \item{observed_df}{Dataframe of the observed outcomes}
#'   \item{po_df}{Dataframe of the potential outcomes}
#'   \item{ir_df}{Dataframe of the outcomes that have been adjusted via
#'   immediate reference (IR) or delta adjustment (Delta) for treatment policy
#'   estimands. The IR outcomes are labelled as ir_data while the delta
#'   adjusted outcomes are labelled as delta_data. The delta adjusted outcomes
#'   will only be available if the correct inputs for delta_adjustment_in are
#'   provided}
#' @examples
#' total_data = 3
#' reference_id = 1
#' threshold = NA
#' timepoints = c(0,24,48,72,96,120,144)
#' IR_display = TRUE
#' delta_adjustment_in = c(0,1)
#'
#' n_patient_ctrl = 120
#' n_patient_expt = 150
#' n_patient_vector = c(n_patient_ctrl, n_patient_expt)
#' n_total = sum(n_patient_vector)
#'
#' mean_control = c(0,0,0,0,0,0,0)
#' mean_treatment = c(0,0.1,0.2,0.4,0.6,0.8,1)
#' mean_list = list(mean_control, mean_treatment)
#'
#' sigma_ar_vec = c(1, 1)
#' pacf_list = list(c(-0.2, 0.4),
#'                  c(-0.2, 0.4))
#'
#' beta_list = list(c(1.25, 1.25),
#'                  c(1.25, 1.25))
#' covariate_df = NA
#'
#' # LoE & EE
#' up_good = "Up"
#' p_loe_max = 0.75
#' z_l_loe = -7
#' z_u_loe = -1
#' p_ee_max = 0.1
#' z_l_ee = 4
#' z_u_ee = 10
#'
#' # Admin & AE
#'
#' p_admin_ctrl = 0.02
#' p_admin_expt = 0.02
#' p_admin = c(p_admin_ctrl, p_admin_expt)
#'
#' prob_ae_ctrl = 0.7
#' prob_ae_expt = 0.9
#' prob_ae = c(prob_ae_ctrl, prob_ae_expt)
#'
#' rate_dc_ae_ctrl = 0.1
#' rate_dc_ae_expt = 0.1
#' rate_dc_ae = c(rate_dc_ae_ctrl, rate_dc_ae_expt)
#'
#' starting_seed_val = 1
#' static_output = TRUE
#'
#' mean_out = plot_means(n_patient_vector = n_patient_vector, timepoints = timepoints,
#' pacf_list = pacf_list, sigma_ar_vec = sigma_ar_vec, mean_list = mean_list,
#' beta_list = beta_list, reference_id = reference_id, seed_val = starting_seed_val,
#' total_data = total_data, threshold = threshold, covariate_df = covariate_df,
#' static_output = static_output)
#'
#'plot_loe_ee (mean_list = mean_list, ref_grp = reference_id,
#'stdev_vec = sigma_ar_vec, p_loe_max = p_loe_max, z_l_loe = z_l_loe,
#'z_u_loe = z_u_loe, p_ee_max = p_ee_max, z_l_ee = z_l_ee, z_u_ee = z_u_ee,
#'up_good = up_good, greyscale = FALSE, static_output = static_output)
#'
#'data_out = data_generator(n_patient_vector = n_patient_vector,
#'p_loe_max = p_loe_max, z_l_loe = z_l_loe, z_u_loe = z_u_loe,
#'p_ee_max = p_ee_max, z_l_ee = z_l_ee, z_u_ee = z_u_ee, timepoints = timepoints,
#'pacf_list = pacf_list, sigma_ar_vec = sigma_ar_vec, mean_list = mean_list,
#'beta_list = beta_list, p_admin = p_admin, rate_dc_ae = rate_dc_ae,
#'prob_ae = prob_ae, seed_val = starting_seed_val, reference_id = reference_id,
#'plot_po = FALSE, up_good = up_good, threshold = threshold,
#'delta_adjustment_in = delta_adjustment_in,
#'covariate_df = covariate_df)
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes
#'
data_generator = function(n_patient_vector,
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
                          seed_val,
                          reference_id,
                          plot_po = FALSE,
                          up_good = "Up",
                          threshold,
                          delta_adjustment_in,
                          covariate_df){

  # Negate the %in% operator in dplyr
  `%notin%`= Negate(`%in%`)

  ### Catch errors for means ###
  if(length(unique(lengths(mean_list))) != 1){
    warning("The lengths of means are not the same.")
    stop()
  }else{
    if(unique(lengths(mean_list)) != length(timepoints)){
      warning("The length of timepoints need to be the as the length of means.")
      stop()
    }

  }

  if(length(unique(lengths(sapply(mean_list,"[[",1)))) != 1){
    warning("The first entry in the means are not all the same.")
    stop()
  }

  set.seed(seed_val)
  total_patients = sum(n_patient_vector)
  n_repeat = length(timepoints) # number of repeat measures, including baseline
  # create a matrix with total number of patients x outcome measurements post baseline
  timepoints_duration = rep_row(diff(timepoints)/timepoints[n_repeat], total_patients)

  # Initialize data structures
  mean_subject_list = po_list = po_after_dc_list = z_list = loe_list = ee_list = admin_list = ae_list = dc_loe_list = dc_ee_list = dc_admin_list = dc_ae_list = dc_final = variance_conditional_list = list()

  ###################################################
  ######## 1.Generate the potential outcomes ########
  ###################################################

  # Check if custom covariate df exists
  if(!(is.data.frame(covariate_df) | is.matrix(covariate_df))){
    # Generate continuous and binary covariates per subject
    covariate_df = data.frame(continuous  = rnorm(n = total_patients, mean = 0, sd = 1),
                              binary = rbinom(n = total_patients, size = 1, prob = 0.5))
  }

  for(i in 1:length(mean_list)){

    # Create the covariance matrix
    if(length(sigma_ar_vec) == 0 | is.na(sum(sigma_ar_vec))){ # If custom covariance matrix provided
      covariance_matrix = as.matrix(pacf_list[[i]])
    }else{
      covariance_matrix = (sigma_ar_vec[i]^2)*pacf_vec_to_acf(pacf_vec = pacf_list[[i]],
                                                              n_repeat = n_repeat)
    }

    # Create the conditional means and covariances given baseline
    Sigma_12 = t(covariance_matrix[1, -1])
    Sigma_22 = covariance_matrix[-1,-1]
    Sigma_21 = as.matrix(covariance_matrix[-1, 1], ncol = 1)
    Sigma_11_inv = as.matrix(1/(covariance_matrix[1,1]))
    variance_conditional_list[[i]] = Sigma_22 - Sigma_21 %*% Sigma_11_inv %*% Sigma_12
  }


  # Generate outcomes at baseline
  baseline_outcome = round(rnorm(n = total_patients, mean = mean_list[[1]][1], sd = 1),2)
  for(i in 1:length(mean_list)){
    # Create a column of NA to include baseline in the po_list
    po_list[[i]] = cbind(baseline_outcome, timepoints_duration)
  }

  for(i in 1:length(mean_list)){
    # if NA or all 0's for betas supplied
    if( is.na(sum(unlist(beta_list))) | sum(((unlist(beta_list))) == 0) == length(unlist(beta_list)) ){
      mean_subject_list[[i]] = cbind(rep_row(mean_list[[i]][1], total_patients),
                                     rep_row(mean_list[[i]][-1], total_patients))

    }else{
      mean_subject_list[[i]] = cbind(rep_row(mean_list[[i]][1], total_patients),
                                     rep_row(mean_list[[i]][-1], total_patients) + rep_col(as.matrix(covariate_df) %*% as.matrix(beta_list[[i]],ncol=1), n_repeat-1))
    }

    # Generate the potential outcomes
    for(patient in 1:total_patients){
      mean_conditional  = mean_subject_list[[i]][patient,-1] + Sigma_21 %*% Sigma_11_inv %*% (po_list[[i]][patient,1] - mean_subject_list[[i]][patient,1])
      po_list[[i]][patient, 2:n_repeat]  =  round((mean_conditional + t(chol(variance_conditional_list[[i]]))%*%rnorm(n_repeat-1)), 2)
    }
  }

  if(plot_po == TRUE){
    po_summary = list()
    n_patient_cumsum = cumsum(n_patient_vector)
    first_patient_reference = max(n_patient_cumsum[reference_id-1]+1, 1)

    for(i in 1:(length(mean_list))){
      first_patient = max((n_patient_cumsum[i-1]+1), 1) # If first subject id is 0, force to be 1

      po_summary[[i]] = unname(colMeans(po_list[[i]][first_patient:n_patient_cumsum[i],]))
    }
    po_summary_out = cbind(do.call(rbind, po_summary), 1:i, seed_val)
    colnames(po_summary_out) = c(paste0("t", timepoints), "arm", "seed")

    return(as.data.frame(po_summary_out))
  }else{
    # Calculate prob of discontinuing due to adverse events
    prob_dc_ae = p_ae_poisson(rate_dc_ae = rate_dc_ae,
                              prob_ae = prob_ae)

    for(i in 1:length(mean_list)){

      # Calculate change from baseline
      z_list[[i]] = po_list[[i]][,-1] - po_list[[i]][,1]

      # Calculate probability of discontinuing due to LoE
      loe_list[[i]] =  p_loe_ee_function(z = z_list[[i]],
                                         p_max = p_loe_max,
                                         z_l = z_l_loe,
                                         p_min = 0,
                                         z_u = z_u_loe,
                                         up_good = (up_good == "Up"))
      # Generate loe discontinuations
      dc_loe_list[[i]] = rbinom(n = length(loe_list[[i]]), size = 1, prob = loe_list[[i]])
      dim(dc_loe_list[[i]]) = dim(timepoints_duration)
      dc_loe_list[[i]] = t(apply(dc_loe_list[[i]], 1, cumsum))
      dc_loe_list[[i]] = ifelse(dc_loe_list[[i]] > 0, 1, 0)

      # Calculate probability of discontinuing due to EE
      ee_list[[i]] = p_loe_ee_function(z = z_list[[i]],
                                       p_max = p_ee_max,
                                       z_l = z_l_ee,
                                       p_min = 0,
                                       z_u = z_u_ee,
                                       up_good = (up_good == "Down"))
      # Generate ee discontinuations
      dc_ee_list[[i]] = rbinom(n = length(ee_list[[i]]), size = 1, prob = ee_list[[i]])
      dim(dc_ee_list[[i]]) = dim(timepoints_duration)
      dc_ee_list[[i]] = t(apply(dc_ee_list[[i]], 1, cumsum))
      dc_ee_list[[i]] = ifelse(dc_ee_list[[i]] > 0, 1, 0)

      # Calculate probability of discontinuing due to admin
      admin_list[[i]] = matrix(p_admin[i], nrow = total_patients, ncol = ncol(ee_list[[i]]))
      # Generate admin discontinuations
      dc_admin_list[[i]] = rbinom(n = length(admin_list[[i]]), size = 1, prob = admin_list[[i]])
      dim(dc_admin_list[[i]]) = dim(timepoints_duration)
      dc_admin_list[[i]] = t(apply(dc_admin_list[[i]], 1, cumsum))
      dc_admin_list[[i]] = ifelse(dc_admin_list[[i]] > 0, 1, 0)

      # Generate AEs
      ae_list[[i]] = matrix(prob_ae[i], nrow = total_patients, ncol = ncol(ee_list[[i]]))
      ae_list[[i]] = rpois(n=length(timepoints_duration), lambda = -log(1-prob_ae[i])*timepoints_duration)
      dim(ae_list[[i]]) = dim(timepoints_duration)

      # Generate ae discontinuations
      dc_ae_list[[i]] = rbinom(n = length(ae_list[[i]]), size = 1, prob = 1-(1-prob_dc_ae[i])^(ae_list[[i]]) )
      dim(dc_ae_list[[i]]) = dim(timepoints_duration)
      dc_ae_list[[i]] = t(apply(dc_ae_list[[i]], 1, cumsum))
      dc_ae_list[[i]] = ifelse(dc_ae_list[[i]] > 0, 1, 0)

      dc_final[[i]] = dc_loe_list[[i]] + dc_ee_list[[i]] + dc_admin_list[[i]] + dc_ae_list[[i]]
      dc_final[[i]] = ifelse(dc_final[[i]] > 0, 1, 0)

      if(!is.na(threshold)){
        po_list[[i]] = ifelse(po_list[[i]] <= threshold, 0, 1)
      }

      po_after_dc_list[[i]] =  po_list[[i]][,-1] + ifelse(dc_final[[i]] > 0, NA, 0)
    }

    ###################################################
    # 2. Estimate of Estimands relative to reference #
    ###################################################

    pp_mean = s_plus_plus_mean = s_star_plus_mean = full_mean = pp_sd = s_plus_plus_sd = s_star_plus_sd = full_sd = list()
    ir_mean = delta_adj_mean = ir_sd = delta_adj_sd = ir_data = delta_adj_data = list()
    n_patient_cumsum = cumsum(n_patient_vector)

    first_patient_reference = max(n_patient_cumsum[reference_id-1]+1, 1) # If first subject id is 0, force to be 1
    reference_group = po_after_dc_list[[reference_id]][first_patient_reference:n_patient_cumsum[reference_id],]
    po_reference_group = po_list[[reference_id]][first_patient_reference:n_patient_cumsum[reference_id],-1]

    for(i in 1:(length(mean_list))){

      first_patient = max((n_patient_cumsum[i-1]+1), 1) # If first subject id is 0, force to be 1
      current_group = po_after_dc_list[[i]][first_patient:n_patient_cumsum[i],]

      # Immediate Reference
      ir_current = current_group
      ir_current[is.na(current_group)] = po_list[[reference_id]][first_patient:n_patient_cumsum[i],-1][is.na(current_group)]
      ir_mean[[i]] = colMeans(ir_current, na.rm = TRUE) - colMeans(po_reference_group, na.rm = TRUE)
      ir_sd[[i]] = sqrt(apply(ir_current, 2, sd, na.rm = TRUE)^2 + apply(po_reference_group, 2, sd, na.rm = TRUE)^2)
      ir_data[[i]] = ir_current

      # delta adjustment
      if(length(delta_adjustment_in) == 0 | is.na(sum(delta_adjustment_in))){
        delta_adjustment = numeric(length(mean_list))
      }else{
        delta_adjustment = delta_adjustment_in
      }

      delta_adj_current = current_group
      delta_adj_current[is.na(current_group)] = po_list[[i]][first_patient:n_patient_cumsum[i],-1][is.na(current_group)] + abs(delta_adjustment[i])*(ifelse(up_good == "Up", -1, 1))

       # If delta far exceeds the reference outcomes, then set it to the reference outcome
      if(up_good == "Up"){
        id_smaller_ref = (delta_adj_current < ir_current)
      }else{
        id_smaller_ref = (delta_adj_current > ir_current)
      }

      id_smaller_ref = (delta_adj_current < ir_current)
      delta_adj_current[id_smaller_ref] = po_list[[reference_id]][first_patient:n_patient_cumsum[i],-1][id_smaller_ref]

      delta_adj_mean[[i]] = colMeans(delta_adj_current, na.rm = TRUE) - colMeans(po_reference_group, na.rm = TRUE)
      delta_adj_sd[[i]] = sqrt(apply(delta_adj_current, 2, sd, na.rm = TRUE)^2 + apply(po_reference_group, 2, sd, na.rm = TRUE)^2)
      delta_adj_data[[i]] = delta_adj_current

      # Calculate ACE
      pp_mean[[i]] = colMeans(current_group, na.rm = TRUE) - colMeans(reference_group, na.rm = TRUE)
      s_plus_plus_mean[[i]] = colMeans(po_after_dc_list[[i]] - po_after_dc_list[[reference_id]], na.rm = TRUE)
      s_star_plus_mean[[i]] = colMeans(po_after_dc_list[[i]] - po_list[[reference_id]][,-1], na.rm = TRUE)
      full_mean[[i]] = colMeans(po_list[[i]][,-1] - po_list[[reference_id]][,-1], na.rm=TRUE)


      # Calculate variance of ACE
      pp_sd[[i]] = sqrt(apply(current_group, 2, sd, na.rm = TRUE)^2 + apply(reference_group, 2, sd, na.rm = TRUE)^2)
      s_plus_plus_sd[[i]] = apply(po_after_dc_list[[i]] - po_after_dc_list[[reference_id]], 2, sd, na.rm = TRUE)
      s_star_plus_sd[[i]] = apply(po_after_dc_list[[i]] - po_list[[reference_id]][,-1], 2, sd, na.rm = TRUE)
      full_sd[[i]] = apply(po_list[[i]][,-1] - po_list[[reference_id]][,-1], 2, sd, na.rm = TRUE)

    }

    pp_mean = pp_mean[-reference_id]
    s_plus_plus_mean = s_plus_plus_mean[-reference_id]
    s_star_plus_mean = s_star_plus_mean[-reference_id]
    full_mean = full_mean[-reference_id]
    ir_mean = ir_mean[-reference_id]
    delta_adj_mean = delta_adj_mean[-reference_id]

    pp_sd = pp_sd[-reference_id]
    s_plus_plus_sd = s_plus_plus_sd[-reference_id]
    s_star_plus_sd = s_star_plus_sd[-reference_id]
    full_sd = full_sd[-reference_id]
    ir_sd = ir_sd[-reference_id]
    delta_adj_sd = delta_adj_sd[-reference_id]

    dc_mean_list = list(lapply(dc_loe_list, colMeans),
                        lapply(dc_ee_list, colMeans),
                        lapply(dc_admin_list, colMeans),
                        lapply(dc_ae_list, colMeans),
                        lapply(dc_final, colMeans))
    names(dc_mean_list) = c("dc_mean_loe_list", "dc_mean_ee_list", "dc_mean_admin_list", "dc_mean_ae_list", "dc_mean_overall_list")

    # causal_estimand_mean = list(pp_mean, s_plus_plus_mean, s_star_plus_mean, full_mean, ir_mean, delta_adj_mean)
    # names(causal_estimand_mean) = c("PP", "S++", "S*+", "Full", "IR", "Delta")
    #
    # causal_estimand_sd = list(pp_sd, s_plus_plus_sd, s_star_plus_sd, full_sd, ir_sd, delta_adj_sd)
    # names(causal_estimand_sd) = c("PP", "S++", "S*+", "Full", "IR", "Delta")

    causal_estimand_mean = list(pp_mean, s_plus_plus_mean, full_mean, ir_mean, delta_adj_mean)
    names(causal_estimand_mean) = c("PP", "S++",  "Full", "IR", "Delta")

    causal_estimand_sd = list(pp_sd, s_plus_plus_sd, full_sd, ir_sd, delta_adj_sd)
    names(causal_estimand_sd) = c("PP", "S++", "Full", "IR", "Delta")

    #########################################################
    ##### Save Potential Outcomes and Observed Outcomes #####
    #########################################################
    po_df = observed_df = ir_df = delta_adj_df = data.frame()


    for(i in 1:length(mean_list)){
      first_patient = max((n_patient_cumsum[i-1]+1), 1) # If first subject id is 0, force to be 1

      # LOE
      dc_loe_temp_df = dc_loe_list[[i]]  %>%
        as.data.frame()
      colnames(dc_loe_temp_df) = paste0("t_", timepoints[-1])

      dc_loe_temp_df = dc_loe_temp_df %>%
        cbind(seed_val) %>%
        mutate(subject = 1:n(),
               arm = i)  %>%
        pivot_longer(cols = starts_with("t_")) %>%
        rename(timepoints = name,
               seed = seed_val,
               dc_loe = value) %>%
        mutate(timepoints = as.numeric(substr(timepoints, 3, nchar(timepoints))))

      # EE
      dc_ee_temp_df = dc_ee_list[[i]]  %>%
        as.data.frame()
      colnames(dc_ee_temp_df) = paste0("t_", timepoints[-1])

      dc_ee_temp_df = dc_ee_temp_df %>%
        cbind(seed_val) %>%
        mutate(subject = 1:n(),
               arm = i)  %>%
        pivot_longer(cols = starts_with("t_")) %>%
        rename(timepoints = name,
               seed = seed_val,
               dc_ee = value) %>%
        mutate(timepoints = as.numeric(substr(timepoints, 3, nchar(timepoints))))

      # admin
      dc_admin_temp_df = dc_admin_list[[i]]  %>%
        as.data.frame()
      colnames(dc_admin_temp_df) = paste0("t_", timepoints[-1])

      dc_admin_temp_df = dc_admin_temp_df %>%
        cbind(seed_val) %>%
        mutate(subject = 1:n(),
               arm = i)  %>%
        pivot_longer(cols = starts_with("t_")) %>%
        rename(timepoints = name,
               seed = seed_val,
               dc_admin = value) %>%
        mutate(timepoints = as.numeric(substr(timepoints, 3, nchar(timepoints))))

      # admin
      dc_ae_temp_df = dc_ae_list[[i]]  %>%
        as.data.frame()
      colnames(dc_ae_temp_df) = paste0("t_", timepoints[-1])

      dc_ae_temp_df = dc_ae_temp_df %>%
        cbind(seed_val) %>%
        mutate(subject = 1:n(),
               arm = i)  %>%
        pivot_longer(cols = starts_with("t_")) %>%
        rename(timepoints = name,
               seed = seed_val,
               dc_ae = value) %>%
        mutate(timepoints = as.numeric(substr(timepoints, 3, nchar(timepoints))))


      dc_df = dc_loe_temp_df %>%
        mutate(dc_ee = dc_ee_temp_df$dc_ee,
               dc_admin = dc_admin_temp_df$dc_admin,
               dc_ae = dc_ae_temp_df$dc_ae)

      # observed
      observed_temp_df = simulated_data_output(n_patient_cumsum = n_patient_cumsum,
                                               i = i,
                                               first_patient = first_patient,
                                               data_in = cbind(po_list[[i]][first_patient:n_patient_cumsum[i], 1],
                                                               po_after_dc_list[[i]][first_patient:n_patient_cumsum[i],]),
                                               covariate_df = covariate_df,
                                               timepoints = timepoints,
                                               beta_list = beta_list,
                                               seed_val = seed_val,
                                               potential_outcomes = FALSE,
                                               observed_indicator = NA)
      observed_df1 = observed_temp_df %>%
        left_join(dc_df, by=c("seed", "subject", "arm", "timepoints"))
      observed_df = rbind(observed_df, observed_df1)

      # PO
      po_temp_df = simulated_data_output(n_patient_cumsum = n_patient_cumsum,
                                         i = i,
                                         first_patient = first_patient,
                                         data_in = po_list[[i]],
                                         covariate_df = covariate_df,
                                         timepoints = timepoints,
                                         beta_list = beta_list,
                                         seed_val = seed_val,
                                         potential_outcomes = TRUE,
                                         observed_indicator = observed_temp_df %>%
                                           distinct(seed, subject, arm) %>%
                                           mutate(observed = 1))
      po_df1 = cbind(po_temp_df, dc_df[,-c(1:4)])
      po_df = rbind(po_df, po_df1)

      # Treatment policy: IR
      ir_temp_df = simulated_data_output(n_patient_cumsum = n_patient_cumsum,
                                         i = i,
                                         first_patient = first_patient,
                                         data_in = cbind(po_list[[i]][first_patient:n_patient_cumsum[i], 1],
                                                         ir_data[[i]]),
                                         covariate_df = covariate_df,
                                         timepoints = timepoints,
                                         beta_list = beta_list,
                                         seed_val = seed_val,
                                         potential_outcomes = FALSE,
                                         observed_indicator = NA)
      ir_df1 = ir_temp_df %>%
        left_join(dc_df, by=c("seed", "subject", "arm", "timepoints"))
      ir_df = rbind(ir_df, ir_df1)

      # Treatment policy: Delta
      delta_adj_temp_df = simulated_data_output(n_patient_cumsum = n_patient_cumsum,
                                                i = i,
                                                first_patient = first_patient,
                                                data_in = cbind(po_list[[i]][first_patient:n_patient_cumsum[i], 1],
                                                                delta_adj_data[[i]]),
                                                covariate_df = covariate_df,
                                                timepoints = timepoints,
                                                beta_list = beta_list,
                                                seed_val = seed_val,
                                                potential_outcomes = FALSE,
                                                observed_indicator = NA)
      delta_adj_df1 = delta_adj_temp_df %>%
        left_join(dc_df, by=c("seed", "subject", "arm", "timepoints"))
      delta_adj_df = rbind(delta_adj_df, delta_adj_df1)

    }



    observed_df = observed_df %>%
      group_by(seed, subject, arm) %>%
      mutate(flag = cumsum(is.na(aval)) ) %>%
      filter(flag < 2) %>%
      dplyr::select(-flag)
  }

  # Check if delta adjustment was requested
  if(length(delta_adjustment_in) == 0 | is.na(sum(delta_adjustment_in))){
    causal_estimand_mean = causal_estimand_mean[names(causal_estimand_mean) %notin% "Delta"]
    causal_estimand_sd = causal_estimand_sd[names(causal_estimand_sd) %notin% "Delta"]

    return(list("estimand_mean" = causal_estimand_mean,
                "estimand_sd" = causal_estimand_sd,
                "dc_mean_list" = dc_mean_list,
                "observed_df" = observed_df,
                "po_df" = po_df,
                "ir_df" = ir_df))
  }else{
    return(list("estimand_mean" = causal_estimand_mean,
                "estimand_sd" = causal_estimand_sd,
                "dc_mean_list" = dc_mean_list,
                "observed_df" = observed_df,
                "po_df" = po_df,
                "ir_df" = ir_df,
                "delta_adj_df" = delta_adj_df))
  }

}
