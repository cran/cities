#' @title simulated_data_output
#'
#' @description Helper function to combine simulated data
#'
#' @param n_patient_cumsum Vector of number of patients
#' @param i Index for arm
#' @param first_patient Index for first patient of arm
#' @param data_in Simulated data from data_generator()
#' @param covariate_df Matrix or dataframe of covariates.
#' Rows correspond to the total number of subjects. Order matters,
#' For instance, if you want to simulate a trial with 3 arms, each of size 30,50 and 80,
#' then covariate_df would have 30+50+80 rows such that the first 30 rows are
#' covariates for arm 1, the next 50 rows are covariates for arm 2
#' and the last 80 rows are covariates for arm 3.
#' @param timepoints Vector of timepoints (e.g. weeks, days, time indices)
#' @param beta_list List of vectors of beta coefficients per arm.
#' All vectors must have the same length and must be the same
#' as the number of columns for the covariate_df.
#' @param seed_val Current seed value
#' @param potential_outcomes TRUE if data to be combined is for potential
#' outcomes, and FALSE otherwise
#' @param observed_indicator Dataframe containing which subjects/arms/timepoints
#' were observed (necessary for potential outcomes), else default to NA
#'
#' @return Dataframe of for either potential outcomes, observed outcomes,
#' outcomes with immediate reference assumption or delta adjustment assumption
#'
#' @examples
#'n_patient_ctrl = 120
#'n_patient_expt = 150
#'n_patient_vector = c(n_patient_ctrl, n_patient_expt)
#'n_patient_cumsum = cumsum(n_patient_vector)
#'total_patients = sum(n_patient_vector)
#'timepoints = c(0,24,48,72,96,120,144)
#'data_in = matrix(rnorm(length(timepoints)*n_patient_ctrl), ncol = length(timepoints))
#'i = 1
#'first_patient = 1

#'covariate_df = data.frame(continuous  = rnorm(n = total_patients, mean = 0, sd = 1),
#'binary = rbinom(n = total_patients, size = 1, prob = 0.5))
#'beta_list = NA
#'seed_val = 1
#'potential_outcomes = FALSE
#'observed_indicator = NA

#'simulated_data_output(n_patient_cumsum = n_patient_cumsum, i = i,
#'first_patient = first_patient, data_in = data_in, covariate_df = covariate_df,
#'timepoints = timepoints, beta_list = beta_list, seed_val = seed_val,
#'potential_outcomes = FALSE, observed_indicator = NA)
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes
#'
simulated_data_output = function(n_patient_cumsum,
                                 i,
                                 first_patient,
                                 data_in,
                                 covariate_df,
                                 timepoints,
                                 beta_list,
                                 seed_val,
                                 potential_outcomes = FALSE,
                                 observed_indicator = NA){

  if(potential_outcomes == TRUE){
    temp_df = data_in %>%
      as.data.frame()
    colnames(temp_df) = paste0("t_", timepoints)
    temp_df = temp_df %>%
      {
        if(is.na(sum(unlist(beta_list))) | sum(((unlist(beta_list))) == 0) == length(unlist(beta_list))){
          cbind(., seed_val)
        }else{
          cbind(., covariate_df, seed_val)
        }
      } %>%
      mutate(subject = 1:n(),
             arm = i) %>%
      rename(seed = seed_val,
             base = t_0) %>%
      pivot_longer(cols = starts_with("t_")) %>%
      rename(timepoints = name,
             aval = value) %>%
      mutate(timepoints = as.numeric(substr(timepoints, 3, nchar(timepoints)))) %>%
      left_join(observed_indicator, by=c("seed", "subject", "arm")) %>%
      mutate(observed = if_else(is.na(observed), 0, observed))
  }else{
    temp_df = data_in %>%
      as.data.frame()
    colnames(temp_df) = paste0("t_", timepoints)
    temp_df = temp_df %>%
      {
        if(is.na(sum(unlist(beta_list))) | sum(((unlist(beta_list))) == 0) == length(unlist(beta_list))){
          cbind(., seed_val)
        }else{
          cbind(., covariate_df[first_patient:n_patient_cumsum[i], ], seed_val)
        }
      } %>%
      mutate(subject = first_patient:n_patient_cumsum[i],
             arm = i) %>%
      rename(seed = seed_val,
             base = t_0) %>%
      pivot_longer(cols = starts_with("t_")) %>%
      rename(timepoints = name,
             aval = value) %>%
      mutate(timepoints = as.numeric(substr(timepoints, 3, nchar(timepoints))))
  }
  return(temp_df)
}
