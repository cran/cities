#' @title plot_means
#'
#' @description Plots the means of simulation parameters.
#'
#' @param n_patient_vector Vector of number of patients
#' @param timepoints Vector of timepoints (e.g. weeks, days, time indices)
#' @param pacf_list List of pacf vectors
#' @param sigma_ar_vec Vector of variances per arm associated with list of pacf vectors
#' @param mean_list List of vectors of means per arm
#' @param beta_list List of vectors of beta coefficients per arm.
#' All vectors must have the same length and must be the same
#' as the number of columns for the covariate_df.
#' @param reference_id ID for pairwise comparisons, e.g. for three arms,
#' if reference_id=1, then arms 2 and 3 will be compared only to arm 1
#' @param seed_val Starting seed value
#' @param threshold Value to dichotomize continuous outcomes on
#' @param total_data Total number of clinical trials to simulate
#' @param covariate_df Matrix or dataframe of covariates.
#' Rows correspond to the total number of subjects. Order matters,
#' For instance, if you want to simulate a trial with 3 arms, each of size 30,50 and 80,
#' then covariate_df would have 30+50+80 rows such that the first 30 rows are
#' covariates for arm 1, the next 50 rows are covariates for arm 2
#' and the last 80 rows are covariates for arm 3.
#' @param static_output TRUE, if static and FALSE if dynamic plot is requested
#'
#' @return The plot of raw means.
#'
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
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes

plot_means = function(n_patient_vector,
                      timepoints,
                      pacf_list,
                      sigma_ar_vec,
                      mean_list,
                      beta_list,
                      reference_id,
                      seed_val,
                      threshold,
                      total_data,
                      covariate_df,
                      static_output = FALSE){


  po_summary = data.frame()
  pb = txtProgressBar(min = (seed_val)-1,      # Minimum value of the progress bar
                      max = (seed_val+total_data-1), # Maximum value of the progress bar
                      style = 3,    # Progress bar style (also available style = 1 and style = 2)
                      width = 50,   # Progress bar width. Defaults to getOption("width")
                      char = "=")   # Character used to create the bar
  for(seed_val in (seed_val):(seed_val+total_data-1)){
    temp = data_generator(n_patient_vector = n_patient_vector,
                          timepoints = timepoints,
                          pacf_list = pacf_list,
                          sigma_ar_vec = sigma_ar_vec,
                          mean_list = mean_list,
                          beta_list = beta_list,
                          seed_val = seed_val,
                          reference_id = reference_id,
                          plot_po = TRUE,
                          covariate_df = covariate_df)
    po_summary = rbind(po_summary, temp)
    # Sets the progress bar to the current state
    setTxtProgressBar(pb, seed_val)
  }

  plot_means_df = po_summary %>%
    pivot_longer(cols = -c(arm, seed)) %>%
    mutate(time = substr(name, 2, nchar(name)))%>%
    dplyr::select(-name, seed) %>%
    group_by(arm, time) %>%
    summarize(mean = round(mean(value),2),
              se = round(sd(value),2)) %>%
    mutate(arm = as.character(arm),
           time = as.numeric(time)) %>%
    rename(Timepoint = time,
           Arm = arm) %>%
    arrange(Arm, Timepoint)

  out_plot = ggplot(plot_means_df ,
                    aes(x=Timepoint, y=mean, color=Arm)) +
    geom_point() +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                  width=0.05,
                  size=0.4)+
    geom_line(aes(group = Arm),
              lwd=2, alpha = 0.75) +
    theme(plot.title = element_text(hjust=0.5),
          legend.position="bottom") +
    ggtitle(paste0("Plot of Mean Responses") )+
    xlab("Time") +
    ylab("Response")+
    scale_color_colorblind()  +
    scale_x_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))


  if(!is.na(threshold)){
    out_plot = out_plot +
      geom_hline(yintercept = threshold, lwd=0.5, lty=2)
  }

  if(static_output == TRUE){
    print(out_plot)
  }else{
    print(ggplotly(out_plot))
    return(plot_means_df)
  }
}
