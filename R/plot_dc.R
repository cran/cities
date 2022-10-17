#' @title plot_dc
#'
#' @description Plots the discontinuation rates by timepoints
#'
#' @param data_out The output from data_generator_loop()
#' @param total_data Total number of clinical trials to simulate
#' @param timepoints Vector of timepoints (e.g. weeks, days, time indices)
#' @param normal_output TRUE if both plots and numeric values of estimands are
#' requested. FALSE if only plots are requested
#' @param static_output TRUE if static mode requested and FALSE if dynamic plot
#' is requested
#' @param greyscale TRUE if greyscale requested and FALSE for color
#'
#' @return Plot and dataframe of proportion of discontinuations.
#'
#' @examples
#'total_data = 3
#'reference_id = 1
#'threshold = NA
#'timepoints = c(0,24,48,72,96,120,144)
#'IR_display = TRUE
#'delta_adjustment_in = c(0,1)

#'n_patient_ctrl = 120
#'n_patient_expt = 150
#'n_patient_vector = c(n_patient_ctrl, n_patient_expt)
#'n_total = sum(n_patient_vector)

#'mean_control = c(0,0,0,0,0,0,0)
#'mean_treatment = c(0,0.1,0.2,0.4,0.6,0.8,1)
#'mean_list = list(mean_control, mean_treatment)

#'sigma_ar_vec = c(1, 1)
#'pacf_list = list(c(-0.2, 0.4),
#'                 c(-0.2, 0.4))

#'beta_list = list(c(1.25, 1.25),
#'                 c(1.25, 1.25))
#'covariate_df = NA

# LoE & EE
#'up_good = "Up"
#'p_loe_max = 0.75
#'z_l_loe = -7
#'z_u_loe = -1
#'p_ee_max = 0.1
#'z_l_ee = 4
#'z_u_ee = 10

# Admin & AE

#'p_admin_ctrl = 0.02
#'p_admin_expt = 0.02
#'p_admin = c(p_admin_ctrl, p_admin_expt)

#'prob_ae_ctrl = 0.7
#'prob_ae_expt = 0.9
#'prob_ae = c(prob_ae_ctrl, prob_ae_expt)

#'rate_dc_ae_ctrl = 0.1
#'rate_dc_ae_expt = 0.1
#'rate_dc_ae = c(rate_dc_ae_ctrl, rate_dc_ae_expt)

#'starting_seed_val = 1
#'static_output = TRUE
#'data_out = data_generator_loop(n_patient_vector = n_patient_vector,
#'p_loe_max = p_loe_max, z_l_loe = z_l_loe, z_u_loe = z_u_loe,
#'p_ee_max = p_ee_max, z_l_ee = z_l_ee, z_u_ee = z_u_ee, timepoints = timepoints,
#'pacf_list = pacf_list, sigma_ar_vec = sigma_ar_vec, mean_list = mean_list,
#'beta_list = beta_list, p_admin = p_admin, rate_dc_ae = rate_dc_ae,
#'prob_ae = prob_ae, seed_val = starting_seed_val, reference_id = reference_id,
#'plot_po = FALSE, up_good = up_good, threshold = threshold,
#'total_data = total_data, delta_adjustment_in = delta_adjustment_in,
#'covariate_df = covariate_df)
#'
#'estimates_out = plot_estimates(data_out = data_out, total_data = total_data,
#'timepoints = timepoints, reference_id = reference_id, IR_display = IR_display,
#'normal_output = TRUE, static_output = static_output)
#'
#'dc_out = plot_dc(data_out = data_out, total_data = total_data,
#'timepoints = timepoints, static_output = static_output)
#'
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes
#'
plot_dc = function(data_out,
                   total_data,
                   timepoints,
                   normal_output = TRUE,
                   static_output = FALSE,
                   greyscale = FALSE){

  dc_df = data.frame()

  for(j in 1:length(data_out$dc_mean_list)){
    for(k in 1:length(data_out$dc_mean_list[[j]])){
      temp = data_out$dc_mean_list[[j]][[k]] %>%
        as.data.frame() %>%
        pivot_longer(cols = everything()) %>%
        dplyr::select(-name) %>%
        mutate(arm = k,
               Timepoints = rep(timepoints[-1], total_data),
               seed = rep(1:total_data, each = length(timepoints)-1),
               reason = toupper(strsplit(names(data_out$dc_mean_list)[j], "[_]")[[1]][3]))
      dc_df = rbind(dc_df, temp)
    }
  }
  dc_out_df = dc_df %>%
    mutate(arm = ifelse(arm == 1, "Ctrl", paste0("Expt ", arm-1))) %>%
    group_by(Timepoints, arm, reason) %>%
    summarize(sd_value = sd(value),
              value = round(mean(value), 2)) %>%
    rename(Reason = reason,
           Arm = arm,
           Value = value) %>%
    arrange(Arm, Reason, Timepoints)

  p_dc = ggplot(data = dc_out_df,
                aes(x=Reason, y=Value, fill=Arm)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_point(position = position_dodge(width=0.9)) +
    facet_wrap(~ordered(paste0("Time = ", Timepoints),
                        levels = paste0("Time = ",timepoints[-1]))) +
    scale_fill_brewer(palette="Dark2")  +
    xlab("Reason")+
    ylab("Proportion")+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="bottom")

  if(greyscale == TRUE){
    p_dc = p_dc  +
      scale_fill_grey(start=0, end=0.6)
  }

  if(static_output == TRUE){
    print(p_dc)
    return(dc_out_df)
  }else{
    if(normal_output == TRUE){
      print(ggplotly(p_dc))
      return(dc_out_df)
    }else{
      print(ggplotly(p_dc))
    }
  }
}
