#' @title plot_loe_ee
#'
#' @description Plots the lack of efficacy (LoE) and excess efficacy (EE) graphs
#'
#' @param mean_list List of vectors of means per arm
#' @param ref_grp ID for pairwise comparisons, e.g. for three arms,
#' if reference_id=1, then arms 2 and 3 will be compared only to arm 1
#' @param stdev_vec Vector of standard deviations per arm.
#' This is used to adjust the x-axis for display
#' @param p_loe_max The maximum probability of discontinuing due to LoE
#' @param z_l_loe The lower (or left) threshold of the LoE curve
#' @param z_u_loe The upper (or right) threshold of the LoE curve
#' @param p_ee_max The maximum probability of discontinuing due to EE
#' @param z_l_ee The lower (or left) threshold of the EE curve
#' @param z_u_ee The upper (or right) threshold of the EE curve
#' @param up_good "Up" if higher outcome values indicate better responses
#' and "Down" otherwise
#' @param greyscale TRUE for greyscale setting and FALSE for color setting
#' @param static_output TRUE, if static and FALSE if dynamic plot is requested
#'
#' @return The plot for LoE and EE.
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
#'
#'plot_loe_ee (mean_list = mean_list, ref_grp = reference_id,
#'stdev_vec = sigma_ar_vec, p_loe_max = p_loe_max, z_l_loe = z_l_loe,
#'z_u_loe = z_u_loe, p_ee_max = p_ee_max, z_l_ee = z_l_ee, z_u_ee = z_u_ee,
#'up_good = up_good, greyscale = FALSE, static_output = static_output)
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes

plot_loe_ee = function(mean_list,
                       ref_grp,
                       stdev_vec,
                       p_loe_max,
                       z_l_loe,
                       z_u_loe,
                       p_ee_max,
                       z_l_ee,
                       z_u_ee,
                       up_good,
                       greyscale,
                       static_output=FALSE){

  stdev_vec = ifelse(is.na(stdev_vec), 0, stdev_vec)

  val_max1 = max(do.call(rbind, mean_list)) - min(mean_list[[ref_grp]])
  val_min1 = min(do.call(rbind, mean_list)) - min(mean_list[[ref_grp]])

  buffer_value = max(stdev_vec*4, 10)
  val_min = val_min1 - buffer_value
  val_max = val_max1 + buffer_value

  z = seq(from = val_min, to = val_max, by = 0.01)

  p_loe = p_loe_ee_function(z = z,
                            p_max = p_loe_max,
                            z_l = z_l_loe,
                            p_min = 0,
                            z_u = z_u_loe,
                            up_good = (up_good == "Up"))

  p_ee = p_loe_ee_function(z = z,
                           p_max = p_ee_max,
                           z_l = z_l_ee,
                           p_min = 0,
                           z_u = z_u_ee,
                           up_good = (up_good == "Down"))

  p_z = data.frame(x = z,
                   p_loe = p_loe,
                   p_ee = p_ee) %>%
    gather(key=reason, value=p_discontinue, p_loe:p_ee) %>%
    mutate(Probability = round(p_discontinue,3),
           `Baseline Difference` = x,
           reason = if_else(reason == "p_ee", "Excess Efficacy", "Lack of Efficacy"),
           Reason = reason)


  plot_loe_ee_out = (ggplot(data=p_z, aes(x=`Baseline Difference`, y=Probability, color=Reason))+
                       geom_line(lwd=2, alpha = 0.75) +
                       xlab("Change of Outcome from Baseline")+
                       ylab("Probability of Discontinuing")+
                       #ggtitle(paste0("Maximum Probability of Discontinuing Group Probability Plot For ", age, "-Year Old ", ifelse(male == 1, "Male", "Female")))+
                       theme(plot.title = element_text(hjust = 0.5),
                             #legend.title = element_blank(),
                             legend.position="bottom")+
                       ggtitle(paste("Prob of Discontinuing") )+
                       ylim(0,1)+
                       xlim(val_min, val_max)+
                       scale_color_colorblind())

  if(greyscale == TRUE){
    plot_loe_ee_out = plot_loe_ee_out +
      scale_color_grey(start=0, end=0.6)
  }

  if(static_output == TRUE){
    print(plot_loe_ee_out)
  }else{
    print(ggplotly(plot_loe_ee_out))
  }

}
