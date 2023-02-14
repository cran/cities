#' @title plot_estimates
#'
#' @description Plots the estimates of the estimands
#'
#' @param data_out The output from data_generator_loop()
#' @param total_data Total number of clinical trials to simulate
#' @param timepoints Vector of timepoints (e.g. weeks, days, time indices)
#' @param reference_id ID for pairwise comparisons, e.g. for three arms,
#' if reference_id=1, then arms 2 and 3 will be compared only to arm 1
#' @param IR_display TRUE if requested to display Immediate Reference estimand.
#' FALSE otherwise
#' @param delta_display TRUE if requested to display Delta estimand.
#' FALSE otherwise
#' @param normal_output TRUE if both plots and numeric values of estimands are
#' requested. FALSE if only plots are requested
#' @param static_output TRUE if static mode requested and FALSE if dynamic plot
#' is requested
#' @param greyscale TRUE if greyscale requested and FALSE for color
#'
#' @return Plot and dataframe of estimands.
#'
#' @examples
#'total_data = 3
#'reference_id = 1
#'threshold = NA
#'timepoints = c(0,24,48,72,96,120,144)
#'IR_display = TRUE
#'delta_adjustment_in = c(0,1)
#'delta_display = TRUE

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
#'delta_display = delta_display, normal_output = TRUE, static_output = static_output)
#'
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes
#'
plot_estimates = function(data_out,
                          total_data,
                          timepoints,
                          reference_id,
                          IR_display = TRUE,
                          delta_display = TRUE,
                          normal_output = TRUE,
                          static_output = FALSE,
                          greyscale = FALSE){

  # Negate the %in% operator in dplyr
  `%notin%`= Negate(`%in%`)

  arm_length = length(data_out$estimand_mean$PP)+1
  timepoint_df = data.frame(name = as.character(1:(length(timepoints)-1)),
                            Timepoints = timepoints[-1])

  estimands_mean = estimands_se = data.frame()
  estimands_label = names(data_out$estimand_mean)


  if(total_data == 1){
    for(i in 1:length(data_out$estimand_mean)){
      temp_mean = do.call(rbind, data_out$estimand_mean[[i]])
      estimands_mean = rbind(estimands_mean, temp_mean)
    }

    estimands_mean = estimands_mean %>%
      mutate(Arm = rep(c(1:arm_length)[-reference_id], length(estimands_label)),
             Estimand = rep(estimands_label, each=arm_length-1)) %>%
      pivot_longer(cols = -c(Arm,Estimand)) %>%
      mutate(name = (substr(name, 2, nchar(name)))) %>%
      rename(mean = value)

    estimands = estimands_mean %>%
      left_join(timepoint_df, by="name") %>%
      mutate(mean = round(mean, 2),
             Timepoints = as.numeric(Timepoints),
             Arm = paste0("Arm ", Arm))

    Timepoints_in = sort(unique(estimands$Timepoints))
    Estimands_in = estimands %>%
      #filter(Estimand != "S*+") %>%
      pull(Estimand)
    Estimands_in = unique(Estimands_in)


    if(IR_display & delta_display){
      Estimands_in = c("Full", "PP", "S*+", "S++", "IR", "Delta")

    }else if(IR_display & !delta_display){
      Estimands_in = c("Full", "PP", "S*+", "S++", "IR")
      estimands = estimands %>%
        filter(Estimand != "Delta")
    }else if(!IR_display & delta_display){
      Estimands_in = c("Full", "PP", "S*+", "S++", "Delta")
      estimands = estimands %>%
        filter(Estimand != "IR")
    }else{
      print("HERE!")
      Estimands_in = c("Full", "PP", "S*+", "S++")
      estimands = estimands %>%
        filter(Estimand != "IR",
               Estimand != "Delta")
      print(unique(estimands$Estimand))
    }

    p_estimands = ggplot(estimands %>%
                           #filter(Estimand != "S*+") %>%
                           mutate(Timepoints = ordered(factor(paste0("Time = ", Timepoints)),
                                                       levels = paste0("Time = ", Timepoints_in)),
                                  Estimand = ordered(factor(Estimand),
                                                     levels = Estimands_in),
                                  Arm = paste0(Arm, " vs Arm ", reference_id)),
                         aes(x=Estimand, y=mean, color=Estimand, shape=Estimand)) +
      geom_point(position=position_dodge(.1), size=3) +
      facet_grid(Arm ~ Timepoints) +
      #xlab("Time") +
      ylab("Response") +
      #guides(color=guide_legend("")) +
      theme(plot.title = element_text(hjust=0.5),
            legend.position="bottom",
            axis.text.x = element_text(angle=90,hjust=0.95, vjust=0.2))  +
      scale_color_brewer(palette="Dark2")

  }else{
    for(i in 1:length(data_out$estimand_mean)){

      temp_mean = do.call(rbind, lapply(data_out$estimand_mean[[i]], colMeans))
      estimands_mean = rbind(estimands_mean, temp_mean)

      temp_sd = do.call(rbind, lapply(data_out$estimand_mean[[i]], colSD))
      estimands_se = rbind(estimands_se, temp_sd)
    }

    estimands_mean = estimands_mean %>%
      mutate(Arm = rep(c(1:arm_length)[-reference_id], length(estimands_label)),
             Estimand = rep(estimands_label, each=arm_length-1)) %>%
      pivot_longer(cols = -c(Arm,Estimand)) %>%
      mutate(name = (substr(name, 2, nchar(name)))) %>%
      rename(mean = value)

    estimands_se = estimands_se %>%
      mutate(Arm = rep(c(1:arm_length)[-reference_id], length(estimands_label)),
             Estimand = rep(estimands_label, each=arm_length-1)) %>%
      pivot_longer(cols = -c(Arm,Estimand)) %>%
      mutate(name = (substr(name, 2, nchar(name)))) %>%
      rename(se = value)

    estimands = estimands_mean %>%
      left_join(estimands_se, by=c("Estimand", "name", "Arm")) %>%
      left_join(timepoint_df, by="name") %>%
      mutate(mean = round(mean, 2),
             Timepoints = as.numeric(Timepoints),
             Arm = paste0("Arm ", Arm))

    Timepoints_in = sort(unique(estimands$Timepoints))
    Estimands_in = estimands %>%
      #filter(Estimand != "S*+") %>%
      pull(Estimand)
    Estimands_in = unique(Estimands_in)

    if(IR_display & delta_display){
      Estimands_in = c("Full", "PP", "S*+", "S++", "IR", "Delta")

    }else if(IR_display & !delta_display){
      Estimands_in = c("Full", "PP", "S*+", "S++", "IR")
      estimands = estimands %>%
        filter(Estimand != "Delta")
    }else if(!IR_display & delta_display){
      Estimands_in = c("Full", "PP", "S*+", "S++", "Delta")
      estimands = estimands %>%
        filter(Estimand != "IR")
    }else{
      Estimands_in = c("Full", "PP", "S*+", "S++")
      estimands = estimands %>%
        filter(Estimand != "IR",
               Estimand != "Delta")
    }

    p_estimands = ggplot(estimands %>%
                           #filter(Estimand != "S*+") %>%
                           mutate(Timepoints = ordered(factor(paste0("Time = ", Timepoints)),
                                                       levels = paste0("Time = ", Timepoints_in)),
                                  Estimand = ordered(factor(Estimand),
                                                     levels = Estimands_in),
                                  Arm = paste0(Arm, " vs Arm ", reference_id)),
                         aes(x=Estimand, y=mean, color=Estimand, shape=Estimand)) +
      geom_point(position=position_dodge(.1), size=3) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                    width=0.1,
                    linewidth=0.5,
                    position=position_dodge(.1))+
      # geom_line(aes(group = Estimand),
      #           lwd=1.5, alpha=0.75,
      #           position=position_dodge(.1)) +
      #ggtitle("Causal Estimands") +
      #geom_boxplot(aes(fill=Estimand)) +
      facet_grid(Arm ~ Timepoints) +
      #xlab("Time") +
      ylab("Response") +
      xlab("")+
      theme_bw() +
      #guides(color=guide_legend("")) +
      theme(axis.text.x = element_text(angle=90, hjust=0.95, vjust=0.2),
            plot.margin = unit(c(0,5,1,1),
                               units="points"),
            legend.position="bottom") +
      scale_color_brewer(palette="Dark2")

  }

  if(greyscale == TRUE){
    p_estimands = p_estimands +
      scale_color_grey(start=0, end=0.6)
  }

  if(static_output == TRUE){
    print(p_estimands)
    return(estimands)
  }else{
    if(normal_output == TRUE){
      print(ggplotly(p_estimands)%>%
              layout(legend = list(orientation = "h", x = 0, y = -0.15)))
      return(estimands)
    }else{
      print(ggplotly(p_estimands)%>%
              layout(legend = list(orientation = "h", x = 0, y = -0.15)))
    }
  }

}
