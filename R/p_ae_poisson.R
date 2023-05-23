#' @title p_ae_poisson
#'
#' @description Helper function that returns probability of discontinuing due
#' to adverse events (AE)
#'
#' @param rate_dc_ae Probability of observing at least one AE
#' @param prob_ae Proportion of discontinuation due to AE
#'
#' @return Probabilities of discontinuing due to AE.
#' @examples
#' p_ae_poisson(c(0.9, 0.8), c(0.1, 0.1))
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes
#'

# Returns probability of discontinuing due to AE via Poisson process
p_ae_poisson = function(rate_dc_ae,
                        prob_ae){
  prob_dc_ae = (log(1-rate_dc_ae))/log(1-prob_ae)
  return(prob_dc_ae)
}
