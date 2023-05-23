#' @title p_loe_ee_function
#'
#' @description Helper function that returns probability of discontinuing due
#' to lack of efficacy (LoE) or excess efficacy (EE) via a piecewise linear
#' function
#'
#' @param z Vector of numeric values, i.e. change from baseline values
#' @param p_max Maximum probability of discontinuing
#' @param z_l The lower (or left) threshold of the piecewise linear function
#' @param p_min Maximum probability of discontinuing (set to 0)
#' @param z_u The upper (or right) threshold of the piecewise linear function
#' @param up_good TRUE if higher outcome values indicate better responses
#' @return Probabilities of discontinuing due to LoE or EE.
#'
#' @examples
#' line_parameters(1,2,4,2)
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes
#'

p_loe_ee_function = function(z,
                             p_max,
                             z_l,
                             p_min = 0,
                             z_u,
                             up_good = TRUE){
  if (up_good == TRUE){
    parameters = line_parameters(x1 = z_l,
                                 y1 = p_max,
                                 x2 = z_u,
                                 y2 = p_min)
    m = parameters[1]
    b = parameters[2]

    p = m*z + b
    set_ceiling = (z <= z_l)
    set_floor = (z > z_u)

    p[set_ceiling] = p_max
    p[set_floor] = p_min
    return(p)
  }else{
    parameters = line_parameters(x1 = z_l,
                                 y1 = p_min,
                                 x2 = z_u,
                                 y2 = p_max)
    m = parameters[1]
    b = parameters[2]

    p = m*z + b
    set_floor = (z <= z_l)
    set_ceiling = (z > z_u)

    p[set_floor] = p_min
    p[set_ceiling] = p_max
    return(p)
  }
}
