#' @title line_parameters
#'
#' @description Helper function that returns slope and intercept for line
#' equation using two points in the cartesian plot: (x1, x2) and (y1, y2)
#'
#' @param x1 first value of the point (x1, x2) in the cartesian plot
#' @param y1 first value of the point (y1, y2) in the cartesian plot
#' @param x2 second value of the point (x1, x2) in the cartesian plot
#' @param y2 second value of the point (y1, y2) in the cartesian plot
#' @return Vector of slope and intercept for equation of line.
#'
#' @examples
#' line_parameters(1,2,4,2)
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes
#'
line_parameters = function(x1, y1, x2, y2){
  m = (y2-y1)/(x2-x1)
  b = y2-m*x2
  return(c(m,b))
}
