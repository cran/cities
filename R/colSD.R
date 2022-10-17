#' @title colSD
#'
#' @description Helper function to calculate standard deviation of matrix by columns
#'
#' @param data_in matrix of numeric values
#'
#' @return Vector of standard deivations of columns of data_in.
#'
#' @examples
#'set.seed(1)
#'colSD(matrix(rnorm(100), ncol=5))
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes
#'
colSD = function(data_in){
  return(apply(data_in, 2, sd))
}
