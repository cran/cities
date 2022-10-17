#' @title rep_rcol
#'
#' @description Helper function to repeat a matrix by column
#'
#' @param x vector to repeat
#' @param n number of repetions
#'
#' @return matrix with vector x repeated n-times by columns.
#'
#' @examples
#'set.seed(1)
#'rep_col(rnorm(5), 5)
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes
#'

rep_col = function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
