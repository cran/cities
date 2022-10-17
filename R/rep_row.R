#' @title rep_row
#'
#' @description Helper function to repeat a matrix by row
#'
#' @param x vector to repeat
#' @param n number of repetions
#'
#' @return Matrix with vector x repeated n-times by rows.
#'
#' @examples
#'set.seed(1)
#'rep_row(rnorm(5), 5)
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes
#'

rep_row = function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
