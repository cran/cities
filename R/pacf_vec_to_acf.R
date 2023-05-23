#' @title pacf_vec_to_acf
#'
#' @description Generate correlation matrix from partial autocorrelations
#'
#' @param pacf_vec Vector of partial autocorrelations
#' @param n_repeat number of repeat measures (must be longer than length of pacf_vec)
#' @return Correlation matrix from partial autocorrelations.
#'
#' @examples
#' pacf_vec_to_acf(c(0.5, -0.1), 5)
#' @export
#' @import dplyr ggplot2 plotly tidyr ggthemes
#'
pacf_vec_to_acf = function(pacf_vec,    # vector of partial autocorrelations
                           n_repeat){   # number of repeat measures
  pacf_lag = length(pacf_vec)
  pacf_matrix_padded = c(1, pacf_vec, rep(0,n_repeat - pacf_lag-1))
  # Generate the PACF toeplitz matrix
  pacf_mat = toeplitz(pacf_matrix_padded)
  # Empty ACF matrix
  acf_mat = matrix(0,
                   nrow=n_repeat,
                   ncol=n_repeat)
  # Specify diagonal of ACF is 1
  diag(acf_mat) = 1
  # Based on (Harry Joe, 2006)
  for (k in 1:(n_repeat-1)){
    acf_mat_temp = matrix(0, nrow=n_repeat, ncol=n_repeat)
    for (j in 1:(n_repeat-k)){
      if (k == 1){
        acf_mat_temp[j,j+k] = pacf_mat[j,j+k]
      }
      else{
        r1_index1 = j
        r1_index2 = (j+1):(j+k-1)
        r1_index = cbind(r1_index1, r1_index2)
        r1_vec = as.matrix(acf_mat[r1_index])

        r3_index1 = j+k
        r3_index2 = (j+1):(j+k-1)
        r3_index = cbind(r3_index1, r3_index2)
        r3_vec = as.matrix(acf_mat[r3_index])

        r2_index_temp = (j+1):(j+k-1)
        r2_index = min(r2_index_temp):max(r2_index_temp)
        r2_mat = acf_mat[r2_index,r2_index]

        D_jk = (sqrt(1-t(r1_vec) %*% solve(r2_mat) %*% r1_vec)) * (sqrt(1-t(r3_vec) %*% solve(r2_mat) %*% r3_vec))
        acf_mat_temp[j,j+k] = t(r1_vec) %*% solve(r2_mat) %*% r3_vec + D_jk*pacf_mat[j,j+k]
      }
    }
    acf_mat = acf_mat + acf_mat_temp + t(acf_mat_temp)
  }
  return(acf_mat)
}
