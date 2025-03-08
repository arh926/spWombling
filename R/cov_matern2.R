#' Mat\'ern(\eqn{\nu=5/2}) covariance kernel
#'
#' @param Delta distance matrix
#' @param phi spatial range
#' @param sig2 spatial variance
#' @keywords cov_matern1
#######################
### Matérn nu = 5/2 ###
#######################
cov_matern2 <- function(Delta = NULL,
                        phi = NULL,
                        sig2 = NULL){
  # covariance
  Sigma = sig2 * (1 + phi * sqrt(5) * Delta + phi^2 * 5 * Delta^2/3) * exp(-phi * sqrt(5) * Delta)
  
  return(Sigma = Sigma)
}
