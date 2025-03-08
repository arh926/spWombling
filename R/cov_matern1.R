#' Mat\'ern(\eqn{\nu=3/2}) covariance kernel
#'
#' @param Delta distance matrix
#' @param phi spatial range
#' @param sig2 spatial variance
#' @keywords cov_matern1
#######################
### Matérn nu = 3/2 ###
#######################
cov_matern1 <- function(Delta = NULL,
                        phi = NULL,
                        sig2 = NULL){
  # covariance
  Sigma = sig2 * (1 + phi * sqrt(3) * Delta) * exp(-phi * sqrt(3) * Delta)
  
  return(Sigma = Sigma)
}
