#' Gaussian covariance kernel
#'
#' @param Delta distance matrix
#' @param phi spatial range
#' @param sig2 spatial variance
#' @keywords cov_gaussian
################
### Gaussian ###
################
cov_gaussian <- function(Delta = NULL,
                        phi = NULL,
                        sig2 = NULL){
  # covariance
  Sigma = sig2 * exp(-phi * Delta^2)
  
  return(Sigma = Sigma)
}
