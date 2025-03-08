#' Mat\'ern(\eqn{\nu}) covariance kernel
#'
#' @param Delta distance matrix
#' @param phi spatial range
#' @param sig2 spatial variance
#' @param nu fractal parameter
#' @keywords cov_matern1
#################
### Matérn nu ###
#################
cov_matern <- function(Delta = NULL,
                       phi = NULL,
                       sig2 = NULL,
                       nu = NULL){
  # covariance
  Sigma = sig2 * 2^(1 - nu)/gamma(nu) * (sqrt(2 * nu) * phi * Delta)^nu *
    besselK(sqrt(2 * nu) * phi * Delta, nu = nu, expon.scaled = TRUE)/
    exp(sqrt(2 * nu) * phi * Delta)
  
  return(Sigma = Sigma)
}
