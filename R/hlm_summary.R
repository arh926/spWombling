#' Model Summary for spatial hierarchical Bayes models
#'
#' Provides a model summary for fitted hlmBayes_sp models
#' 
#' @param sig2 MCMC posterior samples for \eqn{\sigma^2}
#' @param tau2 MCMC posterior samples for \eqn{\tau^2}
#' @param phi MCMC posterior samples for \eqn{\phi}
#' @param beta MCMC posterior samples for \eqn{\beta}
#' @param z MCMC posterior samples for \eqn{Z}
#' @keywords hlm_summary
#' @import stats
#' @importFrom coda as.mcmc
#' @importFrom coda is.mcmc
#' @importFrom coda HPDinterval
#' @export
#' @examples 
##################################################################
### Hierarchical Bayesian spatial model: point referenced data ###
##################################################################
hlm_summary <- function(sig2 = NULL,
                        tau2 = NULL,
                        phi = NULL,
                        beta = NULL,
                        z = NULL){
  
  if(is.null(nrow(beta))){
   param_est =  rbind(sig2 = c(median(sig2), HPDinterval(sig2)),
                      tau2 = c(median(tau2), HPDinterval(tau2)),
                      phi = c(median(phi), HPDinterval(phi)),
                      beta = c(median(beta), HPDinterval(beta)))
   colnames(param_est) = c("median", "lower.hpd", "upper.hpd")
  }else{
    p = ncol(beta)
    param_est =  rbind(sig2 = c(median(sig2), HPDinterval(sig2)),
                       tau2 = c(median(tau2), HPDinterval(tau2)),
                       phi = c(median(phi), HPDinterval(phi)),
                       t(apply(beta, 2, function(x) c(median(x), HPDinterval(as.mcmc(x))))))
    colnames(param_est) = c("median", "lower.hpd", "upper.hpd")
  }
  
  z_est = data.frame(t(apply(z, 2, function(x) c(median(x), HPDinterval(as.mcmc(x))))))
  colnames(z_est) = c("median", "lower.hpd", "upper.hpd")
  z_est$sig = apply(z_est, 1, function(x){
    if(x[2] > 0) return (1)
    if(x[3] < 0) return (-1)
    else return(0)
  })
  
  return(list(param = param_est,
              z = z_est))
}
