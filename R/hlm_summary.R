#' Model Summary for spatial hierarchical Bayes models
#'
#' Provides a model summary for fitted hlmBayes_sp models
#' 
#' @param chain MCMC posterior samples
#' @param niter number of MCMC iterations
#' @param nburn number of burn-in samples
#' @param thin thinning interval
#' @param latent_inf if true returns inference for spatial effects \eqn{Z}
#' @keywords 
#' @import stats
#' @importFrom coda as.mcmc
#' @importFrom coda is.mcmc
#' @importFrom coda HPDinterval
#' @export
#' @examples 
##################################################################
### Hierarchical Bayesian spatial model: point referenced data ###
##################################################################
hlm_summary <- function(chain = NULL,
                        nburn = NULL,
                        niter = NULL,
                        thin = NULL,
                        latent_inf = T){
  if(is.null(chain)) stop("Chain not supplied!")
  if(is.null(nburn)) nburn <- length(chain)/2
  if(is.null(thin)) thin <- 1
  
  if(latent_inf){
    # names.par <- names(unlist(lapply(chain$parameters, function(x) if(is.mcmc(x)) 1)))
    summary.out.pars <- do.call(rbind,lapply(chain$parameters, function(x){
      if(is.mcmc(x)){
        if(is.null(dim(x))){
          c(median(x[seq((nburn+1),niter,by=thin)]),
            HPDinterval(as.mcmc(x[seq((nburn+1),niter,by=thin)])))
        }else{
          cbind.data.frame(apply(x,2,median), HPDinterval(as.mcmc(x[seq((nburn+1),niter,by=thin),])))
        }
      }
    }))
    
    summary.out.latent <- do.call(rbind,lapply(chain$latent.effect, function(x) cbind.data.frame(apply(x,2,median),
                                                                                                 HPDinterval(as.mcmc(x[seq((nburn+1),niter,by=thin),])))))
    colnames(summary.out.pars) <- colnames(summary.out.latent) <- c("median","lower.hpd","upper.hpd")
    
    return(list(summary.pars=summary.out.pars,
                summary.latent=summary.out.latent))
  }else{
    summary.out.pars <- do.call(rbind,lapply(chain$parameters, function(x){
      if(is.mcmc(x)){
        if(is.null(dim(x))){
          c(median(x[seq((nburn+1),niter,by=thin)]),
            HPDinterval(as.mcmc(x[seq((nburn+1),niter,by=thin)])))
        }else{
          cbind.data.frame(apply(x,2,median), HPDinterval(as.mcmc(x[seq((nburn+1),niter,by=thin),])))
        }
      }
    }))
    colnames(summary.out.pars) <- c("median","lower.hpd","upper.hpd")
    # if(verbose) print(summary.out.pars)
    return(list(summary.pars=summary.out.pars))
  }
}
