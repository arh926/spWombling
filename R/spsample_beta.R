#' Posterior samples for spatio-temporal effects and regression coefficients
#'
#' Posterior sampling algorithm for generating samples of spatio-temporal random effects and
#' regression coefficients using a Gibbs sampling algorithm.
#'
#' @param y response
#' @param X design matrix
#' @param coords spatial co-ordinates
#' @param t temporal coordinates
#' @param phis posterior samples of the spatial range parameter
#' @param phit posterior samples of the temporal range parameter
#' @param sig2 posterior samples of the spatio-temporal variance parameter (partial sill)
#' @param tau2 posterior samples of the error variance (nugget)
#' @param cov.type type of covariance kernel being used. Choices include Exponential, Mat\'ern(\eqn{\nu=3/2})), Mat\'ern(\eqn{\nu=5/2}), Gaussian
#' @param silent logical argument for print-statements
#' @keywords spsample_beta
#' @import parallel
#' @export
spsample <- function(y = NULL,
                     X = NULL,
                     coords = NULL,
                     phi = NULL,
                     sig2 = NULL,
                     tau2 = NULL,
                     cov.type = c("matern1", "matern2" ,"gaussian"),
                     silent = TRUE){
  if(is.null(cov.type)) stop("Please provide a covariance kernel!")
  N = length(y)
  Delta = as.matrix(dist(coords))
  XtX = t(X) %*% X
  
  niter = length(phi)
  
  ncores <- detectCores() - 1
  samp.list <- split(1:niter, ceiling(seq_along(1:niter)/(niter/ncores)))
  parallel.index <- 1:ncores
  
  z.beta.list = mclapply(parallel.index, function(x){
    id.x = samp.list[[x]]
    n.x = length(id.x)
    
    z.mcmc = matrix(0, ncol = N, nrow = n.x)
    beta0.mcmc = rep(0, n.x)
    beta.mcmc = matrix(0, nrow = n.x, ncol = ncol(X)) # X should not have intercept
    
    sig2.x <- sig2[id.x]
    phi.x <- phi[id.x]
    tau2.x <- tau2[id.x]
    
    for(i in 1:n.x){
      if(cov.type == "matern1"){
        R.Z = cov_matern1(Delta = Delta,
                          phi = phi.x[i],
                          sig2 = 1)
      }
      if(cov.type == "matern2"){
        R.Z = cov_matern2(Delta = Delta,
                          phi = phi.x[i],
                          sig2 = 1)
      }
      if(cov.type == "gaussian"){
        R.Z = cov_gaussian(Delta = Delta,
                           phi = phi.x[i],
                           sig2 = 1)
      }
      
      ############################################
      # Error Handling:: in case R.Z is singular #
      ############################################
      R.inv.Z = try(chol2inv(chol(R.Z)), silent  = TRUE)
      if("try-error" %in% class(R.inv.Z)) R.inv.Z = chol2inv(chol(R.Z + 1e-4 * diag(N)))
      
      ##############
      # Covariates #
      ##############
      post_sd_beta <- chol2inv(chol(1e-4 * diag(ncol(X)) + XtX/tau2.x[i]))
      post_mean_beta <- post_sd_beta %*% crossprod(X, y - mean(y))/tau2.x[i]
      beta.mcmc[i,] <- beta <- as.vector(crossprod(chol(post_sd_beta), rnorm(ncol(X))) + post_mean_beta)
      
      ######################
      # Gibbs Update for Z #
      ######################
      Sig.Z.in = R.inv.Z/sig2.x[i] + diag(N)/tau2.x[i]
      chol.Sig.Z.in = try(chol(Sig.Z.in), silent = TRUE)
      if("try-error" %in% class(chol.Sig.Z.in)){
        cat("Bad Iteration::", i, "\n")
        if(i == 1) z.mcmc[i,] = 0
        else if(i == 2) z.mcmc[i,] = z.mcmc[(i - 1),]
        else if(i > 0) z.mcmc[i,] = apply(z.mcmc[1:(i - 1),], 2, median)
      }else{
        Sig.Z = chol2inv(chol.Sig.Z.in)
        mu.Z = crossprod(Sig.Z, y - X %*% beta)/tau2.x[i]
        z.mcmc[i,] = z = as.vector(crossprod(chol(Sig.Z), rnorm(N)) + mu.Z)
      }
      
      #######################################
      # Intercept = mean of spatial effects #
      #######################################
      beta0.mcmc[i] = beta0 = mean(z)
    }
    return(list(z = z.mcmc,
                beta0 = beta0.mcmc,
                beta = beta.mcmc))
  }, mc.cores = ncores)
  
  z.mcmc = do.call(rbind, lapply(z.beta.list, function(x) x$z))
  beta0.mcmc = unlist(lapply(z.beta.list, function(x) x$beta0))
  beta.mcmc = unlist(lapply(z.beta.list, function(x) x$beta))
  
  return(list(z = z.mcmc,
              beta0 = beta0.mcmc,
              beta = beta.mcmc))
}
