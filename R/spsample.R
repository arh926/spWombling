#' Posterior samples for spatial effects and regression coefficients
#'
#' Posterior sampling algorithm for generating samples of spatial random effects and
#' regression coefficients using a Gibbs sampling algorithm.
#'
#' @param y response
#' @param coords spatial co-ordinates
#' @param phi posterior samples of the spatial range parameter
#' @param sig2 posterior samples of the spatial variance parameter (partial sill)
#' @param tau2 posterior samples of the error variance (nugget)
#' @param cov.type type of covariance kernel being used. Choices include Exponential, Mat\'ern(\eqn{\nu=3/2})), Mat\'ern(\eqn{\nu=5/2}), Gaussian
#' @param silent logical argument for print-statements
#' @keywords spsample
#' @import parallel
#' @export
spsample <- function(y = NULL,
                     coords = NULL,
                     phi = NULL,
                     sig2 = NULL,
                     tau2 = NULL,
                     cov.type = c("matern1", "matern2" ,"gaussian"),
                     silent = TRUE){
  if(is.null(cov.type)) stop("Please provide a covariance kernel!")
  N = length(y)
  Delta = as.matrix(dist(coords))
  
  niter = length(phi)
  
  ncores <- detectCores() - 1
  samp.list <- split(1:niter, ceiling(seq_along(1:niter)/(niter/ncores)))
  parallel.index <- 1:ncores
  
  z.beta.list = mclapply(parallel.index, function(x){
    id.x = samp.list[[x]]
    n.x = length(id.x)
    
    z.mcmc = matrix(0, ncol = N, nrow = n.x)
    beta.mcmc = rep(0, n.x)
    
    sig2.x = sig2[id.x]
    phi.x = phi[id.x]
    tau2.x = tau2[id.x]
    
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
        mu.Z = crossprod(Sig.Z, y)/tau2.x[i]
        z.mcmc[i,] = z = as.vector(crossprod(chol(Sig.Z), rnorm(N)) + mu.Z)
      }
      
      #######################################
      # Intercept = mean of spatial effects #
      #######################################
      beta.mcmc[i] = beta = mean(z)
    }
    return(list(z = z.mcmc,
                beta = beta.mcmc))
  }, mc.cores = ncores)
  
  z.mcmc = do.call(rbind, lapply(z.beta.list, function(x) x$z))
  beta.mcmc = unlist(lapply(z.beta.list, function(x) x$beta))
  
  return(list(z = z.mcmc,
              beta = beta.mcmc))
}