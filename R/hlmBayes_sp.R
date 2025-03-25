#' A spatial hierarchical Bayes Markov Chain Monte Carlo sampler
#'
#' Fits a univariate Gaussian spatial regression model: \eqn{Y(s)=x(s)^T\beta+Z(s)+\epsilon}. Parameters not listed are optional.
#'
#' @param coords coordinates for observed process (order \eqn{L} x \eqn{2})
#' @param y observed response (order \eqn{L} x  \eqn{1})
#' @param a.sigma shape parameter for inverse-gamma prior on \eqn{\sigma^2}
#' @param b.sigma scale parameter for inverse-gamma prior on \eqn{\sigma^2}
#' @param a.tau shape parameter for inverse-gamma prior on \eqn{\tau^2}
#' @param b.tau scale parameter for inverse-gamma prior on \eqn{\tau^2}
#' @param lower.phi lower limit for uniform prior on \eqn{\phi}
#' @param upper.phi upper limit for uniform prior on \eqn{\phi}
#' @param verbose if true prints output for batches
#' @param niter number of MCMC iterations
#' @param nburn number of burn-in samples
#' @param report batch length
#' @param trgtFn.compute compute posterior
#' @param digits rounding digits
#' @param cov.type covariance type (three available choices: Gaussian, Mat\'ern(\eqn{\nu=3/2})), Mat\'ern(\eqn{\nu=5/2})
#' @keywords spWombling
#' @import stats 
#' @importFrom MASS mvrnorm 
#' @export
###################################
### Collapsed MH Sampler for SP ###
###################################
hlmBayes_sp <- function(y = NULL,
                        coords = NULL,
                        niter = NULL, nburn = NULL, report = NULL,
                        a.sigma = NULL, b.sigma = NULL,
                        a.tau = NULL, b.tau = NULL,
                        lower.phi = NULL, upper.phi = NULL,
                        cov.type = NULL,
                        verbose = TRUE,
                        trgtFn.compute = FALSE,
                        digits = 3){
  # compute distance matrix
  Delta = as.matrix(dist(coords))
  
  # some housekeeping
  N = length(y)
  # MCMC parameters
  if(is.null(niter)){
    niter = 5e3
    nburn = niter/2
    report = 1e2
  }
  
  # learning rates
  e.sig2 = e.tau2 = e.phi = 1e-1
  # covariance type not specified -> fit Matern(\nu)
  if(is.null(cov.type) | !(cov.type %in% c("matern1", "matern2", "gaussian"))){
    e.nu = 1e-2
    lower.nu = 1e-1
    upper.nu = 1e1
  }
  
  # Default hyper-parameter settings  
  if(is.null(a.sigma)) a.sigma = 2
  if(is.null(a.tau)) a.tau = 2
  if(is.null(b.sigma)) b.sigma = 1
  if(is.null(b.tau)) b.tau = 1
  if(is.null(lower.phi)) lower.phi = 1e-3
  if(is.null(upper.phi)) upper.phi = 30
  
  # initial values
  phi = sig2 = tau2 = 1
  if(is.null(cov.type) | !(cov.type %in% c("matern1", "matern2", "gaussian"))) nu = 1
  else{
    if(cov.type == "matern1") nu = 1.5
    if(cov.type == "matern2") nu = 2.5
    if(cov.type == "gaussian") nu = Inf
  }
  
  # batch acceptance probabilities
  accept.p = rep(0, 3)
  if(is.null(cov.type)) accept.p = rep(0, 4)
  
  # initialize storage
  res_sig2 = res_phi = res_tau2 = rep(0, niter)
  accept_m = c()
  if(is.null(cov.type)) res_nu = rep(0, niter)
  if(trgtFn.compute) trgtFn = rep(0, niter)
  
  if(is.null(cov.type) | !(cov.type %in% c("matern1", "matern2", "gaussian"))){
    R = cov_matern(Delta = Delta, phi = phi, sig2 = 1)
  }else{
    if(cov.type == "gaussian"){
      R = cov_gaussian(Delta = Delta, phi = phi, sig2 = 1)
    } 
    if(cov.type == "matern1"){
      R = cov_matern1(Delta = Delta, phi = phi, sig2 = 1)
    } 
    if(cov.type == "matern2"){
     R = cov_matern2(Delta = Delta, phi = phi, sig2 = 1)
    } 
  }
  
  chol.Sigma = chol(sig2 * R + tau2 * diag(N))
  inv.Sigma = chol2inv(chol.Sigma)
  
  if(trgtFn.compute){
    trgtFn.init = -(a.sigma + 1) * log(sig2) - b.sigma/sig2 - 
      -(a.tau + 1) * log(tau2) - b.tau/tau2 - 
      sum(log(diag(chol.Sigma)))/4 - t(y) %*% inv.Sigma %*% y/2
  }
  
  for(i in 1:niter){
    #################
    # update sigma2 #
    #################
    lsig2.draw = log(sig2) + e.sig2 * rnorm(1) # log-normal proposal
    
    chol.Sigma.draw = chol(exp(lsig2.draw) * R + tau2 * diag(N))
    inv.Sigma.draw = chol2inv(chol.Sigma.draw)
    
    r = -b.sigma * (exp(-lsig2.draw) - 1/sig2) - (a.sigma + 1) * (lsig2.draw - log(sig2)) - 
      crossprod(t(crossprod(y, (inv.Sigma.draw - inv.Sigma))), y)/2 - 
      sum(log(diag(chol.Sigma.draw)/diag(chol.Sigma)))
    
    accept.prob = min(r, 0)
    
    if(log(runif(1)) < accept.prob){
      sig2 = exp(lsig2.draw)
      res_sig2[i] = sig2
      
      chol.Sigma = chol.Sigma.draw
      inv.Sigma = inv.Sigma.draw
      
      accept.p[1] = accept.p[1] + 1
    }else{
      res_sig2[i] = sig2
    }
    
    ###############
    # update tau2 #
    ###############
    ltau2.draw = log(tau2) + e.tau2 * rnorm(1) # log-normal proposal
    
    chol.Sigma.draw = chol(sig2 * R + exp(ltau2.draw) * diag(N))
    inv.Sigma.draw = chol2inv(chol.Sigma.draw)
    
    r = -b.tau * (exp(-ltau2.draw) - 1/tau2) - (a.tau + 1) * (ltau2.draw - log(tau2)) - 
      crossprod(t(crossprod(y, (inv.Sigma.draw - inv.Sigma))), y)/2 - 
      sum(log(diag(chol.Sigma.draw)/diag(chol.Sigma)))
    
    accept.prob = min(r, 0)
    
    if(log(runif(1)) < accept.prob){
      tau2 = exp(ltau2.draw)
      res_tau2[i] = tau2
      
      chol.Sigma = chol.Sigma.draw
      inv.Sigma = inv.Sigma.draw
      
      accept.p[2] = accept.p[2] + 1
    }else{
      res_tau2[i] = tau2
    }
    
    
    lphi.draw =  log(phi) + e.phi * rnorm(1)
    
    if((exp(lphi.draw) < lower.phi) | (exp(lphi.draw) > upper.phi)) accept.prob = -Inf
    else{
      
      if(is.null(cov.type) | !(cov.type %in% c("matern1", "matern2", "gaussian"))){
        R.draw = cov_matern(Delta = Delta, phi = exp(lphi.draw), sig2 = 1)
      }else{
        if(cov.type == "gaussian"){
          R.draw = cov_gaussian(Delta = Delta, phi = exp(lphi.draw), sig2 = 1)
        } 
        if(cov.type == "matern1"){
          R.draw = cov_matern1(Delta = Delta, phi = exp(lphi.draw), sig2 = 1)
        } 
        if(cov.type == "matern2"){
          R.draw = cov_matern2(Delta = Delta, phi = exp(lphi.draw), sig2 = 1)
        } 
      }
      
      chol.Sigma.draw = chol(sig2 * R.draw + tau2 * diag(N))
      inv.Sigma.draw = chol2inv(chol.Sigma.draw)
      
      r = - (lphi.draw - log(phi)) - crossprod(t(crossprod(y, (inv.Sigma.draw - inv.Sigma))), y)/2 - 
        sum(log(diag(chol.Sigma.draw)/diag(chol.Sigma)))
      
      accept.prob = min(r, 0)
      
    }
    
    
    if(log(runif(1)) < accept.prob){
      phi = exp(lphi.draw)
      res_phi[i] = phi
      
      R = R.draw
      chol.Sigma = chol.Sigma.draw
      inv.Sigma = inv.Sigma.draw
      
      accept.p[3] = accept.p[3] + 1
    }else{
      res_phi[i] <- phi
    }
    
    #############
    # update nu #
    #############
    if(is.null(cov.type) | !(cov.type %in% c("matern1", "matern2", "gaussian"))){
      nu.draw = nu + e.nu * rnorm(1) 
      
      R.draw = cov_matern(Delta = Delta, phi = phi, sig2 = 1)
      
      chol.Sigma.draw = chol(sig2 * R.draw + tau2 * diag(N))
      inv.Sigma.draw = chol2inv(chol.Sigma.draw)
      
      r = - crossprod(t(crossprod(y, (inv.Sigma.draw - inv.Sigma))), y)/2 - 
        sum(log(diag(chol.Sigma.draw)/diag(chol.Sigma)))
      
      accept.prob = min(r + ifelse(((nu.draw > lower.nu) & (nu.draw < upper.nu)), 0, -Inf), 0)
      
      if(log(runif(1)) < accept.prob){
        res_nu[i] = nu = nu.draw
        
        R = R.draw
        chol.Sigma = chol.Sigma.draw
        inv.Sigma = inv.Sigma.draw
        
        accept.p[4] = accept.p[4] + 1
      }else{
        res_nu[i] = nu
      }
    }
    
    if(trgtFn.compute){
      trgtFn[i] = -(a.sigma + 1) * log(sig2) - b.sigma/sig2 - 
        -(a.tau + 1) * log(tau2) - b.tau/tau2 - 
        sum(log(diag(chol.Sigma)))/4 - t(y) %*% inv.Sigma %*% y/2
    }
    
    ############@@@@@@@@@@@@@@@##############
    # Adaptive Scaling of Proposal Variance #
    #     MH: optimal scales (33%)          #
    ############@@@@@@@@@@@@@@@##############
    if(i %% report == 0){
      accept.p = accept.p/report
      accept_m = rbind(accept_m, accept.p)
      if(verbose){
        if(is.null(cov.type) | !(cov.type %in% c("matern1", "matern2", "gaussian"))){
          if(i > nburn){
            cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:", round(c(e.sig2, e.tau2, e.phi, e.nu), digits), "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
                "Sigma2::", round(median(res_sig2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_sig2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Tau2::", round(median(res_tau2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_tau2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phi::", round(median(res_phi[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_phi[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phi[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Nu::", round(median(res_nu[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_nu[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_nu[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
          }else{
            cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:",round(c(e.sig2, e.tau2, e.phi, e.nu), digits), "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
                "Sigma2::", round(median(res_sig2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_sig2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Tau2::", round(median(res_tau2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_tau2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phi::", round(median(res_phi[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_phi[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phi[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Nu::", round(median(res_nu[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_nu[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_nu[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
          }
        }else{
          if(i > nburn){
            cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:", round(c(e.sig2, e.tau2, e.phi), digits), "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
                "Sigma2::", round(median(res_sig2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_sig2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Tau2::", round(median(res_tau2[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_tau2[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phi::", round(median(res_phi[nburn:i]), digits = (digits - 1)), "(", round(quantile(res_phi[nburn:i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phi[nburn:i], probs = 0.975), digits = (digits - 1)), ")", "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
          }else{
            cat("Iteration::", i, "Acceptance:", accept.p, "Tuning:",round(c(e.sig2, e.tau2, e.phi), digits), "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n",
                "Sigma2::", round(median(res_sig2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_sig2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_sig2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Tau2::", round(median(res_tau2[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_tau2[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_tau2[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\t",
                "Phi::", round(median(res_phi[(i - report + 1):i]), digits = (digits - 1)), "(", round(quantile(res_phi[(i - report + 1):i], probs = 0.025), digits = (digits - 1)), ",", round(quantile(res_phi[(i - report + 1):i], probs = 0.975), digits = (digits - 1)), ")", "\n",
                "=+++++++++++++++++++++++++++++++++++++++++++++++++=", "\n")
          }
        }
      }
      ####################
      # Proposal Scaling #
      ####################
      step = sapply(accept.p, function(x){
        out = 1
        x = max(0.17, min(x, 0.75))
        if(x > 0.50) out = x/0.50
        else if(x < 0.25) out = x/0.25
        out
      })
      e.sig2 = e.sig2 * step[1]
      e.tau2 = e.tau2 * step[2]
      e.phi = e.phi * step[3]
      if(is.null(cov.type) | !(cov.type %in% c("matern1", "matern2", "gaussian"))) e.nu = e.nu * step[4]
      
      accept.p = rep(0, 3)
      if(is.null(cov.type) | !(cov.type %in% c("matern1", "matern2", "gaussian"))) accept.p = rep(0, 4)
    }
  }
  
  if(is.null(cov.type) | !(cov.type %in% c("matern1", "matern2", "gaussian"))){
    if(trgtFn.compute){
      return(list(y = y,
                  coords = coords,
                  trgt_fn = c(trgtFn.init, trgtFn),
                  nu = res_nu,
                  sig2 = res_sig2[(nburn + 1):niter],
                  tau2 = res_tau2[(nburn + 1):niter],
                  phi = res_phi[(nburn + 1):niter]))
    }else{
      return(list(y = y,
                  coords = coords,
                  nu = res_nu,
                  sig2 = res_sig2[(nburn + 1):niter],
                  tau2 = res_tau2[(nburn + 1):niter],
                  phi = res_phi[(nburn + 1):niter]))
    }
  }else{
    if(trgtFn.compute){
      return(list(y = y,
                  coords = coords,
                  trgt_fn = c(trgtFn.init, trgtFn),
                  nu = nu,
                  sig2 = res_sig2[(nburn + 1):niter],
                  tau2 = res_tau2[(nburn + 1):niter],
                  phi = res_phi[(nburn + 1):niter]))
    }else{
      return(list(y = y,
                  coords = coords,
                  nu = nu,
                  sig2 = res_sig2[(nburn + 1):niter],
                  tau2 = res_tau2[(nburn + 1):niter],
                  phi = res_phi[(nburn + 1):niter]))
    }
  }
}