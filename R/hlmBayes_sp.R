#' A spatial hierarchical Bayes Markov Chain Monte Carlo sampler
#'
#' Fits a univariate Gaussian spatial regression model: \eqn{Y(s)=x(s)^T\beta+Z(s)+\epsilon}. Parameters not listed are optional.
#'
#' @param coords coordinates for observed process (order \eqn{L} x \eqn{2})
#' @param y observed response (order \eqn{L} x  \eqn{1})
#' @param X a matrix of covariates (order \eqn{L} x  \eqn{p})
#' @param z_init starting values of spatial effects (order \eqn{L} x  \eqn{1})
#' @param D (use if replication at co-ordinate level) index for observations
#' @param phis_init starting value for \eqn{\phi}
#' @param lower_phis lower bound for uniform prior on \eqn{\phi}
#' @param upper_phis upper bound for uniform prior on \eqn{\phi}
#' @param sigma2_init starting value for \eqn{\phi}
#' @param shape_sigma lower bound for uniform prior on \eqn{\phi}
#' @param scale_sigma upper bound for uniform prior on \eqn{\phi}
#' @param tau2_init starting value for \eqn{\phi}
#' @param shape_tau lower bound for uniform prior on \eqn{\phi}
#' @param scale_tau upper bound for uniform prior on \eqn{\phi}
#' @param beta_init starting value for \eqn{\phi}
#' @param mean_beta lower bound for uniform prior on \eqn{\phi}
#' @param prec_beta upper bound for uniform prior on \eqn{\phi}
#' @param spatial.rand.effect if false perform linear regression
#' @param replication refer to parameter D above
#' @param verbose if true prints output for batches
#' @param nu_est (if cov.type not specified turns to true) smoothness parameter for Mat\'ern
#' @param stepnu_init tuning parameter for \eqn{\nu}
#' @param steps_init tuning parameter for \eqn{\phi}
#' @param nu_init starting value for \eqn{\nu}
#' @param lower_nu lower bound for uniform prior on \eqn{\nu}
#' @param upper_nu upper bound for uniform prior on \eqn{\nu}
#' @param niter number of MCMC iterations
#' @param nburn number of burn-in samples
#' @param report batch length
#' @param cov.type covariance type (three available choices: Gaussian, Mat\'ern(\eqn{\nu=3/2})), Mat\'ern(\eqn{\nu=5/2})
#' @keywords spWombling
#' @import stats Matrix
#' @importFrom MASS mvrnorm 
#' @export
##################################################################
### Hierarchical Bayesian spatial model: point referenced data ###
##################################################################
# Heirarchical Bayes 
hlmBayes_sp <- function(coords = NULL,
                        y = NULL, X = NULL,
                        z_init = NULL,
                        D = NULL,  
                        phis_init = NULL, lower_phis = NULL, upper_phis = NULL, # spatial correlation
                        sigma2_init = NULL, shape_sigma = NULL, scale_sigma = NULL, # overall variance
                        tau2_init = NULL, shape_tau = NULL, scale_tau = NULL, # microscopic variance
                        beta_init = NULL, mean_beta = NULL, prec_beta = NULL, # regression parameters
                        steps_init = NULL, niter = NULL, nburn=NULL, report = NULL,
                        cov.type = c("exponential","gaussian","matern1","matern2"),
                        spatial.rand.effect = T,
                        replication = F,
                        verbose = T,
                        nu_est = F,
                        stepnu_init = NULL,
                        nu_init = NULL, lower_nu = NULL, upper_nu = NULL){
  if(is.null(z_init)) L <- N <- length(y) else L <- length(z_init); N <- length(y)
  
  # some housekeeping
  Delta <- as.matrix(dist(coords))
  d.factor <- 1e-10*diag(L)
  if(is.null(nburn)) nburn <- niter/2
  if(is.null(steps_init)) steps_init <- 1
  if(is.null(X)){
    X <- matrix(1,nrow=L); p <- 1
    XtX <- crossprod(X,X)
  }else{
    XtX <- crossprod(X,X)
  }
  if(is.null(D)){
    D <- 1:L
    DtD <- diag(L)
  }else{
    DtD <- diag(table(D))
  }
  if(is.null(z_init)) z_init <- rep(0,L)
  if(is.null(lower_phis)) lower_phis <- 3/max(Delta)
  if(is.null(upper_phis)) upper_phis <- 300
  if(is.null(phis_init)) phis_init <- runif(1,lower_phis,lower_phis+0.01)
  if(is.null(sigma2_init)) sigma2_init <- 5
  if(is.null(shape_sigma)) shape_sigma <- 2
  if(is.null(scale_sigma)) scale_sigma <- 1
  if(is.null(tau2_init)) tau2_init <- 1
  if(is.null(shape_tau)) shape_tau <- 2
  if(is.null(scale_tau)) scale_tau <- 1
  if(is.null(beta_init)) beta_init <- rep(0,ncol(X))
  if(is.null(mean_beta)) mean_beta <- rep(0,ncol(X))
  if(is.null(prec_beta)) prec_beta <- 1e-6*diag(ncol(X))
  if(is.null(cov.type)){
    nu_est <- T
    lower_nu <- 1
    upper_nu <- 3
  } 
  #browser()
  p <- ncol(X)
  # initialize storage
  res_phis <- res_sigma2 <- res_tau2 <- res_nu <- rep(NA,(niter+1))
  res_beta <- matrix(NA, nrow=(niter+1),ncol=p)
  res_z <- matrix(NA, nrow=(niter+1),ncol=L)
  
  # starting values
  res_phis[1] <- phis <- phis_init
  res_sigma2[1] <- sigma2 <- sigma2_init
  res_tau2[1] <- tau2 <- tau2_init
  res_beta[1,] <- beta <- beta_init
  res_z[1,] <- z <- z_init
  if(nu_est) res_nu[1] <- nu <- nu_init
  
  # MCMC Updates for Parameters 
  accepts_vec <- acceptnu_vec <- c(0)
  accepts  <- acceptnu <- 0
  steps <- steps_init
  if(nu_est) stepnu <- stepnu_init
  if(nu_est){
    R <- Matrix(2^(1-nu)/gamma(nu)*(sqrt(2*nu)*phis*Delta)^(nu)*besselK(as.matrix(sqrt(2*nu)*phis*Delta),nu = nu))
    diag(R) <- 1
    chol.R <- chol(R)
    R.inv <- chol2inv(chol.R)
    E <- diag(chol.R)^2
  }else{
    if(cov.type=="exponential"){
      R <- exp(-phis*Delta)+d.factor
      chol.R <- chol(R)
      R.inv <- chol2inv(chol.R)
      E <- diag(chol.R)^2
    } 
    if(cov.type=="gaussian"){
      R <- exp(-phis*Delta^2)+d.factor
      chol.R <- chol(R)
      R.inv <- chol2inv(chol.R)
      E <- diag(chol.R)^2
    } 
    if(cov.type=="matern1"){
      R <- (1+phis*sqrt(3)*Delta)*exp(-phis*sqrt(3)*Delta)+d.factor
      chol.R <- chol(R)
      R.inv <- chol2inv(chol.R)
      E <- diag(chol.R)^2
    } 
    if(cov.type=="matern2"){
      R <- (1+phis*sqrt(5)*Delta+phis^2*5*Delta^2/3)*exp(-phis*sqrt(5)*Delta)+d.factor
      chol.R <- chol(R)
      R.inv <- chol2inv(chol.R)
      E <- diag(chol.R)^2
    } 
  }
  
  # if(verbose) cat("Iteration ",1,"\n")
  
  for(i in 2:(niter+1)){
    # update beta
    post_sd_beta <- chol2inv(chol(prec_beta+XtX/tau2))
    post_mean_beta <- post_sd_beta%*%(prec_beta%*%mean_beta+crossprod(X,y-z[D])/tau2)
    res_beta[i,] <- beta <- as.vector(mvrnorm(1,post_mean_beta,post_sd_beta))
    
    # update tau
    post_shape_tau <- shape_tau+N/2
    post_rate_tau <- 1/scale_tau+sum((y-crossprod(t(X),beta)-z[D])^2)/2
    res_tau2[i] <- tau2 <- 1/rgamma(1,shape = post_shape_tau,rate = post_rate_tau)
    
    if(spatial.rand.effect){
      # update Z
      Sig.Z.in <- R.inv/sigma2+DtD/tau2
      Sig.Z <- chol2inv(chol(Sig.Z.in))
      if(replication) mu.Z <- crossprod(Sig.Z,aggregate(as.vector(y-crossprod(t(X),beta)),list(D),sum)[,2])/tau2
      else mu.Z <- crossprod(Sig.Z,as.vector(y-crossprod(t(X),beta)))/tau2
      Z <- crossprod(chol(Sig.Z),rnorm(L)) + mu.Z
      res_z[i,] <- z <- as.vector(t(Z)) - mean(Z) # change here
      
      # update sigma
      post_shape_sigma <- shape_sigma+L/2
      post_rate_sigma <- 1/scale_sigma+crossprod(t(crossprod(z,R.inv)),z)/2
      res_sigma2[i] <- sigma2 <- 1/rgamma(1,shape = post_shape_sigma,rate = as.numeric(post_rate_sigma))
      
      # metropolis update phis
      # proposal choices: runif(1,-1,1)
      phis.draw <- rnorm(1)*steps+phis 
      if(phis.draw<0)  accept.prob <- -Inf
      else{
        if(nu_est){
          R.draw_s <- Matrix(2^(1-nu)/gamma(nu)*(sqrt(2*nu)*phis.draw*Delta)^(nu)*besselK(as.matrix(sqrt(2*nu)*phis.draw*Delta),nu = nu))
          diag(R.draw_s) <- 1
          chol.R.draw_s <- chol(R.draw_s)
          R.draw.inv_s <- chol2inv(chol.R.draw_s)
          E.draw_s <- diag(chol.R.draw_s)^2
        }else{
          if(cov.type=="exponential"){
            R.draw_s <- exp(-phis.draw*Delta)+d.factor
            chol.R.draw_s <- chol(R.draw_s)
            R.draw.inv_s <- chol2inv(chol.R.draw_s)
            E.draw_s <- diag(chol.R.draw_s)^2
          } 
          if(cov.type=="gaussian"){
            R.draw_s <- exp(-phis.draw*Delta^2)+d.factor
            chol.R.draw_s <- chol(R.draw_s)
            R.draw.inv_s <- chol2inv(chol.R.draw_s)
            E.draw_s <- diag(chol.R.draw_s)^2
          } 
          if(cov.type=="matern1"){
            R.draw_s <- (1+phis.draw*sqrt(3)*Delta)*exp(-phis.draw*sqrt(3)*Delta)+d.factor
            chol.R.draw_s <- chol(R.draw_s)
            R.draw.inv_s <- chol2inv(chol.R.draw_s)
            E.draw_s <- diag(chol.R.draw_s)^2
          } 
          if(cov.type=="matern2"){
            R.draw_s <- (1+phis.draw*sqrt(5)*Delta+phis.draw^2*5*Delta^2/3)*exp(-phis.draw*sqrt(5)*Delta)+d.factor
            chol.R.draw_s <- chol(R.draw_s)
            R.draw.inv_s <- chol2inv(chol.R.draw_s)
            E.draw_s <- diag(chol.R.draw_s)^2
          }  
        }
        
        ra <- sum(log(sqrt(E/E.draw_s)))
        rb <- -crossprod(t(crossprod(z,(R.draw.inv_s-R.inv))),z)/(2*sigma2)
        accept.prob <- min(ra+rb+ifelse((phis.draw>lower_phis & phis.draw<upper_phis),0,-Inf),0)
      }
      
      
      if(log(runif(1)) < accept.prob){
        res_phis[i] <- phis <- phis.draw
        R <- R.draw_s; R.inv <- R.draw.inv_s; E <- E.draw_s
        accepts <- accepts+1
      }else{
        res_phis[i] <- phis
      }
      
      # metropolis update nu
      # proposal choices: runif(1,-1,1)
      if(nu_est){
        nu.draw <- rnorm(1)*stepnu+nu 
        if(nu.draw<0)  accept.prob <- -Inf
        else{
          R.draw_nu <- Matrix(2^(1-nu.draw)/gamma(nu.draw)*(sqrt(2*nu.draw)*phis*Delta)^(nu.draw)*besselK(as.matrix(sqrt(2*nu.draw)*phis*Delta),nu = nu.draw))
          diag(R.draw_nu) <- 1
          chol.R.draw_nu <- chol(R.draw_nu)
          R.draw.inv_nu <- chol2inv(chol.R.draw_nu)
          E.draw_nu <- diag(chol.R.draw_nu)^2
          
          ra <- sum(log(sqrt(E/E.draw_nu)))
          rb <- -crossprod(t(crossprod(z,(R.draw.inv_nu-R.inv))),z)/(2*sigma2)
          accept.prob <- min(ra+rb+ifelse((nu.draw>lower_nu & nu.draw<upper_nu),0,-Inf),0)
        }
        
        if(log(runif(1)) < accept.prob){
          res_nu[i] <- nu <- nu.draw
          R <- R.draw_nu; R.inv <- R.draw.inv_nu; E <- E.draw_nu
          acceptnu <- acceptnu+1
        }else{
          res_nu[i] <- nu
        }
      }
    }else{
      z <- rep(0,L)
    }
    # Adaptive MHS: defaults at around 40%
    if(i%%report == 0){
      accepts <- accepts/report; accepts_vec <- c(accepts_vec,accepts)
      accepts <- max(0.1667, min(accepts,0.75))
      if(accepts>0.5) steps <- steps*accepts/0.5
      else if(accepts<0.25) steps <- steps*accepts/0.25
      
      if(nu_est){
        acceptnu <- acceptnu/report; acceptnu_vec <- c(acceptnu_vec,acceptnu)
        acceptnu <- max(0.1667, min(acceptnu,0.75))
        if(acceptnu>0.5) stepnu <- stepnu*acceptnu/0.5
        else if(accepts<0.25) stepnu <- stepnu*acceptnu/0.25
      }
      if(verbose){
        cat("Iteration ",i,"\n")
        if(i <= nburn){
          if(nu_est){
            cat("-------------------------","\n",
                "phi.s:","\t",round(median(res_phis[1:i]),3),"\n",
                "nu:","\t",round(median(res_nu[1:i]),3),"\n",
                "sigma.2:","\t",round(median(res_sigma2[1:i]),3),"\n",
                "tau.2:","\t",round(median(res_tau2[1:i]),3),"\n",
                "Acceptance Rate (Batch):","\t",accepts*100,"%","\t",acceptnu*100,"%","\n",
                "Overall Acceptance Rate:", median(accepts_vec[1:(i/report)])*100,"%",median(acceptnu_vec[1:(i/report)])*100,"%","\n",
                "-------------------------","\n")
          }else{
            cat("-------------------------","\n",
                "phi.s:","\t",round(median(res_phis[1:i]),3),"\n",
                "sigma.2:","\t",round(median(res_sigma2[1:i]),3),"\n",
                "tau.2:","\t",round(median(res_tau2[1:i]),3),"\n",
                "Acceptance Rate (Batch):","\t",accepts*100,"%","\n",
                "Overall Acceptance Rate:", median(accepts_vec[1:(i/report)])*100,"%","\n",
                "-------------------------","\n")
          }
        }else{
          if(nu_est){
            cat("-------------------------","\n",
                "phi.s:","\t",round(median(res_phis[nburn:i]),3),"\n",
                "nu:","\t",round(median(res_nu[nburn:i]),3),"\n",
                "sigma.2:","\t",round(median(res_sigma2[nburn:i]),3),"\n",
                "tau.2:","\t",round(median(res_tau2[nburn:i]),3),"\n",
                "Acceptance Rate (Batch):","\t",accepts*100,"%","\t",acceptnu*100,"%","\n",
                "Overall Acceptance Rate:", median(accepts_vec[(nburn/report):(i/report)])*100,"%", median(acceptnu_vec[(nburn/report):(i/report)])*100,"%","\n",
                "-------------------------","\n")
          }else{
            cat("-------------------------","\n",
                "phi.s:","\t",round(median(res_phis[nburn:i]),3),"\n",
                "sigma.2:","\t",round(median(res_sigma2[nburn:i]),3),"\n",
                "tau.2:","\t",round(median(res_tau2[nburn:i]),3),"\n",
                "Acceptance Rate (Batch):","\t",accepts*100,"%","\n",
                "Overall Acceptance Rate:", median(accepts_vec[(nburn/report):(i/report)])*100,"%","\n",
                "-------------------------","\n")
          }
        }
        accepts <- acceptnu <- 0
      }
    }
  }
  if(spatial.rand.effect){
    if(nu_est){
      return(list(parameters = list(post_phis=as.mcmc(res_phis),
                                    post_sigma2=as.mcmc(res_sigma2),
                                    post_tau2=as.mcmc(res_tau2),
                                    post_beta=as.mcmc(res_beta),
                                    post_nu=as.mcmc(res_nu)),
                  latent.effect = list(post_z=as.mcmc(res_z)),
                  ar_s=accepts_vec))
    }else{
      return(list(parameters = list(post_phis=as.mcmc(res_phis),
                                    post_sigma2=as.mcmc(res_sigma2),
                                    post_tau2=as.mcmc(res_tau2),
                                    post_beta=as.mcmc(res_beta)),
                  latent.effect = list(post_z=as.mcmc(res_z)),
                  ar_s=accepts_vec))
    }
  }else{
    return(list(post_tau2=as.mcmc(res_tau2),
                post_beta=as.mcmc(res_beta)))
  }
}