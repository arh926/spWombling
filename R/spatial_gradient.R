#' Gradient and Curvature Assessment
#'
#' performs gradient and curvature assessment on the estimated surface from a hierarchcial Bayesian spatial model. Parameters not listed are optional
#'
#' @param coords coordinates for observed process (order \eqn{L} x \eqn{2})
#' @param chain the posterior samples from the MCMC fit
#' @param cov.type covariance type (three available choices: Gaussian, Mat\'ern(\eqn{\nu=3/2})), Mat\'ern(\eqn{\nu=5/2})
#' @param grid.points coordinates for grid over observed process (order \eqn{n_G} x \eqn{2})
#' @param samples index of samples from the chain to be used
#' @param nbatch number of batches
#' @param niter number of MCMC iterations
#' @param nburn number of burn-in samples
#' @param return.mcmc if true returns mcmc samples of gardients
#' @import stats
#' @importFrom coda as.mcmc
#' @importFrom coda HPDinterval
#' @importFrom MASS mvrnorm
#' @importFrom Matrix Matrix
#' @importFrom Matrix forceSymmetric
#' @import parallel doParallel
#' @keywords 
#' @export
#' @examples


spatial_gradient <- function(coords = NULL,
                             chain = NULL,
                             cov.type = c("gaussian","matern1","matern2"),
                             grid.points = NULL, # on which gradient is to be calculated
                             samples = NULL,
                             niter = NULL,
                             nbatch = 200,
                             nburn = NULL,
                             return.mcmc = TRUE){
  sysinfo = Sys.info()
  Delta <- as.matrix(dist(coords))
  d.factor <- 1e-10 * diag(nrow(Delta))
  if(is.null(niter)) niter <- length(chain$parameters$post_phis)
  if(is.null(nburn)) nburn <- niter/2
  if(is.null(samples)) samples <- (nburn + 1):niter
  
  # parallelizing gradient calculation on MCMC chain
  numCores <- detectCores()
  ncores <- numCores-1
  parallel.index <- 1:ncores
  samp.list <- split(samples, ceiling(seq_along(samples)/(length(samples)/ncores)))
  
  
  dist.s0 <- sapply(1:nrow(grid.points),function(y) apply(coords,1,function(x) sqrt(sum((x-grid.points[y,])^2)) ))
  delta.s0 <- sapply(1:nrow(grid.points), function(y) t(apply(coords,1,function(x) x-grid.points[y,])), simplify = F)
  
  if(sysinfo["sysname"] == "Windows"){
    
    cl <- makeCluster(ncores)
    registerDoParallel(cores = ncores)
    if(cov.type == "gaussian"){
      results.grad <- foreach(x = parallel.index) %do% {
        samp.x <- samp.list[[x]]
        post_phi_thin <- chain$parameters$post_phis[samp.x]
        post_sigma2_thin <- chain$parameters$post_sigma2[samp.x]
        # post_beta_thin <- chain$parameters$post_beta[samp.x]
        post_z_thin <- chain$latent.effect$post_z[samp.x,]
        
        
        mcmc.grad <- list()
        for(i.mcmc in 1:length(post_phi_thin)){
          
          phi.grad.est <- post_phi_thin[i.mcmc]
          sig2.grad.est <- post_sigma2_thin[i.mcmc]
          z.grad.est <- post_z_thin[i.mcmc,]
          # beta.grad.est <- post_beta_thin[i.mcmc]
          
          Sig.Z.grad.est <- sig2.grad.est * exp(-phi.grad.est * Delta^2) + d.factor
          s.grad.in <- chol2inv(chol(Sig.Z.grad.est))
          
          grad.est <- matrix(NA, nrow = nrow(grid.points),ncol = 5)
          
          for(i in 1:nrow(grid.points)){
            # gaussian covariance
            nabla.K <- cbind(-2 * sig2.grad.est * phi.grad.est * exp(-phi.grad.est * dist.s0[,i]^2) * delta.s0[[i]], # gradient
                             -2 * sig2.grad.est * phi.grad.est * exp(-phi.grad.est * dist.s0[,i]^2) * (1 - 2 * phi.grad.est * delta.s0[[i]][, 1]^2), # -ve curvature-11
                             4 * phi.grad.est^2 * sig2.grad.est * exp(-phi.grad.est * dist.s0[,i]^2) * delta.s0[[i]][,1] * delta.s0[[i]][, 2], # -ve curvature-12
                             -2 * sig2.grad.est * phi.grad.est * exp(-phi.grad.est * dist.s0[,i]^2) * (1 - 2 * phi.grad.est * delta.s0[[i]][, 2]^2)) #-ve curvature-22
            V.0 <- sig2.grad.est * diag(c(2 * phi.grad.est,
                                          2 * phi.grad.est,
                                          12 * phi.grad.est^2,
                                          4 * phi.grad.est^2,
                                          12 * phi.grad.est^2))
            nabla.K.t <- t(cbind(2 * sig2.grad.est * phi.grad.est * exp(-phi.grad.est * dist.s0[, i]^2) * delta.s0[[i]], # gradient
                                 -2 * sig2.grad.est * phi.grad.est * exp(-phi.grad.est * dist.s0[, i]^2) * (1 - 2 * phi.grad.est * delta.s0[[i]][, 1]^2), # -ve curvature-11
                                 4 * phi.grad.est^2 * sig2.grad.est * exp(-phi.grad.est * dist.s0[, i]^2) * delta.s0[[i]][,1] * delta.s0[[i]][, 2], # -ve curvature-12
                                 -2 * sig2.grad.est * phi.grad.est * exp(-phi.grad.est * dist.s0[, i]^2) * (1 - 2 * phi.grad.est * delta.s0[[i]][, 2]^2))) #-ve curvature-22
            
            tmp <- t(crossprod(t(nabla.K.t), s.grad.in))
            mean.grad <- crossprod(tmp, z.grad.est)
            var.grad <- forceSymmetric(V.0 + crossprod(tmp, nabla.K))
            
            grad.est[i,] <- as.numeric(mvrnorm(1, mean.grad, var.grad))
          }
          mcmc.grad[[i.mcmc]] <- grad.est
        }
        return(mcmc.grad)
      }
    }else if(cov.type == "matern2"){
      results.grad <- foreach(x = parallel.index) %do% {
        samp.x <- samp.list[[x]]
        post_phi_thin <- chain$parameters$post_phis[samp.x]
        post_sigma2_thin <- chain$parameters$post_sigma2[samp.x]
        # post_beta_thin <- chain$parameters$post_beta[samp.x]
        post_z_thin <- chain$latent.effect$post_z[samp.x,]
        #chain$post_beta[samp.x]
        
        
        mcmc.grad <- list()
        for(i.mcmc in 1:length(post_phi_thin)){
          
          phi.grad.est <- post_phi_thin[i.mcmc]
          sig2.grad.est <- post_sigma2_thin[i.mcmc]
          z.grad.est <- post_z_thin[i.mcmc,]
          
          
          Sig.Z.grad.est <- (1 + sqrt(5) * phi.grad.est * Delta + 5 * phi.grad.est^2 * Delta^2/3) * exp(-sqrt(5) * phi.grad.est * Delta) + d.factor
          s.grad.in <- chol2inv(chol(Sig.Z.grad.est))
          
          grad.est <- matrix(NA, nrow = nrow(grid.points), ncol = 5)
          for(i in 1:nrow(grid.points)){
            #s0 <- grid.points[i,]
            #dist.s0 <- apply(coords,1,function(x) sqrt(sum((x-s0)^2)) )
            #delta.s0 <- t(apply(coords,1,function(x) x-s0 ))
            
            # matern2 covariance
            nabla.K <- 5 * cbind(-phi.grad.est^2 * (1 + sqrt(5) * phi.grad.est * dist.s0[,i]) * exp(-sqrt(5) * phi.grad.est * dist.s0[,i]) * delta.s0[[i]],
                                 -phi.grad.est^2 * exp(-sqrt(5) * phi.grad.est * dist.s0[,i]) * (1 + sqrt(5) * phi.grad.est*dist.s0[,i] - 5 * phi.grad.est^2 * delta.s0[[i]][,1]^2),
                                 5 * phi.grad.est^4 * exp(-sqrt(5) * phi.grad.est * dist.s0[,i]) * delta.s0[[i]][,1] * delta.s0[[i]][,2],
                                 -phi.grad.est^2 * exp(-sqrt(5) * phi.grad.est * dist.s0[,i]) * (1 + sqrt(5) * phi.grad.est * dist.s0[,i] - 5 * phi.grad.est^2 * delta.s0[[i]][,2]^2))/3
            V.0 <- 5 * diag(c(phi.grad.est^2/3,
                              phi.grad.est^2/3,
                              5 * phi.grad.est^4,
                              5 * phi.grad.est^4/3,
                              5 * phi.grad.est^4))
            nabla.K.t <- t(5 * cbind(phi.grad.est^2 * (1 + sqrt(5) * phi.grad.est * dist.s0[,i]) * exp(-sqrt(5) * phi.grad.est * dist.s0[,i]) * delta.s0[[i]],
                                     - phi.grad.est^2 * exp(-sqrt(5) * phi.grad.est * dist.s0[,i]) * (1 + sqrt(5) * phi.grad.est * dist.s0[,i] - 5 * phi.grad.est^2 * delta.s0[[i]][,1]^2),
                                     5 * phi.grad.est^4 * exp(-sqrt(5) * phi.grad.est * dist.s0[,i]) * delta.s0[[i]][,1]*delta.s0[[i]][,2],
                                     - phi.grad.est^2 * exp(-sqrt(5) * phi.grad.est * dist.s0[,i]) * (1 + sqrt(5) * phi.grad.est * dist.s0[,i]- 5 * phi.grad.est^2 * delta.s0[[i]][,2]^2))/3)
            
            tmp <- t(crossprod(t(nabla.K.t),s.grad.in))
            mean.grad <- crossprod(tmp,z.grad.est)
            var.grad <- forceSymmetric(V.0 + crossprod(tmp,nabla.K))
            
            grad.est[i,] <- as.numeric(sqrt(sig2.grad.est) * chol(var.grad) %*% rnorm(5) + mean.grad)
          }
          mcmc.grad[[i.mcmc]] <- grad.est
        }
        return(mcmc.grad)
      }
    }else if(cov.type == "matern1"){
      results.grad <- foreach(x = parallel.index) %do% {
        samp.x <- samp.list[[x]]
        post_phi_thin <- chain$parameters$post_phis[samp.x]
        post_sigma2_thin <- chain$parameters$post_sigma2[samp.x]
        # post_beta_thin <- chain$parameters$post_beta[samp.x]
        post_z_thin <- chain$latent.effect$post_z[samp.x,]
        
        
        mcmc.grad <- list()
        for(i.mcmc in 1:length(post_phi_thin)){
          
          phi.grad.est <- post_phi_thin[i.mcmc]
          sig2.grad.est <- post_sigma2_thin[i.mcmc]
          z.grad.est <- post_z_thin[i.mcmc,]
          # beta.grad.est <- post_beta_thin[i.mcmc]
          
          Sig.Z.grad.est <- (1+phi.grad.est*sqrt(3)*Delta)*exp(-phi.grad.est*sqrt(3)*Delta)+d.factor
          s.grad.in <- chol2inv(chol(Sig.Z.grad.est))
          
          grad.est <- matrix(NA, nrow=nrow(grid.points),ncol=2)
          for(i in 1:nrow(grid.points)){
            # matern1 covariance
            nabla.K <- -3*phi.grad.est^2*exp(-phi.grad.est*sqrt(3)*dist.s0[,i])*delta.s0[[i]]
            V.0 <- 3*phi.grad.est^2*diag(2)
            nabla.K.t <- -t(nabla.K)
            
            tmp <- t(crossprod(t(nabla.K.t),s.grad.in))
            mean.grad <- crossprod(tmp,z.grad.est)
            var.grad <- forceSymmetric(V.0-crossprod(tmp,nabla.K))
            
            grad.est[i,] <- as.numeric(sqrt(sig2.grad.est)*chol(var.grad)%*%rnorm(2)+mean.grad)
          }
          
          mcmc.grad[[i.mcmc]] <- grad.est
        }
        return(mcmc.grad)
      }
    }else{stop("Enter valid covariance! Choices are:: gaussian, matern1, matern2!")}
    
    stopCluster(cl)
    
  }else{
    if(cov.type == "gaussian"){
      results.grad <- mclapply(parallel.index, function(x){
        samp.x <- samp.list[[x]]
        post_phi_thin <- chain$parameters$post_phis[samp.x]
        post_sigma2_thin <- chain$parameters$post_sigma2[samp.x]
        # post_beta_thin <- chain$parameters$post_beta[samp.x]
        post_z_thin <- chain$latent.effect$post_z[samp.x,]
        
        
        mcmc.grad <- list()
        for(i.mcmc in 1:length(post_phi_thin)){
          
          phi.grad.est <- post_phi_thin[i.mcmc]
          sig2.grad.est <- post_sigma2_thin[i.mcmc]
          z.grad.est <- post_z_thin[i.mcmc,]
          # beta.grad.est <- post_beta_thin[i.mcmc]
          
          Sig.Z.grad.est <- sig2.grad.est*exp(-phi.grad.est*Delta^2)+d.factor
          s.grad.in <- chol2inv(chol(Sig.Z.grad.est))
          
          grad.est <- matrix(NA, nrow=nrow(grid.points),ncol=5)
          
          for(i in 1:nrow(grid.points)){
            # gaussian covariance
            nabla.K <- cbind(-2*sig2.grad.est*phi.grad.est*exp(-phi.grad.est*dist.s0[,i]^2)*delta.s0[[i]], # gradient
                             -2*sig2.grad.est*phi.grad.est*exp(-phi.grad.est*dist.s0[,i]^2)*(1-2*phi.grad.est*delta.s0[[i]][,1]^2), # -ve curvature-11
                             4*phi.grad.est^2*sig2.grad.est*exp(-phi.grad.est*dist.s0[,i]^2)*delta.s0[[i]][,1]*delta.s0[[i]][,2], # -ve curvature-12
                             -2*sig2.grad.est*phi.grad.est*exp(-phi.grad.est*dist.s0[,i]^2)*(1-2*phi.grad.est*delta.s0[[i]][,2]^2)) #-ve curvature-22
            V.0 <- sig2.grad.est*diag(c(2*phi.grad.est,
                                        2*phi.grad.est,
                                        12*phi.grad.est^2,
                                        4*phi.grad.est^2,
                                        12*phi.grad.est^2))
            nabla.K.t <- t(cbind(2*sig2.grad.est*phi.grad.est*exp(-phi.grad.est*dist.s0[,i]^2)*delta.s0[[i]], # gradient
                                 -2*sig2.grad.est*phi.grad.est*exp(-phi.grad.est*dist.s0[,i]^2)*(1-2*phi.grad.est*delta.s0[[i]][,1]^2), # -ve curvature-11
                                 4*phi.grad.est^2*sig2.grad.est*exp(-phi.grad.est*dist.s0[,i]^2)*delta.s0[[i]][,1]*delta.s0[[i]][,2], # -ve curvature-12
                                 -2*sig2.grad.est*phi.grad.est*exp(-phi.grad.est*dist.s0[,i]^2)*(1-2*phi.grad.est*delta.s0[[i]][,2]^2))) #-ve curvature-22
            
            tmp <- t(crossprod(t(nabla.K.t),s.grad.in))
            mean.grad <- crossprod(tmp, z.grad.est)
            var.grad <- forceSymmetric(V.0+crossprod(tmp,nabla.K))
            
            grad.est[i,] <- as.numeric(mvrnorm(1,mean.grad,var.grad))
          }
          mcmc.grad[[i.mcmc]] <- grad.est
        }
        return(mcmc.grad)
      }, mc.cores = numCores)
    }else if(cov.type=="matern1"){
      results.grad <- mclapply(parallel.index, function(x){
        samp.x <- samp.list[[x]]
        post_phi_thin <- chain$parameters$post_phis[samp.x]
        post_sigma2_thin <- chain$parameters$post_sigma2[samp.x]
        # post_beta_thin <- chain$parameters$post_beta[samp.x]
        post_z_thin <- chain$latent.effect$post_z[samp.x,]
        
        
        mcmc.grad <- list()
        for(i.mcmc in 1:length(post_phi_thin)){
          
          phi.grad.est <- post_phi_thin[i.mcmc]
          sig2.grad.est <- post_sigma2_thin[i.mcmc]
          z.grad.est <- post_z_thin[i.mcmc,]
          # beta.grad.est <- post_beta_thin[i.mcmc]
          
          Sig.Z.grad.est <- (1+phi.grad.est*sqrt(3)*Delta)*exp(-phi.grad.est*sqrt(3)*Delta)+d.factor
          s.grad.in <- chol2inv(chol(Sig.Z.grad.est))
          
          grad.est <- matrix(NA, nrow=nrow(grid.points),ncol=2)
          for(i in 1:nrow(grid.points)){
            # matern1 covariance
            nabla.K <- -3*phi.grad.est^2*exp(-phi.grad.est*sqrt(3)*dist.s0[,i])*delta.s0[[i]]
            V.0 <- 3*phi.grad.est^2*diag(2)
            nabla.K.t <- -t(nabla.K)
            
            tmp <- t(crossprod(t(nabla.K.t),s.grad.in))
            mean.grad <- crossprod(tmp,z.grad.est)
            var.grad <- forceSymmetric(V.0-crossprod(tmp,nabla.K))
            
            grad.est[i,] <- as.numeric(sqrt(sig2.grad.est)*chol(var.grad)%*%rnorm(2)+mean.grad)
          }
          
          mcmc.grad[[i.mcmc]] <- grad.est
        }
        return(mcmc.grad)
      }, mc.cores = numCores)
    }else if(cov.type=="matern2"){
      results.grad <- mclapply(parallel.index, function(x){
        samp.x <- samp.list[[x]]
        post_phi_thin <- chain$parameters$post_phis[samp.x]
        post_sigma2_thin <- chain$parameters$post_sigma2[samp.x]
        # post_beta_thin <- chain$parameters$post_beta[samp.x]
        post_z_thin <- chain$latent.effect$post_z[samp.x,]
        #chain$post_beta[samp.x]
        
        
        mcmc.grad <- list()
        for(i.mcmc in 1:length(post_phi_thin)){
          
          phi.grad.est <- post_phi_thin[i.mcmc]
          sig2.grad.est <- post_sigma2_thin[i.mcmc]
          z.grad.est <- post_z_thin[i.mcmc,]
          
          
          Sig.Z.grad.est <- (1+sqrt(5)*phi.grad.est*Delta+5*phi.grad.est^2*Delta^2/3)*exp(-sqrt(5)*phi.grad.est*Delta)+d.factor
          s.grad.in <- chol2inv(chol(Sig.Z.grad.est))
          
          grad.est <- matrix(NA, nrow=nrow(grid.points),ncol=5)
          for(i in 1:nrow(grid.points)){
            #s0 <- grid.points[i,]
            #dist.s0 <- apply(coords,1,function(x) sqrt(sum((x-s0)^2)) )
            #delta.s0 <- t(apply(coords,1,function(x) x-s0 ))
            
            # matern2 covariance
            nabla.K <- 5*cbind(-phi.grad.est^2 * (1 + sqrt(5) * phi.grad.est * dist.s0[,i]) * exp(-sqrt(5) * phi.grad.est * dist.s0[,i]) * delta.s0[[i]],
                               -phi.grad.est^2 * exp(-sqrt(5)*phi.grad.est*dist.s0[,i])*(1+sqrt(5)*phi.grad.est*dist.s0[,i]-5*phi.grad.est^2*delta.s0[[i]][,1]^2),
                               5*phi.grad.est^4 * exp(-sqrt(5)*phi.grad.est*dist.s0[,i])*delta.s0[[i]][,1]*delta.s0[[i]][,2],
                               -phi.grad.est^2 * exp(-sqrt(5)*phi.grad.est*dist.s0[,i])*(1+sqrt(5)*phi.grad.est*dist.s0[,i]-5*phi.grad.est^2*delta.s0[[i]][,2]^2))/3
            V.0 <- 5*diag(c(phi.grad.est^2/3,
                            phi.grad.est^2/3,
                            5*phi.grad.est^4,
                            5*phi.grad.est^4/3,
                            5*phi.grad.est^4))
            nabla.K.t <- t(5*cbind(phi.grad.est^2*(1+sqrt(5)*phi.grad.est*dist.s0[,i])*exp(-sqrt(5)*phi.grad.est*dist.s0[,i])*delta.s0[[i]],
                                   -phi.grad.est^2*exp(-sqrt(5)*phi.grad.est*dist.s0[,i])*(1+sqrt(5)*phi.grad.est*dist.s0[,i]-5*phi.grad.est^2*delta.s0[[i]][,1]^2),
                                   5*phi.grad.est^4*exp(-sqrt(5)*phi.grad.est*dist.s0[,i])*delta.s0[[i]][,1]*delta.s0[[i]][,2],
                                   -phi.grad.est^2*exp(-sqrt(5)*phi.grad.est*dist.s0[,i])*(1+sqrt(5)*phi.grad.est*dist.s0[,i]-5*phi.grad.est^2*delta.s0[[i]][,2]^2))/3)
            
            tmp <- t(crossprod(t(nabla.K.t),s.grad.in))
            mean.grad <- crossprod(tmp,z.grad.est)
            var.grad <- forceSymmetric(V.0+crossprod(tmp,nabla.K))
            
            grad.est[i,] <- as.numeric(sqrt(sig2.grad.est)*chol(var.grad)%*%rnorm(5)+mean.grad)
          }
          mcmc.grad[[i.mcmc]] <- grad.est
        }
        return(mcmc.grad)
      }, mc.cores = numCores)
    }else print("Error:: Enter valid covariance")
    
  }
  
  
  if(cov.type=="matern1"){
    grad.s1.temp <- grad.s2.temp <- c()
    for(i in 1:length(results.grad)){
      grad.s1.temp <- rbind(grad.s1.temp,
                            do.call(rbind,lapply(results.grad[[i]],function(x) x[,1])))
      grad.s2.temp <- rbind(grad.s2.temp,
                            do.call(rbind,lapply(results.grad[[i]],function(x) x[,2])))
    }
    # compute median gradients
    ngrad <- samples-(samples[1]-1)
    slist <- split(ngrad, ceiling(seq_along(ngrad)/(length(ngrad)/nbatch)))
    grad.s1.temp1 <- t(sapply(slist, function(x) apply(grad.s1.temp[x,],2,median)))
    grad.s2.temp1 <- t(sapply(slist, function(x) apply(grad.s2.temp[x,],2,median)))
    
    grad.s1.mcmc <- as.mcmc(grad.s1.temp1)
    grad.s2.mcmc <- as.mcmc(grad.s2.temp1)
    
    grad.s1.est <- apply(grad.s1.mcmc,2,median)
    grad.s2.est <- apply(grad.s2.mcmc,2,median)
    
    # compute HPD for median gradients
    grad.s1.hpd <- round(cbind(grad.s1.est, HPDinterval(grad.s1.mcmc)),2)
    grad.s2.hpd <- round(cbind(grad.s2.est, HPDinterval(grad.s2.mcmc)),2)
    
    if(return.mcmc){
      return(list(grad1.est=grad.s1.hpd,
                  grad2.est=grad.s2.hpd,
                  grad.s1.mcmc=grad.s1.temp,
                  grad.s2.mcmc=grad.s2.temp))
    }else{
      return(list(grad1.est=grad.s1.hpd,
                  grad2.est=grad.s2.hpd))
    }
  }else{
    grad.s1.temp <- grad.s2.temp <- grad.s11.temp <- grad.s12.temp <- grad.s22.temp <- c()
    for(i in 1:length(results.grad)){
      grad.s1.temp <- rbind(grad.s1.temp,
                            do.call(rbind,lapply(results.grad[[i]],function(x) x[,1])))
      grad.s2.temp <- rbind(grad.s2.temp,
                            do.call(rbind,lapply(results.grad[[i]],function(x) x[,2])))
      grad.s11.temp <- rbind(grad.s11.temp,
                             do.call(rbind,lapply(results.grad[[i]],function(x) x[,3])))
      grad.s12.temp <- rbind(grad.s12.temp,
                             do.call(rbind,lapply(results.grad[[i]],function(x) x[,4])))
      grad.s22.temp <- rbind(grad.s22.temp,
                             do.call(rbind,lapply(results.grad[[i]],function(x) x[,5])))
    }
    ngrad <- samples-(samples[1]-1)
    slist <- split(ngrad, ceiling(seq_along(ngrad)/(length(ngrad)/nbatch)))
    grad.s1.temp1 <- t(sapply(slist, function(x) apply(grad.s1.temp[x,],2,median)))
    grad.s2.temp1 <- t(sapply(slist, function(x) apply(grad.s2.temp[x,],2,median)))
    grad.s11.temp1 <- t(sapply(slist, function(x) apply(grad.s11.temp[x,],2,median)))
    grad.s12.temp1 <- t(sapply(slist, function(x) apply(grad.s12.temp[x,],2,median)))
    grad.s22.temp1 <- t(sapply(slist, function(x) apply(grad.s22.temp[x,],2,median)))
    
    grad.s1.mcmc <- as.mcmc(grad.s1.temp1)
    grad.s2.mcmc <- as.mcmc(grad.s2.temp1)
    grad.s11.mcmc <- as.mcmc(grad.s11.temp1)
    grad.s12.mcmc <- as.mcmc(grad.s12.temp1)
    grad.s22.mcmc <- as.mcmc(grad.s22.temp1)
    
    
    grad.s1.est <- apply(grad.s1.mcmc,2,median)
    grad.s2.est <- apply(grad.s2.mcmc,2,median)
    grad.s11.est <- apply(grad.s11.mcmc,2,median)
    grad.s12.est <- apply(grad.s12.mcmc,2,median)
    grad.s22.est <- apply(grad.s22.mcmc,2,median)
    
    # gradient-x
    grad.s1.hpd <- data.frame(round(cbind(grad.s1.est, HPDinterval(grad.s1.mcmc)),4))
    # gradient-y
    grad.s2.hpd <- data.frame(round(cbind(grad.s2.est, HPDinterval(grad.s2.mcmc)),4))
    # curvature-entries
    grad.s11.hpd <- data.frame(round(cbind(grad.s11.est, HPDinterval(grad.s11.mcmc)),4))
    grad.s12.hpd <- data.frame(round(cbind(grad.s12.est, HPDinterval(grad.s12.mcmc)),4))
    grad.s22.hpd <- data.frame(round(cbind(grad.s22.est, HPDinterval(grad.s22.mcmc)),4))
    
    
    
    
    if(return.mcmc){
      return(list(grad1.est=grad.s1.hpd, #gradient
                  grad2.est=grad.s2.hpd, #gradient
                  grad11.est=grad.s11.hpd,
                  grad12.est=grad.s12.hpd,
                  grad22.est=grad.s22.hpd,
                  grad.s1.mcmc=grad.s1.temp,
                  grad.s2.mcmc=grad.s2.temp,
                  grad.s11.mcmc=grad.s11.temp,
                  grad.s12.mcmc=grad.s12.temp,
                  grad.s22.mcmc=grad.s22.temp))
    }else{
      return(list(grad1.est=grad.s1.hpd, #gradient
                  grad2.est=grad.s2.hpd, #gradient
                  grad11.est=grad.s11.hpd,
                  grad12.est=grad.s12.hpd,
                  grad22.est=grad.s22.hpd))
    }
  }
}
