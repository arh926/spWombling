#' Curvilinear Bayesian Wombling
#'
#' performs Bayesian wombling on curves within the spatial surface.
#'
#' @param coords coordinates for observed process (order \eqn{L} x \eqn{2})
#' @param model the posterior samples from the MCMC fit
#' @param cov.type covariance type (three available choices: Gaussian, Mat\'ern(\eqn{\nu=3/2})), Mat\'ern(\eqn{\nu=5/2})
#' @param nbatch number of batches
#' @param curve coordinates for the curve
#' @param type specify type or rule of quadrature
#' @param approx if TRUE curve will be approximated
#' @param verbose if TRUE will print status
#' @param rule.len specifies accuracy for quadrature
#' @import stats
#' @importFrom coda as.mcmc
#' @importFrom coda is.mcmc
#' @importFrom coda HPDinterval
#' @importFrom MASS mvrnorm
#' @importFrom Matrix Matrix
#' @keywords 
#' @export
#' @examples

######################################################################
### Bayesian Wombling for Curves effected in parallel for segments ###
######################################################################
# if(c("coda") %in% rownames(installed.packages()) == FALSE) {install.packages("coda")}
# if(c("Matrix") %in% rownames(installed.packages()) == FALSE) {install.packages("Matrix")}
# ifelse("coda" %in% (.packages()),"## Required package:: coda already attached",require(coda))
# ifelse("Matrix" %in% (.packages()),"## Required package:: Matrix already attached",require(Matrix))

bayes_cwomb <- function(coords = NULL,
                        model = NULL,
                        cov.type = c("gaussian","matern1","matern2"),
                        nbatch = NULL,
                        curve = NULL,
                        type = c("rectilinear","riemann.sum"),
                        approx = NULL,
                        verbose = T,
                        rule.len = 10){
  len <- rule.len; Delta <- as.matrix(dist(coords))
  # quadrature rule
  rule.uv <-  seq(0,1,length.out=len)
  rule.mv <- expand.grid(rule.uv,rule.uv)
  # t and u values
  tval <- sapply(1:(nrow(curve)-1), function(x) sqrt(sum((curve[(x+1),]-curve[x,])^2)))
  uval <- t(sapply(1:(nrow(curve)-1), function(x) (curve[(x+1),]-curve[x,])))/tval
  
  if(cov.type=="matern1"){
    if(type=="rectilinear"){
      if(is.null(approx)){
        
        # calculation of the wombling measure
        ntimes <- 1:length(model$phis)
        womb_measure1 <- matrix(NA,nrow=(nrow(curve)-1),ncol=ntimes)
        for(i.mcmc in 1:ntimes){
          # MCMC parameters
          sig2.est <- model$sig2[i.mcmc]
          phi.est <- model$phi[i.mcmc]
          z.est <- model$z[i.mcmc,]
          Sig.Z.grad.est <- cov_matern1(Delta = Delta, sig2 = sig2.est, phi = phi.est)
          s.grad.in <- chol2inv(chol(Sig.Z.grad.est))
          
          for(i in 1:(nrow(curve)-1)){
            s0 <- curve[i,]
            u <- uval[i,]
            ################################################
            # Line-Integral Process: Cross-Covariance term #
            ################################################
            # Univariate Quadrature: Area 
            gamma.cov <- t(apply(coords,1,function(xy){
              sj <- xy
              grad.li <- sapply(rule.uv*tval[i], function(x) Gamma.t.1.m1(x,u=u,s0=s0,sj=sj,phi.est=phi.est))
              sum(grad.li[-len])*(tval[i]/len)
            }))
            tgamma.cov <- t(apply(coords,1,function(xy){
              sj <- xy
              grad.li <- sapply(rule.uv*tval[i], function(x) Gamma.1.m1(x,u=u,s0=s0,sj=sj,phi.est=phi.est))
              sum(grad.li[-len])*(tval[i]/len)
            }))
            
            ###################################################
            # Line-Integral Process: Variance-Covariance term #
            ###################################################
            # Bivariate Quadrature: Volume
            k.g.11 <- sum(apply(rule.mv*tval[i],1,function(x) -K.11.m1(t=x,u=u,phi.est=phi.est))[-len^2])*tval[i]^2/len^2
            k.gamma <- matrix(k.g.11,nrow=1,ncol=1,byrow=T)
            
            tmp <- t(crossprod(tgamma.cov,s.grad.in))
            mean.womb <- crossprod(tmp,z.est)
            var.womb <-  k.gamma-crossprod(tmp,gamma.cov)
            
            wm <- mvrnorm(1,mean.womb,sig2.est*var.womb)
            
            womb_measure1[i,i.mcmc] <- wm[1]
          }
          if(verbose) cat("Iteration",i.mcmc,"\n")
        }
        slist <- split(1:ntimes, ceiling(seq_along(1:ntimes)/(ntimes/nbatch)))
        w1.batch <- do.call(rbind,lapply(slist,function(x) apply(womb_measure1[,x],1,median)))
        w1.est <- cbind.data.frame(apply(w1.batch,2,median),t(apply(w1.batch,2,function(x) HPDinterval(as.mcmc(x)))))
        colnames(w1.est) <- c("median","lower","upper")
        w1.est$sig <- apply(w1.est,1,function(x){
          if(x[2]>0 & x[3]>0) return (1)
          if(x[2]<0 & x[3]<0) return (-1)
          else return(0)
        })
        
        womb.measure.inf <- matrix(c(sum(w1.est[,"median"])/sum(tval),median(rowSums(w1.batch)/sum(tval)), HPDinterval(as.mcmc(rowSums(w1.batch)/sum(tval)))),nrow=1)
        colnames(womb.measure.inf) <- c("line-segment-level","batch-median","batch-lhpd","batch-uhpd")
        rownames(womb.measure.inf) <- c("Gamma1(C)")
        return(list(tval=tval,
                    uval=uval,
                    womb1.mcmc = womb_measure1,
                    womb.grad.inf=w1.est,
                    womb.measure.inf=womb.measure.inf))
      }
    }else if(type=="riemann.sum"){
      u.perp.mat <- t(apply(uval,1,function(x){
        u.perp <- rev(x)
        u.perp[2] <- -u.perp[2]
        u.perp
      }))

      gradient_est <- spatial_gradient(coords = coords,
                                       grid.points = curve,
                                       model = model,
                                       cov.type = cov.type,
                                       nbatch = nbatch,
                                       return.mcmc = TRUE)
      grad.s1.curve <- gradient_est$grad.s1.mcmc
      grad.s2.curve <- gradient_est$grad.s2.mcmc
      
      estval <- estval1 <- c()
      for(j in 1:nrow(grad.s1.curve)){
        estval <- rbind(estval,grad.s1.curve[j,-nrow(curve)]*u.perp.mat[,1]+grad.s2.curve[j,-nrow(curve)]*u.perp.mat[,2])
      }
      slist <- split(1:nrow(grad.s1.curve), ceiling(seq_along(1:nrow(grad.s1.curve))/(length(1:nrow(grad.s1.curve))/nbatch)))
      estval.1 <- do.call(rbind,lapply(slist, function(x) apply(estval[x,],2,median)))
      
      # line-segment level
      estval.inf <- data.frame(unlist(round(t(apply(estval.1,2,function(x){
        x <- as.mcmc(x)
        c(median(x),HPDinterval(x))
      })),3)))
      estval.inf$sig <- apply(estval.inf,1,function(x){
        if(x[2]>0) return (1)
        if(x[3]<0) return (-1)
        else return(0)
      })
      
      w1.c <- c(sum(estval.inf[,1]*tval)/sum(tval),median(apply(estval,1,function(x) sum(x*tval)/sum(tval))), HPDinterval(as.mcmc(apply(estval,1,function(x) sum(x*tval)/sum(tval)))))
      womb.measure.inf <- matrix(w1.c, nrow=1)
      colnames(womb.measure.inf) <- c("line-segment-level","batch-median","batch-lhpd","batch-uhpd")
      rownames(womb.measure.inf) <- c("Gamma1(C)")
      return(list(tval=tval,
                  uval=uval,
                  womb1.mcmc = estval,
                  womb.grad.inf=estval.inf,
                  womb.measure.inf=womb.measure.inf))
    }
  }else if(cov.type=="matern2"){
    if(type=="rectilinear"){
      if(is.null(approx)){
        
        # calculation of the wombling measure
        ntimes <- 1:length(model$phis)
        womb_measure1 <- womb_measure2 <- matrix(NA,nrow=(nrow(curve)-1),ncol=ntimes)
        for(i.mcmc in 1:ntimes){
          # MCMC parameters
          sig2.est <- model$sig2[i.mcmc]
          phi.est <- model$phi[i.mcmc]
          z.est <- model$z[i.mcmc,]
          Sig.Z.grad.est <- cov_matern2(Delta = Delta, sig2 = sig2.est, phi = phi.est)
          s.grad.in <- chol2inv(chol(Sig.Z.grad.est))
          
          for(i in 1:(nrow(curve)-1)){
            s0 <- curve[i,]
            u <- uval[i,]
            ################################################
            # Line-Integral Process: Cross-Covariance term #
            ################################################
            # Univariate Quadrature: Area 
            gamma.cov <- t(apply(coords,1,function(xy){
              sj <- xy
              grad.li <- sapply(rule.uv*tval[i], function(x) Gamma.t.1.m2(x,u=u,s0=s0,sj=sj,phi.est=phi.est))
              curv.li <- sapply(rule.uv*tval[i], function(x) Gamma.t.2.m2(x,u=u,s0=s0,sj=sj,phi.est=phi.est))
              c(sum(grad.li[-len])*(tval[i]/len),
                sum(curv.li[-len])*(tval[i]/len))
            }))
            tgamma.cov <- t(apply(coords,1,function(xy){
              sj <- xy
              grad.li <- sapply(rule.uv*tval[i], function(x) Gamma.1.m2(x,u=u,s0=s0,sj=sj,phi.est=phi.est))
              curv.li <- sapply(rule.uv*tval[i], function(x) Gamma.2.m2(x,u=u,s0=s0,sj=sj,phi.est=phi.est))
              c(sum(grad.li[-len])*(tval[i]/len),
                sum(curv.li[-len])*(tval[i]/len))
            }))
            
            ###################################################
            # Line-Integral Process: Variance-Covariance term #
            ###################################################
            # Bivariate Quadrature: Volume
            k.g.11 <- sum(apply(rule.mv*tval[i],1,function(x) -K.11.m2(t=x,u=u,phi.est=phi.est))[-len^2])*tval[i]^2/len^2
            k.g.12 <- sum(apply(rule.mv*tval[i],1,function(x) K.12.m2(t=x,u=u,phi.est=phi.est))[-len^2])*tval[i]^2/len^2
            k.g.22 <- sum(apply(rule.mv*tval[i],1,function(x) K.22.m2(t=x,u=u,phi.est=phi.est))[-len^2])*tval[i]^2/len^2
            k.g.21 <- -k.g.12
            k.gamma <- matrix(c(k.g.11,k.g.12,k.g.21,k.g.22),nrow=2,ncol=2,byrow=T)
            
            tmp <- t(crossprod(tgamma.cov,s.grad.in))
            mean.womb <- crossprod(tmp,z.est)
            var.womb <-  k.gamma-crossprod(tmp,gamma.cov)
            
            wm <- mvrnorm(1,mean.womb,sig2.est*var.womb)
            
            womb_measure1[i,i.mcmc] <- wm[1]
            womb_measure2[i,i.mcmc] <- wm[2]
          }
          if(verbose) cat("Iteration",i.mcmc,"\n")
        }
        slist <- split(1:ntimes, ceiling(seq_along(1:ntimes)/(ntimes/nbatch)))
        w1.batch <- do.call(rbind,lapply(slist,function(x) apply(womb_measure1[,x],1,median)))
        w2.batch <- do.call(rbind,lapply(slist,function(x) apply(womb_measure2[,x],1,median)))
        w1.est <- cbind.data.frame(apply(w1.batch,2,median),t(apply(w1.batch,2,function(x) HPDinterval(as.mcmc(x)))))
        w2.est <- cbind.data.frame(apply(w2.batch,2,median),t(apply(w2.batch,2,function(x) HPDinterval(as.mcmc(x)))))
        colnames(w1.est) <- colnames(w2.est) <- c("median","lower","upper")
        w1.est$sig <- apply(w1.est,1,function(x){
          if(x[2]>0 & x[3]>0) return (1)
          if(x[2]<0 & x[3]<0) return (-1)
          else return(0)
        })
        w2.est$sig <- apply(w2.est,1,function(x){
          if(x[2]>0 & x[3]>0) return (1)
          if(x[2]<0 & x[3]<0) return (-1)
          else return(0)
        })
        
        womb.measure.inf <- rbind(c(sum(w1.est[,"median"])/sum(tval),median(rowSums(w1.batch)/sum(tval)), HPDinterval(as.mcmc(rowSums(w1.batch)/sum(tval)))),
                                  c(sum(w2.est[,"median"])/sum(tval),median(rowSums(w2.batch)/sum(tval)), HPDinterval(as.mcmc(rowSums(w2.batch)/sum(tval)))))
        colnames(womb.measure.inf) <- c("line-segment-level","batch-median","batch-lhpd","batch-uhpd")
        rownames(womb.measure.inf) <- c("Gamma1(C)","Gamma2(C)")
        return(list(tval=tval,
                    uval=uval,
                    womb1.mcmc = womb_measure1,
                    womb2.mcmc = womb_measure2,
                    womb.grad.inf=w1.est,
                    womb.curv.inf=w2.est,
                    womb.measure.inf=womb.measure.inf))
      }
    }else if(type=="riemann.sum"){
      u.perp.mat <- t(apply(uval,1,function(x){
        u.perp <- rev(x)
        u.perp[2] <- -u.perp[2]
        u.perp
      }))
      
      gradient_est <- spatial_gradient(coords=coords,
                                       grid.points = curve,
                                       model = model,
                                       cov.type = cov.type,
                                       nbatch = nbatch,
                                       return.mcmc = TRUE)
      grad.s1.curve <- gradient_est$grad.s1.mcmc
      grad.s2.curve <- gradient_est$grad.s2.mcmc
      grad.s11.curve <- gradient_est$grad.s11.mcmc
      grad.s12.curve <- gradient_est$grad.s12.mcmc
      grad.s22.curve <- gradient_est$grad.s22.mcmc
      
      estval <- estval1 <- c()
      for(j in 1:nrow(grad.s1.curve)){
        estval <- rbind(estval,grad.s1.curve[j,-nrow(curve)]*u.perp.mat[,1]+grad.s2.curve[j,-nrow(curve)]*u.perp.mat[,2])
        estval1 <- rbind(estval1,grad.s11.curve[j,-nrow(curve)]*u.perp.mat[,1]^2+grad.s22.curve[j,-nrow(curve)]*u.perp.mat[,2]^2+2*grad.s12.curve[j,-nrow(curve)]*u.perp.mat[,1]*u.perp.mat[,2])
      }
      slist <- split(1:nrow(grad.s1.curve), ceiling(seq_along(1:nrow(grad.s1.curve))/(length(1:nrow(grad.s1.curve))/nbatch)))
      estval.1 <- do.call(rbind,lapply(slist, function(x) apply(estval[x,],2,median)))
      estval1.1 <- do.call(rbind,lapply(slist, function(x) apply(estval1[x,],2,median)))
      
      # line-segment level
      estval.inf <- data.frame(unlist(round(t(apply(estval.1,2,function(x){
        x <- as.mcmc(x)
        c(median(x),HPDinterval(x))
      })),3)))
      estval.inf$sig <- apply(estval.inf,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      
      estval1.inf <- data.frame(unlist(round(t(apply(estval1.1,2,function(x){
        x <- as.mcmc(x)
        c(median(x),HPDinterval(x))
      })),3)))
      estval1.inf$sig <- apply(estval1.inf,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      
      w1.c <- c(sum(estval.inf[,1]*tval)/sum(tval),median(apply(estval,1,function(x) sum(x*tval)/sum(tval))), HPDinterval(as.mcmc(apply(estval,1,function(x) sum(x*tval)/sum(tval)))))
      w2.c <- c(sum(estval1.inf[,1]*tval)/sum(tval),median(apply(estval1,1,function(x) sum(x*tval)/sum(tval))), HPDinterval(as.mcmc(apply(estval1,1,function(x) sum(x*tval)/sum(tval)))))
      womb.measure.inf <- rbind(w1.c, w2.c)
      colnames(womb.measure.inf) <- c("line-segment-level","batch-median","batch-lhpd","batch-uhpd")
      rownames(womb.measure.inf) <- c("Gamma1(C)","Gamma2(C)")
      return(list(tval=tval,
                  uval=uval,
                  womb1.mcmc = estval,
                  womb2.mcmc = estval1,
                  womb.grad.inf=estval.inf,
                  womb.curv.inf=estval1.inf,
                  womb.measure.inf=womb.measure.inf))
    }
  }else if(cov.type=="gaussian"){
    if(type=="rectilinear"){
      if(is.null(approx)){
        
        # calculation of the wombling measure
        ntimes <- 1:length(model$phi)
        womb_measure1 <- womb_measure2 <- matrix(NA,nrow=(nrow(curve)-1),ncol=ntimes)
        for(i.mcmc in 1:ntimes){
          # MCMC parameters
          sig2.est <- model$sig2[i.mcmc]
          phi.est <- model$phi[i.mcmc]
          z.est <- model$z[i.mcmc,]
          Sig.Z.grad.est <- cov_gaussian(Delta = Delta, sig2 = sig2.est, phi = phi.est)
          s.grad.in <- chol2inv(chol(Sig.Z.grad.est))
          
          for(i in 1:(nrow(curve)-1)){
            s0 <- curve[i,]
            u <- uval[i,]
            ################################################
            # Line-Integral Process: Cross-Covariance term #
            ################################################
            # Univariate Quadrature: Area 
            gamma.cov <- t(apply(coords,1,function(xy){
              sj <- xy
              grad.li <- sapply(rule.uv*tval[i], function(x) Gamma.t.1.g(x,u=u,s0=s0,sj=sj,phi.est=phi.est))
              curv.li <- sapply(rule.uv*tval[i], function(x) Gamma.t.2.g(x,u=u,s0=s0,sj=sj,phi.est=phi.est))
              c(sum(grad.li[-len])*(tval[i]/len),
                sum(curv.li[-len])*(tval[i]/len))
            }))
            tgamma.cov <- t(apply(coords,1,function(xy){
              sj <- xy
              grad.li <- sapply(rule.uv*tval[i], function(x) Gamma.1.g(x,u=u,s0=s0,sj=sj,phi.est=phi.est))
              curv.li <- sapply(rule.uv*tval[i], function(x) Gamma.2.g(x,u=u,s0=s0,sj=sj,phi.est=phi.est))
              c(sum(grad.li[-len])*(tval[i]/len),
                sum(curv.li[-len])*(tval[i]/len))
            }))
            
            ###################################################
            # Line-Integral Process: Variance-Covariance term #
            ###################################################
            # Bivariate Quadrature: Volume
            k.g.11 <- sum(apply(rule.mv*tval[i],1,function(x) -K.11.g(t=x,u=u,phi.est=phi.est))[-len^2])*tval[i]^2/len^2
            k.g.12 <- sum(apply(rule.mv*tval[i],1,function(x) K.12.g(t=x,u=u,phi.est=phi.est))[-len^2])*tval[i]^2/len^2
            k.g.22 <- sum(apply(rule.mv*tval[i],1,function(x) K.22.g(t=x,u=u,phi.est=phi.est))[-len^2])*tval[i]^2/len^2
            k.g.21 <- -k.g.12
            k.gamma <- matrix(c(k.g.11,k.g.12,k.g.21,k.g.22),nrow=2,ncol=2,byrow=T)
            
            tmp <- t(crossprod(tgamma.cov,s.grad.in))
            mean.womb <- crossprod(tmp,z.est)
            var.womb <-  k.gamma-crossprod(tmp,gamma.cov)
            
            wm <- mvrnorm(1,mean.womb,sig2.est*var.womb)
            
            womb_measure1[i,i.mcmc] <- wm[1]
            womb_measure2[i,i.mcmc] <- wm[2]
          }
          if(verbose) cat("Iteration",i.mcmc,"\n")
        }
        slist <- split(1:ntimes, ceiling(seq_along(1:ntimes)/(ntimes/nbatch)))
        w1.batch <- do.call(rbind,lapply(slist,function(x) apply(womb_measure1[,x],1,median)))
        w2.batch <- do.call(rbind,lapply(slist,function(x) apply(womb_measure2[,x],1,median)))
        w1.est <- cbind.data.frame(apply(w1.batch,2,median),t(apply(w1.batch,2,function(x) HPDinterval(as.mcmc(x)))))
        w2.est <- cbind.data.frame(apply(w2.batch,2,median),t(apply(w2.batch,2,function(x) HPDinterval(as.mcmc(x)))))
        colnames(w1.est) <- colnames(w2.est) <- c("median","lower","upper")
        w1.est$sig <- apply(w1.est,1,function(x){
          if(x[2]>0 & x[3]>0) return (1)
          if(x[2]<0 & x[3]<0) return (-1)
          else return(0)
        })
        w2.est$sig <- apply(w2.est,1,function(x){
          if(x[2]>0 & x[3]>0) return (1)
          if(x[2]<0 & x[3]<0) return (-1)
          else return(0)
        })
        
        womb.measure.inf <- rbind(c(sum(w1.est[,"median"])/sum(tval),median(rowSums(w1.batch)/sum(tval)), HPDinterval(as.mcmc(rowSums(w1.batch)/sum(tval)))),
                                  c(sum(w2.est[,"median"])/sum(tval),median(rowSums(w2.batch)/sum(tval)), HPDinterval(as.mcmc(rowSums(w2.batch)/sum(tval)))))
        colnames(womb.measure.inf) <- c("line-segment-level","batch-median","batch-lhpd","batch-uhpd")
        rownames(womb.measure.inf) <- c("Gamma1(C)","Gamma2(C)")
        return(list(tval=tval,
                    uval=uval,
                    womb1.mcmc = womb_measure1,
                    womb2.mcmc = womb_measure2,
                    womb.grad.inf=w1.est,
                    womb.curv.inf=w2.est,
                    womb.measure.inf=womb.measure.inf))
      }
    }else if(type=="riemann.sum"){
      u.perp.mat <- t(apply(uval,1,function(x){
        u.perp <- rev(x)
        u.perp[2] <- -u.perp[2]
        u.perp
      }))
      
      gradient_est <- spatial_gradient(coords=coords,
                                       grid.points = curve,
                                       model = model,
                                       cov.type = cov.type,
                                       nbatch = nbatch,
                                       return.mcmc = TRUE)
      grad.s1.curve <- gradient_est$grad.s1.mcmc
      grad.s2.curve <- gradient_est$grad.s2.mcmc
      grad.s11.curve <- gradient_est$grad.s11.mcmc
      grad.s12.curve <- gradient_est$grad.s12.mcmc
      grad.s22.curve <- gradient_est$grad.s22.mcmc
      
      estval <- estval1 <- c()
      for(j in 1:nrow(grad.s1.curve)){
        estval <- rbind(estval,grad.s1.curve[j,-nrow(curve)]*u.perp.mat[,1]+grad.s2.curve[j,-nrow(curve)]*u.perp.mat[,2])
        estval1 <- rbind(estval1,grad.s11.curve[j,-nrow(curve)]*u.perp.mat[,1]^2+grad.s22.curve[j,-nrow(curve)]*u.perp.mat[,2]^2+2*grad.s12.curve[j,-nrow(curve)]*u.perp.mat[,1]*u.perp.mat[,2])
      }
      slist <- split(1:nrow(grad.s1.curve), ceiling(seq_along(1:nrow(grad.s1.curve))/(length(1:nrow(grad.s1.curve))/nbatch)))
      estval.1 <- do.call(rbind,lapply(slist, function(x) apply(estval[x,],2,median)))
      estval1.1 <- do.call(rbind,lapply(slist, function(x) apply(estval1[x,],2,median)))
      
      # line-segment level
      estval.inf <- data.frame(unlist(round(t(apply(estval.1,2,function(x){
        x <- as.mcmc(x)
        c(median(x),HPDinterval(x))
      })),3)))
      estval.inf$sig <- apply(estval.inf,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      
      estval1.inf <- data.frame(unlist(round(t(apply(estval1.1,2,function(x){
        x <- as.mcmc(x)
        c(median(x),HPDinterval(x))
      })),3)))
      estval1.inf$sig <- apply(estval1.inf,1,function(x){
        if(x[2]>0 & x[3]>0) return (1)
        if(x[2]<0 & x[3]<0) return (-1)
        else return(0)
      })
      
      w1.c <- c(sum(estval.inf[,1]*tval)/sum(tval),median(apply(estval,1,function(x) sum(x*tval)/sum(tval))), HPDinterval(as.mcmc(apply(estval,1,function(x) sum(x*tval)/sum(tval)))))
      w2.c <- c(sum(estval1.inf[,1]*tval)/sum(tval),median(apply(estval1,1,function(x) sum(x*tval)/sum(tval))), HPDinterval(as.mcmc(apply(estval1,1,function(x) sum(x*tval)/sum(tval)))))
      womb.measure.inf <- rbind(w1.c, w2.c)
      colnames(womb.measure.inf) <- c("line-segment-level","batch-median","batch-lhpd","batch-uhpd")
      rownames(womb.measure.inf) <- c("Gamma1(C)","Gamma2(C)")
      return(list(tval=tval,
                  uval=uval,
                  womb1.mcmc = estval,
                  womb2.mcmc = estval1,
                  womb.grad.inf=estval.inf,
                  womb.curv.inf=estval1.inf,
                  womb.measure.inf=womb.measure.inf))
    }
  }else{
    stop("Enter a valid covariance function. Choices are: matern1, matern2, gaussian")
  }
}

####################################
# Covariance - Terms:: Gaussian    # 
####################################
Gamma.t.1.g <- function(t,u,s0,sj, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.j.t <- -t*u+(sj-s0)
  
  # womb-measure-grad
  t.0 <- -2*phi.est
  t.1 <- phi.est*sum((delta.j.t)^2)
  sum(t.0*exp(-t.1)*delta.j.t*u.perp)
}

Gamma.t.2.g <- function(t,u,s0,sj, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.j.t <- -t*u+(sj-s0)
  
  # womb-measure-curvature
  t.0 <- -2*phi.est
  t.1 <- phi.est*sum((delta.j.t)^2)
  
  k11 <- t.0*exp(-t.1)*(1-2*phi.est*delta.j.t[1]^2)
  k22 <- t.0*exp(-t.1)*(1-2*phi.est*delta.j.t[2]^2)
  k12 <- -2*phi.est*t.0*exp(-t.1)*delta.j.t[1]*delta.j.t[2]
  (k11*u.perp[1]^2+2*k12*u.perp[1]*u.perp[2]+k22*u.perp[2]^2)
}
Gamma.1.g <- function(t,u,s0,sj, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.j.t <- t*u+(s0-sj)
  
  # womb-measure-grad
  t.0 <- -2*phi.est
  t.1 <- phi.est*sum((delta.j.t)^2)
  sum(t.0*exp(-t.1)*delta.j.t*u.perp)
}

Gamma.2.g <- function(t,u,s0,sj, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.j.t <- t*u+(s0-sj)
  
  # womb-measure-curvature
  t.0 <- -2*phi.est
  t.1 <- phi.est*sum((delta.j.t)^2)
  
  k11 <- t.0*exp(-t.1)*(1-2*phi.est*delta.j.t[1]^2)
  k22 <- t.0*exp(-t.1)*(1-2*phi.est*delta.j.t[2]^2)
  k12 <- -2*phi.est*t.0*exp(-t.1)*delta.j.t[1]*delta.j.t[2]
  (k11*u.perp[1]^2+2*k12*u.perp[1]*u.perp[2]+k22*u.perp[2]^2)
}

# Covariance terms
K.11.g <- function(t,u, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.t2.t1 <- (t[2]-t[1])*u
  
  # womb-measure-curvature
  t.0 <- -2*phi.est
  t.1 <- phi.est*sum((delta.t2.t1)^2)
  
  k11 <- t.0*exp(-t.1)*(1-2*phi.est*delta.t2.t1[1]^2)
  k22 <- t.0*exp(-t.1)*(1-2*phi.est*delta.t2.t1[2]^2)
  k12 <- -2*phi.est*t.0*exp(-t.1)*delta.t2.t1[1]*delta.t2.t1[2]
  (k11*u.perp[1]^2+2*k12*u.perp[1]*u.perp[2]+k22*u.perp[2]^2)
}

K.12.g <- function(t,u, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  c.u <- c(u.perp[1]^2,2*u.perp[1]*u.perp[2],u.perp[2]^2)
  delta.t2.t1 <- (t[2]-t[1])*u; norm.delta <- sqrt(sum((delta.t2.t1)^2))
  if(delta.t2.t1[1]==0 & delta.t2.t1[2]==0) return(0)
  else{
    t.1 <- phi.est*norm.delta^2
    k111 <- 4*phi.est^2*exp(-t.1)*(3-2*phi.est*delta.t2.t1[1]^2)*delta.t2.t1[1]
    k112 <- 4*phi.est^2*exp(-t.1)*(1-2*phi.est*delta.t2.t1[1]^2)*delta.t2.t1[2]
    k122 <- 4*phi.est^2*exp(-t.1)*(1-2*phi.est*delta.t2.t1[2]^2)*delta.t2.t1[1]
    k211 <- k112; k212 <- k122
    k222 <- 4*phi.est^2*exp(-t.1)*(3-2*phi.est*delta.t2.t1[2]^2)*delta.t2.t1[2]
    nabla3.K <- matrix(c(k111,k112,k122,k211,k212,k222), nrow=2,ncol=3,byrow=T)
    crossprod(t(crossprod(u.perp,nabla3.K)),c.u)
  }
}

K.22.g <- function(t,u, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  c.u <- c(u.perp[1]^2,2*u.perp[1]*u.perp[2],u.perp[2]^2)
  delta.t2.t1 <- (t[2]-t[1])*u; norm.delta <- sqrt(sum((delta.t2.t1)^2))
  if(delta.t2.t1[1]==0 & delta.t2.t1[2]==0){
    nabla4.K <- 4*diag(c(phi.est^2,
                         phi.est^2/3,
                         phi.est^2))
    return(crossprod(t(crossprod(c.u,nabla4.K)),c.u))
  }else{
    t.1 <- phi.est*norm.delta^2
    k1111 <- 4*phi.est^2*exp(-t.1)*(3-12*phi.est*delta.t2.t1[1]^2+4*phi.est^2*delta.t2.t1[1]^4)
    k1212 <- 4*phi.est^2*exp(-t.1)*(1-2*delta.t2.t1[1]^2)*(1-2*delta.t2.t1[2]^2)
    k2222 <- 4*phi.est^2*exp(-t.1)*(3-12*phi.est*delta.t2.t1[2]^2+4*phi.est^2*delta.t2.t1[2]^4)
    k1112 <- -8*phi.est^3*exp(-t.1)*(3-2*phi.est*delta.t2.t1[1]^2)*delta.t2.t1[1]*delta.t2.t1[2]
    k1122 <- k1211 <- k2211 <- k1212
    k1222 <- -8*phi.est^3*exp(-t.1)*(3-2*phi.est*delta.t2.t1[2]^2)*delta.t2.t1[2]*delta.t2.t1[1]
    k2212 <- k1222
    nabla4.K <- matrix(c(k1111,k1112,k1122,
                         k1211,k1212,k1222,
                         k2211,k2212,k2222), nrow=3,ncol=3,byrow=T)
    return(crossprod(t(crossprod(c.u,nabla4.K)),c.u))
  }
}


####################################
# Covariance - Terms:: Matern(3/2) # 
####################################
Gamma.t.1.m1 <- function(t,u,s0,sj, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.j.t <- -t*u+(sj-s0)
  
  # womb-measure-grad
  t.0 <- -3*phi.est^2
  t.1 <- sqrt(3)*phi.est*sqrt(sum((delta.j.t)^2))
  sum(t.0*exp(-t.1)*delta.j.t*u.perp)
}
Gamma.1.m1 <- function(t,u,s0,sj, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.j.t <- t*u+(s0-sj)
  
  # womb-measure-grad
  t.0 <- -3*phi.est^2
  t.1 <- sqrt(3)*phi.est*sqrt(sum((delta.j.t)^2))
  sum(t.0*exp(-t.1)*delta.j.t*u.perp)
}
# Covariance terms
K.11.m1 <- function(t,u, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.t2.t1 <- (t[2]-t[1])*u
  
  # womb-measure-curvature
  t.0 <- -3*phi.est^2
  t.1 <- sqrt(3)*phi.est*sqrt(sum((delta.t2.t1)^2))
  
  k11 <- t.0*exp(-t.1)*(1-sqrt(3)*phi.est*delta.t2.t1[1]^2/sqrt(sum((delta.t2.t1)^2)))
  k22 <- t.0*exp(-t.1)*(1-sqrt(3)*phi.est*delta.t2.t1[2]^2/sqrt(sum((delta.t2.t1)^2)))
  k12 <- -sqrt(3)*phi.est*t.0*exp(-t.1)*delta.t2.t1[1]*delta.t2.t1[2]/sqrt(sum((delta.t2.t1)^2))
  (k11*u.perp[1]^2+2*k12*u.perp[1]*u.perp[2]+k22*u.perp[2]^2)
}


####################################
# Covariance - Terms:: Matern(5/2) # 
####################################
# Cross-covariance terms
Gamma.t.1.m2 <- function(t,u,s0,sj, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.j.t <- -t*u+(sj-s0)
  
  # womb-measure-grad
  t.0 <- -5*phi.est^2
  t.1 <- sqrt(5)*phi.est*sqrt(sum((delta.j.t)^2))
  sum(t.0*exp(-t.1)*(1+t.1)*delta.j.t/3*u.perp)
}

Gamma.t.2.m2 <- function(t,u,s0,sj, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.j.t <- -t*u+(sj-s0)
  
  # womb-measure-curvature
  t.0 <- -5*phi.est^2
  t.1 <- sqrt(5)*phi.est*sqrt(sum((delta.j.t)^2))
  
  k11 <- t.0*exp(-t.1)*(1+t.1-5*phi.est^2*delta.j.t[1]^2)/3
  k22 <- t.0*exp(-t.1)*(1+t.1-5*phi.est^2*delta.j.t[2]^2)/3
  k12 <- -5*phi.est^2*t.0*exp(-t.1)*delta.j.t[1]*delta.j.t[2]/3
  (k11*u.perp[1]^2+2*k12*u.perp[1]*u.perp[2]+k22*u.perp[2]^2)
}
Gamma.1.m2 <- function(t,u,s0,sj, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.j.t <- t*u+(s0-sj)
  
  # womb-measure-grad
  t.0 <- -5*phi.est^2
  t.1 <- sqrt(5)*phi.est*sqrt(sum((delta.j.t)^2))
  sum(t.0*exp(-t.1)*(1+t.1)*delta.j.t/3*u.perp)
}

Gamma.2.m2 <- function(t,u,s0,sj, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.j.t <- t*u+(s0-sj)
  
  # womb-measure-grad
  t.0 <- -5*phi.est^2
  t.1 <- sqrt(5)*phi.est*sqrt(sum((delta.j.t)^2))
  # womb-measure-curvature
  k11 <- t.0*exp(-t.1)*(1+t.1-5*phi.est^2*delta.j.t[1]^2)/3
  k22 <- t.0*exp(-t.1)*(1+t.1-5*phi.est^2*delta.j.t[2]^2)/3
  k12 <- -5*phi.est^2*t.0*exp(-t.1)*delta.j.t[1]*delta.j.t[2]/3
  (k11*u.perp[1]^2+2*k12*u.perp[1]*u.perp[2]+k22*u.perp[2]^2)
}

# Covariance terms
K.11.m2 <- function(t,u, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  delta.t2.t1 <- (t[2]-t[1])*u
  # womb-measure-grad
  t.0 <- -5*phi.est^2
  t.1 <- sqrt(5)*phi.est*sqrt(sum((delta.t2.t1)^2))
  # womb-measure-curvature
  k11 <- t.0*exp(-t.1)*(1+t.1-5*phi.est^2*delta.t2.t1[1]^2)/3
  k22 <- t.0*exp(-t.1)*(1+t.1-5*phi.est^2*delta.t2.t1[2]^2)/3
  k12 <- -5*phi.est^2*t.0*exp(-t.1)*delta.t2.t1[1]*delta.t2.t1[2]/3
  (k11*u.perp[1]^2+2*k12*u.perp[1]*u.perp[2]+k22*u.perp[2]^2)
}

K.12.m2 <- function(t,u, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  c.u <- c(u.perp[1]^2,2*u.perp[1]*u.perp[2],u.perp[2]^2)
  delta.t2.t1 <- (t[2]-t[1])*u; norm.delta <- sqrt(sum((delta.t2.t1)^2))
  if(delta.t2.t1[1]==0 & delta.t2.t1[2]==0) return(0)
  else{
    t.1 <- sqrt(5)*phi.est*norm.delta
    k111 <- 25*phi.est^4*exp(-t.1)*(3-sqrt(5)*phi.est*delta.t2.t1[1]^2/norm.delta)*delta.t2.t1[1]/3
    k112 <- 25*phi.est^4*exp(-t.1)*(1-sqrt(5)*phi.est*delta.t2.t1[1]^2/norm.delta)*delta.t2.t1[2]/3
    k122 <- 25*phi.est^4*exp(-t.1)*(1-sqrt(5)*phi.est*delta.t2.t1[2]^2/norm.delta)*delta.t2.t1[1]/3
    k211 <- k112; k212 <- k122
    k222 <- 25*phi.est^4*exp(-t.1)*(3-sqrt(5)*phi.est*delta.t2.t1[2]^2/norm.delta)*delta.t2.t1[2]/3
    nabla3.K <- matrix(c(k111,k112,k122,k211,k212,k222), nrow=2,ncol=3,byrow=T)
    return(crossprod(t(crossprod(u.perp,nabla3.K)),c.u))
  }
}

K.22.m2 <- function(t,u, phi.est){
  u.perp <- rev(u); u.perp[2] <- -u.perp[2]
  c.u <- c(u.perp[1]^2,2*u.perp[1]*u.perp[2],u.perp[2]^2)
  delta.t2.t1 <- (t[2]-t[1])*u; norm.delta <- sqrt(sum((delta.t2.t1)^2))
  if(delta.t2.t1[1]==0 & delta.t2.t1[2]==0){
    nabla4.K <- 25*diag(c(phi.est^4,
                          phi.est^4/3,
                          phi.est^4))
    return(crossprod(t(crossprod(c.u,nabla4.K)),c.u))
  }else{
    t.1 <- sqrt(5)*phi.est*norm.delta
    k1111 <- 25*phi.est^4*exp(-t.1)*(3-6*sqrt(5)*phi.est*delta.t2.t1[1]^2/norm.delta+sqrt(5)*phi.est*delta.t2.t1[1]^4/norm.delta^3+5*phi.est^2*delta.t2.t1[1]^4/norm.delta^2)/3
    k1212 <- 25*phi.est^4*exp(-t.1)*((1-sqrt(5)*delta.t2.t1[1]^2/norm.delta)*(1-sqrt(5)*delta.t2.t1[2]^2/norm.delta)+sqrt(5)*phi.est*delta.t2.t1[1]^2*delta.t2.t1[2]^2/norm.delta^3)/3
    k2222 <- 25*phi.est^4*exp(-t.1)*(3-6*sqrt(5)*phi.est*delta.t2.t1[2]^2/norm.delta+sqrt(5)*phi.est*delta.t2.t1[2]^4/norm.delta^3-5*phi.est^2*delta.t2.t1[2]^4/norm.delta^2)/3
    k1112 <- 25*sqrt(5)*phi.est^5*exp(-t.1)*(delta.t2.t1[1]/norm.delta^3-(3-sqrt(5)*phi.est*delta.t2.t1[1]^2/norm.delta)/norm.delta)*delta.t2.t1[1]^2*delta.t2.t1[2]/3
    k1122 <- k1211 <- k2211 <- k1212
    k1222 <- 25*sqrt(5)*phi.est^5*exp(-t.1)*(delta.t2.t1[2]/norm.delta^3-(3-sqrt(5)*phi.est*delta.t2.t1[2]^2/norm.delta)/norm.delta)*delta.t2.t1[2]^2*delta.t2.t1[1]/3
    k2212 <- k1222
    nabla4.K <- matrix(c(k1111,k1112,k1122,
                         k1211,k1212,k1222,
                         k2211,k2212,k2222), nrow=3,ncol=3,byrow=T)
    return(crossprod(t(crossprod(c.u,nabla4.K)),c.u))
  }
}

# # True Value
# TrueGamma.t.1.m2 <- function(t,u,s0){
#   u.perp <- rev(u); u.perp[2] <- -u.perp[2]
#   t.pt <- s0+t*u
#   30*pi*cos(3*pi*(t.pt[1]))*u.perp[1]-30*pi*sin(3*pi*(t.pt[2]))*u.perp[2]
# }
# 
# TrueGamma.t.2.m2 <- function(t,u,s0){
#   u.perp <- rev(u); u.perp[2] <- -u.perp[2]
#   t.pt <- s0+t*u
# 
#   h11 <- -90*pi^2*sin(3*pi*(t.pt[1]))
#   h22 <- -90*pi^2*cos(3*pi*(t.pt[2]))
#   h12 <- 0
#   h11*u.perp[1]^2+2*h12*u.perp[1]*u.perp[2]+h22*u.perp[2]^2
# }
# 
# true_val <- matrix(NA,nr=nrow(curve)-1,nc=2)
# for(i in 1:(nrow(curve)-1)){
#   s0 <- curve[i,]
#   u <- uval[i,]
#   grad.li <- sapply(rule.uv*tval[i], function(x) TrueGamma.t.1.m2(x,u=u,s0=s0))
#   curv.li <- sapply(rule.uv*tval[i], function(x) TrueGamma.t.2.m2(x,u=u,s0=s0))
# 
#   true_val[i,]<- c(sum(grad.li[-length(rule.uv)]*(tval[i]/length(rule.uv))),
#                    sum(curv.li[-length(rule.uv)]*(tval[i]/length(rule.uv))))
# }

# colSums(true_val)/sum(tval)
# 
# # Checks
# # sapply(1:nrow(true_val), function(x){
# #   ifelse(true_val[x,1]<w1.est[x,3] & true_val[x,1]>w1.est[x,2],1,0)+ifelse(true_val[x,2]<w2.est[x,3] & true_val[x,2]>w2.est[x,2],1,0)
# # })
# 
# mat <- matrix(c(1:3), nr=1,nc=3, byrow=T)
# layout(mat,
#        widths = c(5,5,5),
#        heights = c(3,3,3))
# 
# sp_plot(11,"Spectral",cbind(coords,y), contour.plot = T,raster.surf = F, legend = F)
# # col.grad <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))
# # col.segment <- col.grad(nrow(w1.est))[order(w1.est[,1])]
# # for(i in 1:(nrow(curve)-1)){
# #   if(w1.est$sig[i]%in%c(1,-1)) lines(rbind(curve[i,],curve[(i+1),]), col=col.segment[i],lwd=4)
# #   else lines(rbind(curve[i,],curve[(i+1),]))
# # }
# # sp_plot(11,"Spectral",cbind(coords,y), contour.plot = T,raster.surf = F, legend = F)
# # col.segment <- col.grad(nrow(true_val))[order(true_val[,1])]
# # for(i in 1:(nrow(curve)-1)){
# #   lines(rbind(curve[i,],curve[(i+1),]), col=col.segment[i],lwd=4)
# # }
# col.segment <- sapply(w1.est$sig, function(x){
#   if(x==0) "white"
#   else if(x==1) "green"
#   else if(x==-1) "cyan"
# })
# for(i in 1:(nrow(curve)-1)) lines(rbind(curve[i,],curve[(i+1),]), col=col.segment[i],lwd=4)
# 
# sp_plot(11,"Spectral",cbind(coords,y), contour.plot = T,raster.surf = F, legend = F)
# col.segment <- sapply(w2.est$sig, function(x){
#   if(x==0) "white"
#   else if(x==1) "green"
#   else if(x==-1) "cyan"
# })
# for(i in 1:(nrow(curve)-1)) lines(rbind(curve[i,],curve[(i+1),]), col=col.segment[i],lwd=4)
