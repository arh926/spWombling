# spWombling


Illustration of curvilinear Bayesian Wombling on spatial data.

<p align="center">
  <img width="600" height="550" src="https://user-images.githubusercontent.com/73150479/143766893-6adbde31-34d4-4f68-aea2-71c69b563455.png">
<p>
  
Reference to the paper titled, Curvature Processes: Directional Concavity in Gaussian random fields. (Aritra Halder, Sudipto Banerjee, Dipak K. Dey)


Code for performing curvilinear Bayesian wombling on data

## Contents

1. Load the data and separate (a) co-ordinates (b) response (c) covariates
2. Fit a spatial Bayesian hierarchical model to the data
3. Perform gradient, curvature and posterior surface analysis
4. Locate or Annotate curves of interest
5. Perform rectilinear wombling at line segments

Demonstration with synthetic data.

### Load the data and separate (a) co-ordinates (b) response (c) covariates
```
if(!require(devtools)) install.packages("devtools")
devtools::install_github('arh926/spWombling')
require(spWombling)
# loading additional packages
require(coda)
require(sp)
require(Matrix)

N <- 100
tau <- 1
# synthetic location
coords <- matrix(runif(2*N),nr=N,nc=2) 
# synthetic pattern 
y <- rnorm(N,10*(sin(3*pi*coords[,1])+cos(3*pi*coords[,2])),tau)
# y <- rnorm(N,10*(sin(3*pi*coords[,1])*cos(3*pi*coords[,2])),tau)
y=y # response
X=matrix(1,nr=N) # intercept
```
![pattern](https://user-images.githubusercontent.com/73150479/143763604-311b0763-7148-45c3-9979-d732bedc01bf.jpg)

### Fit a spatial Bayesian hierarchical model to the data

```
niter <- 1e4
nburn <- niter/2
report <- 1e2
mc_sp <- hlmBayes_sp(coords = coords,
                     y = y,
                     X = X,
                     niter = niter,
                     nburn = nburn,
                     report = report,
                     cov.type = "matern2",
                     steps_init = 1,
                     verbose = T)
# inference from model
model_summary <- hlm_summary(chain = mc_sp,
                             nburn = nburn,
                             niter = niter,
                             thin = 1)
coef <- model_summary$summary.pars; round(coef,4)

z <- model_summary$summary.latent
```

### Perform gradient, curvature and posterior surface analysis

```
grid.points <- as.matrix(expand.grid(seq(0,1,by=0.05)[-c(1,21)],
                                     seq(0,1,by=0.05)[-c(1,21)]),nc=2)
samples <- (nburn+1):niter
gradient_est <- spatial_gradient(coords=coords,
                                 grid.points = grid.points,
                                 samples = samples,
                                 chain = mc_sp,
                                 cov.type = "matern2",
                                 niter = niter,
                                 nburn = nburn,
                                 nbatch = 100,
                                 return.mcmc = T)
# Gradient and Curvature Assessment
grad.s1.hpd <- data.frame(gradient_est$grad1.est)
grad.s1.hpd$signif <- apply(grad.s1.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  

grad.s2.hpd <- data.frame(gradient_est$grad2.est)
grad.s2.hpd$signif <- apply(grad.s2.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  

grad.s11.hpd <- data.frame(gradient_est$grad11.est)
grad.s11.hpd$signif <- apply(grad.s11.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  
grad.s12.hpd <- data.frame(gradient_est$grad12.est)
grad.s12.hpd$signif <- apply(grad.s12.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  
grad.s22.hpd <- data.frame(gradient_est$grad22.est)
grad.s22.hpd$signif <- apply(grad.s22.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  

# Posterior Surface Analysis
det_est <- eigen_est_1 <- eigen_est_2 <- laplace_est <- div_est <- matrix(NA,nr=niter,nc=nrow(grid.points))
for(i in 1:nburn){
  for(j in 1:nrow(grid.points)){
    hess.mat <- matrix(c(gradient_est$grad.s11.mcmc[i,j],
                         gradient_est$grad.s12.mcmc[i,j],
                         gradient_est$grad.s12.mcmc[i,j],
                         gradient_est$grad.s22.mcmc[i,j]),nr=2,nc=2,byrow = T)
    det_est[i,j] <- det(hess.mat)
    eigen.hess.mat <- eigen(hess.mat)$values
    eigen_est_1[i,j] <- eigen.hess.mat[1]; eigen_est_2[i,j] <- eigen.hess.mat[2]
    laplace_est[i,j] <- sum(diag(hess.mat))
    div_est[i,j] <- gradient_est$grad.s1.mcmc[i,j]+gradient_est$grad.s2.mcmc[i,j]
  }
}
slist <- split(1:nburn, ceiling(seq_along(1:nburn)/(length(1:nburn)/100)))
# Determinant Surface
det_est_val <- cbind.data.frame(apply(do.call(rbind,lapply(slist, function(x) apply(det_est[x,],2,median))),2,median),
                                HPDinterval(as.mcmc(do.call(rbind,lapply(slist, function(x) apply(det_est[x,],2,median))))))
det_est_val$sig <- apply(det_est_val,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
# Eigen Value Surface
eigen_est_val1 <- cbind.data.frame(apply(do.call(rbind,lapply(slist, function(x) apply(eigen_est_1[x,],2,median))),2,median),
                                   HPDinterval(as.mcmc(do.call(rbind,lapply(slist, function(x) apply(eigen_est_1[x,],2,median))))))
eigen_est_val1$sig <- apply(eigen_est_val1,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
eigen_est_val2 <- cbind.data.frame(apply(do.call(rbind,lapply(slist, function(x) apply(eigen_est_2[x,],2,median))),2,median),
                                   HPDinterval(as.mcmc(do.call(rbind,lapply(slist, function(x) apply(eigen_est_2[x,],2,median))))))
eigen_est_val2$sig <- apply(eigen_est_val2,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
# Laplace Operator Surface
laplace_est_val <- cbind.data.frame(apply(do.call(rbind,lapply(slist, function(x) apply(laplace_est[x,],2,median))),2,median),
                                    HPDinterval(as.mcmc(do.call(rbind,lapply(slist, function(x) apply(laplace_est[x,],2,median))))))
laplace_est_val$sig <- apply(laplace_est_val,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
# Divergence Operator
div_est_val <- cbind.data.frame(apply(do.call(rbind,lapply(slist, function(x) apply(div_est[x,],2,median))),2,median),
                                HPDinterval(as.mcmc(do.call(rbind,lapply(slist, function(x) apply(div_est[x,],2,median))))))
div_est_val$sig <- apply(div_est_val,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
```
#### True and Estimated Gradients
![pat1-grad](https://user-images.githubusercontent.com/73150479/143763736-66adb1db-e453-42ea-8742-defdef444fad.jpg)
#### True and Estimated Curvature
![pat1-curv](https://user-images.githubusercontent.com/73150479/143764182-21d77245-2005-4e9a-8586-3ce19076fd7e.jpg)


### Locate or Annotate curves of interest

```
#####################################
# locate curves:: level sets/contours
#####################################
mat <- matrix(c(1:2), nr=1,nc=2, byrow=T)
layout(mat,
       widths = c(5,2),
       heights = c(3,3))

rastr.obj <- sp_plot(11,"Spectral",cbind(coords,y), contour.plot = T,raster.surf = T, legend = F)
x <- raster::rasterToContour(rastr.obj,nlevels=20)
# Test for wombling
x.levels <- as.numeric(as.character(x$level))
level.choice <- c(x.levels[2],x.levels[length(x.levels)-1])
subset.points.1 <- subset(x,level==level.choice[1])
subset.points.2 <- subset(x,level==level.choice[2])
x <- raster::rasterToContour(rastr.obj,nlevels=10)
subset.points.3 <- subset(x,level==15)
# curve-1
lines(subset.points.1@lines[[1]]@Lines[[1]]@coords,lwd=2.5)
# curve-2
lines(subset.points.2@lines[[1]]@Lines[[3]]@coords,lwd=2.5)
# curve 3
lines(subset.points.3@lines[[1]]@Lines[[3]]@coords,lwd=2.5)

#####################################
# locate curves:: annotate curves
#####################################

pts.hull <- locator(20)
subset.points.4 <- bezierCurve(pts.hull$x,pts.hull$y,100)
# curve 4
lines(subset.points.4, type="l",lwd=2.5)
points(pts.hull, cex=0.5, col="white")
points(pts.hull, pch="+", cex=0.5)

legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(y),max(y),l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)
```
### Perform rectilinear wombling at line segments
```
curve <- subset.points.1@lines[[1]]@Lines[[1]]@coords
womb.measure <- bayes_cwomb(coords = coords,
                            chain = mc_sp,
                            cov.type = "matern2",
                            niter = niter,
                            nburn = nburn,
                            nbatch = 100,
                            curve = curve,
                            type = "rectilinear")
str(womb.measure)
womb.measure$womb.measure.inf
```
#### Curvilinear Gradient
![post-grad-1](https://user-images.githubusercontent.com/73150479/143763716-3931ea11-d8f6-4819-9828-57140a04f8bc.jpg)
#### Curvilinear Curvature
![post-grad-2](https://user-images.githubusercontent.com/73150479/143763720-e4ce5560-db24-42a9-95e4-e3fc52d51719.jpg)

## Authors

| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Aritra Halder (maintainer)| gxk9jg@virginia.edu   | Research Assistant Professor, Biocomplexity Institute, UVA|                         
| Sudipto Banerjee | sudipto@ucla.edu   | Professor and Chair, Department of Biostatistics,  UCLA |
| Dipak K. Dey | dipak.dey@uconn.edu   | Board of Trustees Distinguished Professor, Department of Statistics,  UCONN |
<!--- --->

