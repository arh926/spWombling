% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_gradient.R
\name{spatial_gradient}
\alias{spatial_gradient}
\title{Gradient and Curvature Assessment}
\usage{
spatial_gradient(
  coords = NULL,
  model = NULL,
  cov.type = c("matern1", "matern2", "gaussian", ),
  grid.points = NULL,
  nbatch = 200,
  return.mcmc = TRUE
)
}
\arguments{
\item{coords}{coordinates for observed process (order \eqn{L} x \eqn{2})}

\item{model}{the posterior samples from the MCMC fit}

\item{cov.type}{covariance type (three available choices: Gaussian,  Mat\'ern(\eqn{\nu=3/2})), Mat\'ern(\eqn{\nu=5/2})}

\item{grid.points}{coordinates for grid over observed process (order \eqn{n_G} x \eqn{2})}

\item{nbatch}{number of batches}

\item{return.mcmc}{if true returns MCMC samples of gradients}
}
\description{
performs gradient and curvature assessment on the estimated surface from a hierarchical Bayesian spatial model. Parameters not listed are optional
}
