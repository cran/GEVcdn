\name{GEVcdn-package}
\alias{GEVcdn-package}
\alias{GEVcdn}
\docType{package}
\title{
GEV Conditional Density Estimation Network
}
\description{
Parameters in a Generalized Extreme Value (GEV) distribution are specified as a
function of covariates using a conditional density estimation network (CDN),
which is a probabilistic variant of the multilayer perceptron neural network.
If the covariate is time or is dependent on time, then the GEV CDN model can be
used to perform nonlinear, nonstationary GEV analysis of hydrological or
climatological time series. Owing to the flexibility of the neural network
architecture, the model is capable of representing a wide range of nonstationary
relationships, including those involving interactions between covariates.
Model parameters are estimated by generalized maximum likelihood, an approach
that is tailored to the estimation of GEV parameters from geophysical time
series.
}
\details{
Procedures for fitting GEV CDN models are provided by the functions
\code{\link{gevcdn.fit}} and \code{\link{gevcdn.bag}}. Once a model has been
developed, \code{\link{gevcdn.evaluate}} is used to evaluate the GEV
distribution parameters as a function of covariates. Confidence intervals for
GEV parameters and specified quantiles can be estimated using
\code{\link{gevcdn.bootstrap}}. All other functions are used internally and
should not normally need to be called directly by the user.

Note: the GEV distribution functions \code{dgev}, \code{pgev}, 
\code{qgev}, and \code{rgev} are from the archived 
\code{VGAM 1.0-6} package. The convention for the sign of the shape parameter
is opposite to that used in hydrology and thus differs from Cannon (2010).
}

\references{
Cannon, A.J., 2010. A flexible nonlinear modelling framework for
nonstationary generalized extreme value analysis in hydroclimatology.
Hydrological Processes, 24: 673-685. DOI: 10.1002/hyp.7506

Cannon, A.J., 2011. GEVcdn: an R package for nonstationary extreme
value analysis by generalized extreme value conditional density
estimation network, Computers & Geosciences, 37: 1532-1533.
DOI: 10.1016/j.cageo.2011.03.005

Yee, T.W., 2018. VGAM: Vector Generalized Linear and Additive Models. 
R package version 1.0-6.
}

\keyword{ package }
