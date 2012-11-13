\name{gevcdn.cost}
\alias{gevcdn.cost}
\title{
Cost function for GEV CDN model fitting
}
\description{
The generalized maximum likelihood (GML) cost function used for GEV CDN model
fitting (Martins and Stedinger, 2000). Calculates the negative of the logarithm
of the GML, which includes a shifted beta distribution prior for the GEV shape
parameter. A normal distribution prior can also be set for the magnitude of the
input-hidden layer weights, thus leading to weight penalty regularization.
}

\usage{
gevcdn.cost(weights, x, y, n.hidden, Th, fixed, scale.min, beta.p,
            beta.q, sd.norm)
}

\arguments{
  \item{weights}{
      weight vector of length returned by \code{\link{gevcdn.initialize}}.
}
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
}
  \item{y}{
      column matrix of target values with number of rows equal to the number of samples.
}
  \item{n.hidden}{
      number of hidden nodes in the GEV CDN model.
}
  \item{Th}{
      hidden layer transfer function; defaults to \code{\link{gevcdn.logistic}}.
}
  \item{fixed}{
      vector indicating GEV parameters to be held constant; elements chosen from \code{c("location", "scale", "shape")}
}
  \item{scale.min}{
      minimum allowable value for the GEV scale parameter.
}
  \item{beta.p}{
      \code{shape1} parameter for shifted beta distribution prior for GEV shape parameter.
}
  \item{beta.q}{
      \code{shape2} parameter for shifted beta distribution prior for GEV shape parameter.
}
  \item{sd.norm}{
      \code{sd} parameter for normal distribution prior for the magnitude of input-hidden layer weights; equivalent to weight penalty regularization.
}
}

\references{

Martins, E.S. and J.R. Stedinger, 2000. Generalized maximum-likelihood
generalized extreme-value quantile estimators for hydrologic data. Water
Resources Research, 36: 737-744. DOI: 10.1029/1999WR900330

}

\seealso{
    \code{\link{gevcdn.fit}}, \code{\link{gevcdn.bag}}, \code{\link{dgev}},
    \code{\link{optim}}
}
