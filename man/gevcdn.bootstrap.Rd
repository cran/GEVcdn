\name{gevcdn.bootstrap}
\alias{gevcdn.bootstrap}

\title{
Bootstrapped confidence intervals for GEV CDN parameters and quantiles
}
\description{
Used to assist in the calculation of bootstrapped confidence intervals for
GEV location, scale, and shape parameters, as well as for specified quantiles.
Residual and parametric bootstrap estimates are supported.
}
\details{
Note: the \code{\link[boot:boot]{boot}} package provides a more
comprehensive suit of methods for estimating bootstrap-based confidence
intervals.
}
\usage{
gevcdn.bootstrap(n.bootstrap, x, y, iter.max = 1000, n.hidden = 2,
                 Th = gevcdn.logistic, fixed = NULL,
                 init.range = c(-0.25, 0.25),
                 scale.min = .Machine$double.eps,
                 beta.p = 3.3, beta.q = 2, sd.norm = Inf,
                 n.trials = 5, method = c("BFGS", "Nelder-Mead"),
                 boot.method = c("residual", "parametric"),
                 init.from.prev = TRUE, max.fails = 100,
                 probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9,
                 0.95, 0.99), ...)
}

\arguments{
  \item{n.bootstrap}{
      number of bootstrap trials used to calculate confidence intervals.
}
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
}
  \item{y}{
      column matrix of target values with number of rows equal to the number of samples.
}
  \item{iter.max}{
      maximum number of iterations of optimization algorithm.
}
  \item{n.hidden}{
      number of hidden nodes in each GEV CDN ensemble member.
}
  \item{Th}{
      hidden layer transfer function; defaults to \code{\link{gevcdn.logistic}}.
}
  \item{fixed}{
      vector indicating GEV parameters to be held constant; elements chosen from \code{c("location", "scale", "shape")}
}
  \item{init.range}{
      range for random weights on [\code{min(init.range)}, \code{max(init.range)}]
}
  \item{scale.min}{
      minimum allowable value for the GEV scale parameter.
}
  \item{beta.p}{
      \code{shape1} parameter for shifted beta distribution prior for GEV shape parameter
}
  \item{beta.q}{
      \code{shape2} parameter for shifted beta distribution prior for GEV shape parameter
}
  \item{sd.norm}{
      \code{sd} parameter for normal distribution prior for the magnitude of input-hidden layer weights; equivalent to weight penalty regularization.
}
  \item{n.trials}{
      number of repeated trials used to avoid shallow local minima during optimization
}
  \item{method}{
      optimization method for \code{\link{optim}} function; must be chosen from \code{c("BFGS", "Nelder-Mead")}.
}
  \item{boot.method}{
      bootstrap method; must be chosen from \code{c("residual", "parametric")}.
}
  \item{init.from.prev}{
      logical value indicating whether or not optimization runs should be initialized from the final weights of the previous run.
}
  \item{max.fails}{
      maximum number of repeated exceptions allowed during optimization
}
  \item{probs}{
      vector of cumulative probabilities for which quantile confidence intervals are desired.
}
  \item{\dots}{
      additional arguments passed to the \code{control} list of \code{\link{optim}} function
}
}

\value{
  a list consisting of
    \item{weights.bootstrap}{a list of length \code{n.bootstrap} consisting of weights following the format returned by \code{\link{gevcdn.fit}}}
    \item{parms.bootstrap}{a list of length \code{n.bootstrap} consisting of GEV parameters following the format returned by \code{\link{gevcdn.evaluate}}}
    \item{location.bootstrap}{a matrix of GEV location parameters with number of rows equal to that of \code{x} and number of columns equal to \code{n.bootstrap}}
    \item{scale.bootstrap}{a matrix of GEV scale parameters with number of rows equal to that of \code{x} and number of columns equal to \code{n.bootstrap}}
    \item{shape.bootstrap}{a matrix of GEV shape parameters with number of rows equal to that of \code{x} and number of columns equal to \code{n.bootstrap}}
    \item{quantiles.bootstrap}{a list of length \code{n.bootstrap} with each consisting of a matrix with number of rows equal to that of \code{x} and columns corresponding to \code{probs}}
}

\references{

Cannon, A.J., 2010. A flexible nonlinear modelling framework for
nonstationary generalized extreme value analysis in hydroclimatology.
Hydrological Processes, 24: 673-685. DOI: 10.1002/hyp.7506

}

\seealso{
    \code{\link{gevcdn.cost}}, \code{\link{gevcdn.evaluate}},
    \code{\link{gevcdn.fit}}, \code{\link{gevcdn.bag}},
    \code{\link{dgev}}, \code{\link{optim}}
}

\examples{
## Generate synthetic data

x <- as.matrix(seq(0.1, 1, length = 50))

loc <- x^2
scl <- x/2
shp <- seq(-0.1, 0.3, length = length(x))

set.seed(100)
y <- as.matrix(rgev(length(x), location = loc, scale = scl,
               shape = shp))

# Not run:
## Fit 30 bootstrapped models
#
# CI <- gevcdn.bootstrap(n.bootstrap = 30, x = x, y = y,
#                       iter.max = 100, n.hidden = 2,
#                       Th = gevcdn.logistic, n.trials = 1,
#                       boot.method = "residual",
#                       probs = c(0.1, 0.5, 0.9))
#
## Plot data and percentile confidence intervals for GEV parameters
#
# par(mfrow = c(2, 2))
# matplot(x, y, type = "b", pch = 19, col = "red", xlab = "x",
#        ylab = "y", main = "gevcdn.bootstrap")
# matplot(x, cbind(loc, t(apply(CI$location.bootstrap, 1, quantile,
#         p = c(0.025, 0.975)))), type = c("l", "b", "b"), pch = 20,
#         lwd = 3, col = c("black", rep("green", 2)), xlab = "x",
#         ylab = "location", main = "location CI")
# matplot(x, cbind(scl, t(apply(CI$scale.bootstrap, 1, quantile,
#         p = c(0.025, 0.975)))), type = c("l", "b", "b"), pch = 20,
#         lwd = 3, col = c("black", rep("orange", 2)), xlab = "x",
#         ylab = "scale", main = "scale CI")
# matplot(x, cbind(shp, t(apply(CI$shape.bootstrap, 1, quantile,
#         p = c(0.025, 0.975)))), type = c("l", "b", "b"), pch = 20,
#         lwd = 3, col = c("black", rep("brown", 2)), xlab = "x",
#         ylab = "shape", main = "shape CI")
}
