\name{gevcdn.fit}
\alias{gevcdn.fit}

\title{
Fit a GEV CDN model
}
\description{
Fit a GEV CDN model via nonlinear optimization of the generalized maximum
likelihood cost function.
}
\details{
Fit a nonstationary GEV CDN model (Cannon, 2010) by minimizing a cost function
based on the generalized maximum likelihood of Martins and Stedinger (2000).
The hidden layer transfer function \code{Th} should be set to
\code{\link{gevcdn.logistic}} for a nonlinear model and to
\code{\link{gevcdn.identity}} for a linear model. In the nonlinear case, the
number of hidden nodes \code{n.hidden} controls the overall complexity of the
model. GEV parameters can be held constant (i.e., stationary) via the
\code{fixed} argument. The form of the shifted beta distribution prior for the
GEV shape parameter is controlled by the \code{beta.p} and \code{beta.q}
arguments. Weight penalty regularization for the magnitude of the input-hidden
layer weights can be applied by setting \code{sd.norm} to a value less than
\code{Inf}.

Values of the Akaike information criterion (AIC), Akaike information criterion
with small sample size correction (AICc), and Bayesian information criterion
(BIC) are calculated to assist in model selection. It is possible for such
criteria to fail in the face of overfitting, for example with a nonlinear model
and \code{n.hidden} set too high, as the GEV distribution may converge on one or
more samples. This can usually be diagnosed by inspecting the scale parameter
for near zero values. In this case, one can apply a weight penalty
(via \code{sd.norm}), although this rules out the use of AIC/AICc/BIC for model
selection as the effective number of model parameters will no longer equal the
number of weights in the GEV CDN model. Alternatively, a lower threshold (as a
proportion of the range of \code{y}) for the scale parameter can be set
via \code{scale.min}. Finally, bootstrap aggregation is available via
\code{\link{gevcdn.bag}} as a second method for fitting GEV CDN models.

Note: values of \code{x} and \code{y} need not be standardized or rescaled
by the user. All variables are automatically scaled to range between 0 and 1
prior to fitting and parameters are automatically rescaled by
\code{\link{gevcdn.evaluate}}.
}

\usage{
gevcdn.fit(x, y, iter.max = 1000, n.hidden = 2,
           Th = gevcdn.logistic, fixed = NULL,
           init.range = c(-0.25, 0.25),
           scale.min = .Machine$double.eps, beta.p = 3.3,
           beta.q = 2, sd.norm = Inf, n.trials = 5,
           method = c("BFGS", "Nelder-Mead"), max.fails = 100, ...)
}

\arguments{
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
      number of hidden nodes in the GEV CDN model.
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
      \code{shape1} parameter for shifted beta distribution prior for GEV shape parameter.
}
  \item{beta.q}{
      \code{shape2} parameter for shifted beta distribution prior for GEV shape parameter.
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
  \item{max.fails}{
      maximum number of repeated exceptions allowed during optimization
}
  \item{\dots}{
      additional arguments passed to the \code{control} list of \code{\link{optim}} function
}
}

\value{
  a list consisting of
    \item{W1}{input-hidden layer weights}
    \item{W2}{hidden-output layer weights}
  Attributes indicating the minimum/maximum values of \code{x} and \code{y}; the
  values of \code{Th}, \code{fixed}, \code{scale.min}; the negative of the
  logarithm of the generalized maximum likelihood \code{GML}, the negative of
  the logarithm of the likelihood \code{NLL}, the value of the penalty term
  \code{penalty}, the Bayesian information criterion \code{BIC}, the Akaike
  information criterion with (\code{AICc}) and without (\code{AIC}) small sample
  size correction; and the number of model parameters \code{k} are also
  returned.
}

\references{

Cannon, A.J., 2010. A flexible nonlinear modelling framework for
nonstationary generalized extreme value analysis in hydroclimatology.
Hydrological Processes, 24: 673-685. DOI: 10.1002/hyp.7506

Martins, E.S. and J.R. Stedinger, 2000. Generalized maximum-likelihood
generalized extreme-value quantile estimators for hydrologic data. Water
Resources Research, 36:737-744. DOI: 10.1029/1999WR900330

}

\seealso{
    \code{\link{gevcdn.cost}}, \code{\link{gevcdn.evaluate}},
    \code{\link{gevcdn.bag}}, \code{\link{gevcdn.bootstrap}},
    \code{\link{dgev}}, \code{\link{optim}}
}
\examples{
## Generate synthetic data, quantiles

x <- as.matrix(seq(0.1, 1, length = 50))

loc <- x^2
scl <- x/2
shp <- seq(-0.1, 0.3, length = length(x))

set.seed(100)
y <- as.matrix(rgev(length(x), location = loc, scale = scl,
               shape = shp))
q <- sapply(c(0.1, 0.5, 0.9), qgev, location = loc, scale = scl,
            shape = shp)

## Define a hierarchy of models of increasing complexity

models <- list()
# Stationary model
models[[1]] <- list(Th = gevcdn.identity,
                    fixed = c("location", "scale", "shape"))
# Linear model
models[[2]] <- list(Th = gevcdn.identity)
# Nonlinear, 1 hidden node
models[[3]] <- list(n.hidden = 1, Th = gevcdn.logistic)
# Nonlinear, 2 hidden nodes
models[[4]] <- list(n.hidden = 2, Th = gevcdn.logistic)

## Fit models

weights.models <- list()
for(i in seq_along(models)){
    weights.models[[i]] <- gevcdn.fit(x = x, y = y, n.trials = 1,
                                      n.hidden = models[[i]]$n.hidden,
                                      Th = models[[i]]$Th,
                                      fixed = models[[i]]$fixed)
}

## Select model with minimum AICc

models.AICc <- sapply(weights.models, attr, which = "AICc")
weights.best <- weights.models[[which.min(models.AICc)]]
parms.best <- gevcdn.evaluate(x, weights.best)

## 10th, 50th, and 90th percentiles

q.best <- sapply(c(0.1, 0.5, 0.9), qgev,
                 location = parms.best[,"location"],
                 scale = parms.best[,"scale"],
                 shape = parms.best[,"shape"])

## Plot data and quantiles

matplot(x, cbind(y, q, q.best), type = c("b", rep("l", 6)),
        lty = c(1, rep(c(1, 2, 1), 2)),
        lwd = c(1, rep(c(3, 2, 3), 2)),
        col = c("red", rep("orange", 3), rep("blue", 3)),
        pch = 19, xlab = "x", ylab = "y", main = "gevcdn.fit")
}
