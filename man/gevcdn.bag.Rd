\name{gevcdn.bag}
\alias{gevcdn.bag}

\title{
Fit an ensemble of GEV CDN models via bagging
}
\description{
Used to fit an ensemble of GEV CDN models using bootstrap aggregation
(bagging) and, optionally, early stopping.
}
\details{
Bootstrap aggregation (bagging) (Breiman, 1996) is used as an alternative to
\code{\link{gevcdn.fit}} for fitting an ensemble of nonstationary GEV CDN models
(Cannon, 2010). Each ensemble member is trained on bootstrapped \code{x} and
\code{y} sample pairs. As an added check on overfitting, early stopping,
whereby training is stopped prior to convergence of the optimization algorithm,
can be turned on. In this case, the "best" iteration for stopping optimization
is chosen based on model performance on out-of-bag samples. Additional details on
the bagging/early stopping procedure can be found in Cannon and Whitfield
(2001, 2002). Unlike \code{\link{gevcdn.fit}}, where models of increasing complexity
are fitted and the one that satisfies some model selection criterion is chosen for
final use, model selection in \code{gevcdn.bag} is implicit. Ensemble averaging and,
optionally, early stopping are used to limit model complexity. One need only
set \code{n.hidden} to avoid underfitting.

Note: values of \code{x} and \code{y} need not be standardized or rescaled by
the user. All variables are automatically scaled to range between 0 and 1 prior
to fitting and parameters are automatically rescaled by
\code{\link{gevcdn.evaluate}}.
}

\usage{
gevcdn.bag(x, y, iter.max = 1000, iter.step = 10,
           n.bootstrap = 30, n.hidden = 3, Th = gevcdn.logistic,
           fixed = NULL, init.range = c(-0.25, 0.25),
           scale.min = .Machine$double.eps, beta.p = 3.3,
           beta.q = 2, sd.norm = Inf,
           method = c("BFGS", "Nelder-Mead"), max.fails = 100,
           silent = TRUE, ...)
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
  \item{iter.step}{
      number of iterations between which the value of the cost function is calculated on the out-of-bootstrap samples; used during stopped training.
}
  \item{n.bootstrap}{
      number of ensemble members trained during bootstrap aggregation (bagging).
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
  \item{method}{
      optimization method for \code{\link{optim}} function; must be chosen from \code{c("BFGS", "Nelder-Mead")}.
}
  \item{max.fails}{
      maximum number of repeated exceptions allowed during optimization.
}
  \item{silent}{
      logical value indicating whether or not cost function information should be reported every \code{iter.step} iterations.
}
  \item{\dots}{
      additional arguments passed to the \code{control} list of \code{\link{optim}} function.
}
}

\value{
  a list of length \code{n.bootstrap} with elements consisting of
    \item{W1}{input-hidden layer weights}
    \item{W2}{hidden-output layer weights}
  Attributes indicating the minimum/maximum values of \code{x} and \code{y}; the
  values of \code{Th}, \code{fixed}, \code{scale.min}; a logical value
  \code{stopped.training} indicating whether or not early stopping was used;
   and the minimum value of the cost function on the out-of-bag samples
  \code{cost.valid} are also returned.
}

\references{

Breiman, L., 1996. Bagging predictors. Machine Learning, 24(2): 123-140.
DOI: 10.1007/BF00058655

Cannon, A.J., 2010. A flexible nonlinear modelling framework for
nonstationary generalized extreme value analysis in hydroclimatology.
Hydrological Processes, 24: 673-685. DOI: 10.1002/hyp.7506

Cannon, A.J. and P.H. Whitfield, 2001. Modeling transient pH depressions in
coastal streams of British Columbia using neural networks, 37(1): 73-89.
DOI: 10.1111/j.1752-1688.2001.tb05476.x

Cannon, A.J. and P.H. Whitfield, 2002. Downscaling recent streamflow conditions
in British Columbia, Canada using ensemble neural network models. Journal of
Hydrology, 259: 136-151. DOI: 10.1016/S0022-1694(01)00581-9

}

\seealso{
    \code{\link{gevcdn.cost}}, \code{\link{gevcdn.evaluate}},
    \code{\link{gevcdn.fit}}, \code{\link{gevcdn.bootstrap}},
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

\dontrun{
## Fit ensemble of models with early stopping turned on

weights.on <- gevcdn.bag(x = x, y = y, iter.max = 100,
                         iter.step = 10, n.bootstrap = 10,
                         n.hidden = 2)

parms.on <- lapply(weights.on, gevcdn.evaluate, x = x)

## 10th, 50th, and 90th percentiles

q.10.on <- q.50.on <- q.90.on <- matrix(NA, ncol=length(parms.on),
                                        nrow=nrow(x))
for(i in seq_along(parms.on)){
    q.10.on[,i] <- qgev(p = 0.1,
                        location = parms.on[[i]][,"location"],
                        scale = parms.on[[i]][,"scale"],
                        shape = parms.on[[i]][,"shape"])
    q.50.on[,i] <- qgev(p = 0.5,
                        location = parms.on[[i]][,"location"],
                        scale = parms.on[[i]][,"scale"],
                        shape = parms.on[[i]][,"shape"])
    q.90.on[,i] <- qgev(p = 0.9,
                        location = parms.on[[i]][,"location"],
                        scale = parms.on[[i]][,"scale"],
                        shape = parms.on[[i]][,"shape"])
}

## Plot data and quantiles

matplot(cbind(y, q, rowMeans(q.10.on), rowMeans(q.50.on),
        rowMeans(q.90.on)), type = c("b", rep("l", 6)),
        lty = c(1, rep(c(1, 2, 1), 2)),
        lwd = c(1, rep(c(3, 2, 3), 2)),
        col = c("red", rep("orange", 3), rep("blue", 3)),
        pch = 19, xlab = "x", ylab = "y",
        main = "gevcdn.bag (early stopping on)")
}
}
