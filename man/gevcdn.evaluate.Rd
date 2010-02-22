\name{gevcdn.evaluate}
\alias{gevcdn.evaluate}

\title{
Evaluate parameters from trained GEV CDN model
}
\description{
Evaluate a trained GEV CDN model, resulting in a matrix with columns
corresponding to the location, scale, and shape parameters of the GEV
distribution.
}
\usage{
gevcdn.evaluate(x, weights)
}

\arguments{
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
}
  \item{weights}{
      list containing GEV CDN input-hidden and hidden-output layer weight matrices from \code{\link{gevcdn.fit}} or \code{\link{gevcdn.bag}}.
}
}

\value{
  a matrix with number of rows equal to that of \code{x} and columns
  corresponding to the GEV location, scale, and shape parameters.
}

\references{
Cannon, A.J., 2010. A flexible nonlinear modelling framework for
nonstationary generalized extreme value analysis in hydroclimatology.
Hydrological Processes, 24: 673-685. DOI: 10.1002/hyp.7506
}

\seealso{
    \code{\link{gevcdn.fit}}, \code{\link{gevcdn.bag}}, \code{\link{dgev}}
}
\examples{
# Generate synthetic data

x <- as.matrix(1:50)
y <- as.matrix(rgev(length(x), location = 0, scale = 1, shape = 0.2))

## Fit stationary model

weights <- gevcdn.fit(x = x, y = y, Th = gevcdn.identity,
                      fixed = c("location", "scale", "shape"))

## Evaluate GEV parameters

parms <- gevcdn.evaluate(x, weights)
cat("GEV parameters", parms[1,], "\n")
}
