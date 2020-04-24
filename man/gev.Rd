\name{gev}
\alias{dgev}
\alias{pgev}
\alias{qgev}
\alias{rgev}
\alias{dgumbel}
\alias{pgumbel}
\alias{qgumbel}
\alias{rgumbel}
\alias{is.Numeric}
\title{
Generalized extreme value distribution
}
\description{
Generalized extreme value density, distribution, quantile, and
random variate functions.
}
\usage{
dgev(x, location = 0, scale = 1, shape = 0, log = FALSE,
    tolshape0 = sqrt(.Machine$double.eps))
pgev(q, location = 0, scale = 1, shape = 0, lower.tail = TRUE,
     log.p = FALSE)
qgev(p, location = 0, scale = 1, shape = 0, lower.tail = TRUE,
     log.p = FALSE)
rgev(n, location = 0, scale = 1, shape = 0)

dgumbel(x, location = 0, scale = 1, log = FALSE)
pgumbel(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
qgumbel(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
rgumbel(n, location = 0, scale = 1)

is.Numeric(x, length.arg = Inf, integer.valued = FALSE,
    positive = FALSE)
}
\arguments{
    \item{x, q}{vector of quantiles.}
    \item{p}{vector of probabilities.}
    \item{n}{number of observations.}
    \item{location}{the location parameter mu.}
    \item{scale}{the (positive) scale parameter sigma. Must consist of positive values.}
    \item{shape}{the shape parameter xi.}
    \item{log}{Logical.}
    \item{lower.tail, log.p}{Logical.}
    \item{tolshape0}{Positive numeric. Threshold/tolerance value for resting
          whether xi is zero. If the absolute value of the estimate of
          xi is less than this value then it will be assumed zero and a
          Gumbel distribution will be used.}
    \item{length.arg, integer.valued, positive}{Arguments of internally called is.Numeric function.}
}
