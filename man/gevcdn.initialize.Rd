\name{gevcdn.initialize}
\alias{gevcdn.initialize}
\title{
Initialize GEV CDN weight vector
}
\description{
Random initialization of the weight vector used during fitting of the GEV CDN model.
}
\usage{
gevcdn.initialize(x, n.hidden, init.range)
}

\arguments{
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
}
  \item{n.hidden}{
      number of hidden nodes in the GEV CDN model.
}
  \item{init.range}{
      range for random weights on [\code{min(init.range)}, \code{max(init.range)}]
}
}

\seealso{
    \code{\link{gevcdn.reshape}}
}
