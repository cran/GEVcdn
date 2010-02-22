\name{gevcdn.reshape}
\alias{gevcdn.reshape}
\title{
Reshape a GEV CDN weight vector
}
\description{
Reshapes a weight vector used during fitting of the GEV CDN model into input-hidden and hidden-output layer weight matrices.
}
\usage{
gevcdn.reshape(x, weights, n.hidden)
}

\arguments{
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
}
  \item{weights}{
      weight vector of length returned by \code{\link{gevcdn.initialize}}.
}
  \item{n.hidden}{
      number of hidden nodes in the GEV CDN model.
}
}

\seealso{
    \code{\link{gevcdn.initialize}}
}
