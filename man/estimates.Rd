\name{estimates}
\alias{estimates}
\alias{estimates.hlme}
\alias{estimates.lcmm}
\alias{estimates.Jointlcmm}
\alias{estimates.multlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum likelihood estimates
}
\description{
This function provides the vector of maximum likelihood estimates of a model estimated with \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}.
}
\usage{
\method{estimates}{hlme}(x,cholesky=TRUE)
\method{estimates}{lcmm}(x,cholesky=TRUE)
\method{estimates}{Jointlcmm}(x,cholesky=TRUE)
\method{estimates}{multlcmm}(x,cholesky=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}
}
  \item{cholesky}{
optional logical indicating if the parameters of variance-covariance of the random effets should be displayed instead of their cholesky transformations used in the estimation process.  
}
}

\value{
a vector with all estimates of the model.
}

\author{
Cecile Proust-Lima, Viviane Philipps
}
\seealso{
\code{\link{VarCov}}, \code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{multlcmm}}, \code{\link{Jointlcmm}}
}

