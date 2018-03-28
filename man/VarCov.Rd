\name{VarCov}
\alias{VarCov}
\alias{VarCov.hlme}
\alias{VarCov.lcmm}
\alias{VarCov.Jointlcmm}
\alias{VarCov.multlcmm}

\title{
Variance-covariance of the estimates
}
\description{
This function provides the variance-covariance matrix of the
estimates. vcov is an alias for it.
}
\usage{
\method{VarCov}{hlme}(x)
\method{VarCov}{lcmm}(x)
\method{VarCov}{Jointlcmm}(x)
\method{VarCov}{multlcmm}(x)
}

\arguments{
  \item{x}{
an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}
}
}

\value{
a matrix containing the variance-covariance of the estimates. For the parameters of the matrix of variance-covariance of the random effects, the Cholesky transformed parameters are considered so that VarCov provides the covariance matrix of function \code{estimates} with cholesky=TRUE.
}

\author{
Cecile Proust-Lima, Viviane Philipps
}
\seealso{
\code{\link{estimates}}, \code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{multlcmm}}, \code{\link{Jointlcmm}}
}

