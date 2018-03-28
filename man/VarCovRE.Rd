\name{VarCovRE}
\alias{VarCovRE}
\alias{VarCovRE.hlme}
\alias{VarCovRE.lcmm}
\alias{VarCovRE.Jointlcmm}
\alias{VarCovRE.multlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates, standard errors and Wald test for the parameters of the variance-covariance matrix of the random effects.
}
\description{
Fromm the Cholesky transformed parameters, this function provides estimates, standard errors and Wald test for the parameters of the variance-covariance matrix of the random effects.
}
\usage{
\method{VarCovRE}{hlme}(Mod)
\method{VarCovRE}{lcmm}(Mod)
\method{VarCovRE}{Jointlcmm}(Mod)
\method{VarCovRE}{multlcmm}(Mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mod}{
an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}
}
}

\value{
a matrix containing the estimates of the parameters of the variance-covariance matrix of the random effects, their standard errors, and, for the covariance parameters,
the Wald statistic and the associated p-value.
}

\author{
Cecile Proust-Lima, Lionelle Nkam and Viviane Philipps
}