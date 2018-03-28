\name{summary.lcmm}
\alias{summary.hlme}
\alias{summary.lcmm}
\alias{summary.Jointlcmm}
\alias{summary.multlcmm}
\alias{summary.epoce}
\alias{summary.Diffepoce}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of a \code{hlme}, \code{lcmm}, \code{Jointlcmm}, \code{multlcmm}, \code{epoce} or \code{Diffepoce} objects
}
\description{
The function provides a summary of \code{hlme}, \code{lcmm}, \code{multlcmm} and \code{Jointlcmm} estimations, or \code{epoce} and \code{Diffepoce} computations.
}
\usage{
\method{summary}{hlme}(object,\dots)
\method{summary}{lcmm}(object,\dots)
\method{summary}{Jointlcmm}(object,\dots)
\method{summary}{epoce}(object,\dots)
\method{summary}{Diffepoce}(object,\dots)
\method{summary}{multlcmm}(object,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object inheriting from classes \code{hlme}, \code{lcmm}, \code{multlcmm} for fitted latent class mixed-effects, or class \code{Jointlcmm} for a Joint latent class mixed model or \code{epoce} or \code{Diffepoce} for predictive accuracy computations.
}
  \item{\dots}{
further arguments to be passed to or from other methods.
          They are ignored in this function.
}
}

\value{
For  \code{epoce} or \code{Diffepoce} objects, returns NULL.
For \code{hlme}, \code{lcmm}, \code{Jointlcmm} or \code{multlcmm} returns also
a matrix containing the fixed effect estimates in the longitudinal model, their standard errors, Wald statistics and p-values
}

\author{
Cecile Proust-Lima, Viviane Philipps, Amadou Diakite and Benoit Liquet
}
\seealso{
\code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{multlcmm}},
\code{\link{Jointlcmm}}, \code{epoce}, \code{Diffepoce}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
