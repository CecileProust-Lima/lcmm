% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.lcmm.R
\name{summary.lcmm}
\alias{summary.lcmm}
\alias{summary.hlme}
\alias{summary.Jointlcmm}
\alias{summary.multlcmm}
\alias{summary.epoce}
\alias{summary.Diffepoce}
\alias{summary.mpjlcmm}
\alias{summary.externSurv}
\alias{summary.externX}
\title{Summary of a \code{hlme}, \code{lcmm}, \code{Jointlcmm}, \code{multlcmm},
\code{mpjlcmm}, \code{externSurv}, \code{externX}
\code{epoce} or \code{Diffepoce} objects}
\usage{
\method{summary}{lcmm}(object, ...)
}
\arguments{
\item{object}{an object inheriting from classes \code{hlme}, \code{lcmm},
\code{multlcmm} for fitted latent class mixed-effects, or class
\code{Jointlcmm}, \code{mpjlcmm} for a Joint latent class mixed model or \code{epoce} or
\code{Diffepoce} for predictive accuracy computations or \code{externSurv}, \code{externX}
for secondary regression models.}

\item{\dots}{further arguments to be passed to or from other methods.  They
are ignored in this function.}
}
\value{
For \code{epoce} or \code{Diffepoce} objects, returns NULL. For
\code{hlme}, \code{lcmm}, \code{Jointlcmm} or \code{multlcmm} returns also a
matrix containing the fixed effect estimates in the longitudinal model,
their standard errors, Wald statistics and p-values
}
\description{
The function provides a summary of \code{hlme}, \code{lcmm}, \code{multlcmm}
and \code{Jointlcmm} estimations, or \code{epoce} and \code{Diffepoce}
computations.
}
\seealso{
\code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{multlcmm}},
\code{\link{Jointlcmm}}, \code{epoce}, \code{Diffepoce}
}
\author{
Cecile Proust-Lima, Viviane Philipps, Amadou Diakite and Benoit
Liquet
}
\keyword{print}
