\name{print.lcmm}
\alias{print.hlme}
\alias{print.lcmm}
\alias{print.Jointlcmm}
\alias{print.multlcmm}
\alias{print.epoce}
\alias{print.Diffepoce}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Brief summary of a \code{hlme}, \code{lcmm}, \code{Jointlcmm},\code{multlcmm}, \code{epoce} or \code{Diffepoce} objects
}
\description{
The function provides a brief summary of \code{hlme}, \code{lcmm},\code{multlcmm} or \code{Jointlcmm} estimations, and \code{epoce} or \code{Diffepoce} computations.
}
\usage{
\method{print}{hlme}(x,\dots)
\method{print}{Jointlcmm}(x,\dots)
\method{print}{epoce}(x,\dots)
\method{print}{lcmm}(x,\dots)
\method{print}{Diffepoce}(x,\dots)
\method{print}{multlcmm}(x,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from classes \code{hlme}, \code{lcmm} for fitted latent class mixed-effects, or class \code{Jointlcmm}for a Joint latent class mixed model or \code{epoce} for predictive accuracy computations.
}
  \item{\dots}{
further arguments to be passed to or from other methods.
          They are ignored in this function.
}
}
\author{
Cecile Proust-Lima, Viviane Philipps, Amadou Diakite and Benoit Liquet
}
\seealso{
  \code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{Jointlcmm}}, \code{epoce}, \code{Diffepoce}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
