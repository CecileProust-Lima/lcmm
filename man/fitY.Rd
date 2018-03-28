\name{fitY}
\alias{fitY}
\alias{fitY.lcmm}
\alias{fitY.multlcmm}
\alias{fitY.Jointlcmm}
\title{
Marginal predictions of the longitudinal outcome(s) in their
natural scale from \code{lcmm}, \code{Jointlcmm} or \code{multlcmm} objects 
}
\description{
The function computes the marginal predictions of the longitudinal outcome(s) in their
natural scale on the individual data used for the estimation  from \code{lcmm}, \code{Jointlcmm} or \code{multlcmm} objects.
}
\usage{

\method{fitY}{lcmm}(x)
\method{fitY}{multlcmm}(x)
\method{fitY}{Jointlcmm}(x)
}

\arguments{
  \item{x}{
an object inheriting from classes \code{lcmm} or \code{multlcmm}.
}

}

\value{
  For \code{lcmm} and \code{Jointlcmm} objects, returns a matrix with ng+1 columns
  containing the subject identifier and the ng class-specific marginal predicted
  values.
  
  For \code{multlcmm} objects, returns a matrix with ng+2 columns
  containing the subject identifier, the outcome indicator and the ng
  class-specific predicted values.
  
}
\author{
Cecile Proust-Lima, Viviane Philipps
}
\seealso{
\code{\link{predictY}}, \code{\link{plot.lcmm}}
}

