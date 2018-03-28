\name{cuminc}
\alias{cuminc}
\title{Predicted cumulative incidence of event according to a profile of covariates
}
\description{
This function computes the predicted cumulative incidence of each cause of event according to a profile of covariates from a joint latent class model. Confidence bands can be computed by a Monte-Carlo method.
}
\usage{
cuminc(x, time, draws = FALSE, ndraws = 2000, ...)
}

\arguments{
  \item{x}{
an object inheriting from class \code{Jointlcmm}
}
  \item{time}{
a vector of times at which the cumulative incidence is calculated
}
\item{draws}{
  optional boolean specifying whether a Monte Carlo approximation of the posterior distribution of the cumulative incidence is computed and the median, 2.5\% and 97.5\% percentiles are given. Otherwise, the predicted cumulative incidence is computed at the point estimate. By default, draws=FALSE.
}
  \item{ndraws}{
if draws=TRUE, ndraws specifies the number of draws that should be generated to approximate the posterior distribution of the predicted cumulative incidence. By default, ndraws=2000.
}
  \item{\dots}{
further arguments, in particular values of the covariates specified in
  the survival part of the joint model. 
}
}

\value{
  An object of class \code{cuminc} containing as many matrices as profiles defined by the covariates
  values. Each of these matrices contains the event-specific cumulative incidences in
  each latent class at the different times specified.
}

\author{
Viviane Philipps and Cecile Proust-Lima
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Jointlcmm}},\code{\link{plot.Jointlcmm}},\code{\link{plot.cuminc}}
}

