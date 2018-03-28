\name{data_lcmm}
\alias{data_lcmm}
\docType{data}
\title{
Simulated dataset for lcmm and Jointlcmm functions
}
\description{
The data were simulated from a joint latent class mixed model with 3 classes. 
Repeated data of 3 longitudinal outcomes (Ydep1, Ydep2, Ydep3) and censored time of event (Tevent, Event) with delayed entry (Tentry) were simulated for a total of 300 subjects. 
The three latent classes were predicted by the continuous covariate X3. 
In each latent class, the longitudinal outcome Ydep1 followed a linear mixed model including intercept, time and squared time both with correlated random-effects and class-specific fixed effects. In addition, the binary covariate X1 and its interaction with time X1:Time had a common impact (over classes) on the Ydep1 trajectory.  The longitudinal ordinal outcomes Ydep2 and Ydep3 were generated from Ydep1 using threshold models with respectively 30 and 10 thresholds. 
In each latent class, the time of event followed a class-specific Weibull hazard with a common proportional effect of the binary covariate X2. Both time of entry Tentry and time of censoring had a uniform distribution
}
\usage{data_lcmm}
\format{
  A data frame with 1678 observations over 300 different subjects and 22 variables.
  \describe{
    \item{\code{ID}}{subject identification number}
    \item{\code{Ydep1}}{longitudinal continuous outcome}
    \item{\code{Ydep2}}{longitudinal ordinal outcome with 31 levels}
    \item{\code{Ydep3}}{longitudinal ordinal outcome with 11 levels}
    \item{\code{Tentry}}{delayed entry for the time-to-event}
    \item{\code{Tevent}}{observed time-to-event: either censoring time or time of event}
    \item{\code{Event}}{indicator that Tevent is the time of event}
    \item{\code{Time}}{time of measurement}
    \item{\code{X1}}{binary covariate}
    \item{\code{X2}}{binary covariate}
    \item{\code{X3}}{continuous covariate}
    \item{\code{X4}}{categorical covariate}
  }
}
\seealso{
\code{\link{Jointlcmm}}, \code{\link{lcmm}}, \code{\link{hlme}}   
}

\keyword{datasets}






