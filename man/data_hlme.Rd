\name{data_hlme}
\alias{data_hlme}
\docType{data}
\title{
Simulated dataset for hlme function
}
\description{
The data were simulated from a 3-latent class linear mixed model. Repeated data for 100 subjects were simulated. The three latent classes are predicted by X2 and X3. In each latent class, Y follows a linear mixed model including intercept and time both with correlated random-effects and class-specific fixed effects. In addition, X1 and X1*time have a common impact over classes on the Y trajectory.  
}
\usage{data_hlme}
\format{
  A data frame with 326 observations on the following 9 variables.
  \describe{
    \item{\code{ID}}{subject identification number}
    \item{\code{Y}}{longitudinal outcome}
    \item{\code{Time}}{time of measurement}
    \item{\code{X1}}{binary covariate}
    \item{\code{X2}}{binary covariate}
    \item{\code{X3}}{binary covariate}
  }
}
\seealso{
\code{\link{hlme}}, \code{\link{postprob}}, \code{\link{summary.lcmm}}, \code{\link{plot.predict}}    
}


\keyword{datasets}
