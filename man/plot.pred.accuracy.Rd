\name{plot.pred.accuracy}
\alias{plot.epoce}
\alias{plot.Diffepoce}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots 
}
\description{
This function displays plots related to predictive accuracy functions: \code{epoce} and \code{Diffepoce}.


These functions do not apply for the moment with multiple causes of event (competing risks).
}
\usage{
\method{plot}{epoce}(x,\dots)
\method{plot}{Diffepoce}(x,\dots)
}
\arguments{
  \item{x}{an object inheriting from classes \code{epoce} or \code{Diffepoce}
}
  \item{\dots}{
other parameters to be passed through to plotting functions
}
}
\value{
Returns plots related to \code{epoce} and \code{Diffepoce}
}
\details{
For \code{epoce} objects, the function displays the EPOCE estimate (either MPOL or CVPOL) according to the time of prediction. 
For \code{Diffepoce} objects, \code{plot} displays the difference in EPOCE estimates (either MPOL or CVPOL) and its 95\% tracking interval between two joint latent class models
}
\author{
Cecile Proust-Lima and Viviane Philipps
}



\seealso{
\code{\link{epoce}},\code{\link{Diffepoce}}
}


\examples{
\dontrun{
# estimation of the joint latent class model
m3 <- Jointlcmm(fixed= Ydep1~Time*X1,mixture=~Time,random=~Time,
classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~X1+mixture(X2),
hazard="3-quant-splines",hazardtype="PH",ng=3,data=data_lcmm,
B=c(0.7667, 0.4020, -0.8243, -0.2726, 0.0000, 0.0000, 0.0000, 0.3020,
-0.6212, 2.6247, 5.3139, -0.0255, 1.3595, 0.8172, -11.6867, 10.1668,
10.2355, 11.5137, -2.6209, -0.4328, -0.6062, 1.4718, -0.0378, 0.8505,
0.0366, 0.2634, 1.4981))
# predictive accuracy of the model evaluated with EPOCE
VecTime <- c(1,3,5,7,9,11,13,15)
cvpl <- epoce(m3,var.time="Time",pred.times=VecTime)
summary(cvpl)
plot(cvpl,bty="l",ylim=c(0,2))
}
}



