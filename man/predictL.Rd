\name{predictL}
\alias{predictL.lcmm}
\alias{predictL}
\alias{predictL.multlcmm}
\alias{predictL.Jointlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Class-specific marginal predictions in the latent process scale for
\code{lcmm}, \code{Jointlcmm} and \code{multlcmm} objects
}
\description{
This function provides a matrix containing the class-specific predicted trajectories computed in the latent process scale, that is the latent process underlying the curvilinear outcome(s), for a profile of covariates specified by the user. This function applies only to \code{lcmm} and \code{multlcmm} objects. The function \code{plot.predict} provides directly the plot of these class-specific predicted trajectories. The function \code{predictY} provides the class-specific predicted trajectories computed in the natural scale of the outcome(s).
} 
\usage{
\method{predictL}{lcmm}(x,newdata,var.time,na.action=1,confint=FALSE,...)
\method{predictL}{multlcmm}(x,newdata,var.time,na.action=1,confint=FALSE,...)
\method{predictL}{Jointlcmm}(x,newdata,var.time,na.action=1,confint=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from class \code{lcmm},\code{multlcmm} or
\code{Jointlcmm} representing a (joint) (latent class) mixed model involving a latent process and estimated link function(s).
}
  \item{newdata}{
data frame containing the data from which predictions are computed. The
data frame should include at least all the covariates listed in
x$Xnames2. Names in the data frame should be exactly x$Xnames2 that are the names of covariates specified in \code{lcmm} or \code{multlcmm} calls.
}
  \item{var.time}{
A character string containing the name of the variable that corresponds to time in the data frame (x axis in the plot).
}
 \item{na.action}{
Integer indicating how NAs are managed. The default is 1 for
'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
'na.pass' or 'na.exclude' are not implemented in the current version.
}
\item{confint}{
logical indicating if confidence should be provided. Default to FALSE.
  }
  \item{\dots}{
further arguments to be passed to or from other methods.
          They are ignored in this function.
}
}
\value{
  An object of class \code{predictL} with values :

  - \code{pred} : a matrix containing the class-specific
  predicted values in the latent process scale, the lower and the upper
  limits of the confidence intervals (if calculated).

   - \code{times} : the \code{var.time} variable from \code{newdata}
}
\author{
Cecile Proust-Lima, Viviane Philipps
}
\seealso{
\code{\link{plot.predict}}, \code{\link{predictY}}, \code{\link{lcmm}}   
}
\examples{
#### Prediction from a 2-class model with a Splines link function
\dontrun{
## fitted model
m<-lcmm(Ydep2~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
subject='ID',ng=2,data=data_lcmm,link="splines",B=c(
-0.175,      -0.191,       0.654,      -0.443, 
-0.345,      -1.780,       0.913,       0.016, 
 0.389,       0.028,       0.083,      -7.349, 
 0.722,       0.770,       1.376,       1.653, 
 1.640,       1.285))
summary(m)
## predictions for times from 0 to 5 for X1=0
newdata<-data.frame(Time=seq(0,5,length=100),
X1=rep(0,100),X2=rep(0,100),X3=rep(0,100))
predictL(m,newdata,var.time="Time")
## predictions for times from 0 to 5 for X1=1
newdata$X1 <- 1
predictY(m,newdata,var.time="Time")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
