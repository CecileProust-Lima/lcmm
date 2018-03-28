\name{predictY}
\alias{predictY}
\alias{predictY.lcmm}
\alias{predictY.hlme}
\alias{predictY.Jointlcmm}
\alias{predictY.multlcmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Marginal predictions (possibly class-specific) of a \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm} object in the natural scale of the longitudinal outcome(s) for a specified profile of covariates.
}
\description{
For \code{hlme} and \code{Jointlcmm} objects, the function computes the predicted values of the longitudinal marker in each latent class for a specified profile of covariates. 
For \code{lcmm} and \code{multlcmm} objects, the function computes predicted values in the natural scale of the outcomes for a specified profile of covariates. For linear and threshold links, the predicted values are computed analytically. For splines and Beta links, a Gauss-Hermite or Monte-Carlo integration are used to numerically compute the predictions. In addition, for any type of link function, confidence bands (and median) can be computed by a Monte Carlo approximation of the posterior distribution of the predicted values.
}
\usage{
\method{predictY}{lcmm}(x,newdata,var.time,methInteg=0,nsim=20,draws=FALSE,
ndraws=2000,na.action=1,...)
\method{predictY}{hlme}(x,newdata,var.time,draws=FALSE,na.action=1,...)
\method{predictY}{Jointlcmm}(x,newdata,var.time,methInteg=0,nsim=20,draws=FALSE,
ndraws=2000,na.action=1,...)
\method{predictY}{multlcmm}(x,newdata,var.time,methInteg=0,nsim=20,draws=FALSE,
ndraws=2000,na.action=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting from class \code{lcmm}, \code{hlme}, \code{Jointlcmm} or \code{multlcmm} representing a general latent class mixed model.
}
  \item{newdata}{
data frame containing the data from which predictions are computed. The data frame should include at least all the covariates listed in x$Xnames2. Names in the data frame should be exactly x$Xnames2 that are the names of covariates specified in \code{lcmm}, \code{hlme}, \code{Jointlcmm} or \code{multlcmm} calls.
}
  \item{var.time}{
A character string containing the name of the variable that corresponds to time in the data frame (x axis in the plot).
}
  \item{methInteg}{
optional integer specifying the type of numerical integration required only for predictions with splines or Beta link functions. Value 0 (by default) specifies a Gauss-Hermite integration which is very rapid but neglects the correlation between the predicted values (in presence of random-effects). Value 1 refers to a Monte-Carlo integration which is slower but correctly account for the correlation between the predicted values.
}
  \item{nsim}{For a \code{lcmm}, \code{multlcmm} or \code{Jointlcmm} object only; optional number of points used in the numerical integration with splines or Beta link functions. For methInteg=0, nsim should be chosen among the following values: 5, 7, 9, 15, 20, 30, 40 or 50 (nsim=20 by default). If methInteg=1, nsim should be relatively important (more than 200).
}
  \item{draws}{optional boolean specifying whether median and confidence bands of the predicted values should be computed (TRUE) - whatever the type of link function. For a \code{lcmm}, \code{multlcmm} or \code{Jointlcmm} object, a Monte Carlo approximation of the posterior distribution of the predicted values is computed and the median, 2.5\% and 97.5\% percentiles are given. Otherwise, the predicted values are computed at the point estimate. By default, draws=FALSE.  
}
  \item{ndraws}{For a \code{lcmm}, \code{multlcmm} or \code{Jointlcmm} object only; if draws=TRUE, ndraws specifies the number of draws that should be generated to approximate the posterior distribution of the predicted values. By default, ndraws=2000.
}
 \item{na.action}{
Integer indicating how NAs are managed. The default is 1 for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as 'na.pass' or 'na.exclude' are not implemented in the current version.}
  \item{\dots}{
further arguments to be passed to or from other methods.
          They are ignored in this function.
}
}
\value{
  An object of class \code{predictY} with values :
  
  - \code{pred} : a matrix  with the same rows (number and order) as in
  newdata.
  
For \code{hlme}  objects and  \code{lcmm} or  \code{Jointlcmm} with \code{draws=FALSE}, returns a matrix with ng columns corresponding to the ng class-specific vectors of predicted values computed at the point estimate 

For objects of class \code{lcmm} or \code{Jointlcmm} with \code{draws=TRUE}, returns a matrix with ng*3 columns representing the ng class-specific 50\%, 2.5\% and 97.5\% percentiles of the approximated posterior distribution of the class-specific predicted values.

For objects of class \code{multlcmm} with \code{draws=FALSE}, returns a matrix with ng+1 columns: the first column indicates the name of the outcome which is predicted and the ng subsequent columns correspond to the ng class-specific vectors of predicted values computed at the point estimate 

For objects of class \code{multlcmm} with \code{draws=TRUE}, returns a matrix with ng*3+1 columns: the first column indicates the name of the outcome which is predicted and the ng*3 subsequent columns correspond to the ng class-specific 50\%, 2.5\% and 97.5\% percentiles of the approximated posterior distribution of the class-specific predicted values.

- \code{times} : the \code{var.time} variable from \code{newdata}
}
\author{
Cecile Proust-Lima, Viviane Philipps
}
\seealso{
\code{\link{lcmm}}, \code{\link{multlcmm}},  \code{\link{hlme}}, \code{\link{Jointlcmm}}
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
pred0 <- predictY(m,newdata,var.time="Time")
head(pred0)
## Option draws=TRUE to compute a MonteCarlo 
# approximation of the predicted value distribution 
# (quite long with ndraws=2000 by default)
\dontrun{
pred0MC <- predictY(m,newdata,draws=TRUE,var.time="Time")
}
## predictions for times from 0 to 5 for X1=1
newdata$X1 <- 1
pred1 <- predictY(m,newdata,var.time="Time")
## Option draws=TRUE to compute a MonteCarlo 
# approximation of the predicted value distribution 
# (quite long with ndraws=2000 by default)
\dontrun{
pred1MC <- predictY(m,newdata,draws=TRUE,var.time="Time")
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
