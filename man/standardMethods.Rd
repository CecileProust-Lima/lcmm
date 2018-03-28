\name{methods}
\alias{coef.hlme}
\alias{coef.lcmm}
\alias{coef.Jointlcmm}
\alias{coef.multlcmm}  
\alias{vcov.hlme}
\alias{vcov.lcmm}
\alias{vcov.Jointlcmm}
\alias{vcov.multlcmm}  
\alias{fixef}
\alias{fixef.hlme}      
\alias{fixef.lcmm}      
\alias{fixef.Jointlcmm} 
\alias{fixef.multlcmm}  
\alias{ranef}
\alias{ranef.hlme}      
\alias{ranef.lcmm}      
\alias{ranef.Jointlcmm} 
\alias{ranef.multlcmm}
\alias{fitted.hlme}      
\alias{fitted.lcmm}      
\alias{fitted.Jointlcmm} 
\alias{fitted.multlcmm}  
\alias{residuals.hlme}      
\alias{residuals.lcmm}      
\alias{residuals.Jointlcmm} 
\alias{residuals.multlcmm}  

\title{
  Standard methods for estimated models
}
\description{
 coef, vcov, fixef, ranef, fitted and residuals methods for
estimated hlme, lcmm, mutlcmm and Jointlcmm models.
}
\usage{

\method{coef}{hlme}(object,...)
\method{coef}{lcmm}(object,...)
\method{coef}{Jointlcmm}(object,...)
\method{coef}{multlcmm}(object,...)

\method{vcov}{hlme}(object,...)
\method{vcov}{lcmm}(object,...)
\method{vcov}{Jointlcmm}(object,...)
\method{vcov}{multlcmm}(object,...)

\method{fixef}{hlme}(object,...)
\method{fixef}{lcmm}(object,...)
\method{fixef}{Jointlcmm}(object,...)
\method{fixef}{multlcmm}(object,...)

\method{ranef}{hlme}(object,...)
\method{ranef}{lcmm}(object,...)
\method{ranef}{Jointlcmm}(object,...)
\method{ranef}{multlcmm}(object,...)

\method{fitted}{hlme}(object,...)
\method{fitted}{lcmm}(object,...)
\method{fitted}{Jointlcmm}(object,...)
\method{fitted}{multlcmm}(object,...)

\method{residuals}{hlme}(object,...)
\method{residuals}{lcmm}(object,...)
\method{residuals}{Jointlcmm}(object,...)
\method{residuals}{multlcmm}(object,...)

}

\arguments{
  \item{object}{
an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}
}
\item{...}{
 other arguments. There are ignored in these functions. }
}

\value{
  For \code{coef}, the vector of the estimates.

  For \code{vcov}, the variance-covariance matrix of the estimates.

  For \code{fixef} :
   - for \code{hlme}, \code{lcmm} and \code{multlcmm} objects, a list
  containing the fixed effects estimates in the
  class-membership model and in the longitudinal model.
  - for \code{Jointlcmm} objects,  a list containing the fixed effects
  estimates in the class-membership model, the survival model and in the
  longitudinal model.
  
  For \code{ranef}, a matrix (nrow=number of subjects, ncol=number of
  covariates with random effect) containing the individual random
  effects. 

  For \code{fitted}, a vector containing the subject-specific
  predictions extracted from \code{object}.

  For \code{residuals}, a vector containing the subject-specific
  residuals extracted from \code{object}.
}

\author{
Cecile Proust-Lima, Viviane Philipps
}


