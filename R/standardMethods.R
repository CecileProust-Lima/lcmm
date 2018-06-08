
#' Standard methods for estimated models
#' 
#' coef, vcov, fixef, ranef, fitted and residuals methods for estimated hlme,
#' lcmm, mutlcmm and Jointlcmm models.
#' 
#' 
#' 
#' 
#'
#' 
#' @aliases coef.hlme coef.lcmm coef.Jointlcmm coef.multlcmm vcov.hlme
#' vcov.lcmm vcov.Jointlcmm vcov.multlcmm fixef fixef.hlme fixef.lcmm
#' fixef.Jointlcmm fixef.multlcmm ranef ranef.hlme ranef.lcmm ranef.Jointlcmm
#' ranef.multlcmm fitted.hlme fitted.lcmm fitted.Jointlcmm fitted.multlcmm
#' residuals.hlme residuals.lcmm residuals.Jointlcmm residuals.multlcmm
#' @param object an object of class \code{hlme}, \code{lcmm}, \code{multlcmm}
#' or \code{Jointlcmm}
#' @param ...  other arguments. There are ignored in these functions.
#' @return For \code{coef}, the vector of the estimates.
#' 
#' For \code{vcov}, the variance-covariance matrix of the estimates.
#' 
#' For \code{fixef} : - for \code{hlme}, \code{lcmm} and \code{multlcmm}
#' objects, a list containing the fixed effects estimates in the
#' class-membership model and in the longitudinal model.  - for
#' \code{Jointlcmm} objects, a list containing the fixed effects estimates in
#' the class-membership model, the survival model and in the longitudinal
#' model.
#' 
#' For \code{ranef}, a matrix (nrow=number of subjects, ncol=number of
#' covariates with random effect) containing the individual random effects.
#' 
#' For \code{fitted}, a vector containing the subject-specific predictions
#' extracted from \code{object}.
#' 
#' For \code{residuals}, a vector containing the subject-specific residuals
#' extracted from \code{object}.
#' @author Cecile Proust-Lima, Viviane Philipps
#' 
#' @name StandardMethods
#' 

#' 





# coef = estimates
#' @export
coef.hlme <- function(object,...)
    {
        estimates.hlme(object)
    }
#' @export
coef.lcmm <- function(object,...)
    {
        estimates.lcmm(object)
    }
#' @export
coef.multlcmm <- function(object,...)
    {
        estimates.multlcmm(object)
    }
#' @export
coef.Jointlcmm <- function(object,...)
    {
        estimates.Jointlcmm(object)
    }

# variance-covariance matrix of the estimates
#' @export
vcov.hlme <- function(object,...)
    {
        VarCov.hlme(object)
    }
#' @export
vcov.lcmm <- function(object,...)
    {
        VarCov.lcmm(object)
    }
#' @export
vcov.multlcmm <- function(object,...)
    {
        VarCov.multlcmm(object)
    }
#' @export
vcov.Jointlcmm <- function(object,...)
    {
        VarCov.Jointlcmm(object)
    }



# fixed effects 
#' @export
fixef.hlme <- function(object,...)
    {
        res <- list(NA,NA)
        if(object$ng>1) res[[1]] <- as.vector(object$best[1:object$N[1]])
        if(object$N[2]>1) res[[2]] <- as.vector(object$best[object$N[1]+1:object$N[2]])
        return(res)
    }

#' @export
fixef.lcmm <- function(object,...)
    {
        res <- list(NA,NA)
        if(object$ng>1) res[[1]] <- as.vector(object$best[1:object$N[1]])
        if(object$N[2]>1) res[[2]] <- as.vector(object$best[object$N[1]+1:object$N[2]])
        return(res)
    }

#' @export
fixef.multlcmm <- function(object,...)
    {
        res <- list(NA,NA)
        if(object$ng>1) res[[1]] <- as.vector(object$best[1:object$N[1]])
        if(object$N[3]>1) res[[2]] <- as.vector(object$best[(object$N[1]+1):object$N[3]])
        return(res)
    }

#' @export
fixef.Jointlcmm <- function(object,...)
    {
        res <- list(NA,NA,NA)
        if(object$ng>1) res[[1]] <- as.vector(object$best[1:object$N[1]])
        if(object$N[3]>1) res[[2]] <- as.vector(object$best[object$N[1]+object$N[2]+1:object$N[3]])
        if(object$N[4]>1) res[[3]] <- as.vector(object$best[object$N[1]+object$N[2]+object$N[3]+1:object$N[4]])
        return(res)
    }

#' @export
fixef <- function(object,...) UseMethod("fixef")


# random effects
#' @export
ranef.hlme <- function(object,...)
    {
        if(object$N[3]>0) return(object$predRE[,-1,drop=FALSE])
        else return(NA)
    }

#' @export
ranef.lcmm <- function(object,...)
    {
        if(object$N[3]>0) return(object$predRE[,-1,drop=FALSE])
        else return(NA)
    }

#' @export
ranef.multlcmm <- function(object,...)
    {
        if(object$N[4]>0) return(object$predRE[,-1,drop=FALSE])
        else return(NA)
    }

#' @export
ranef.Jointlcmm <- function(object,...)
    {
        if(object$N[5]>0) return(object$predRE[,-1,drop=FALSE])
        else return(NA)
    }

#' @export
ranef <- function(object,...) UseMethod("ranef")




#  subject-specific predictions
#' @export
fitted.hlme <- function(object,...)
    {
        return(object$pred[,"pred_ss"])
    }
#' @export
fitted.lcmm <- function(object,...)
    {
        return(object$pred[,"pred_ss"])
    }
#' @export
fitted.multlcmm <- function(object,...)
    {
        return(object$pred[,"pred_ss"])
    }
#' @export
fitted.Jointlcmm <- function(object,...)
    {
        return(object$pred[,"pred_ss"])
    }



# subject-specific  residuals
#' @export
residuals.hlme <- function(object,...)
    {
        return(object$pred[,"obs"]-object$pred[,"pred_ss"])
    } 
#' @export
residuals.lcmm <- function(object,...)
    {
        return(object$pred[,"obs"]-object$pred[,"pred_ss"])
    } 
#' @export
residuals.multlcmm <- function(object,...)
    {
        return(object$pred[,"obs"]-object$pred[,"pred_ss"])
    } 
#' @export
residuals.Jointlcmm <- function(object,...)
    {
        return(object$pred[,"obs"]-object$pred[,"pred_ss"])
} 














