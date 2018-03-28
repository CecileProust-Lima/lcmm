
# coef = estimates
coef.hlme <- function(object,...)
    {
        estimates.hlme(object)
    }
coef.lcmm <- function(object,...)
    {
        estimates.lcmm(object)
    }
coef.multlcmm <- function(object,...)
    {
        estimates.multlcmm(object)
    }
coef.Jointlcmm <- function(object,...)
    {
        estimates.Jointlcmm(object)
    }

# variance-covariance matrix of the estimates
vcov.hlme <- function(object,...)
    {
        VarCov.hlme(object)
    }
vcov.lcmm <- function(object,...)
    {
        VarCov.lcmm(object)
    }
vcov.multlcmm <- function(object,...)
    {
        VarCov.multlcmm(object)
    }
vcov.Jointlcmm <- function(object,...)
    {
        VarCov.Jointlcmm(object)
    }



# fixed effects 
fixef.hlme <- function(object,...)
    {
        res <- list(NA,NA)
        if(object$ng>1) res[[1]] <- as.vector(object$best[1:object$N[1]])
        if(object$N[2]>1) res[[2]] <- as.vector(object$best[object$N[1]+1:object$N[2]])
        return(res)
    }

fixef.lcmm <- function(object,...)
    {
        res <- list(NA,NA)
        if(object$ng>1) res[[1]] <- as.vector(object$best[1:object$N[1]])
        if(object$N[2]>1) res[[2]] <- as.vector(object$best[object$N[1]+1:object$N[2]])
        return(res)
    }

fixef.multlcmm <- function(object,...)
    {
        res <- list(NA,NA)
        if(object$ng>1) res[[1]] <- as.vector(object$best[1:object$N[1]])
        if(object$N[3]>1) res[[2]] <- as.vector(object$best[(object$N[1]+1):object$N[3]])
        return(res)
    }

fixef.Jointlcmm <- function(object,...)
    {
        res <- list(NA,NA,NA)
        if(object$ng>1) res[[1]] <- as.vector(object$best[1:object$N[1]])
        if(object$N[3]>1) res[[2]] <- as.vector(object$best[object$N[1]+object$N[2]+1:object$N[3]])
        if(object$N[4]>1) res[[3]] <- as.vector(object$best[object$N[1]+object$N[2]+object$N[3]+1:object$N[4]])
        return(res)
    }

fixef <- function(object,...) UseMethod("fixef")


# random effects
ranef.hlme <- function(object,...)
    {
        if(object$N[3]>0) return(object$predRE[,-1,drop=FALSE])
        else return(NA)
    }

ranef.lcmm <- function(object,...)
    {
        if(object$N[3]>0) return(object$predRE[,-1,drop=FALSE])
        else return(NA)
    }

ranef.multlcmm <- function(object,...)
    {
        if(object$N[4]>0) return(object$predRE[,-1,drop=FALSE])
        else return(NA)
    }

ranef.Jointlcmm <- function(object,...)
    {
        if(object$N[5]>0) return(object$predRE[,-1,drop=FALSE])
        else return(NA)
    }

ranef <- function(object,...) UseMethod("ranef")




#  subject-specific predictions
fitted.hlme <- function(object,...)
    {
        return(object$pred[,"pred_ss"])
    }
fitted.lcmm <- function(object,...)
    {
        return(object$pred[,"pred_ss"])
    }
fitted.multlcmm <- function(object,...)
    {
        return(object$pred[,"pred_ss"])
    }
fitted.Jointlcmm <- function(object,...)
    {
        return(object$pred[,"pred_ss"])
    }



# subject-specific  residuals
residuals.hlme <- function(object,...)
    {
        return(object$pred[,"obs"]-object$pred[,"pred_ss"])
    } 
residuals.lcmm <- function(object,...)
    {
        return(object$pred[,"obs"]-object$pred[,"pred_ss"])
    } 
residuals.multlcmm <- function(object,...)
    {
        return(object$pred[,"obs"]-object$pred[,"pred_ss"])
    } 
residuals.Jointlcmm <- function(object,...)
    {
        return(object$pred[,"obs"]-object$pred[,"pred_ss"])
    } 
