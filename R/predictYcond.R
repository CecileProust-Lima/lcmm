#' Conditional predictions of a \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}
#' object in the natural scale of the longitudinal outcome(s) for specified
#' latent process values.
#'
#' The function computes the predicted values of the longitudinal markers in their
#' natural scale for specified values of the latent process. For splines and Beta
#' links, a Gauss-Hermite integration is used to numerically compute the predictions.
#' In addition, for any type of link function, confidence bands (and median) can be
#' computed by a Monte Carlo approximation of the posterior distribution of the
#' predicted values.
#' 
#' @param x an object inheriting from class \code{lcmm}, 
#' \code{Jointlcmm} or \code{multlcmm} representing a general latent class
#' mixed model.
#' @param lprocess numeric vector containing the latent process values at which the
#' predictions should be computed.
#' @param condRE_Y for multlcmm objects only, logical indicating if the predictions
#' are conditional to the outcome specific random effects or not. Default to FALSE,
#' the predictions are marginal to these random effects.
#' @param nsim  number of points used in the numerical integration (Monte-Carlo) with
#' splines or Beta link functions. nsim should be relatively important
#' (nsim=200 by default).
#' @param draws optional boolean specifying whether median and confidence bands
#' of the predicted values should be computed (TRUE) - whatever the type of
#' link function. A Monte Carlo approximation of the posterior distribution of the
#' predicted values is computed and the median, 2.5\% and 97.5\% percentiles
#' are given. Otherwise, the predicted values are computed at the point
#' estimate. By default, draws=FALSE.
#' @param ndraws if draws=TRUE, ndraws specifies the number of draws that should be
#' generated to approximate the posterior distribution of the predicted values.
#' By default, ndraws=2000.
#' @param \dots further arguments to be passed to or from other methods.  They
#' are ignored in this function.
#'
#' @return An object of class \code{predictYcond} with values :
#'
#' - \code{pred} : 
#' If draws=FALSE, returns a matrix with 3 columns : the first column indicates the
#' name of the outcome, the second indicates the latent process value and the last
#' is the computed prediction.
#' If draws=TRUE, returns a matrix with 5 columns : the name of the outcome, the
#' latent process value and the 50\%, 2.5\% and 97.5\% percentiles of the approximated
#' posterior distribution of predicted values.
#'
#' - \code{object} : the model from which the predictions are computed.
#'
#' @examples
#' \dontrun{
#' m12 <- lcmm(Ydep2~Time+I(Time^2),random=~Time,subject='ID',ng=1,
#' data=data_lcmm,link="3-equi-splines")
#' predm12 <- predictYcond(m12,lprocess=seq(-8,2,length.out=100),draws=TRUE)
#' plot(predm12)
#' }
#'
#' @author Cecile Proust-Lima, Viviane Philipps
#' @seealso \code{\link{predictY}}, \code{\link{predictlink}}
#' 
#' @export
#' 
predictYcond <- function(x,lprocess,condRE_Y=FALSE,nsim=200,draws=FALSE,ndraws=2000,...)
{
    if((!class(x) %in% c("lcmm","multlcmm","Jointlcmm"))) stop("Use only with lcmm, multlcmm or Jointlcmm objects")
    if((class(x)=="lcmm") & any(x$linktype==3)) stop("This function is not available for ordinal outcome")
    if((class(x)=="Jointlcmm") & any(x$linktype==-1)) stop("This function is not available without any link function")
    
    if(x$conv!=1 & draws==TRUE)
    {
        cat("No confidence interval will be provided since the program did not converge properly \n")
        draws <- FALSE
    }

    condRE_Y <- as.numeric(condRE_Y)
    
    if(x$conv %in% c(1,2,3))
    {
        if(class(x)=="multlcmm")
        {
            ## cas multlcmm:
            debut <- x$N[3]+x$N[4]+x$N[5]+x$N[7] # nef+nvc+nw+ncor
            nalea <- x$N[6]
            ny <- x$N[8]
            nerr <- ny
            Ynames <- x$Ynames
            if(nalea==0) condRE_Y <- 1
        }
        else
        {
            if(class(x)=="lcmm") # lcmm continu
            {
                debut <- sum(x$N[1:4])-1 # nprob+nef+nvc+nw-1
            }
            else #Jointlcmm avec link
            {
                debut <- sum(x$N[1:7])-1  # nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor-1
            }
            
            nalea <- 0
            ny <- 1
            nerr <- 0 # variance erreur fixee a 1
            Ynames <- as.character(x$call$fixed[2])
        }
        
        npm <- length(x$best)
        nbzitr <- rep(2,ny)
        if(class(x)=="multlcmm")
        {
            nbzitr[which(x$linktype==2)] <- x$nbnodes
        }
        else
        {
            nbzitr[which(x$linktype==2)] <- length(x$linknodes)
        }
                
        lambda <- na.omit(lprocess)
        maxmes <- length(lambda)
        lambdatot <- rep(lambda,ny)
        
        Ycond <- rep(0,maxmes*ny)

        if(!draws)
        {
            out <- .Fortran(C_predictcondmult,
                            as.double(lambdatot),
                            as.integer(condRE_Y),
                            as.integer(nalea),
                            as.integer(ny),
                            as.integer(nerr),
                            as.integer(maxmes),
                            as.integer(npm),
                            as.double(x$best),
                            as.integer(debut),
                            as.integer(x$epsY),
                            as.integer(x$linktype),
                            as.integer(nbzitr),
                            as.double(x$linknodes),
                            as.integer(unlist(x$modalites)),
                            as.integer(x$nbmod),
                            as.integer(nsim),
                            Ycond=as.double(Ycond))
            
            out$Ycond[out$Ycond==9999] <- NA
            
            Ycond <- data.frame(rep(Ynames,each=maxmes),rep(lambda,ny),out$Ycond)
            colnames(Ycond) <- c("Yname","Lprocess","Ypred")
        }
        else
        {
            ndraws <- as.integer(ndraws)
            ydraws <- NULL
            
            posfix <- eval(x$call$posfix)
            
            if(ndraws>0)
            {
                Mat <- matrix(0,ncol=npm,nrow=npm)
                Mat[upper.tri(Mat,diag=TRUE)]<- x$V
                if(length(posfix))
                {
                    Mat2 <- Mat[-posfix,-posfix]
                    Chol2 <- chol(Mat2)
                    Chol <- matrix(0,npm,npm)
                    Chol[setdiff(1:npm,posfix),setdiff(1:npm,posfix)] <- Chol2
                    Chol <- t(Chol)
                }
                else
                {
                    Chol <- chol(Mat)
                    Chol <- t(Chol)
                }
            }
            
            best <- x$best
            nef <- x$N[3]
            nvc <- x$N[4]
            if(nvc>0) 
            {
                if(x$idiag==0) best[nef+1:nvc] <- x$cholesky[-1]
                else best[nef+1:nvc] <- x$cholesky[c((1:(nvc+1)*2:(nvc+2))/2)[-1]]
            }
            
            
            for (j in 1:ndraws)
            {  
                bdraw <- rnorm(npm)
                bdraw <- best + Chol %*% bdraw
                
                out <- .Fortran(C_predictcondmult,
                                as.double(lambdatot),
                                as.integer(condRE_Y),
                                as.integer(nalea),
                                as.integer(ny),
                                as.integer(nerr),
                                as.integer(maxmes),
                                as.integer(npm),
                                as.double(bdraw),
                                as.integer(debut),
                                as.integer(x$epsY),
                                as.integer(x$linktype),
                                as.integer(nbzitr),
                                as.double(x$linknodes),
                                as.integer(unlist(x$modalites)),
                                as.integer(x$nbmod),
                                as.integer(nsim),
                                Ycond=as.double(Ycond))
                
                out$Ycond[out$Ycond==9999] <- NA
                ydraws <- cbind(ydraws,out$Ycond)
            }

            f <- function(x) {
                quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
            }
            ydistr <- apply(ydraws,1,FUN=f)
            Ypred_50 <- matrix(ydistr[2,],ncol=1,byrow=FALSE)
            Ypred_2.5 <- matrix(ydistr[1,],ncol=1,byrow=FALSE)
            Ypred_97.5 <- matrix(ydistr[3,],ncol=1,byrow=FALSE)
            
            Ycond <- data.frame(rep(Ynames,each=maxmes),rep(lambda,ny),Ypred_50,Ypred_2.5,Ypred_97.5)
            colnames(Ycond) <- c("Yname","Lprocess","Ypred_50","Ypred_2.5","Ypred_97.5")
        }

        res <- list(pred=Ycond,object=x)
    }
    else
    {
        res <- list(pred=NA,object=x)
    }
    
    class(res) <- "predictYcond"
    return(res)
}
