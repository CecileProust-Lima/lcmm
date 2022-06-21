#' Conditional probabilities and item information given specified latent process values
#' for \code{lcmm} or \code{multlcmm}
#' object with ordinal outcomes.
#'
#' The function computes the conditional probability and information function of
#' each level of each ordinal outcome and the information function at the item level.
#' Confidence bands (and median) can be computed by a Monte Carlo approximation.
#' 
#' @param x an object inheriting from class  \code{lcmm} or \code{multlcmm}, 
#' representing a general (latent class) mixed model.
#' @param lprocess numeric vector containing the latent process values at which the
#' predictions should be computed.
#' @param condRE_Y for multlcmm objects only, logical indicating if the predictions
#' are conditional to the outcome-specific random-effects or not. Default to FALSE=
#' the predictions are marginal to these random effects.
#' @param nsim  number of points used in the numerical integration (Monte-Carlo) with
#' splines or Beta link functions. nsim should be relatively important
#' (nsim=200 by default).
#' @param draws optional boolean specifying whether median and confidence bands
#' of the predicted values should be computed (TRUE). A Monte Carlo approximation 
#' of the posterior distribution of the predicted values is computed and the median, 
#' 2.5\% and 97.5\% percentiles are given. Otherwise, the predicted values are 
#' computed at the point estimate. By default, draws=FALSE.
#' @param ndraws if draws=TRUE, ndraws specifies the number of draws that should be
#' generated to approximate the posterior distribution of the predicted values.
#' By default, ndraws=2000.
#' @param \dots further arguments to be passed to or from other methods.  They
#' are ignored in this function.
#'
#' @return An object of class \code{ItemInfo} with values :
#'
#' - \code{ItemInfo}: 
#' If draws=FALSE, returns a matrix with 3 columns: the first column indicates the
#' name of the outcome, the second indicates the latent process value and the last
#' is the computed Fisher information.
#' If draws=TRUE, returns a matrix with 5 columns: the name of the outcome, the
#' latent process value and the 50\%, 2.5\% and 97.5\% percentiles of the approximated
#' posterior distribution of information.
#' 
#' - \code{LevelInfo}: 
#' If draws=FALSE, returns a matrix with 5 columns: the first column indicates the
#' name of the outcome, the second indicates the outcome's level, the third indicates the
#' latent process value and the two last contain the probability and Fisher information.
#' If draws=TRUE, returns a matrix with 5 columns: the name of the outcome,
#' the outcome's level, the latent process value and the 50\%, 2.5\% and 97.5\%
#' percentiles of the approximated posterior distribution of the probability and information.
#'
#' - \code{object}: the model from which the computations are done.
#' 
#' - \code{IC}: indicator specifying if confidence intervals are computed.
#'
#' @examples
#' \dontrun{
#' ## This is a toy example to illustrate the information functions.
#' ## The binary outcomes are arbitrarily created, please do not
#' ## consider them as relevent indicators.
#' data_lcmm$Yord1 <- as.numeric(data_lcmm$Ydep1>10)
#' data_lcmm$Yord2 <- as.numeric(data_lcmm$Ydep2>25)
#' m <- multlcmm(Yord1+Yord2~Time+I(Time^2),random=~Time,subject='ID',ng=1,
#' data=data_lcmm,link="thresholds")
#' info <- ItemInfo(m,lprocess=seq(-4,4,length.out=100),draws=TRUE)
#' plot(info)
#' par(mfrow=c(1,2))
#' plot(info, which="LevelInfo", outcome="Yord1")
#' plot(info, which="LevelInfo", outcome="Yord2")
#' plot(info, which="LevelProb", outcome="Yord1")
#' plot(info, which="LevelProb", outcome="Yord2")
#' }
#'
#' @author Cecile Proust-Lima, Viviane Philipps
#' 
#' @export
#' 
ItemInfo <- function(x,lprocess,condRE_Y=FALSE,nsim=200,draws=FALSE,ndraws=2000,...)
{
    if(!inherits(x, c("lcmm","multlcmm"))) stop("Use only with lcmm or multlcmm objects")
    if(all(x$linktype!=3)) stop("This function is only available for ordinal outcome")
    
    if(x$conv!=1 & draws==TRUE)
    {
        cat("No confidence interval will be provided since the program did not converge properly \n")
        draws <- FALSE
    }

    condRE_Y <- as.numeric(condRE_Y)
    
    if(x$conv %in% c(1,2,3))
    {
        if(inherits(x,"multlcmm"))
        {
            ## cas multlcmm:
            debut <- x$N[3]+x$N[4]+x$N[5]+x$N[7] # nef+nvc+nw+ncor
            nalea <- x$N[6]
            ny <- x$N[8]
            nerr <- ny
            Ynames <- x$Ynames
            if(nalea==0) condRE_Y <- 1
            modalites <- unlist(x$modalites)
            nbmod <- x$nbmod
        }
        else
        {
            debut <- sum(x$N[1:4])-1 # nprob+nef+nvc+nw-1
            nalea <- 0
            ny <- 1
            nerr <- 0 # variance erreur fixee a 1
            Ynames <- as.character(x$call$fixed[2])
            modalites <- seq(x$linknodes[1],x$linknodes[2])[x$ide==1]
            nbmod <- sum(x$ide)+1
        }
        
        npm <- length(x$best)
        nbzitr <- rep(2,ny)
        if(inherits(x,"multlcmm"))
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
        
        info <- rep(0,maxmes*(2*sum(nbmod)+ny)) # Pkl, IklPkl, Ik

        if(!draws)
        {
            out <- .Fortran(C_iteminfo,
                            as.double(lambdatot),
                            as.integer(condRE_Y),
                            as.integer(nalea),
                            as.integer(ny),
                            as.integer(maxmes),
                            as.integer(npm),
                            as.double(x$best),
                            as.integer(debut),
                            as.integer(nbzitr),
                            as.integer(x$linktype),
                            as.integer(unlist(x$modalites)),
                            as.integer(nbmod),
                            as.integer(nsim),
                            as.integer(2*sum(nbmod)+ny),
                            info=as.double(info))

            out$info[which(!is.finite(out$info))] <- NA

            infolevel <- as.data.frame(matrix(NA, maxmes, 5))
            infoitem <- data.frame(Yname=rep(Ynames[which(x$linktype==3)],each=maxmes),
                                   Lprocess=rep(lambda,ny),
                                   Info=rep(NA, ny*maxmes))
            j <- 0
            jj <- 0
            for(k in 1:ny)
            {
                if(x$linktype[k] != 3) next
                
                infolevel[jj+1:(nbmod[k]*maxmes),1] <- Ynames[k]
                infolevel[jj+1:(nbmod[k]*maxmes),2] <- rep(x$modalites[[k]], each=maxmes)
                infolevel[jj+1:(nbmod[k]*maxmes),3] <- rep(lambda, nbmod[k])
                infolevel[jj+1:(nbmod[k]*maxmes),4] <- out$info[j+1:(nbmod[k]*maxmes)]
                infolevel[jj+1:(nbmod[k]*maxmes),5] <- out$info[j+nbmod[k]*maxmes+1:(nbmod[k]*maxmes)]

                infoitem[(k-1)*maxmes+1:maxmes,3] <- out$info[j+2*nbmod[k]*maxmes+1:maxmes]

                j <- j+(2*x$nbmod[k]+1)*maxmes
                jj <- jj+nbmod[k]*maxmes
            }

            colnames(infolevel) <- c("Yname","Level","Lprocess","Prob","Info")
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
                
                out <- .Fortran(C_iteminfo,
                                as.double(lambdatot),
                                as.integer(condRE_Y),
                                as.integer(nalea),
                                as.integer(ny),
                                as.integer(maxmes),
                                as.integer(npm),
                                as.double(bdraw),
                                as.integer(debut),
                                as.integer(nbzitr),
                                as.integer(x$linktype),
                                as.integer(unlist(x$modalites)),
                                as.integer(nbmod),
                                as.integer(nsim),
                                as.integer(2*sum(nbmod)+ny),
                                info=as.double(info))
                
                out$info[which(!is.finite(out$info))] <- NA
                ydraws <- cbind(ydraws,out$info)
            }

            f <- function(x) {
                quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
            }
            ydistr <- apply(ydraws,1,FUN=f)
            Ypred_50 <- matrix(ydistr[2,],ncol=1,byrow=FALSE)
            Ypred_2.5 <- matrix(ydistr[1,],ncol=1,byrow=FALSE)
            Ypred_97.5 <- matrix(ydistr[3,],ncol=1,byrow=FALSE)
            
            ## info <- data.frame(rep(Ynames,each=maxmes),
            ##                    rep(unlist(x$modalites), each=maxmes),
            ##                    rep(lambda,ny),
            ##                    Ypred_50,Ypred_2.5,Ypred_97.5)
            ## colnames(info) <- c("Yname","Lprocess","Ypred_50","Ypred_2.5","Ypred_97.5")

            
            infolevel <- as.data.frame(matrix(NA, maxmes, 9))
            infoitem <- data.frame(Yname=rep(Ynames[which(x$linktype==3)],each=maxmes),
                                   Lprocess=rep(lambda,ny),
                                   Info_50=rep(NA, ny*maxmes),
                                   Info_2.5=rep(NA, ny*maxmes),
                                   Info_97.5=rep(NA, ny*maxmes))
            j <- 0
            jj <- 0
            for(k in 1:ny)
            {
                if(x$linktype[k] != 3) next
                
                infolevel[jj+1:(nbmod[k]*maxmes),1] <- Ynames[k]
                infolevel[jj+1:(nbmod[k]*maxmes),2] <- rep(x$modalites[[k]], each=maxmes)
                infolevel[jj+1:(nbmod[k]*maxmes),3] <- rep(lambda, nbmod[k])
                infolevel[jj+1:(nbmod[k]*maxmes),4] <- Ypred_50[j+1:(nbmod[k]*maxmes)]
                infolevel[jj+1:(nbmod[k]*maxmes),5] <- Ypred_2.5[j+1:(nbmod[k]*maxmes)]
                infolevel[jj+1:(nbmod[k]*maxmes),6] <- Ypred_97.5[j+1:(nbmod[k]*maxmes)]
                infolevel[jj+1:(nbmod[k]*maxmes),7] <- Ypred_50[j+nbmod[k]*maxmes+1:(nbmod[k]*maxmes)]
                infolevel[jj+1:(nbmod[k]*maxmes),8] <- Ypred_2.5[j+nbmod[k]*maxmes+1:(nbmod[k]*maxmes)]
                infolevel[jj+1:(nbmod[k]*maxmes),9] <- Ypred_97.5[j+nbmod[k]*maxmes+1:(nbmod[k]*maxmes)]

                infoitem[(k-1)*maxmes+1:maxmes,3] <- Ypred_50[j+2*nbmod[k]*maxmes+1:maxmes]
                infoitem[(k-1)*maxmes+1:maxmes,4] <- Ypred_2.5[j+2*nbmod[k]*maxmes+1:maxmes]
                infoitem[(k-1)*maxmes+1:maxmes,5] <- Ypred_97.5[j+2*nbmod[k]*maxmes+1:maxmes]

                j <- j+(2*nbmod[k]+1)*maxmes
                jj <- jj+nbmod[k]*maxmes
            }

            colnames(infolevel) <- c("Yname","Level","Lprocess","Prob_50","Prob_2.5","Prob_97.5","Info_50","Info_2.5","Info_97.5")
        }

        res <- list(ItemInfo=infoitem, LevelInfo=infolevel, object=x, IC=draws)
    }
    else
    {
        res <- list(ItemInfo=NA, LevelInfo=NA, object=x)
    }
    
    class(res) <- "ItemInfo"
    return(res)
}
