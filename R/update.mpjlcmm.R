#' Update the longitudinal submodels
#'
#' This function updates the longitudinal submodels of a mpjlcmm object by injecting
#' the estimated parameters and their variances in each hlme/lcmm/multlcmm model used
#' to define the multi-process joint model.
#' The same (uni-process) models as specified in the mpjlcmm call are returned, with
#' updated outputs for best, V, conv, cholesky, pred, predRE, predRE_Y, pprob.
#' All postfit functions (plots, predictions, ...) can then be called on the updated hlme/lcmm/multlcmm models.
#' See mpjlcmm's help page for examples.
#'
#' @param object an estimated mpjlcmm model
#' @param \dots not used
#'
#' @return A list of hlme/lcmm/multlcmm models. The models appear in the same order as
#' specified in the call to the mpjlcmm function.
#' 
#' @export
update.mpjlcmm <- function(object,...)
{

    # nombre de parametres
    K <- object$K
    ng <- object$ng
    ny <- object$ny
    nprob <- object$N[1]
    nrisqtot <- object$N[2]
    nvarxevt <- object$N[3]
    avt <- 3
    if(object$nbevt>0) avt <- 2+object$nbevt
    nef <- object$Nprm[avt+1:K]
    ncontr <- object$Nprm[avt+K+1:K]
    nvc <- object$Nprm[avt+2*K+1:K]
    nw <- object$Nprm[avt+3*K+1:K]
    ncor <- object$Nprm[avt+4*K+1:K]
    nerr <- object$Nprm[avt+5*K+1:K]
    nalea <- object$Nprm[avt+6*K+1:K]
    ntr <- object$Nprm[avt+7*K+1:sum(ny)]
    ntrtotK <- sapply(1:K,function(k) sum(ntr[sum(ny[1:k])-ny[k]+1:ny[k]]))

    npmtot <- nef+ncontr+nvc+nw+ncor+nerr+nalea+ntrtotK

    posfix <- eval(object$call$posfix)

    ##liste des modeles longitudinaux
    lK <- eval(object$call$longitudinal)
    res <- vector("list",K)

    ## variances des prm
    Vtot <- matrix(0,length(object$best),length(object$best))
    Vtot[upper.tri(Vtot,diag=TRUE)] <- object$V
    
    sumnpm <- 0
    for(k in 1:K)
    {
        ## le k-ieme modele mixte avec les estimations du conjoint:
        mcall <- object$longicall[[k]] #lK[[k]]$call
        fixk <- NULL
        if(ng>1)
        {
            mcall$B <- c(object$best[1:(ng-1)],object$best[nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k]])
            if(length(posfix))
            {
                if(any(posfix %in% c(1:(ng-1))))
                {
                    fixk <- posfix[which(posfix %in% c(1:(ng-1)))]
                }
                if(any(posfix %in% c(nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k])))
                {
                    fixk <- c(fixk, posfix[which(posfix %in% c(nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k]))] - (nrisqtot+nvarxevt+sumnpm))
                }
            }
        }
        else
        {
            mcall$B <- c(object$best[nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k]])
            if(length(posfix))
            {
                if(any(posfix %in% c(nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k])))
                {
                    fixk <- c(fixk, posfix[which(posfix %in% c(nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k]))] - (nrisqtot+nvarxevt+sumnpm))
                }
            }
        }
        mcall$posfix <- fixk
        mcall$maxiter <- 0
        mcall$verbose <- FALSE
        if(length(fixk) == (ng-1+npmtot[k])) m <- lK[[k]] else m <- eval(mcall)

        ## ajouter les variances
        if(ng>1)
        {
            V <- Vtot[c(1:(ng-1),nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k]),c(1:(ng-1),nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k])]
        }
        else
        {
            V <- Vtot[c(nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k]),c(nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k])]
        }
        
        m$V <- V[upper.tri(V,diag=TRUE)]
        m$conv <- object$conv

        ##predictions
        if(inherits(m,"multlcmm")) m$pred[,1:ncol(object$pred)] <- object$pred[which(object$pred[,2] %in% m$Ynames),]
        else m$pred[,1:(ncol(object$pred)-1)] <- object$pred[which(object$pred[,2]==object$Names$Yname[k]),-2]

        m$pprob <- object$pprob
        m$pprobY <- object$pprobY
        

        res[[k]] <- m

        sumnpm <- sumnpm + npmtot[k]
    }

    return(res)
}
