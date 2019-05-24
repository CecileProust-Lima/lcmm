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
    nef <- object$Nprm[3+1:K]
    ncontr <- object$Nprm[3+K+1:K]
    nvc <- object$Nprm[3+2*K+1:K]
    nw <- object$Nprm[3+3*K+1:K]
    ncor <- object$Nprm[3+4*K+1:K]
    nerr <- object$Nprm[3+5*K+1:K]
    nalea <- object$Nprm[3+6*K+1:K]
    ntr <- object$Nprm[3+7*K+1:sum(ny)]
    ntrtotK <- sapply(1:K,function(k) sum(ntr[sum(ny[1:k])-ny[k]+1:ny[k]]))

    npmtot <- nef+ncontr+nvc+nw+ncor+nerr+nalea+ntrtotK

    ##liste des modeles longitudinaux
    lK <- eval(object$call$longitudinal)
    res <- vector("list",K)

    ## variances des prm
    Vtot <- matrix(0,length(object$best),length(object$best))
    Vtot[upper.tri(Vtot,diag=TRUE)] <- object$V
    
    sumnpm <- 0
    sumV <- (nprob+nrisqtot+nvarxevt)*(nprob+nrisqtot+nvarxevt+1)/2
    sumnobs <- 0
    for(k in 1:K)
    {
        ## le k-ieme modele mixte avec les estimations du conjoint:
        mcall <- lK[[k]]$call
        if(ng>1)
        {
            mcall$B <- c(object$best[1:(ng-1)],object$best[nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k]])
        }
        else
        {
            mcall$B <- c(object$best[nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k]])
        }
        mcall$maxiter <- 0
        mcall$verbose <- FALSE
        m <- eval(mcall)

        ## ajouter les variances
        V <- Vtot[c(1:(ng-1),nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k]),c(1:(ng-1),nprob+nrisqtot+nvarxevt+sumnpm+1:npmtot[k])]
        m$V <- V[upper.tri(V,diag=TRUE)]
        m$conv <- object$conv

        ##predictions
        if(class(m)=="multlcmm") m$pred <- object$pred[which(object$pred[,2] %in% m$Ynames),]
        else m$pred <- object$pred[which(object$pred[,2]==object$Names$Yname[k]),-2]

        m$pprob <- object$pprob
        m$pprobY <- object$pprobY
        

        res[[k]] <- m

        sumnpm <- sumnpm + npmtot[k]
        sumV <- sumV + npmtot[k]*(npmtot[k]+1)/2
        sumnobs <- sumnobs + object$N[11+k]
    }

    return(res)
}
