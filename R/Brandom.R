#' @export
Brandom <- function(theta0,v0,w,b0,chol=NULL,mult=0)
{
    b <- rep(NA,length(w)) # initialisation

    ## indices des prm non generes a partir de theta0
    pasrandom <- which((w==0))

    ## w sans les 0
    wr <- w[which(w!=0)]

    ## quels prm de theta0 sont utilises plusieurs fois
    po1 <- duplicated(wr)
    po2 <- duplicated(wr,fromLast=TRUE)
    po <- po1+po2
    #po[pasrandom] <- 0

    ## prm theta0 avec des repetitions (si mixture) et la variance associee
    bb <- theta0[w]
    vbb <- v0[w,w]

    ## mettre les covariance a 0 si les prm sont repetees
    d <- diag(vbb) # sauvegarder la diagonale
    vbb[which(po>0),] <- 0 
    vbb[,which(po>0)] <- 0
    diag(vbb) <- d # remettre la diagonale

    ## tirer aleatoirement
    if(any(vbb!=0))
    {
        ch <- chol(vbb)
        ch <- t(ch)
        br <- bb + ch %*% rnorm(length(bb))
    }
    else
    {
        br <- bb
    }
    
    ## rassembler br et b0
    b[which(w!=0)] <- br
    b[which(w==0)] <- b0

    ## transformer chol en var-cov
    if(!is.null(chol))
    {
        if(is.list(chol))
        {
            nvc <- sapply(chol,length)
        }
        else
        {
            nvc <- length(chol)
            chol <- list(chol)
        }
        
        for(k in 1:length(nvc))
        {
            if(nvc[k]>0)
            {
                if(mult==0) # cas general
                {
                    nea <- (-1+sqrt(1+8*nvc[k]))/2
                    chea <- matrix(0,nrow=nea,ncol=nea)
                    chea[upper.tri(chea,diag=TRUE)] <- b[chol[[k]]]
                    v <- t(chea)%*%chea
                    
                    b[chol[[k]]] <- v[upper.tri(v,diag=TRUE)]            
                }
        
                if(mult==1) # cas multlcmm : mettre 1 pour la premiere variance
                {
                    nea <- (-1+sqrt(1+8*(nvc[k]+1)))/2
                    chea <- matrix(0,nrow=nea,ncol=nea)
                    chea[upper.tri(chea,diag=TRUE)] <- c(1,b[chol[[k]]])
                    v <- t(chea)%*%chea
                    
                    b[chol[[k]]] <- v[upper.tri(v,diag=TRUE)][-1]
                }
            }
        }
    }
    

    return(b)
}
