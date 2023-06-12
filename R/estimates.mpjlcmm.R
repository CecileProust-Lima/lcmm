#' @export
estimates.mpjlcmm <- function(x,cholesky=TRUE)
{
    if(missing(x)) stop("The argument x should be specified")
    if (!inherits(x, "mpjlcmm")) stop("use only with \"mpjlcmm\" objects")
    if(is.na(as.logical(cholesky))) stop("cholesky should be TRUE or FALSE")

    if(x$conv==1 | x$conv==2 | x$conv==3)
    {
        res <- x$best
        cholesky <- as.logical(cholesky)

        if(isTRUE(cholesky) & x$N[6]>0)
        {
            K <- x$K
            ny <- x$ny
            l <- 3
            if(x$nbevt>1) l <- 2+x$nbevt
            nef <- x$Nprm[l+1:K]
            ncontr <- x$Nprm[l+K+1:K]
            nvc <- x$Nprm[l+2*K+1:K]
            nw <- x$Nprm[l+3*K+1:K]
            ncor <- x$Nprm[l+4*K+1:K]
            nerr <- x$Nprm[l+5*K+1:K]
            nalea <- x$Nprm[l+6*K+1:K]
            ntr <- x$Nprm[l+7*K+1:sum(ny)]
            ch <- x$cholesky
            contrainte <- x$contrainte

            tmp <- sum(x$N[1:3]) #nprob+nrisqtot+nvarxevt
            jch <- 0
            sumny <- 0
            for(k in 1:K)
            {
                if(nvc[k]>0)
                {
                    res[tmp + nef[k]+ncontr[k]+1:nvc[k]] <- ch[jch + 1:nvc[k]]
                    jch <- jch + nvc[k]
                }
                
                tmp <- tmp + nef[k] + ncontr[k] + nvc[k] + nw[k] + ncor[k] + nerr[k] + nalea[k] + sum(ntr[sumny+1:ny[k]])
                sumny <- sumny + ny[k]
            }
            
            names(res) <- sub("varcov","cholesky",names(res))
        }
    }
    else
    {
        res <- NA
        cat("Output can not be produced since the program stopped abnormally. \n")
    }
    
    return(res)
}
