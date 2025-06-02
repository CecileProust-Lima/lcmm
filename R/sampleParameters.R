#' Sample models parameters
#'
#' Generates parameters from models estimated with \code{hlme},
#' \code{lcmm}, \code{multlcmm}, \code{Jointlcmm} or \code{mpjlcmm} functions.
#' Parameters are randomly sampled from a Gaussian distribution using estimates
#' and corresponding variances of a model.
#'
#' @param x an object of class \code{hlme}, \code{lcmm}, \code{multlcmm},
#' \code{Jointlcmm}, or \code{mpjlcmm}.
#' @param cholesky optional logical indicating if cholesky parameters should be
#' returned. Default to FALSE, the variance-covariance parameters of the
#' random effects are returned.
#' @return a named vector containing the sampled parameters.
#' @author Viviane Philipps
#' 
#' @examples
#' m <- hlme(Y ~ Time * X1, random = ~ Time, subject = 'ID', data = data_hlme)
#'
#' sampleParameters(m)
#'
#' @export
sampleParameters <- function(x, cholesky = FALSE)
{
    if(!inherits(x, c("hlme", "lcmm", "multlcmm", "Jointlcmm", "mpjlcmm"))) stop('Use only with "hlme", "lcmm", "multlcmm", "Jointlcmm", or "mpjlcmm" objects')

    ## estimated parameters and their variance
    b <- estimates(x, cholesky = TRUE)
    V <- vcov(x)

    ## sample new parameters
    bdraw <- as.vector(mvtnorm::rmvnorm(1, mean = b, sigma = V))
    
    if(cholesky == TRUE)
    {
        names(bdraw) <- names(b)
        return(bdraw)
    }

    ## replacing cholesky parameters by variance-covariance parameters
    nea <- sum(x$idea)
    if(nea > 0)
    {
        if(!inherits(x, "mpjlcmm"))
        {
            
            if(inherits(x, c("hlme", "lcmm")))
            {
                avt <- sum(x$N[1:2])
                nvc <- x$N[3]
            }
            if(inherits(x, "multlcmm"))
            {
                avt <- x$N[3]
                nvc <- x$N[4]
            }
            if(inherits(x, "Jointlcmm"))
            {
                avt <- sum(x$N[1:4])
                nvc <- x$N[5]
            }

            if(nvc > 0)
            {
                if(x$idiag == 0)
                {
                    ch <- matrix(0, nea, nea)
                    if(!inherits(x, "multlcmm"))
                        ch[upper.tri(ch, diag = TRUE)] <- bdraw[avt + 1:nvc]
                    else
                        ch[upper.tri(ch, diag = TRUE)] <- c(1, bdraw[avt + 1:nvc])
                    
                    varcov <- t(ch) %*% ch
                    
                    if(!inherits(x, "multlcmm"))
                        bdraw[avt + 1:nvc] <- varcov[upper.tri(varcov, diag = TRUE)]
                    else
                        bdraw[avt + 1:nvc] <- varcov[upper.tri(varcov, diag = TRUE)][-1]
                }
                else
                {
                    bdraw[avt + 1:nvc] <- bdraw[avt + 1:nvc]^2
                }
            }
        }
        else # mpjlcmm
        {
            K <- x$K
            l <- 3
            if(x$nbevt>1) l <- 2+x$nbevt
            nef <- x$Nprm[l+1:K]
            ncontr <- x$Nprm[l+K+1:K]
            nvc <- x$Nprm[l+2*K+1:K]            
            contrainte <- x$contrainte
            idiag = x$idiag
            nea <- as.numeric(tapply(x$idea, rep(1:x$K, x$nv), sum))
            
            tmp <- sum(x$N[1:3]) 
            for(k in 1:K)
            {
                if(nvc[k]>0)
                {
                    if(idiag[k] == 0)
                    {
                        ch <- matrix(0, nea[k], nea[k])
                        if(contrainte[k] != 2)
                        {
                            ch[upper.tri(ch, diag = TRUE)] <- bdraw[tmp + nef[k] + ncontr[k] + 1:nvc[k]]
                            varcov <- t(ch) %*% ch
                            bdraw[tmp + nef[k] + ncontr[k] + 1:nvc[k]] <- varcov[upper.tri(varcov, diag = TRUE)]
                        }
                        else
                        {
                            ch[upper.tri(ch, diag = TRUE)] <- c(1, bdraw[tmp + nef[k] + ncontr[k] + 1:nvc[k]])
                            varcov <- t(ch) %*% ch
                            bdraw[tmp + nef[k] + ncontr[k] + 1:nvc[k]] <- varcov[upper.tri(varcov, diag = TRUE)][-1]
                        }
                    }
                    else
                    {
                        bdraw[tmp + nef[k] + ncontr[k] + 1:nvc[k]] <- bdraw[tmp + nef[k] + ncontr[k] + 1:nvc[k]]^2
                    }
                }

                tmp <- tmp + x$npmK[k]
            }
        }
    }

    names(bdraw) <- names(x$best)
    return(bdraw) 
}
