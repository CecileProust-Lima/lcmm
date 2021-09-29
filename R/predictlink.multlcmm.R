#' @export
predictlink.multlcmm <- function(x,ndraws=2000,Yvalues,...)
{
    ## verification des arguments
    if(missing(x)) stop("The model should be specified.")
    if(!(inherits(x,"multlcmm"))) stop("To use only with \"multlcmm\" objects")
    if(!missing(Yvalues) & all(x$linktype==3)) warning("With thresholds links, no \"Yvalues\" should be specified. Default values will be used. \n")
    if(x$conv!=1 & ndraws!=0) stop("No confidence intervals can be produced since the program did not converge properly")


    if(x$conv %in% c(1,2,3))
        {  
            ny <- x$N[8]

            if(missing(Yvalues))
                {
                    new.transf <- FALSE
                    Yvalues <- as.vector(x$estimlink[,2*1:ny-1])
                }
            else
            {
                if(!is.matrix(Yvalues)) stop("Yvalues should be a matrix")
                    new.transf <- TRUE
                    na.fail(Yvalues)
                    Yvalues <- apply(Yvalues,2,sort)
                    if(is.null(dim(Yvalues)))
                        {
                            Yvalues <- matrix(Yvalues,ncol=ny)
                        }

                    ## en ordinal, remplacer Yvalues par les modalites
                    nsim <- nrow(Yvalues)
                    nmax <- max(2*sapply(x$modalites, length), nsim)
                    Yvalues2 <- matrix(NA, nrow=nmax, ncol=ny)
                    for(yk in 1:ny)
                    {
                        if(x$linktype[yk] == 3)
                        {
                            Yvalues2[1:(2*length(x$modalites[[yk]])), yk] <- rep(x$modalites[[yk]], each=2)
                            if(2*length(x$modalites[[yk]]) < nmax) Yvalues2[(2*length(x$modalites[[yk]])+1):nmax, yk] <- Yvalues2[2*length(x$modalites[[yk]]), yk]
                        }
                        else
                        {
                            Yvalues2[1:nsim,yk] <- Yvalues[1:nsim,yk]
                            if(nsim < nmax) Yvalues2[(nsim+1):nmax,yk] <- Yvalues2[nsim,yk]
                        }
                    }
                    Yvalues <- Yvalues2
                    
                    ##controler si minY<Yvalues<maxY
                    for(yk in 1:ny)
                        {
                            if(any(Yvalues[,yk]<x$estimlink[1,2*yk-1]) | any(Yvalues[,yk]>x$estimlink[nrow(x$estimlink),2*yk-1])) stop("The values specified in \"Yvalues\" are not in the range of the outcome")
                        }
                    
                    Yvalues <- as.vector(Yvalues)
                }
            
            
            ##preparation des arguments pour Fortran
            npm <- length(x$best)
            best <- x$best
            if(x$N[4]>0)
                {
                    if(x$idiag==0) best[x$N[3]+1:x$N[4]] <- x$cholesky[-1]
                    else best[x$N[3]+1:x$N[4]] <- x$cholesky[c((1:(x$N[4]+1)*2:(x$N[4]+2))/2)[-1]]
                }

            minY <- x$estimlink[1,2*1:ny-1]
            maxY <- x$estimlink[nrow(x$estimlink),2*1:ny-1]
            
            ntrtot <- rep(NA,x$N[8])
            numSPL <- 0
            dimide <- rep(1,ny)
            for (yk in 1:x$N[8])
                {
                    if (x$linktype[yk]==0)
                        {
                            ntrtot[yk] <- 2
                        }
                    if (x$linktype[yk]==1)
                        {
                            ntrtot[yk] <- 4
                        }
                    if (x$linktype[yk]==2)
                        {
                            numSPL <-  numSPL+1
                            ntrtot[yk] <- x$nbnodes[numSPL]+2
                        }
                    if (x$linktype[yk]==3)
                        {
                            ntrtot[yk] <- x$nbmod[yk]-1
                            minY[yk] <- 1
                            maxY[yk] <- x$nbmod[yk]
                            dimide[yk] <- x$nbmod[yk]-1
                        }
                }
            
            imoins <- sum(x$N[3:8])
            zitr <- x$linknodes
            nbzitr <- rep(2,ny)
            nbzitr[which(x$linktype==2)] <- x$nbnodes
            maxnbzitr <- max(nbzitr)
            epsY <- x$epsY
            
            nsim <- length(Yvalues)/ny
            
            ide <- matrix(1,nrow=max(dimide),ncol=ny) 
            
            ndraws <- as.integer(ndraws)

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
            
            
            ## calcul des valeurs trasnformees si necessaire
            if(isTRUE(new.transf) & ndraws==0)
                {
                    resFortran <- rep(0,nsim*ny)
                    
                    out0 <- .Fortran(C_calculustransfo,
                                     as.double(best),
                                     as.integer(npm),
                                     as.integer(ny),
                                     as.integer(x$linktype),
                                     as.integer(ntrtot),
                                     as.integer(imoins),
                                     as.double(zitr),
                                     as.integer(maxnbzitr),
                                     as.double(Yvalues),
                                     as.integer(nsim),
                                     as.double(minY),
                                     as.double(maxY),
                                     as.double(epsY),
                                     as.integer(ide),
                                     as.integer(dimide),
                                     transfo=as.double(resFortran))
                    
                    transfY <- out0$transfo
                }
            else
                {
                    transfY <- as.vector(x$estimlink[,2*1:ny])
                }

            if(ndraws>0)
                {
                    if(x$conv==1)
                        {
                            ## boucle pour b=1,...B :
                            Hydraws <- NULL
                            for (j in 1:ndraws)
                                {
                                    bdraw <- rnorm(npm)
                                    bdraw <- best + Chol %*% bdraw
                                    
                                    resFortran <- rep(0,nsim*ny)
                                    
                                    out <- .Fortran(C_calculustransfo,
                                                    as.double(bdraw),
                                                    as.integer(npm),
                                                    as.integer(ny),
                                                    as.integer(x$linktype),
                                                    as.integer(ntrtot),
                                                    as.integer(imoins),
                                                    as.double(zitr),
                                                    as.integer(maxnbzitr),
                                                    as.double(Yvalues),
                                                    as.integer(nsim),
                                                    as.double(minY),
                                                    as.double(maxY),
                                                    as.double(epsY),
                                                    as.integer(ide),
                                                    as.integer(dimide),
                                                    transfo=as.double(resFortran))
                                    
                                    Hydraws <- cbind(Hydraws,out$transfo)
                                }
                            
                            ## calcul des bornes IC
                            f <- function(x)
                                {
                                    quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
                                }
                            
                            Hydistr <- apply(Hydraws,1,FUN=f)
                            mediane <- as.vector(Hydistr[2,])
                            borne_inf <- as.vector(Hydistr[1,])
                            borne_sup <- as.vector(Hydistr[3,])
                            
                            ## resultat a renvoyer
                            Yname <- rep(x$Ynames,each=nsim)
                            res <- data.frame(Yname=Yname,Yvalues=Yvalues,transfY_50=mediane,transfY_2.5=borne_inf,transfY_97.5=borne_sup)    
                        }
                    
                    if(x$conv==2 | x$conv==3)
                        {
                            ## resultat a renvoyer
                            Yname <- rep(x$Ynames,each=nsim)
                            borne_inf <- rep(NA,length(Yname))
                            borne_sup <- rep(NA,length(Yname))
                            mediane <- rep(NA,length(Yname))
                            res <- data.frame(Yname=Yname,Yvalues=Yvalues,transfY_50=mediane,transfY_2.5=borne_inf,transfY_97.5=borne_sup)   
                        }
                }
            else
                {
                    Yname <- rep(x$Ynames,each=nsim)
                    res <- data.frame(Yname=Yname,Yvalues=Yvalues,transY=transfY)
                }
            
            
        }
    else
        {
            cat("Output can not be produced since the program stopped abnormally.")
            res <- NA 
        }

    
    res.list <- NULL
    res.list$pred <- res
    res.list$object <- x
    class(res.list) <- "predictlink"
    return(res.list)
}




#' Confidence intervals for the estimated link functions from \code{lcmm},
#' \code{Jointlcmm} and \code{multlcmm}
#' 
#' This function provides 95\% confidence intervals around the estimated
#' transformation given in estimlink attribute of \code{lcmm}, \code{Jointlcmm}
#' and \code{multlcmm} objects. It can also be used to evaluate the link
#' functions at other values than those given in attribute \code{estimlink} of
#' \code{lcmm}, \code{Jointlcmm} or \code{multlcmm} object.
#' 
#' 
#' @aliases predictlink.lcmm predictlink.multlcmm predictlink.Jointlcmm
#' predictlink
#' @param x an object inheriting from classes \code{lcmm}, \code{Jointlcmm} or
#' \code{multlcmm}.
#' @param ndraws the number of draws that should be generated to approximate
#' the posterior distribution of the transformed values. By default,
#' ndraws=2000.
#' @param Yvalues a vector (for a \code{lcmm} or \code{Jointlcmm} object) or a
#' matrix (for a \code{multlcmm} object) containing the values at which to
#' compute the transformation(s). Default to the values in \code{x$estimlink}.
#' @param \dots other parameters (ignored)
#' @return An object of class \code{predictlink} with values :
#' 
#' - \code{pred} :
#' 
#' For a \code{lcmm} or \code{Jointlcmm} object, a data frame containing the
#' values at which the transformation is evaluated, the transformed values and
#' the lower and the upper limits of the confidence intervals (if ndraws>0).
#' 
#' For a \code{multlcmm} object, a data frame containing the indicator of the
#' outcome, the values at which the transformations are evaluated,the
#' transformed values and the lower and the upper limits of the confidence
#' intervals (if ndraws>0).
#' 
#' - \code{object} : the object from which the link function is predicted
#' @author Cecile Proust-Lima and Viviane Philipps
#' @seealso \code{\link{lcmm}}, \code{\link{multlcmm}},
#' \code{\link{plot.lcmm}}, \code{\link{plot.predictlink}}
#' @examples
#' 
#'  \dontrun{
#' 
#' ## Univariate mixed model with splines link funciton
#' m14<-lcmm(Ydep2~Time+I(Time^2),random=~Time,subject='ID',ng=1,
#' data=data_lcmm,link="5-manual-splines",intnodes=c(10,20,25),
#' B=c(-0.89255, -0.09715, 0.56335, 0.21967, 0.61937, -7.90261, 0.75149, 
#' -1.22357, 1.55832, 1.75324, 1.33834, 1.0968))
#' 
#' ##Transformed values of several scores and their confidence intervals
#' transf.m14 <- predictlink(m14,ndraws=2000,Yvalues=c(0,1,7:30))
#' plot(transf.m14)
#' 
#' 
#' ## Multivariate mixed model with splines link functions
#' m1 <- multlcmm(Ydep1+Ydep2~1+Time*X2+contrast(X2),random=~1+Time,
#' subject="ID",randomY=TRUE,link=c("4-manual-splines","3-manual-splines"),
#' intnodes=c(8,12,25),data=data_lcmm,
#' B=c(-1.071, -0.192,  0.106, -0.005, -0.193,  1.012,  0.870,  0.881,
#'   0.000,  0.000, -7.520,  1.401,  1.607 , 1.908,  1.431,  1.082,
#'  -7.528,  1.135 , 1.454 , 2.328, 1.052))
#' ##Confidence intervals for the transformed values (given in m1$estimlink)
#' transf.m1 <- predictlink(m1,ndraws=200)
#' plot(transf.m1)
#' }
#' 
#' @export
#' 
predictlink <- function(x,ndraws,Yvalues,...) UseMethod("predictlink")
