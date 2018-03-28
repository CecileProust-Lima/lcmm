

predictlink.multlcmm <- function(x,ndraws=2000,Yvalues,...)
{
    ## verification des arguments
    if(missing(x)) stop("The model should be specified.")
    if(!(inherits(x,"multlcmm"))) stop("To use only with \"multlcmm\" objects")
    #if(!missing(Yvalues) & x$linktype==3) warning("With thresholds links, no \"Yvalues\" should be specified. Default values will be used. \n")
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
                    new.transf <- TRUE
                    na.fail(Yvalues)
                    Yvalues <- apply(Yvalues,2,sort)
                    if(is.null(dim(Yvalues)))
                        {
                            Yvalues <- matrix(Yvalues,ncol=ny)
                        }
                    
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
            
            ntrtot <- rep(NA,x$N[8])
            numSPL <- 0
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
                }
            
            imoins <- sum(x$N[3:8])
            zitr <- x$linknodes
            nbzitr <- rep(2,ny)
            nbzitr[which(x$linktype==2)] <- x$nbnodes
            maxnbzitr <- max(nbzitr)
            epsY <- x$epsY
            minY <- x$estimlink[1,2*1:ny-1]
            maxY <- x$estimlink[nrow(x$estimlink),2*1:ny-1]
            
            nsim <- length(Yvalues)/ny
            
            ide <- matrix(0,nrow=1,ncol=ny) #pas encore de threshold link
            dimide <- rep(1,ny)
            
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


predictlink <- function(x,ndraws,Yvalues,...) UseMethod("predictlink")
