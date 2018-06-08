#' @export
predictlink.Jointlcmm <- function(x,ndraws=2000,Yvalues,...)
{
    ## verification des arguments
    if(missing(x)) stop("The model should be specified.")
    if(!(inherits(x,"Jointlcmm"))) stop("To use only with \"Jointlcmm\" objects")
    if(x$linktype==-1) stop("The model does not define any link function.")
    if(x$conv!=1 & ndraws!=0) stop("No confidence intervals can be produced since the program did not converge properly")
    
    if(x$conv %in% c(1,2,3)) 
        {
            if(missing(Yvalues))
                {
                    new.transf <- FALSE
                    Yvalues <- x$estimlink[,1]
                }
            else
                {
                    new.transf <- TRUE
                    Yvalues <- na.omit(Yvalues)
                    if(any(Yvalues<x$estimlink[1,1]) | any(Yvalues>x$estimlink[nrow(x$estimlink),1])) stop("The values specified in \"Yvalues\" are not in the range of the outcome")
                    Yvalues <- sort(Yvalues)
                }
            
            
            ##preparation des arguments pour Fortran
            npm <- length(x$best)
            best <- x$best
            if(x$idiag==0 & x$N[5]>0) best[sum(x$N[1:4])+1:x$N[5]] <- x$cholesky
            if(x$idiag==1 & x$N[5]>0) best[sum(x$N[1:4])+1:x$N[5]] <- sqrt(best[sum(x$N[1:4])+1:x$N[5]])
            
            if(x$linktype==0) ntrtot <- 2
            if(x$linktype==1) ntrtot <- 4
            if(x$linktype==2) ntrtot <- length(x$linknodes)+2
            
            imoins <- sum(x$N[1:7])
            zitr <- x$linknodes
            maxnbzitr <- ifelse(x$linktype==2,length(x$linknodes),2)
            epsY <- x$epsY
            minY <- x$estimlink[1,1]
            maxY <- x$estimlink[nrow(x$estimlink),1]
            ny <- 1
            nsim <- length(Yvalues)
            
            if(x$linktype==3)
                {
                    ide <- x$ide
                    dimide <- length(ide)
                }
            else
                {
                    ide <- rep(0,1) #pas utilise si pas ordinal
                    dimide <- 1
                }
            
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
                    resFortran <- rep(0,nsim)
                    
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
                    transfY <- x$estimlink[,2]
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
                                    
                                    resFortran <- rep(0,nsim)
                                    
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
                                    quantile(x[!is.na(x)],probs=c(0.5,0.025,0.975))
                                }
                            
                            Hydistr <- apply(Hydraws,1,FUN=f)
                            borne_inf <- as.vector(Hydistr[2,])
                            borne_sup <- as.vector(Hydistr[3,])
                            mediane <- as.vector(Hydistr[1,])
                            
                            ## resultat a renvoyer
                            res <- data.frame(Yvalues=Yvalues,transfY_50=mediane,transfY_2.5=borne_inf,transfY_97.5=borne_sup)
                        }
                    
                    if(x$conv==2 | x$conv==3)
                        {
                            borne_inf <- rep(NA,length(Yvalues))
                            borne_sup <- rep(NA,length(Yvalues))
                            mediane <- rep(NA,length(Yvalues))
                            res <- data.frame(Yvalues=Yvalues,transfY_50=mediane,transfY_2.5=borne_inf,transfY_97.5=borne_sup)
                        }
                }
            else
                {
                    res <- data.frame(Yvalues=Yvalues,transY=transfY)
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



