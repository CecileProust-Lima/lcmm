summary.hlme <- function(object,...){
    x <- object
    if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")

    cat("Heterogenous linear mixed model", "\n")
    cat("     fitted by maximum likelihood method", "\n")

    cl <- x$call
    cl$B <- NULL
    if(is.data.frame(cl$data))
        {
            cl$data <- NULL
            x$call$data <- NULL    
        }
    cat(" \n")
    dput(cl)
    cat(" \n")

    posfix <- eval(cl$posfix)


    cat("Statistical Model:", "\n")
    cat(paste("     Dataset:", as.character(as.expression(x$call$data))),"\n")
    cat(paste("     Number of subjects:", x$ns),"\n")
    cat(paste("     Number of observations:", length(x$pred[,1])),"\n")
    if(length(x$na.action))cat(paste("     Number of observations deleted:",length(x$na.action)),"\n")
    cat(paste("     Number of latent classes:", x$ng), "\n")
    cat(paste("     Number of parameters:", length(x$best))," \n")
    if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")
    cat(" \n")
    cat("Iteration process:", "\n")

    if(x$conv==1) cat("     Convergence criteria satisfied")
    if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
    if(x$conv==4|x$conv==12)
        {
            cat("     The program stopped abnormally. No results can be displayed.\n")
        }
    else
        {

            cat(" \n")
            cat("     Number of iterations: ", x$niter, "\n")
            cat("     Convergence criteria: parameters=", signif(x$gconv[1],2), "\n")
            cat("                         : likelihood=", signif(x$gconv[2],2), "\n") 
            cat("                         : second derivatives=", signif(x$gconv[3],2), "\n")
            cat(" \n")
            cat("Goodness-of-fit statistics:", "\n")
            cat(paste("     maximum log-likelihood:", round(x$loglik,2))," \n")
            cat(paste("     AIC:", round(x$AIC,2))," \n")
            cat(paste("     BIC:", round(x$BIC,2))," \n")
            cat(" \n")
            cat(" \n")

            cat("Maximum Likelihood Estimates:", "\n")
            cat(" \n")

            NPROB <- x$N[1]
            NEF   <- x$N[2]
            NVC   <- x$N[3]
            NW    <- x$N[4]
            ncor <- x$N[5]
            NPM   <- length(x$best)


            se <- rep(1,NPM)
            if (x$conv==1)
                {
                    ##recuperation des indices de V
                    id <- 1:NPM
                    indice <- rep(id*(id+1)/2)
                    se <- sqrt(x$V[indice])
                    if (NVC>0) se[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- 1
                    wald <- x$best/se
                    pwald <- 1-pchisq(wald**2,1)
                    coef <- x$best
                }
            else
                {
                    se <- NA
                    wald <- NA
                    pwald <- NA
                    coef <- x$best

                    sech <- rep(NA,length(coef))
                    waldch <- rep(NA,length(coef))
                    pwaldch <- rep(NA,length(coef))
                }

            ##prendre abs pour les parametres mis au carre
            if(NW>0) coef[NPROB+NEF+NVC+1:NW] <- abs(coef[NPROB+NEF+NVC+1:NW])
            if(ncor>0) coef[NPROB+NEF+NVC+NW+ncor] <- abs(coef[NPROB+NEF+NVC+NW+ncor])
            coef[NPROB+NEF+NVC+NW+ncor+1] <- abs(coef[NPROB+NEF+NVC+NW+ncor+1])


            ## convertir en character
            if(x$conv!=2)
                {
                    coefch <- format(as.numeric(sprintf("%.5f",coef)),nsmall=5,scientific=FALSE)
                    sech <- format(as.numeric(sprintf("%.5f",se)),nsmall=5,scientific=FALSE)
                    waldch <- format(as.numeric(sprintf("%.3f",wald)),nsmall=3,scientific=FALSE)
                    pwaldch <- format(as.numeric(sprintf("%.5f",pwald)),nsmall=5,scientific=FALSE)
                }
            else
                {
                    coefch <- format(as.numeric(sprintf("%.5f",coef)),nsmall=5,scientific=FALSE)
                }

            if(length(posfix))
                {
                    coefch[posfix] <- paste(coefch[posfix],"*",sep="")
                    sech[posfix] <- ""
                    waldch[posfix] <- ""
                    pwaldch[posfix] <- ""
                }
            



            if(NPROB>0)
                {
                    cat("Fixed effects in the class-membership model:\n" )
                    cat("(the class of reference is the last class) \n")

                    tmp <- cbind(coefch[1:NPROB],sech[1:NPROB],waldch[1:NPROB],pwaldch[1:NPROB])
                    maxch <- apply(tmp,2,function(x) max(nchar(x)))
                    if(any(c(1:NPROB) %in% posfix)) maxch[1] <- maxch[1]-1
                    dimnames(tmp) <- list(names(coef)[1:NPROB],
                                          c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                            paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                            paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                            paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
                    
                    cat("\n")
                    print(tmp,quote=FALSE,na.print="")
                    cat("\n")
                }


            cat("Fixed effects in the longitudinal model:\n" )

            tmp <- cbind(coefch[(NPROB+1):(NPROB+NEF)],
                         sech[(NPROB+1):(NPROB+NEF)],
                         waldch[(NPROB+1):(NPROB+NEF)],
                         pwaldch[(NPROB+1):(NPROB+NEF)])
            
            maxch <- apply(tmp,2,function(x) max(nchar(x)))
            if(any(c(NPROB+1:NEF) %in% posfix)) maxch[1] <- maxch[1]-1
            dimnames(tmp) <- list(names(coef)[NPROB+1:NEF],
                                          c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                            paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                            paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                            paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
            
            cat("\n")
            print(tmp,quote=FALSE,na.print="")
            cat("\n")

            tTable <- cbind(round(coef[(NPROB+1):(NPROB+NEF)],5),
                            round(se[(NPROB+1):(NPROB+NEF)],5),
                            round(wald[(NPROB+1):(NPROB+NEF)],3),
                            round(pwald[(NPROB+1):(NPROB+NEF)],5))
            dimnames(tTable) <- list(names(coef)[NPROB+1:NEF], c("coef", "Se", "Wald", "p-value"))

            if(NVC>0)
                {
                    cat("\n")
                    cat("Variance-covariance matrix of the random-effects:\n" )
                    if(x$idiag==1)
                        {
                            if (NVC>1)
                                {
                                    Mat.cov <- diag(coef[(NPROB+NEF+1):(NPROB+NEF+NVC)])
                                }
                            else
                                {
                                    Mat.cov <- matrix(coef[(NPROB+NEF+1)],ncol=1)
                                }
                            colnames(Mat.cov) <- x$Xnames[x$idea0==1]
                            rownames(Mat.cov) <- x$Xnames[x$idea0==1]
                            Mat.cov[lower.tri(Mat.cov)] <- 0
                            Mat.cov[upper.tri(Mat.cov)] <- NA

                            #print(Mat.cov,na.print="")
                            #cat("\n")
                        }


                    if(x$idiag==0)
                        {
                            Mat.cov <- matrix(0,ncol=sum(x$idea0),nrow=sum(x$idea0))
                            colnames(Mat.cov) <- x$Xnames[ x$idea0==1]
                            rownames(Mat.cov) <- x$Xnames[ x$idea0==1]
                            Mat.cov[upper.tri(Mat.cov,diag=TRUE)] <- coef[(NPROB+NEF+1):(NPROB+NEF+NVC)]
                            Mat.cov <- t(Mat.cov)
                            Mat.cov[upper.tri(Mat.cov)] <- NA

                            #print(Mat.cov,na.print="")
                            #cat("\n")
                        }

                    if(any(posfix %in% c(NPROB+NEF+1:NVC)))
                        {
                            Mat.cov <- apply(Mat.cov,2,format,digits=5,nsmall=5)
                            Mat.cov[upper.tri(Mat.cov)] <- ""
                            pf <- sort(intersect(c(NPROB+NEF+1:NVC),posfix))
                            p <- matrix(0,sum(x$idea0),sum(x$idea0))
                            if(x$idiag==FALSE) p[upper.tri(p,diag=TRUE)] <- c(NPROB+NEF+1:NVC)
                            if(x$idiag==TRUE & NVC>1) diag(p) <- c(NPROB+NEF+1:NVC)
                            if(x$idiag==TRUE & NVC==1) p <- matrix(c(NPROB+NEF+1),1,1)
                            Mat.cov[which(t(p) %in% pf)] <- paste(Mat.cov[which(t(p) %in% pf)],"*",sep="")
                            print(Mat.cov,quote=FALSE)
                        }
                    else
                        {
                            prmatrix(round(Mat.cov,5),na.print="")
                        }
                    cat("\n")
                }
            
            std <- NULL
            nom <- NULL
            if((NW>=1) & (x$ng>1))
                {
                    nom <- paste("Proportional coefficient class",c(1:(x$ng-1)),sep="")
                    std <-cbind(coefch[NPROB+NEF+NVC+1:NW],
                                sech[NPROB+NEF+NVC+1:NW])
                }
            if(ncor==2)
                {
                    nom <- c(nom,"AR correlation parameter:","AR standard error:")
                    std <-rbind(std,c(coefch[(NPROB+NEF+NVC+NW+1)],
                                      sech[(NPROB+NEF+NVC+NW+1)]),
                                c(coefch[(NPROB+NEF+NVC+NW+2)],
                                  sech[(NPROB+NEF+NVC+NW+2)]))
                }
            if(ncor==1)
                {
                    nom <- c(nom,"BM standard error:")
                    std <-rbind(std,c(coefch[(NPROB+NEF+NVC+NW+1)],sech[(NPROB+NEF+NVC+NW+1)]))
                }
            std <- rbind(std,c(coefch[NPM],sech[NPM]))
            nom <- c(nom,"Residual standard error:")
            
            rownames(std) <- nom
            maxch <- apply(std,2,function(x) max(nchar(x)))
            if(any(c(NPROB+NEF+NVC+1:(NW+ncor+1)) %in% posfix)) maxch[1] <- maxch[1]-1
            colnames(std) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                               paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""))
            
            print(std,quote=FALSE,na.print="")
            cat("\n")

            if(length(posfix))
                {
                    cat(" * coefficient fixed by the user \n \n")
                }

            return(invisible(tTable))
        }
}


