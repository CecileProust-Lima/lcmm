summary.Jointlcmm <- function(object,...)
{
    x <- object
    if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")

    cat("Joint latent class model for quantitative outcome and competing risks", "\n")
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

    cat(paste("     Number of observations:", x$N[9]),"\n")
    cat(paste("     Number of latent classes:", x$ng), "\n")
    cat(paste("     Number of parameters:", length(x$best))," \n")
    if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")

    nbevt <- length(x$hazard[[1]])
    nprisq <- rep(NA,nbevt)
    nrisq <- rep(NA,nbevt)

    typrisq <- x$hazard[[1]]
    hazardtype <- x$hazard[[2]]
    nz <- x$hazard[[4]] 

    for(ke in 1:nbevt)
        {
            if(typrisq[ke]==1) nprisq[ke] <- nz[ke]-1
            if(typrisq[ke]==2) nprisq[ke] <- 2
            if(typrisq[ke]==3) nprisq[ke] <- nz[ke]+2

            if(hazardtype[ke]=="Common") nrisq[ke] <- nprisq[ke]
            if(hazardtype[ke]=="PH") nrisq[ke] <- nprisq[ke]+x$ng-1
            if(hazardtype[ke]=="Specific") nrisq[ke] <- nprisq[ke]*x$ng

            cat(paste("     Event ",ke,": \n",sep=""))
            cat(paste("        Number of events: ", x$N[9+ke],"\n",sep=""))
            if(x$ng>1)
                {
                    if (hazardtype[ke]=="Specific") cat("        Class-specific hazards and \n")
                    if (hazardtype[ke]=="PH") cat("        Proportional hazards over latent classes and \n")
                    if (hazardtype[ke]=="Common") cat("        Common hazards over classes and \n")
                }
            
            if (typrisq[ke]==2)
                {
                    cat("        Weibull baseline risk function \n")
                }
            if (typrisq[ke]==1)
                {
                    cat("        Piecewise constant baseline risk function with nodes \n")
                    cat("       ",x$hazard[[3]][1:nz[ke],ke]," \n")
                }
            if (typrisq[ke]==3)
                {
                    cat("        M-splines constant baseline risk function with nodes \n")
                    cat("       ",x$hazard[[3]][1:nz[ke],ke]," \n")
                }
            
            
        }

    ntrtot <- x$N[8]
    numSPL <- 0
    if(x$linktype!=-1)
        {
            cat(paste("     Link function for ",x$Names$Yname,": ",sep=""))
            if (x$linktype==0)
                {
                    cat("Linear \n")
                }
            if (x$linktype==1)
                {
                    cat("Standardised Beta CdF \n")
                }
            if (x$linktype==2) 
                {
                    cat("Quadratic I-splines with nodes ", x$linknodes ,"\n")
                }
        }

    cat(" \n")
    cat("Iteration process:", "\n")

    if(x$conv==1) cat("     Convergence criteria satisfied")
    if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
    if(x$conv==3) cat("     Convergence with restrained Hessian matrix")
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
            if(!is.na(x$scoretest[1])&(length(x$hazard[[1]])==1)){
                cat(paste("     Score test statistic for CI assumption: ", round(x$scoretest[1],3)," (p-value=",round((1-pchisq(x$scoretest[1],sum(x$idea))),4),")" ,sep=""))
            }
            if(!is.na(x$scoretest[1])&(length(x$hazard[[1]])>1)){
                cat(paste("     Score test statistic for global CI assumption: ", round(x$scoretest[1],3)," (p-value=",round((1-pchisq(x$scoretest[1],sum(x$idea))),4),")" ,sep=""),"\n")
            }
            if(!is.na(x$scoretest[1])&(length(x$hazard[[1]])>1)){
                cat("     Score test statistic for event-specific CI assumption: \n")
                for (ke in 1:length(x$hazard[[1]])){ 
                    if(!is.na(x$scoretest[1+ke])){
                        cat(paste("           event ",ke,":", round(x$scoretest[1+ke],3)," (p-value=",round((1-pchisq(x$scoretest[1+ke],sum(x$idea))),4),")" ,sep=""),"\n")
                    }
                    else{
                        cat(paste("           event ",ke,": problem in the computation", "\n"))
                    }
                }
            }
            
            cat(" \n")
            cat(" \n")

            cat("Maximum Likelihood Estimates:", "\n")
            cat(" \n")
            
            nprob <- x$N[1]
            nrisqtot <- x$N[2]
            nvarxevt <- x$N[3]
            nef <- x$N[4]
            nvc <- x$N[5]
            nw <- x$N[6]
            ncor <- x$N[7]
            ntrtot <- x$N[8]
            NPM <- length(x$best)

            #nvdepsurv <- length(x$Name$TimeDepVar.name)

            
            se <- rep(NA,NPM)
            if (x$conv==1 | x$conv==3)
                {
                    ##recuperation des indices de V
                    id <- 1:NPM
                    indice <- id*(id+1)/2
                    se <-sqrt(x$V[indice])
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

            if(nw>0) coef[nprob+nrisqtot+nvarxevt+nef+nvc+1:nw] <- abs(coef[nprob+nrisqtot+nvarxevt+nef+nvc+1:nw])
            if(ncor>0) coef[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor] <- abs(coef[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor])
            if(ntrtot==1) coef[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1] <- abs(coef[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1])

            
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

            
            if(nprob>0)
                {
                    cat("Fixed effects in the class-membership model:\n" )
                    cat("(the class of reference is the last class) \n")

                    tmp <- cbind(coefch[1:nprob],sech[1:nprob],waldch[1:nprob],pwaldch[1:nprob])
                    maxch <- apply(tmp,2,function(x) max(nchar(x)))
                    if(any(c(1:nprob) %in% posfix)) maxch[1] <- maxch[1]-1
                    dimnames(tmp) <- list(names(coef)[1:nprob],
                                          c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                            paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                            paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                            paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
                    
                    cat("\n")
                    print(tmp,quote=FALSE,na.print="")
                    cat("\n")
                }


            cat("Parameters in the proportional hazard model:\n" )

            tmp <- cbind(coefch[nprob+1:(nrisqtot+nvarxevt)],
                         sech[nprob+1:(nrisqtot+nvarxevt)],
                         waldch[nprob+1:(nrisqtot+nvarxevt)],
                         pwaldch[nprob+1:(nrisqtot+nvarxevt)])
            maxch <- apply(tmp,2,function(x) max(nchar(x)))
            if(any(c(nprob+1:(nrisqtot+nvarxevt)) %in% posfix)) maxch[1] <- maxch[1]-1
            dimnames(tmp) <- list(names(coef)[nprob+1:(nrisqtot+nvarxevt)],
                                  c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                    paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                    paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                    paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
            cat("\n")
            print(tmp,quote=FALSE,na.print="")
            cat("\n")


            

            cat("Fixed effects in the longitudinal model:\n" )
            
            if(x$linktype!=-1)
                {
                    tmp <- matrix(c(paste(c(rep(" ",max(nchar(coefch[nprob+nrisqtot+nvarxevt+1:nef]))-ifelse(any(c(nprob+nrisqtot+nvarxevt+1:nef) %in% posfix),2,1)),0),collapse=""),"","",""),nrow=1,ncol=4)
                    tTable <- matrix(c(0,NA,NA,NA),nrow=1,ncol=4)
                }
            if(x$linktype==-1)
                {
                    tmp <- NULL
                    tTable <-NULL
                }
            
            if (nef>0)
                {
                    tmp2 <- cbind(coefch[nprob+nrisqtot+nvarxevt+1:nef],
                                  sech[nprob+nrisqtot+nvarxevt+1:nef],
                                  waldch[nprob+nrisqtot+nvarxevt+1:nef],
                                  pwaldch[nprob+nrisqtot+nvarxevt+1:nef])
                    tmp <- rbind(tmp,tmp2)
                    tTable <- rbind(tTable,cbind(round(coef[nprob+nrisqtot+nvarxevt+1:nef],5),
                                                 round(se[nprob+nrisqtot+nvarxevt+1:nef],5),
                                                 round(wald[nprob+nrisqtot+nvarxevt+1:nef],3),
                                                 round(pwald[nprob+nrisqtot+nvarxevt+1:nef],5)))
                }
            
            interc <- "intercept"
            if (x$ng>1)
                {
                    interc <- paste(interc,"class1")
                }
            if(x$linktype!=-1) interc <- paste(interc,"(not estimated)")
            if(x$linktype==-1) interc <- NULL
            
            if(nef>0)
                {
                    maxch <- apply(tmp,2,function(x) max(nchar(x)))
                    if(any(c(nprob+nrisqtot+nvarxevt+1:nef) %in% posfix)) maxch[1] <- maxch[1]-1

                    dimnames(tmp) <- list(c(interc,names(coef)[nprob+nrisqtot+nvarxevt+1:nef]),
                                          c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                            paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                            paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                            paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
                }
            else
                {
                    dimnames(tmp) <- list(interc, c("coef", "Se", "Wald", "p-value"))
                }

            rownames(tTable) <- rownames(tmp)
            colnames(tTable) <-  c("coef", "Se", "Wald", "p-value")
            
            cat("\n")
            print(tmp,quote=FALSE,na.print="")
            cat("\n")

            
            if(nvc>0)
                {    
                    cat("\n")
                    cat("Variance-covariance matrix of the random-effects:\n" )
                    if(x$idiag==1)
                        {
                            Mat.cov <- diag(coef[nprob+nrisqtot+nvarxevt+nef+1:nvc])
                            Mat.cov[lower.tri(Mat.cov)] <- 0
                            Mat.cov[upper.tri(Mat.cov)] <- NA
                            if(nvc==1) Mat.cov <- matrix(coef[nprob+nrisqtot+nvarxevt+nef+1:nvc],1,1)
                        }
                    if(x$idiag==0)
                        {
                            Mat.cov<-matrix(0,ncol=sum(x$idea),nrow=sum(x$idea))
                            Mat.cov[upper.tri(Mat.cov,diag=TRUE)] <- coef[nprob+nrisqtot+nvarxevt+nef+1:nvc]
                            Mat.cov <-t(Mat.cov)
                            Mat.cov[upper.tri(Mat.cov)] <- NA
                        }
                    colnames(Mat.cov) <-x$Names$Xnames[x$idea==1]
                    rownames(Mat.cov) <-x$Names$Xnames[x$idea==1]

                    if(any(posfix %in% c(nprob+nrisqtot+nvarxevt+nef+1:nvc)))
                        {
                            Mat.cov <- apply(Mat.cov,2,format,digits=5,nsmall=5)
                            Mat.cov[upper.tri(Mat.cov)] <- ""
                            pf <- sort(intersect(c(nprob+nrisqtot+nvarxevt+nef+1:nvc),posfix))
                            p <- matrix(0,sum(x$idea),sum(x$idea))
                            if(x$idiag==FALSE) p[upper.tri(p,diag=TRUE)] <- c(nprob+nrisqtot+nvarxevt+nef+1:nvc)
                            if(x$idiag==TRUE & nvc>1) diag(p) <- c(nprob+nrisqtot+nvarxevt+nef+1:nvc)
                            if(x$idiag==TRUE & nvc==1) p <- matrix(c(nprob+nrisqtot+nvarxevt+nef+1),1,1)
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
            if(nw>=1) 
                {
                    nom <- paste("Proportional coefficient class",c(1:(x$ng-1)),sep="")
                    std <-cbind(coefch[nprob+nrisqtot+nvarxevt+nef+nvc+1:nw],
                                sech[nprob+nrisqtot+nvarxevt+nef+nvc+1:nw])
                }
            if(ncor==2)
                {
                    nom <- c(nom,"AR correlation parameter:","AR standard error:")
                    std <-rbind(std,c(coefch[nprob+nrisqtot+nvarxevt+nef+nvc+nw+1],
                                      sech[nprob+nrisqtot+nvarxevt+nef+nvc+nw+1]),
                                c(coefch[nprob+nrisqtot+nvarxevt+nef+nvc+nw+2],
                                  sech[nprob+nrisqtot+nvarxevt+nef+nvc+nw+2]))
                }
            if(ncor==1) 
                {
                    nom <- c(nom,"BM standard error:")
                    std <-rbind(std,c(coefch[nprob+nrisqtot+nvarxevt+nef+nvc+nw+1],
                                      sech[nprob+nrisqtot+nvarxevt+nef+nvc+nw+1]))
                }
            
            if (!is.null(std)) 
                {
                    rownames(std) <- nom
                    maxch <- apply(std,2,function(x) max(nchar(x)))
                    if(any(c(nprob+nrisqtot+nvarxevt+nef+nvc+1:(nw+ncor)) %in% posfix)) maxch[1] <- maxch[1]-1
                    colnames(std) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                       paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""))

                    print(std,quote=FALSE,na.print="")
                    cat("\n")
                }

            if(x$linktype==-1)
                {
                    tmp <- cbind(coefch[NPM],sech[NPM])
                    rownames(tmp) <- "Residual standard error"
                    maxch <- apply(tmp,2,function(x) max(nchar(x)))
                    if(c(NPM) %in% posfix) maxch[1] <- maxch[1]-1
                    colnames(tmp) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                       paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""))
                    
                    print(tmp,quote=FALSE,na.print="")
                    cat("\n")
                }
            else
                {
                    cat("Residual standard error (not estimated) = 1\n")
                    cat("\n")
                    
                    
                    cat("Parameters of the link function:\n" )
                    
                    tmp <- cbind(coefch[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1):NPM],
                                 sech[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1):NPM],
                                 waldch[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1):NPM],
                                 pwaldch[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1):NPM])
                    rownames(tmp) <- names(x$best[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1):NPM])
                    maxch <- apply(tmp,2,function(x) max(nchar(x)))
                    if(any(c((nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1):NPM) %in% posfix)) maxch[1] <- maxch[1]-1
                    colnames(tmp) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                       paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                       paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                       paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep=""))
                    cat("\n")
                    print(tmp,quote=FALSE,na.print="")
                    cat("\n")
                }

            if(length(posfix))
                {
                    cat(" * coefficient fixed by the user \n \n")
                }
            
            return(invisible(tTable))
        }
}
