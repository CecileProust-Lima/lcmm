#' @export
#'
summary.irt <- function(object,...)
{
    x <- object

    cat("IRT model with shared random effects", "\n")
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

    cat(paste("     Number of observations:", x$N[13]),"\n")
    cat(paste("     Number of parameters:", length(x$best))," \n")
    if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")


    nbevt <- x$N[14]
    if(nbevt>0)
    {
        nprisq <- rep(NA,nbevt)
        nrisq <- rep(NA,nbevt)
        for(ke in 1:nbevt)
        {
            if(x$typrisq[ke]==1) nprisq[ke] <- x$nz[ke]-1
            if(x$typrisq[ke]==2) nprisq[ke] <- 2
            if(x$typrisq[ke]==3) nprisq[ke] <- x$nz[ke]+2

            
            cat(paste("     Event",ke,": \n"))
            cat(paste("        Number of events: ", x$nevent[ke],"\n",sep=""))
            
            if (x$typrisq[ke]==2)
            {
                cat("        Weibull baseline risk function \n")
            }
            if (x$typrisq[ke]==1)
            {
                cat("        Piecewise constant baseline risk function with nodes \n")
                cat("        ",x$hazardnodes[1:x$nz[ke],ke]," \n")
            }
            if (x$typrisq[ke]==3)
            {
                cat("        M-splines constant baseline risk function with nodes \n")
                cat("        ",x$hazardnodes[1:x$nz[ke],ke]," \n")
            }
                        
        }
    }

    
    ny <- x$N[12]
    ntr <- rep(NA,ny)
    numSPL <- 0
    cat("     Link functions: ")
    for (yk in 1:ny)
        {
            if (x$linktype[yk]==0)
                {
                    ntr[yk] <- 2
                    if (yk>1) cat("                     ")
                    cat("Linear for",x$Names$Ynames[yk]," \n")
                }
            if (x$linktype[yk]==1)
                {
                    ntr[yk] <- 4
                    if (yk>1) cat("                     ")
                    cat("Standardised Beta CdF for",x$Names$Ynames[yk]," \n")
                }
            if (x$linktype[yk]==2) 
                {
                    numSPL <- numSPL+1
                    ntr[yk] <- x$nbnodes[numSPL]+2
                    if (yk>1) cat("                     ")
                    cat("Quadratic I-splines with nodes", x$linknodes[1:x$nbnodes[numSPL],yk],"for",x$Names$Ynames[yk], "\n")
                }
            if (x$linktype[yk]==3) 
                {
                    ntr[yk] <- x$nbmod[yk]-1
                    if (yk>1) cat("                     ")
                    cat("Thresholds for",x$Names$Ynames[yk], "\n")
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
            cat(" \n")


            cat("Maximum Likelihood Estimates:", "\n")
            cat(" \n")

            nprob <- x$N[1]
            nrisqtot <- x$N[2]
            nvarxevt <- x$N[3]
            nasso <- x$N[4]
            nef <- x$N[5]
            ncontr <- x$N[6]
            nvc <- x$N[7]
            nw <- x$N[8]
            ncor <- x$N[9]
            ntrtot <- x$N[10]
            nalea <- x$N[11]
            ny <- x$N[12]
            nbevt <- x$N[13]
            NPM <- length(x$best)
            
            
            se <- rep(NA,NPM)
            if (x$conv==1 | x$conv==3)
                {
                    ##recuperation des indices de V
                    id <- 1:NPM
                    indice <- id*(id+1)/2
                    se <- sqrt(x$V[indice])
                    se[which(is.na(se))] <- 1
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

            if(ncor>0) coef[nrisqtot+nvarxevt+nef+ncontr+nvc+ncor] <- abs(coef[nrisqtot+nvarxevt+nef+ncontr+nvc+ncor])
            coef[nrisqtot+nvarxevt+nef+ncontr+nvc+ncor+ntrtot+nalea+1:ny] <- abs(coef[nrisqtot+nvarxevt+nef+ncontr+nvc+ncor+ntrtot+nalea+1:ny])
            if(nalea>0) coef[nrisqtot+nvarxevt+nef+ncontr+nvc+ncor+ntrtot+1:ny] <- abs(coef[nrisqtot+nvarxevt+nef+ncontr+nvc+ncor+ntrtot+1:ny])

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

            ## fct pr determiner la longueur max d'une chaine de caracteres
            ## (avec gestion des NA)
            maxchar <- function(x)
                {
                    xx <- na.omit(x)
                    if(length(xx))
                        {
                            res <- max(nchar(xx))
                        }
                    else
                        {
                            res <- 2
                        }
                    return(res)
                }

#browser()

            if(nbevt>0)
            {
                cat("\n")
                cat("Parameters in the proportional hazard model:\n" )
                cat("\n")
                
                tmp <- cbind(coefch[c(1:(nrisqtot+nvarxevt), nrisqtot+nvarxevt+1:nasso)],
                             sech[c(1:(nrisqtot+nvarxevt), nrisqtot+nvarxevt+1:nasso)],
                             waldch[c(1:(nrisqtot+nvarxevt), nrisqtot+nvarxevt+1:nasso)],
                             pwaldch[c(1:(nrisqtot+nvarxevt), nrisqtot+nvarxevt+1:nasso)])
                maxch <- apply(tmp,2,maxchar)
                if(any(c(1:(nrisqtot+nvarxevt), nrisqtot+nvarxevt+1:nasso) %in% posfix)) maxch[1] <- maxch[1]-1
                dimnames(tmp) <- list(names(coef)[1:(nrisqtot+nvarxevt+nasso)],
                                      c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                        paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                        paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                        paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
                cat("\n")
                print(tmp,quote=FALSE,na.print="")
                cat("\n")
            }



            cat("Fixed effects in the longitudinal model:\n" )

            tmp <- matrix(c(0,NA,NA,NA),nrow=1,ncol=4)
            
            if (nef>0)
                {
                    tmp2 <- cbind(round(coef[nrisqtot+nvarxevt+nasso+1:nef],5),round(se[nrisqtot+nvarxevt+nasso+1:nef],5),round(wald[nrisqtot+nvarxevt+nasso+1:nef],3),round(pwald[nrisqtot+nvarxevt+nasso+1:nef],5))
                    tmp <- rbind(tmp,tmp2)
                }
            interc <- "intercept (not estimated)"

            if(nef>0) dimnames(tmp) <- list(c(interc,names(coef)[nrisqtot+nvarxevt+nasso+1:nef]), c("coef", "Se", "Wald", "p-value"))
            else dimnames(tmp) <- list(interc, c("coef", "Se", "Wald", "p-value"))
            cat("\n")
            
            if(ncontr>0)
                {
                    indice2 <- 1:NPM*(1:NPM+1)/2
                    nom.contr <- x$Xnames[as.logical(x$idcontr)]
                    for (i in 1:sum(x$idcontr))
                        {
                            ##matrice de variance pour test et se du dernier coef
                            indtmp <- indice2[(nrisqtot+nvarxevt+nasso+nef+(i-1)*(ny-1)+1):(nrisqtot+nvarxevt+nasso+nef+i*(ny-1))]
                            indtmp <- cbind(indtmp-0:(length(indtmp)-1),indtmp)
                            indV <- NULL
                            for (j in 1:dim(indtmp)[1])
                                {
                                    indV <- c(indV,seq(indtmp[j,1],indtmp[j,2]))
                                }
                            Vcontr <- matrix(0,ny-1,ny-1)
                            Vcontr[upper.tri(Vcontr,diag=TRUE)] <- x$V[indV]
                            Vcontr <- t(Vcontr)
                            Vcontr[upper.tri(Vcontr)] <- Vcontr[lower.tri(Vcontr)]
                            
                            vect.gamma <- coef[(nrisqtot+nvarxevt+nasso+nef+(i-1)*(ny-1)+1):(nrisqtot+nvarxevt+nasso+nef+i*(ny-1))]
                            if(any(c((nrisqtot+nvarxevt+nasso+nef+(i-1)*(ny-1)+1):(nrisqtot+nvarxevt+nasso+nef+i*(ny-1))) %in% posfix))
                                {
                                    wald.contr <- NA
                                    p.wald.contr <- NA
                                }
                            else
                                {
                                    wald.contr <- t(vect.gamma) %*% solve(Vcontr,vect.gamma)
                                    p.wald.contr <- 1-pchisq(wald.contr,ny-1)
                                }
                            
                            tmp2 <- cbind(round(vect.gamma,5),
                                          round(se[(nrisqtot+nvarxevt+nasso+nef+(i-1)*(ny-1)+1):(nrisqtot+nvarxevt+nasso+nef+i*(ny-1))],5),
                                          round(wald[(nrisqtot+nvarxevt+nasso+nef+(i-1)*(ny-1)+1):(nrisqtot+nvarxevt+nasso+nef+i*(ny-1))],3),
                                          round(pwald[(nrisqtot+nvarxevt+nasso+nef+(i-1)*(ny-1)+1):(nrisqtot+nvarxevt+nasso+nef+i*(ny-1))],5))
                            tmp2 <- rbind(rep(NA,4),tmp2)
                            
                            if(x$conv %in% c(1,3))
                                {
                                    pp <- -sum(na.omit(tmp2[,1]))/sqrt(sum(Vcontr))
                                    tmp2 <- rbind(tmp2,c(round(-sum(na.omit(tmp2[,1])),5),round(sqrt(sum(Vcontr)),5),round(pp,3),round(1-pchisq(pp*pp,1),5)))
                                    if(is.na(p.wald.contr)) rownames(tmp2) <- c(paste("Contrasts on ",nom.contr[i],sep=""),x$Ynames)
                                    else
                                        {
                                            if(round(p.wald.contr,5)!=0) rownames(tmp2) <- c(paste("Contrasts on ",nom.contr[i]," (p=",round(p.wald.contr,5),")",sep=""),x$Ynames)
                                            if(round(p.wald.contr,5)==0) rownames(tmp2) <- c(paste("Contrasts on ",nom.contr[i]," (p<0.00001)",sep=""),x$Ynames)
                                        }
                                }
                            if(x$conv==2)
                                {
                                    tmp2 <- rbind(tmp2,c(-sum(na.omit(tmp2[,1])),NA,NA,NA))
                                    rownames(tmp2) <- c(paste("Contrasts on ",nom.contr[i],sep=""),x$Ynames)
                                }
                            rownames(tmp2)[nrow(tmp2)] <- paste(rownames(tmp2)[nrow(tmp2)],"**",sep="")
                            if(!is.finite(tmp2[nrow(tmp2),3])) tmp2[nrow(tmp2),2:4] <- NA
                            tmp <- rbind(tmp,tmp2)
                        }
                }
            
            tTable <- tmp

            if((nef>0) & any(c(nrisqtot+nvarxevt+nasso+1:nef) %in% posfix))
                {      
                    col1 <- rep(NA,length(tmp[,1]))
                    col1[which(!is.na(tmp[,1]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmp[,1]))),nsmall=5,scientific=FALSE)
                    col2 <- rep(NA,length(tmp[,2]))
                    col2[which(!is.na(tmp[,2]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmp[,2]))),nsmall=5,scientific=FALSE)
                    col3 <- rep(NA,length(tmp[,3]))
                    col3[which(!is.na(tmp[,3]))] <- format(as.numeric(sprintf("%.3f",na.omit(tmp[,3]))),nsmall=3,scientific=FALSE)
                    col4 <- rep(NA,length(tmp[,4]))
                    col4[which(!is.na(tmp[,4]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmp[,4]))),nsmall=5,scientific=FALSE)

                    pf <- sort(intersect(c(nrisqtot+nvarxervt+1:(nef+ncontr)),posfix))
                    p <- rep(0,length(tmp[,1]))
                    p[which(rownames(tmp) %in% c(x$Names$Xnames,x$Names$Ynames[-ny]))] <- c(nrisqtot+nvarxervt+1:(nef+ncontr))
                    col1[which(p %in% pf)] <- paste(col1[which(p %in% pf)],"*",sep="")
                    col2[which(p %in% pf)] <- NA
                    col3[which(p %in% pf)] <- NA
                    col4[which(p %in% pf)] <- NA

                    tmp <- cbind(col1,col2,col3,col4)
                    rownames(tmp) <- rownames(tTable)
                    maxch <- apply(tmp,2,maxchar)
                    maxch[1] <- maxch[1]-1

                    colnames(tmp) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                       paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                       paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                       paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep=""))
                    
                    print(tmp,quote=FALSE,na.print="")
                    cat("\n")
                }
            else
                {
                    prmatrix(round(tmp,5),na.print="")
                    cat("\n")
                }
            

            
            cat("\n")
            cat("Variance-covariance matrix of the random-effects:\n" )
            cat("(the variance of the first random effect is not estimated)\n")
            if(x$idiag==1)
                {
                    if (nvc>0) 
                        {
                            Mat.cov <- diag(c(1,coef[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc]))
                        }
                    else
                        {
                            Mat.cov <- matrix(1,ncol=1)
                        }
                    Mat.cov[lower.tri(Mat.cov)] <- 0
                    Mat.cov[upper.tri(Mat.cov)] <- NA
                }
            if(x$idiag==0)
                {
                    Mat.cov<-matrix(0,ncol=sum(x$idea),nrow=sum(x$idea))
                    if(nvc>0) 
                        {
                            Mat.cov[upper.tri(Mat.cov,diag=TRUE)]<-c(1,coef[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc])
                            Mat.cov <-t(Mat.cov)
                            Mat.cov[upper.tri(Mat.cov)] <- NA
                        }
                    else Mat.cov[1,1] <- 1
                }
            colnames(Mat.cov) <-x$Names$Xnames[x$idea==1]
            rownames(Mat.cov) <-x$Names$Xnames[x$idea==1]
            
            if(any(posfix %in% c(nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc)))
                {
                    Mat.cov <- apply(Mat.cov,2,format,digits=5,nsmall=5)
                    Mat.cov[upper.tri(Mat.cov)] <- ""
                    pf <- sort(intersect(c(nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc),posfix))
                    p <- matrix(0,sum(x$idea),sum(x$idea))
                    if(x$idiag==FALSE) p[upper.tri(p,diag=TRUE)] <- c(0,nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc)
                    if(x$idiag==TRUE) diag(p) <- c(0,nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc)
                    Mat.cov[which(t(p) %in% pf)] <- paste(Mat.cov[which(t(p) %in% pf)],"*",sep="")
                    print(Mat.cov,quote=FALSE)
                }
            else
                {
                    prmatrix(round(Mat.cov,5),na.print="")
                }
            cat("\n")
            
            std <- NULL
            nom <- NULL
            if(ncor==2)
                {
                    nom <- c(nom,"AR correlation parameter:","AR standard error:")
                    std <-rbind(std,c(coefch[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+1],sech[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+1]),c(coefch[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+2],sech[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+2]))
                }
            if(ncor==1) 
                {
                    nom <- c(nom,"BM standard error:")
                    std <-rbind(std,c(coefch[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+1],sech[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+1]))
                }
            if (!is.null(std)) 
                {
                    rownames(std) <- nom
                    maxch <- apply(std,2,maxchar)
                    if(any(c(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+1:ncor) %in% posfix)) maxch[1] <- maxch[1]-1
                    colnames(std) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                       paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""))

                    print(std,quote=FALSE,na.print="")
                    cat("\n")
                }

            
            std.err <- NULL
            nom <- NULL
            std.err <- rbind(std.err,coefch[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+1:ny])
            nom <- c(nom, "Residual standard error:")
            if(nalea>0)
                {
                    std.err <- rbind(std.err,coefch[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+1:nalea])
                    nom <- c(nom, "Standard error of the random effect:")
                }  

            rownames(std.err) <- nom
            maxch <- apply(std.err,2,maxchar)
            if(any(c(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+1:(nalea+ny)) %in% posfix))
                {
                    if(nalea>0)
                        {
                            maxch[union(grep("*",std.err[1,]),grep("*",std.err[2,]))] <- maxch[union(grep("*",std.err[1,]),grep("*",std.err[2,]))]-1
                        }
                    else
                        {
                            maxch[grep("*",std.err[1,])] <- maxch[grep("*",std.err[1,])]-1
                        }
                }
            colnames(std.err) <- sapply(1:ny,function(k) paste(paste(rep(" ",max(0,maxch[k]-maxchar(x$Ynames[k]))),collapse=""),x$Ynames[k],sep=""))
            print(std.err,quote=FALSE,na.print="")
            
            cat("\n")
            
            cat("Parameters of the link functions:\n" )
            
            tmp <- cbind(coefch[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntrtot],
                         sech[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntrtot],
                         waldch[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntrtot],
                         pwaldch[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntrtot])
            tmp.rownames <- NULL
            for (yk in 1:ny)
                {
                    tmp.rownames <- c(tmp.rownames, paste(rep(x$Names$Ynames[yk],ntr[yk]),names(coef[(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:yk])-ntr[yk]+1):(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:yk]))]),sep="-"))
                }
            rownames(tmp) <- tmp.rownames
            maxch <- apply(tmp,2,maxchar)
            if(any(c(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntrtot) %in% posfix)) maxch[1] <- maxch[1]-1
            colnames(tmp) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                               paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                               paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                               paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep=""))
            cat("\n")
            print(tmp,quote=FALSE,na.print="")
            cat("\n")

            if(length(posfix))
                {
                    cat(" *  coefficient fixed by the user \n \n")
                }
            if(ncontr>0)
                {
                    cat(" ** coefficient not estimated but obtained from the others as minus the sum of them \n \n")
                }
            
            return(invisible(tTable))
        }
}

