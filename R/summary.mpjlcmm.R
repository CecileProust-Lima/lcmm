#' @export
#'
summary.mpjlcmm <- function(object,...)
{
    x <- object
    cat("Multivariate joint latent class model for quantitative outcome and competing risks", "\n")
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
    nbevt <- x$nbevt
    K <- x$K
    
    cat("Statistical Model:", "\n")
    cat(paste("     Dataset:", x$call$data),"\n")
    cat(paste("     Number of subjects:", x$ns),"\n")
    cat(paste("     Number of longitudinal models:", x$K),"\n")
    cat(paste("     Number of observations:", paste(x$N[11+1:K],collapse=" ")),"\n")
    cat(paste("     Number of latent classes:", x$ng), "\n")
    cat(paste("     Number of parameters:", length(x$best))," \n")
    if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")

    if(nbevt>0)
    {
        nprisq <- rep(NA,nbevt)
        nrisq <- rep(NA,nbevt)
        for(ke in 1:nbevt)
        {
            if(x$typrisq[ke]==1) nprisq[ke] <- x$nz[ke]-1
            if(x$typrisq[ke]==2) nprisq[ke] <- 2
            if(x$typrisq[ke]==3) nprisq[ke] <- x$nz[ke]+2

            nrisq[ke] <- x$Nprm[1+ke]
            
            cat(paste("     Event",ke,": \n"))
            cat(paste("        Number of events: ", x$N[11+K+ke],"\n",sep=""))
            if(x$ng>1)
            {
                if (x$hazardtype[ke]=="Specific") cat("        Class-specific hazards and \n")
                if (x$hazardtype[ke]=="PH") cat("        Proportional hazards over latent classes and \n")
                if (x$hazardtype[ke]=="Common") cat("        Common hazards over classes and \n")
            }
            
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

    ny <- x$ny
    Ynames <- x$Names$Yname
    if(any(x$linktype!=-1))
    {
        cat("     Link functions: ")
        for (yk in 1:sum(x$ny))
        {
            if (x$linktype[yk]==0)
            {
                if (yk>1) cat("                     ")
                cat("Linear for",Ynames[yk]," \n")
            }
            if (x$linktype[yk]==1)
            {
                if (yk>1) cat("                     ")
                cat("Standardised Beta CdF for",Ynames[yk]," \n")
            }
            if (x$linktype[yk]==2)
            {
                if (yk>1) cat("                     ")
                cat("Quadratic I-splines with nodes", x$linknodes[1:x$nbzitr[yk],yk]," for ",Ynames[yk], "\n")
            }
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
        
        cat(" \n")
        cat(" \n")
        cat("Maximum Likelihood Estimates:", "\n")
        cat(" \n")
        
        nprob <- x$N[1]
        nrisqtot <- x$N[2]
        nvarxevt <- x$N[3]
        nef <- x$Nprm[3+1:K]
        ncontr <- x$Nprm[3+K+1:K]
        nvc <- x$Nprm[3+2*K+1:K]
        nw <- x$Nprm[3+3*K+1:K]
        ncor <- x$Nprm[3+4*K+1:K]
        nerr <- x$Nprm[3+5*K+1:K]
        nalea <- x$Nprm[3+6*K+1:K]
        ntr <- x$Nprm[3+7*K+1:sum(x$ny)]
        NPM <- length(x$best)

        ## shorten names if > 20 characters
        names_best <- names(x$best)
        if(any(sapply(names_best, nchar)>20))
        {
            islong <- which(sapply(names_best, nchar)>20)
            split_names_best <- strsplit(names_best, split=":", fixed=TRUE)
            short_names_best <- lapply(split_names_best, gsub, pattern="\\(.*\\)", replacement="(...)")
            new_names <- lapply(short_names_best, paste, collapse=":")
            names_best[islong] <- unlist(new_names)[islong]
            names_best[nprob+1:nrisqtot] <- names(x$best)[nprob+1:nrisqtot]
            names(x$best) <- names_best
            
            islong <- which(sapply(x$Names$Xnames, nchar)>20)
            if(length(islong))
            {
                x$Names$Xnames[islong] <- sapply(x$Names$Xnames[islong], gsub, pattern="\\(.*\\)", replacement="(...)")
            }
        }
        

        se <- rep(NA,NPM)
        if (x$conv==1 | x$conv==3)
        {
            ##recuperation des indices de V
            id <- 1:NPM
            indice <- id*(id+1)/2
            se <- sqrt(x$V[indice])
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
        
        if(nprob>0)
        {
            cat("Fixed effects in the class-membership model:\n" )
            cat("(the class of reference is the last class) \n")

            tmp <- cbind(coefch[1:nprob],sech[1:nprob],waldch[1:nprob],pwaldch[1:nprob])
            maxch <- apply(tmp,2,maxchar)
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


        if(nbevt>0)
        {
            cat("\n")
            cat("Parameters in the proportional hazard model:\n" )
            cat("\n")

            tmp <- cbind(coefch[nprob+1:(nrisqtot+nvarxevt)],
                         sech[nprob+1:(nrisqtot+nvarxevt)],
                         waldch[nprob+1:(nrisqtot+nvarxevt)],
                         pwaldch[nprob+1:(nrisqtot+nvarxevt)])
            maxch <- apply(tmp,2,maxchar)
            if(any(c(nprob+1:(nrisqtot+nvarxevt)) %in% posfix)) maxch[1] <- maxch[1]-1
            dimnames(tmp) <- list(names(coef)[nprob+1:(nrisqtot+nvarxevt)],
                                  c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                    paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                    paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                    paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
            cat("\n")
            print(tmp,quote=FALSE,na.print="")
            cat("\n")
        }


        tTable <- vector("list",K)
        
        sumny <- 0
        sumnpm <- 0
        sumnv <- 0
        for (k in 1:K)
        {
            sumntr <- 0
            
            cat("\n")
            cat("Longitudinal model for",paste(Ynames[sumny+1:ny[k]],collapse="/"),":\n" )
            cat("\n")
            
            cat("Fixed effects in the longitudinal model:\n" )            
            if (x$ng>1)
            {
                interc <- "intercept class1"
            }
            else
            {
                interc <- "intercept"
            }
            if(x$contrainte!=0)
            {
                tmp <- matrix(c(0,NA,NA,NA),nrow=1,ncol=4)
                interc <- paste(interc, "(not estimated)")
            }
            else
            {
                interc <- NULL
                tmp <- NULL
            }
            
            if(nef[k]>0)
            {
                tmp2 <- cbind(round(coef[nprob+nrisqtot+nvarxevt+sumnpm+1:nef[k]],5),
                              round(se[nprob+nrisqtot+nvarxevt+sumnpm+1:nef[k]],5),
                              round(wald[nprob+nrisqtot+nvarxevt+sumnpm+1:nef[k]],3),
                              round(pwald[nprob+nrisqtot+nvarxevt+sumnpm+1:nef[k]],5))
                tmp <- rbind(tmp,tmp2)
                dimnames(tmp) <- list(c(interc,names(coef)[nprob+nrisqtot+nvarxevt+sumnpm+1:nef[k]]), c("coef", "Se", "Wald", "p-value"))
            }
            else
            {
                dimnames(tmp) <- list(interc, c("coef", "Se", "Wald", "p-value"))
            }

            if(ncontr[k]>0)
            {
                indice2 <- 1:NPM*(1:NPM+1)/2
                nom.contr <- x$Names$Xnames[as.logical(x$idcontr[sumnv+1:x$nv[k]])]
                for (i in 1:sum(x$idcontr[sumnv+1:x$nv[k]]))
                {
                    ##matrice de variance pour test et se du dernier coef
                    indtmp <- indice2[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+((i-1)*(ny[k]-1)+1):(i*(ny[k]-1))]
                    indtmp <- cbind(indtmp-0:(length(indtmp)-1),indtmp)
                    indV <- NULL
                    for (j in 1:dim(indtmp)[1])
                    {
                        indV <- c(indV,seq(indtmp[j,1],indtmp[j,2]))
                    }
                    Vcontr <- matrix(0,ny[k]-1,ny[k]-1)
                    Vcontr[upper.tri(Vcontr,diag=TRUE)] <- x$V[indV]
                    Vcontr <- t(Vcontr)
                    Vcontr[upper.tri(Vcontr)] <- Vcontr[lower.tri(Vcontr)]
                    
                    vect.gamma <- coef[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+((i-1)*(ny[k]-1)+1):(i*(ny[k]-1))]
                    if(any(c(nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+((i-1)*(ny[k]-1)+1):(i*(ny[k]-1))) %in% posfix))
                    {
                        wald.contr <- NA
                        p.wald.contr <- NA
                    }
                    else
                    {
                        wald.contr <- t(vect.gamma) %*% solve(Vcontr,vect.gamma)
                        p.wald.contr <- 1-pchisq(wald.contr,ny[k]-1)
                    }
                    
                    tmp2 <- cbind(round(vect.gamma,5),
                                  round(se[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+((i-1)*(ny[k]-1)+1):(i*(ny[k]-1))],5),
                                  round(wald[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+((i-1)*(ny[k]-1)+1):(i*(ny[k]-1))],3),
                                  round(pwald[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+((i-1)*(ny[k]-1)+1):(i*(ny[k]-1))],5))
                    tmp2 <- rbind(rep(NA,4),tmp2)
                    
                    if(x$conv %in% c(1,3))
                    {
                        pp <- -sum(na.omit(tmp2[,1]))/sqrt(sum(Vcontr))
                        tmp2 <- rbind(tmp2,c(round(-sum(na.omit(tmp2[,1])),5),round(sqrt(sum(Vcontr)),5),round(pp,3),round(1-pchisq(pp*pp,1),5)))
                        if(is.na(p.wald.contr)) rownames(tmp2) <- c(paste("Contrasts on ",nom.contr[i],sep=""),Ynames[sumny+1:ny[k]])
                        else
                        {
                            if(round(p.wald.contr,5)!=0) rownames(tmp2) <- c(paste("Contrasts on ",nom.contr[i]," (p=",round(p.wald.contr,5),")",sep=""),Ynames[sumny+1:ny[k]])
                            if(round(p.wald.contr,5)==0) rownames(tmp2) <- c(paste("Contrasts on ",nom.contr[i]," (p<0.00001)",sep=""),Ynames[sumny+1:ny[k]])
                        }
                    }
                    if(x$conv==2)
                    {
                        tmp2 <- rbind(tmp2,c(-sum(na.omit(tmp2[,1])),NA,NA,NA))
                        rownames(tmp2) <- c(paste("Contrasts on ",nom.contr[i],sep=""),Ynames[sumny+1:ny[k]])
                    }
                    rownames(tmp2)[nrow(tmp2)] <- paste(rownames(tmp2)[nrow(tmp2)],"**",sep="")
                    if(!is.finite(tmp2[nrow(tmp2),3])) tmp2[nrow(tmp2),2:4] <- NA
                    tmp <- rbind(tmp,tmp2)
                }
            }
            
            tTable[[k]] <- tmp

            if(nef[k]>0 & any(c(nprob+nrisqtot+nvarxevt+sumnpm+1:(nef[k]+ncontr[k])) %in% posfix))
            {      
                col1 <- rep(NA,length(tmp[,1]))
                col1[which(!is.na(tmp[,1]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmp[,1]))),nsmall=5,scientific=FALSE)
                col2 <- rep(NA,length(tmp[,2]))
                col2[which(!is.na(tmp[,2]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmp[,2]))),nsmall=5,scientific=FALSE)
                col3 <- rep(NA,length(tmp[,3]))
                col3[which(!is.na(tmp[,3]))] <- format(as.numeric(sprintf("%.3f",na.omit(tmp[,3]))),nsmall=3,scientific=FALSE)
                col4 <- rep(NA,length(tmp[,4]))
                col4[which(!is.na(tmp[,4]))] <- format(as.numeric(sprintf("%.5f",na.omit(tmp[,4]))),nsmall=5,scientific=FALSE)

                pf <- sort(intersect(c(nprob+nrisqtot+nvarxevt+sumnpm+1:(nef[k]+ncontr[k])),posfix))
                p <- rep(0,length(tmp[,1]))
                a0 <- 1:nef[k]
                if(x$contrainte!=0) a0 <- c(0,1:nef[k])
                a1 <- rep(c(NA,1:(ny[k]-1),NA),sum(x$idcontr))
                a2 <- rep(nef[k]+cumsum(c(0:(sum(x$idcontr)-1))),each=ny[k]+1)
                a <- c(a0,a1+a2)
                p[which(a>0)] <- c(nprob+nrisqtot+nvarxevt+sumnpm+1:(nef[k]+ncontr[k]))
                #p[which(rownames(tmp) %in% c(x$Names$Xnames,Ynames[sumny+1:(ny[k]-1)]))] <- c(nprob+nrisqtot+nvarxevt+sumnpm+1:(nef[k]+ncontr[k]))
                col1[which(p %in% pf)] <- paste(col1[which(p %in% pf)],"*",sep="")
                col2[which(p %in% pf)] <- NA
                col3[which(p %in% pf)] <- NA
                col4[which(p %in% pf)] <- NA

                tmp <- cbind(col1,col2,col3,col4)
                rownames(tmp) <- rownames(tTable[[k]])
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
            
            if(sum(x$idea[sumnv+1:x$nv[k]])>0) cat("Variance-covariance matrix of the random-effects:\n" )
            if(x$contrainte==2)
            {
                cat("(the variance of the first random effect is not estimated)\n")
                if(nvc[k]==0)
                {
                    Mat.cov <- matrix(1,nrow=1,ncol=1)
                    colnames(Mat.cov) <- x$Names$Xnames[which(x$idea[sumnv+1:x$nv[k]]==1)]
                    rownames(Mat.cov) <- x$Names$Xnames[which(x$idea[sumnv+1:x$nv[k]]==1)]
                    prmatrix(Mat.cov)
                }
            }
            if(nvc[k]>0)
            {
                if(x$idiag[k]==1)
                {
                    if(x$contrainte==2)
                    {
                        Mat.cov <- diag(c(1,coef[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]]))
                    }
                    else
                    {
                        Mat.cov <- diag(coef[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]],nrow=nvc[k],ncol=nvc[k])
                    }
                    Mat.cov[lower.tri(Mat.cov)] <- 0
                    Mat.cov[upper.tri(Mat.cov)] <- NA
                }
                if(x$idiag[k]==0)
                {
                    Mat.cov<-matrix(0,ncol=sum(x$idea[sumnv+1:x$nv[k]]),nrow=sum(x$idea[sumnv+1:x$nv[k]]))
                    if(x$contrainte==2)
                    {
                        Mat.cov[upper.tri(Mat.cov,diag=TRUE)] <- c(1,coef[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]])
                    }
                    else
                    {
                        Mat.cov[upper.tri(Mat.cov,diag=TRUE)] <- coef[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]]
                    }
                    Mat.cov <-t(Mat.cov)
                    Mat.cov[upper.tri(Mat.cov)] <- NA
                }
                colnames(Mat.cov) <- x$Names$Xnames[which(x$idea[sumnv+1:x$nv[k]]==1)]
                rownames(Mat.cov) <- x$Names$Xnames[which(x$idea[sumnv+1:x$nv[k]]==1)]
            
            
                if(any(posfix %in% c(nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k])))
                {
                    Mat.cov <- apply(Mat.cov,2,format,digits=5,nsmall=5)
                    Mat.cov[upper.tri(Mat.cov)] <- ""
                    pf <- sort(intersect(c(nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]),posfix))
                    p <- matrix(0,sum(x$idea[sumnv+1:x$nv[k]]),sum(x$idea[sumnv+1:x$nv[k]]))
                    if(x$idiag[k]==FALSE)
                    {
                        if(x$contrainte==2)
                        {
                            p[upper.tri(p,diag=TRUE)] <- c(0,nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k])
                        }
                        else
                        {
                            p[upper.tri(p,diag=TRUE)] <- nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]
                        }
                    }
                    if(x$idiag[k]==TRUE)
                    {
                        if(x$contrainte==2)
                        {
                            diag(p) <- c(0,nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k])
                        }
                        else
                        {
                            diag(p) <- nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]
                        }
                    }
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
            if(nw[k]>=1) 
            {
                nom <- paste("Proportional coefficient class",c(1:(x$ng-1)),sep="")
                std <- cbind(coefch[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+1:nw[k]],
                             sech[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+1:nw[k]])
            }
            if(ncor[k]==2)
            {
                nom <- c(nom,"AR correlation parameter:","AR standard error:")
                std <- rbind(std,c(coefch[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+1],
                                   sech[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+1]),
                             c(coefch[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+2],
                               sech[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+2]))
            }
            if(ncor[k]==1) 
            {
                nom <- c(nom,"BM standard error:")
                std <- rbind(std,c(coefch[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+1],
                                   sech[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+1]))
            }
            if (!is.null(std)) 
            {
                rownames(std) <- nom
                maxch <- apply(std,2,maxchar)
                if(any(c(nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+1:(nw[k]+ncor[k])) %in% posfix)) maxch[1] <- maxch[1]-1
                colnames(std) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                   paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""))

                print(std,quote=FALSE,na.print="")
                cat("\n")
            }

            if(x$contrainte==1)
            {
                cat("Residual standard error: 1 (not estimated)\n")
            }
            else
            {
                std.err <- NULL
                nom <- NULL
                std.err <- rbind(std.err,coefch[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+1:nerr[k]])
                nom <- c(nom, "Residual standard error:")
                if(nalea[k]>0)
                {
                    std.err <- rbind(std.err,coefch[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+1:nalea[k]])
                    nom <- c(nom, "Standard error of the random effect:")
                }  

                rownames(std.err) <- nom
                maxch <- apply(std.err,2,maxchar)
                if(any(c(nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+1:(nerr[k]+nalea[k])) %in% posfix))
                {
                    if(nalea[k]>0)
                    {
                        maxch[union(grep("*",std.err[1,]),grep("*",std.err[2,]))] <- maxch[union(grep("*",std.err[1,]),grep("*",std.err[2,]))]-1
                    }
                    else
                    {
                        maxch[grep("*",std.err[1,])] <- maxch[grep("*",std.err[1,])]-1
                    }
                }
                colnames(std.err) <- sapply(1:ny[k],function(k) paste(paste(rep(" ",max(0,maxch[k]-maxchar(Ynames[sumny+k]))),collapse=""),Ynames[sumny+k],sep=""))
                print(std.err,quote=FALSE,na.print="")
            }
            cat("\n")

            if(any(ntr[sumny+1:ny[k]]>0))
            {
                cat("Parameters of the link functions:\n" )
                
                tmp <- cbind(coefch[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+1:sum(ntr[sumny+1:ny[k]])],
                             sech[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+1:sum(ntr[sumny+1:ny[k]])],
                             waldch[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+1:sum(ntr[sumny+1:ny[k]])],
                             pwaldch[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+1:sum(ntr[sumny+1:ny[k]])])
                tmp.rownames <- NULL
                for (yk in 1:ny[k])
                {
                    tmp.rownames <- c(tmp.rownames, paste(rep(Ynames[sumny+yk],ntr[sumny+yk]), names(coef[(nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+sum(ntr[sumny+1:yk])-ntr[sumny+yk]+1):(nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+sum(ntr[sumny+1:yk]))]),sep="-"))
                }
                rownames(tmp) <- tmp.rownames
                maxch <- apply(tmp,2,maxchar)
                if(any(c(nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+1:sum(ntr[sumny+1:ny[k]])) %in% posfix)) maxch[1] <- maxch[1]-1
                colnames(tmp) <- c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                                   paste(paste(rep(" ",max(maxch[2]-2,0)),collapse=""),"Se",sep=""),
                                   paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                                   paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep=""))
                cat("\n")
                print(tmp,quote=FALSE,na.print="")
                cat("\n") 
            }
            
            sumnpm <- sumnpm + nef[k] + ncontr[k] + nvc[k] + nw[k] + ncor[k] + nerr[k] + nalea[k] + sum(ntr[sumny+1:ny[k]])
            sumny <- sumny + x$ny[k]
            sumnv <- sumnv + x$nv[k]
        }

        if(length(posfix))
        {
            cat(" *  coefficient fixed by the user \n \n")
        }
        if(any(ncontr>0))
        {
            cat(" ** coefficient not estimated but obtained from the others as minus the sum of them \n \n")
        }
        
        return(invisible(tTable))
    }
}


