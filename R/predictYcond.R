predictYcond <- function(x,lambda,nsim=200,draws=FALSE,ndraws=2000,...)
    {
        if(x$conv!=1 & draws==TRUE)
            {
                cat("No confidence interval will be provided since the program did not converge properly \n")
                draws <- FALSE
            }
        
        if(x$conv %in% c(1,2))
            {
                ## cas multlcmm:
                debut <- x$N[3]+x$N[4]+x$N[5]+x$N[7] # nef+nvc+nw+ncor
                nalea <- x$N[6]
                ny <- x$N[8]
                
                npm <- length(x$best)
                nbzitr <- rep(2,ny)
                nbzitr[which(x$linktype==2)] <- x$nbnodes


                lambda <- na.omit(lambda)
                maxmes <- length(lambda)
                lambdatot <- rep(lambda,ny)
                
                Ycond <- rep(0,maxmes*ny)

                if(!draws)
                    {
                        out <- .Fortran(C_predictcondmult,
                                        as.double(lambdatot),
                                        as.integer(nalea),
                                        as.integer(ny),
                                        as.integer(maxmes),
                                        as.integer(npm),
                                        as.double(x$best),
                                        as.integer(debut),
                                        as.integer(x$epsY),
                                        as.integer(x$linktype),
                                        as.integer(nbzitr),
                                        as.double(x$linknodes),
                                        as.integer(nsim),
                                        Ycond=as.double(Ycond))
                        
                        out$Ycond[out$Ycond==9999] <- NA
                        
                        Ycond <- data.frame(rep(x$Ynames,each=maxmes),rep(lambda,ny),out$Ycond)
                        colnames(Ycond) <- c("Yname","Lambda","Ypred")
                    }
                else
                    {
                        ndraws <- as.integer(ndraws)
                        ydraws <- NULL
                        
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
                        
                        best <- x$best
                        nef <- x$N[3]
                        nvc <- x$N[4]
                        if(nvc>0) 
                            {
                                if(x$idiag==0) best[nef+1:nvc] <- x$cholesky[-1]
                                else best[nef+1:nvc] <- x$cholesky[c((1:(nvc+1)*2:(nvc+2))/2)[-1]]
                            }
                        
                        
                        for (j in 1:ndraws)
                            {  
                                bdraw <- rnorm(npm)
                                bdraw <- best + Chol %*% bdraw
                                
                                out <- .Fortran(C_predictcondmult,
                                                as.double(lambdatot),
                                                as.integer(nalea),
                                                as.integer(ny),
                                                as.integer(maxmes),
                                                as.integer(npm),
                                                as.double(bdraw),
                                                as.integer(debut),
                                                as.integer(x$epsY),
                                                as.integer(x$linktype),
                                                as.integer(nbzitr),
                                                as.double(x$linknodes),
                                                as.integer(nsim),
                                                Ycond=as.double(Ycond))
                                
                                out$Ycond[out$Ycond==9999] <- NA
                                ydraws <- cbind(ydraws,out$Ycond)
                            }

                        f <- function(x) {
                            quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
                        }
                        ydistr <- apply(ydraws,1,FUN=f)
                        Ypred_50 <- matrix(ydistr[2,],ncol=1,byrow=FALSE)
                        Ypred_2.5 <- matrix(ydistr[1,],ncol=1,byrow=FALSE)
                        Ypred_97.5 <- matrix(ydistr[3,],ncol=1,byrow=FALSE)
                        
                        Ycond <- data.frame(rep(x$Ynames,each=maxmes),rep(lambda,ny),Ypred_50,Ypred_2.5,Ypred_97.5)
                        colnames(Ycond) <- c("Yname","Lambda","Ypred_50","Ypred_2.5","Ypred_97.5")
                    }

                res <- list(pred=Ycond,object=x)
            }
        else
            {
                res <- list(pred=NA,object=x)
            }
        
        class(res) <- "predictYcond"
        return(res)
    }
