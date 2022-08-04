#' @name predictY
#' @rdname predictY
#' @export
#'
predictY.Jointlcmm <- function(x,newdata,var.time,methInteg=0,nsim=20,draws=FALSE,ndraws=2000,na.action=1,...)
{
    if(missing(newdata)) stop("The argument newdata should be specified")
    if(missing(x)) stop("The argument x should be specified")
    if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
    Xnames2 <- setdiff(x$Names$Xnames2,x$Names$Tnames)
    if(!all(Xnames2 %in% c(colnames(newdata),"intercept")))  stop(paste("newdata should at least include the following covariates: ",paste(Xnames2,collapse=" "), sep="\n"))
    if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
    if (!(methInteg %in% c(0,1))) stop("The integration method must be either 0 for Gauss-Hermite or 1 for Monte-Carlo")
    if ((methInteg==0)&(!(nsim %in% c(5,7,9,15,20,30,40,50)))) stop("For Gauss-Hermite integration method, 'nsim' should be either 5,7,9,15,20,30,40 or 50")
#    if(missing(var.time)) stop("missing argument 'var.time'")
#    if(!(var.time %in% colnames(newdata))) stop("'var.time' should be included in newdata")


    if(x$conv==1 | x$conv==2 | x$conv==3)
        {

            if(x$conv==2 & draws==TRUE)
                {
                    cat("No confidence interval will be provided since the program did not converge properly \n")
                    draws <- FALSE
                }
 
            if(x$conv==3 & draws==TRUE)
                {
                    cat("No confidence interval will be provided since the program did not converge properly \n")
                    draws <- FALSE
                }

            
            if (x$Names$Xnames2[1]!="intercept")
                {
                    newdata1 <- newdata[,x$Names$Xnames2,drop=FALSE]
                    colnames(newdata1) <- x$Names$Xnames2
                    newdata1 <- data.frame(newdata1)
                }
            else
                {
                    newdata1 <- cbind(rep(1,length(newdata[,1])),newdata[,x$Names$Xnames2[-1]])
                    colnames(newdata1) <- c("intercept",x$Names$Xnames2[-1])
                    newdata1 <- data.frame(newdata1)
                }


            X1 <- NULL
            X2 <- NULL
            b1 <- NULL
            b2 <- NULL

            call_fixed <- x$call$fixed[3]
            if(is.null(x$call$random)) {call_random <- -1} else call_random <- x$call$random[2]
            if(is.null(x$call$classmb)) {call_classmb <- -1} else call_classmb <- x$call$classmb[2]
            if(is.null(x$call$mixture)) {call_mixture <- -1} else call_mixture <- x$call$mixture[2]
            if(is.null(x$call$survival)) {call_survival <- -1} else call_survival <- x$call$survival[3]

            nbevt <- length(x$hazard[[1]])
            call_survival <- gsub("mixture","",call_survival)
            for(k in 1:nbevt)
                {
                    call_survival <- gsub(paste("cause",k,sep=""),"",call_survival)
                }
            call_survival <- gsub("cause","",call_survival)
            call_survival <- call(call_survival) 

            if(!(na.action %in% c(1,2))) stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

            if(na.action==1)
                {
                    na.action=na.omit
                }
            else
                {
                    na.action=na.fail
                }




### pour les facteurs
    ## ##donnees de l estimation
    ## if(!is.null(x$data))
    ## {
    ##     olddata <- x$data
    ## }
    ## else
    ## {
    ##     olddata <- eval(x$call$data)
    ## }

    ##         ##cas ou une variable du dataset est un facteur
    ##         for(v in Xnames2[-1])
    ##             {
    ##                 if (is.factor(olddata[,v]))     
    ##                     {
    ##                         mod <- levels(olddata[,v])
    ##                         if (!(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
    ##                         newdata1[,v] <- factor(newdata1[,v], levels=mod)
    ##                     }
    ##             }

            ## transform to factor is the variable appears in levels$levelsdata
            for(v in colnames(newdata1))
            {
                if(v %in% names(x$levels$levelsdata))
                {
                    if(!is.null(x$levels$levelsdata[[v]]))
                    {
                        newdata1[,v] <- factor(newdata1[,v], levels=x$levels$levelsdata[[v]])
                    }
                }
            }

            
            ##cas ou on a factor() dans l'appel
            ## z <- all.names(call_fixed)
            ## ind_factor <- which(z=="factor")
            ## if(length(ind_factor))
            ##     {
            ##         nom.factor <- z[ind_factor+1]  
            ##         for (v in nom.factor)
            ##             {
            ##                 mod <- levels(as.factor(olddata[,v]))
            ##                 if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
            ##                 newdata1[,v] <- factor(newdata1[,v], levels=mod)
            ##             }
            ##     }
            call_fixed <- gsub("factor","",call_fixed)

            ## z <- all.names(call_random)
            ## ind_factor <- which(z=="factor")
            ## if(length(ind_factor))
            ##     {
            ##         nom.factor <- z[ind_factor+1]
            ##         for (v in nom.factor)
            ##             {
            ##                 mod <- levels(as.factor(olddata[,v]))
            ##                 if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
            ##                 newdata1[,v] <- factor(newdata1[,v], levels=mod)
            ##             }
            ##     }
            call_random <- gsub("factor","",call_random)
            
            ## z <- all.names(call_classmb)
            ## ind_factor <- which(z=="factor")
            ## if(length(ind_factor))
            ##     {
            ##         nom.factor <- z[ind_factor+1]
            ##         for (v in nom.factor)
            ##             {
            ##                 mod <- levels(as.factor(olddata[,v]))
            ##                 if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
            ##                 newdata1[,v] <- factor(newdata1[,v], levels=mod)
            ##             }
            ##     }
            call_classmb <- gsub("factor","",call_classmb)
            
            ## z <- all.names(call_mixture)
            ## ind_factor <- which(z=="factor")
            ## if(length(ind_factor))
            ##     {
            ##         nom.factor <- z[ind_factor+1]
            ##         for (v in nom.factor)
            ##             {
            ##                 mod <- levels(as.factor(olddata[,v]))
            ##                 if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
            ##                 newdata1[,v] <- factor(newdata1[,v], levels=mod)
            ##             }
            ##     }
            call_mixture <- gsub("factor","",call_mixture) 

            ## z <- all.names(call_survival)
            ## ind_factor <- which(z=="factor")
            ## if(length(ind_factor))
            ##     {
            ##         nom.factor <- z[ind_factor+1]  
            ##         for (v in nom.factor)
            ##             {
            ##                 mod <- levels(as.factor(olddata[,v]))
            ##                 if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
            ##                 newdata1[,v] <- factor(newdata1[,v], levels=mod)
            ##             }
            ##     }
            call_survival <- gsub("factor","",call_survival)
            
call_mixture <- formula(paste("~",call_mixture,sep=""))
call_random <- formula(paste("~",call_random,sep=""))
call_classmb <- formula(paste("~",call_classmb,sep=""))
call_survival <- formula(paste("~",call_survival,sep=""))   
            
### Traitement des donnees manquantes
            ## fixed
            ##mcall <- x$call[c(1,match(c("data"),names(x$call),0))]
            mcall <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
            mcall$na.action <- na.action
            mcall$data <- newdata1


            ## fixed
            m <- mcall
            m$formula <- formula(paste("~",call_fixed,sep=""))
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.fixed <- attr(m,"na.action")

            ## mixture
            if((length(attr(terms(call_mixture),"term.labels"))+attr(terms(call_mixture),"intercept"))>0)
                {
                    id.X_mixture <- 1
                    m <- mcall
                    m$formula <- call_mixture
                    m[[1]] <- as.name("model.frame")	
                    m <- eval(m, sys.parent()) 
                    na.mixture <- attr(m,"na.action")
                }
            else
                {
                    id.X_mixture <- 0
                    na.mixture <- NULL
                }

            ## random
            if((length(attr(terms(call_random),"term.labels"))+attr(terms(call_random),"intercept"))>0)
                {
                    id.X_random <- 1
                    m <- mcall
                    m$formula <- call_random
                    m[[1]] <- as.name("model.frame")	
                    m <- eval(m, sys.parent()) 
                    na.random <- attr(m,"na.action")
                }
            else
                {
                    id.X_random <- 0
                    na.random <- NULL
                }

            ## classmb
            if((length(attr(terms(call_classmb),"term.labels"))+attr(terms(call_classmb),"intercept"))>0)
                { 
                    id.X_classmb <- 1
                    m <- mcall	
                    m$formula <- call_classmb
                    m[[1]] <- as.name("model.frame")	
                    m <- eval(m, sys.parent()) 
                    na.classmb <- attr(m,"na.action")
                }
            else
                {
                    id.X_classmb <- 0
                    na.classmb <- NULL
                }

            ##survival
            if((length(attr(terms(call_survival),"term.labels"))+attr(terms(call_survival),"intercept"))>0)
                {
                    id.X_survival <- 1
                    m <- mcall
                    m$formula <- call_survival
                    m[[1]] <- as.name("model.frame")	
                    m <- eval(m, sys.parent()) 
                    na.survival <- attr(m,"na.action")
                }
            else
                {
                    id.X_survival <- 0
                    na.survival <- NULL
                }

            ##cor
            na.cor <- NULL
            if(x$N[7]>0)
                {
                    z <- which(x$idcor==1)
                    var.cor <- newdata1[,x$Names$Xnames[z]]
                    na.cor <- which(is.na(var.cor))
                }

            ##var.time
            if(!missing( var.time))
                {
                    if(!(var.time %in% colnames(newdata))) stop("'var.time' should be included in newdata")
                    if(var.time %in% colnames(newdata1))
                        {
                            times <- newdata1[,var.time,drop=FALSE]
                        }
                    else
                        {
                            times <- newdata[,var.time,drop=FALSE]
                        }
                }
            else
                {
                    times <- newdata[,1,drop=FALSE]
                }

            
            ## Table sans donnees manquante: newdata1
            na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.survival,na.cor))
            if(length(na.action))
                {
                    newdata1 <- newdata1[-na.action,,drop=FALSE]
                    times <- times[-na.action,,drop=FALSE]
                }
            

            ## create one data frame for each formula (useful with factors)
            newdata1fixed <- newdata1
            for(v in colnames(newdata1fixed))
            {
                if(v %in% names(x$levels$levelsfixed))
                {
                    if(!is.null(x$levels$levelsfixed[[v]]))
                    {
                        newdata1fixed[,v] <- factor(newdata1fixed[,v], levels=x$levels$levelsfixed[[v]])
                        if(any(is.na(newdata1fixed[,v]))) stop(paste("Wrong factor level in variable",v))
                    }
                }
            }
            newdata1mixture <- newdata1
            for(v in colnames(newdata1mixture))
            {
                if(v %in% names(x$levels$levelsmixture))
                {
                    if(!is.null(x$levels$levelsmixture[[v]]))
                    {
                        newdata1mixture[,v] <- factor(newdata1mixture[,v], levels=x$levels$levelsmixture[[v]])
                        if(any(is.na(newdata1mixture[,v]))) stop(paste("Wrong factor level in variable",v))
                    }
                }
            }
            newdata1random <- newdata1
            for(v in colnames(newdata1random))
            {
                if(v %in% names(x$levels$levelsrandom))
                {
                    if(!is.null(x$levels$levelsrandom[[v]]))
                    {
                        newdata1random[,v] <- factor(newdata1random[,v], levels=x$levels$levelsrandom[[v]])
                        if(any(is.na(newdata1random[,v]))) stop(paste("Wrong factor level in variable",v))
                    }
                }
            }
            newdata1classmb <- newdata1
            for(v in colnames(newdata1classmb))
            {
                if(v %in% names(x$levels$levelsclassmb))
                {
                    if(!is.null(x$levels$levelsclassmb[[v]]))
                    {
                        newdata1classmb[,v] <- factor(newdata1classmb[,v], levels=x$levels$levelsclassmb[[v]])
                        if(any(is.na(newdata1classmb[,v]))) stop(paste("Wrong factor level in variable",v))
                    }
                }
            }
            newdata1surv <- newdata1
            for(v in colnames(newdata1surv))
            {
                if(v %in% names(x$levels$levelssurv))
                {
                    if(!is.null(x$levels$levelssurv[[v]]))
                    {
                        newdata1surv[,v] <- factor(newdata1surv[,v], levels=x$levels$levelssurv[[v]])
                        if(any(is.na(newdata1surv[,v]))) stop(paste("Wrong factor level in variable",v))
                    }
                }
            }
            
            ## Construction de nouvelles var explicatives sur la nouvelle table

            ## intercept (car plus forcement dans fixed si link=NULL)
            X_intercept <- model.matrix(~1,data=newdata1)
            colnames(X_intercept) <- "intercept"
            
            ## fixed            
            X_fixed <- model.matrix(formula(paste("~",call_fixed,sep="")),data=newdata1fixed)
            if(colnames(X_fixed)[1]=="(Intercept)")
                {
                    X_fixed <- X_fixed[,-1,drop=FALSE]
                }
            
            ## mixture
            if(id.X_mixture==1)
                {
                    X_mixture <- model.matrix(call_mixture,data=newdata1mixture)	
                    if(colnames(X_mixture)[1]=="(Intercept)")
                        {
                            colnames(X_mixture)[1] <- "intercept"
                        }
                }
            
            ## random
            if(id.X_random==1)
                {
                    X_random <- model.matrix(call_random,data=newdata1random)	
                    if(colnames(X_random)[1]=="(Intercept)")
                        {
                            colnames(X_random)[1] <- "intercept"
                        }
                }
            
            ## classmb
            if(id.X_classmb==1)
                { 
                    X_classmb <- model.matrix(call_classmb,data=newdata1classmb)
                    if(colnames(X_classmb)[1]=="(Intercept)")
                        {
                            colnames(X_classmb)[1] <- "intercept"
                        }
                }

            ## survival
            if(id.X_survival==1)
                { 
                    X_survival <- model.matrix(call_survival,data=newdata1surv)
                    if(colnames(X_survival)[1]=="(Intercept)")
                    {
                        colnames(X_survival)[1] <- "intercept"
                    }
                }	
            
            ##cor
            if(x$N[7]>0)  #on reprend la variable de temps de cor
                {
                    z <- which(x$idcor==1)
                    var.cor <- newdata1[,x$Names$Xnames[z]]
                }
            
            
            ## Construction des var expli
            newdata1 <- cbind(X_intercept,X_fixed)
            colX <- c("intercept",colnames(X_fixed))

            if(id.X_mixture == 1)
                {
                    for(i in 1:length(colnames(X_mixture)))
                        {
                            if((colnames(X_mixture)[i] %in% colnames(newdata1))==FALSE)
                                {
                                    newdata1 <- cbind(newdata1,X_mixture[,i])
                                    colnames(newdata1) <- c(colX,colnames(X_mixture)[i])
                                    colX <- colnames(newdata1)
                                }
                        }
                }
            if(id.X_random == 1)
                {
                    for(i in 1:length(colnames(X_random)))
                        {
                            if((colnames(X_random)[i] %in% colnames(newdata1))==FALSE)
                                {
                                    newdata1 <- cbind(newdata1,X_random[,i])
                                    colnames(newdata1) <- c(colX,colnames(X_random)[i])
                                    colX <- colnames(newdata1)
                                }	 
                        }
                }
            if(id.X_classmb == 1)
                {
                    for(i in 1:length(colnames(X_classmb)))
                        {
                            if((colnames(X_classmb)[i] %in% colnames(newdata1))==FALSE)
                                {
                                    newdata1 <- cbind(newdata1,X_classmb[,i])
                                    colnames(newdata1) <- c(colX,colnames(X_classmb)[i])
                                    colX <- colnames(newdata1)
                                }	
                        }
                }
            if(id.X_survival == 1)
                {
                    for(i in 1:length(colnames(X_survival)))
                        {
                            if((colnames(X_survival)[i] %in% colnames(newdata1))==FALSE)
                                {
                                    newdata1 <- cbind(newdata1,X_survival[,i])
                                    colnames(newdata1) <- c(colX,colnames(X_survival)[i])
                                    colX <- colnames(newdata1)
                                }
                        }
                }

            if(x$N[7]>0)
                {
                    idspecif <- matrix(x$idspecif,nbevt,length(x$idg),byrow=TRUE)
                    if(x$idea[z]==0 & x$idprob[z]==0 & x$idg[z]==0 & x$idcom[z]==0 & all(idspecif[,z]==0))
                        {
                            newdata1 <- cbind(newdata1,var.cor)
                            colnames(newdata1) <- c(colX,x$Names$Xnames[z])
                            colX <- colnames(newdata1)
                        }
                    
                }



            ## calcul des predictions

            ## sans transfo
            if(x$linktype==-1)
                {
                    placeV <- list() #places pour les variances des beta
                    placeV$commun <- NA
                    for(i in 1:x$ng)
                        {
                            placeV[paste("class",i,sep="")] <- NA
                        }
                    
                    kk<-0
                    for(k in 1:length(x$idg))
                        {
                            if(x$idg[k]==1)
                                {
                                    X1 <- cbind(X1,newdata1[,k])
                                    place <- sum(x$N[1:3])+kk
                                    b1 <- c(b1,x$best[place+1])
                                    placeV$commun <- c(placeV$commun,place+1)  
                                    kk <- kk+1
                                }
                            
                            if(x$idg[k]==2)
                                {
                                    X2 <- cbind(X2,newdata1[,k])
                                    place1 <- sum(x$N[1:3])+kk+1
                                    place2 <- sum(x$N[1:3])+kk+x$ng
                                    b2 <- rbind(b2,x$best[place1:place2])
                                    for(i in 1:x$ng)
                                        {
                                            placeV[[paste("class",i,sep="")]] <- c(placeV[[paste("class",i,sep="")]],sum(x$N[1:3])+kk+i)
                                        }
                                    kk <- kk+x$ng
                                }
                        }

                    Ypred<-matrix(0,length(newdata1[,1]),x$ng)
                    for(g in 1:x$ng)
                        {
                            if(length(b1) != 0)
                                {
                                    Ypred[,g]<- X1 %*% b1 
                                }
                            if(length(b2) != 0)
                                {
                                    Ypred[,g]<- Ypred[,g] + X2 %*% b2[,g]
                                }
                        }


                    colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
                    if (x$ng==1) colnames(Ypred) <- "Ypred"
                    res <- Ypred

                    if(isTRUE(draws))
                        {
                            ##IC pour les predictions 
                            lower <- matrix(0,nrow(Ypred),ncol(Ypred))  
                            upper <- matrix(0,nrow(Ypred),ncol(Ypred))
                            colnames(lower) <- paste("lower.class",1:x$ng,sep="")
                            colnames(upper) <- paste("upper.class",1:x$ng,sep="")
                            
                            ##extraction de Var(beta)
                            Vbeta <- matrix(0,x$N[4],x$N[4])
                            npm <- length(x$best)
                            indice <- 1:npm * (1:npm+1) /2
                            indtmp <- indice[sum(x$N[1:3])+1:x$N[4]]
                            indtmp <- cbind(indtmp-0:(length(indtmp)-1),indtmp)

                            indV <- NULL
                            for(i in 1:nrow(indtmp))
                                {
                                    indV <- c(indV,seq(indtmp[i,1],indtmp[i,2]))
                                }
                            
                            Vbeta[upper.tri(Vbeta, diag=TRUE)] <- x$V[indV]
                            Vbeta <- t(Vbeta)
                            Vbeta[upper.tri(Vbeta, diag=TRUE)] <- x$V[indV]

                            if(x$ng==1)
                                {                                   
                                    varY <- apply(X1,1,function(x) matrix(x,nrow=1) %*% Vbeta %*% matrix(x,ncol=1))

                                    lower[,1] <- Ypred[,1] - 1.96 * sqrt(varY)
                                    upper[,1] <- Ypred[,1] + 1.96 * sqrt(varY)

                                }
                            else
                                {
                                    for(g in 1:x$ng)
                                        {
                                            ind <- na.omit(c(placeV[["commun"]],placeV[[paste("class",g,sep="")]]))
                                            
                                            X12 <- cbind(X1,X2)    
                                            
                                            X12 <- X12[,order(ind),drop=FALSE]
                                            
                                            Vbeta_g <- Vbeta[sort(ind)-sum(x$N[1:3]),sort(ind)-sum(x$N[1:3])]

                                            varY <- apply(X12,1, function(x) matrix(x,nrow=1) %*% Vbeta_g %*% matrix(x,ncol=1))

                                            lower[,g] <- Ypred[,g] -1.96 * sqrt(varY)
                                            upper[,g] <- Ypred[,g] +1.96 * sqrt(varY)
                                        }
                                    
                                }

                            res <- cbind(res,lower,upper)
                        }                            
                    
                }
            
            
            ## cas du lien lineaire
            if(x$linktype==0)
                {
                    placeV <- list() #places pour les variances des beta
                    placeV$commun <- NA
                    for(i in 1:x$ng)
                        {
                            placeV[paste("class",i,sep="")] <- NA
                        }
                    
                    kk<-0
                    for(k in 1:length(x$idg))
                        {
                            if(x$idg[k]==1)
                                {
                                    X1 <- cbind(X1,newdata1[,k])
                                    if (k==1) b1 <- c(b1,0)
                                    if (k>1)
                                        {
                                            place <- sum(x$N[1:3])+kk
                                            b1 <- c(b1,x$best[place+1])
                                            placeV$commun <- c(placeV$commun,place+1)  
                                            kk <- kk+1
                                        }
                                }
                            
                            if(x$idg[k]==2)
                                {
                                    X2 <- cbind(X2,newdata1[,k])
                                    if (k==1)
                                        {
                                            place1 <- sum(x$N[1:3])+kk+1
                                            place2 <- sum(x$N[1:3])+kk+x$ng-1
                                            b2 <- rbind(b2,c(0,x$best[place1:place2]))
                                            for(i in 2:x$ng)
                                                {
                                                    placeV[[paste("class",i,sep="")]] <- c(placeV[[paste("class",i,sep="")]],sum(x$N[1:3])+kk+i-1)
                                                } 
                                            kk <- kk+x$ng-1
                                        }
                                    if (k>1)
                                        {
                                            place1 <- sum(x$N[1:3])+kk+1
                                            place2 <- sum(x$N[1:3])+kk+x$ng
                                            b2 <- rbind(b2,x$best[place1:place2])
                                            for(i in 1:x$ng)
                                                {
                                                    placeV[[paste("class",i,sep="")]] <- c(placeV[[paste("class",i,sep="")]],sum(x$N[1:3])+kk+i)
                                                }
                                            kk <- kk+x$ng
                                        }
                                }
                        }

                    Ypred<-matrix(0,length(newdata1[,1]),x$ng)
                    for(g in 1:x$ng)
                        {
                            if(length(b1) != 0)
                                {
                                    Ypred[,g]<- X1 %*% b1 
                                }
                            if(length(b2) != 0)
                                {
                                    Ypred[,g]<- Ypred[,g] + X2 %*% b2[,g]
                                }

                            Ypred[,g] <- Ypred[,g]*abs(x$best[length(x$best)])+x$best[length(x$best)-1]
                        }


                    colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
                    if (x$ng==1) colnames(Ypred) <- "Ypred"
                    res <- Ypred

                    if(isTRUE(draws))
                        {
                            ##IC pour les predictions 
                            lower <- matrix(0,nrow(Ypred),ncol(Ypred))  
                            upper <- matrix(0,nrow(Ypred),ncol(Ypred))
                            colnames(lower) <- paste("lower.class",1:x$ng,sep="")
                            colnames(upper) <- paste("upper.class",1:x$ng,sep="")
                            
                            ##extraction de Var(beta)
                            Vbeta <- matrix(0,x$N[4],x$N[4])
                            npm <- length(x$best)
                            indice <- 1:npm * (1:npm+1) /2
                            indtmp <- indice[sum(x$N[1:3])+1:x$N[4]]
                            indtmp <- cbind(indtmp-0:(length(indtmp)-1),indtmp)

                            indV <- NULL
                            for(i in 1:nrow(indtmp))
                                {
                                    indV <- c(indV,seq(indtmp[i,1],indtmp[i,2]))
                                }
                            
                            Vbeta[upper.tri(Vbeta, diag=TRUE)] <- x$V[indV]
                            Vbeta <- t(Vbeta)
                            Vbeta[upper.tri(Vbeta, diag=TRUE)] <- x$V[indV]

                            ## Var(link)
                            Vlink <- matrix(0,2,2)
                            Vlink[1,1] <- x$V[length(x$V)-npm]
                            Vlink[2,2] <- x$V[length(x$V)]
                            Vlink[1,2] <- x$V[length(x$V)-1]
                            Vlink[2,1] <- Vlink[1,2]

                            ## Cov(beta,link)
                            covbl <- matrix(0,x$N[4],2)
                            covbl[,1] <- x$V[(npm-2)*(npm-1)/2+sum(x$N[1:3])+1:x$N[4]]
                            covbl[,2] <- x$V[npm*(npm+1)/2-sum(x$N[4:8])+1:x$N[4]]

                            
                            if(x$ng==1)
                                {
                                    varbl <- cbind(Vbeta,covbl)
                                    varbl <- rbind(varbl,cbind(t(covbl),Vlink))
                                    
                                    X <- cbind(X1[,-1,drop=FALSE],1,1)
                                    varY <- apply(X,1,function(x) matrix(x,nrow=1) %*% varbl %*% matrix(x,ncol=1))

                                    lower[,1] <- Ypred[,1] - 1.96 * sqrt(varY)
                                    upper[,1] <- Ypred[,1] + 1.96 * sqrt(varY)

                                }
                            else
                                {
                                    for(g in 1:x$ng)
                                        {
                                            ind <- na.omit(c(placeV[["commun"]],placeV[[paste("class",g,sep="")]]))
                                            
                                            if(g==1)
                                                {
                                                    if(x$idg[1]==1)
                                                        {
                                                            X12 <- X12 <- cbind(X1[,-1,drop=FALSE],X2)
                                                        }
                                                    
                                                    if(x$idg[1]==2)
                                                        {
                                                            X12 <- X12 <- cbind(X1,X2[,-1,drop=FALSE])
                                                        }
                                                }
                                            else
                                                {
                                                    X12 <- cbind(X1,X2)    
                                                }
                                            
                                            X12 <- X12[,order(ind),drop=FALSE]

                                            covbl_g <- covbl[sort(ind)-sum(x$N[1:3]),,drop=FALSE]
                                            Vbeta_g <- Vbeta[sort(ind)-sum(x$N[1:3]),sort(ind)-sum(x$N[1:3])]

                                            varbl_g <- cbind(Vbeta_g,covbl_g)
                                            varbl_g <- rbind(varbl_g,cbind(t(covbl_g),Vlink))

                                            varY <- apply(cbind(X12,1,1),1, function(x) matrix(x,nrow=1) %*% varbl_g %*% matrix(x,ncol=1))

                                            lower[,g] <- Ypred[,g] -1.96 * sqrt(varY)
                                            upper[,g] <- Ypred[,g] +1.96 * sqrt(varY)
                                        }
                                    
                                }

                            res <- cbind(res,lower,upper)
                            
                        }
                }


            if(x$linktype %in% c(1,2))
                {
                    ## arg pour fortran
                    maxmes <- nrow(newdata1)
                    nv <- length(x$Names$Xnames)
                    npm <- length(x$best)
                    best <- x$best
                    if(x$idiag==0 & x$N[5]>0) best[sum(x$N[1:4])+1:x$N[5]] <- x$cholesky
                    if(x$idiag==1 & x$N[5]>0) best[sum(x$N[1:4])+1:x$N[5]] <- sqrt(best[sum(x$N[1:4])+1:x$N[5]])
                    nbzitr <- length(x$linknodes)
                    Ydiscret <- 0
                    Ymarg <- rep(0,maxmes*x$ng)

                    ## ordonner best comme ds lcmm
                    b1 <- NULL
                    if(x$N[1]>0) b1 <- c(b1,best[1:x$N[1]])
                    if(x$N[4]>0) b1 <- c(b1,best[sum(x$N[1:3])+1:x$N[4]])
                    if(x$N[5]>0) b1 <- c(b1,best[sum(x$N[1:4])+1:x$N[5]])
                    if(x$N[6]>0) b1 <- c(b1,best[sum(x$N[1:5])+1:x$N[6]])
                    if(x$N[8]>0) b1 <- c(b1,best[sum(x$N[1:7])+1:x$N[8]])
                    if(x$N[7]>0) b1 <- c(b1,best[sum(x$N[1:6])+1:x$N[7]])

                    idprob1 <- x$idprob
                    idprob1[1] <- 0
                    npm1 <- length(b1)

                    
                    if(!isTRUE(draws))
                        {
                            out <- .Fortran(C_predictcont,
                                            as.double(newdata1),
                                            as.integer(idprob1),
                                            as.integer(x$idea),
                                            as.integer(x$idg),
                                            as.integer(x$idcor),
                                            as.integer(x$ng),
                                            as.integer(x$N[7]),
                                            as.integer(nv),
                                            as.integer(maxmes),
                                            as.integer(x$idiag),
                                            as.integer(x$N[6]),
                                            as.integer(npm1),
                                            as.double(b1),
                                            as.double(x$epsY),
                                            as.integer(x$linktype),
                                            as.integer(nbzitr),
                                            as.double(x$linknodes),
                                            as.integer(nsim),
                                            as.integer(methInteg),
                                            as.integer(Ydiscret),
                                            y=as.double(Ymarg))
                            
                            out$y[which(out$y==9999)] <- NA
                            res <- matrix(out$y,ncol=x$ng,byrow=FALSE)

                            if (x$ng==1) colnames(res) <- "Ypred"
                            if (x$ng>1) colnames(res) <- paste("Ypred_class",1:x$ng,sep="")
                        }
                    else
                        {
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
                            
                            doOneDraw <- function()
                                {
                                    bdraw <- rnorm(npm)
                                    bdraw <- best + Chol %*% bdraw
                                    
                                    bdraw1 <- NULL
                                    if(x$N[1]>0) bdraw1 <- c(bdraw1,bdraw[1:x$N[1]])
                                    if(x$N[4]>0) bdraw1 <- c(bdraw1,bdraw[sum(x$N[1:3])+1:x$N[4]])
                                    if(x$N[5]>0) bdraw1 <- c(bdraw1,bdraw[sum(x$N[1:4])+1:x$N[5]])
                                    if(x$N[6]>0) bdraw1 <- c(bdraw1,bdraw[sum(x$N[1:5])+1:x$N[6]])
                                    if(x$N[8]>0) bdraw1 <- c(bdraw1,bdraw[sum(x$N[1:7])+1:x$N[8]])
                                    if(x$N[7]>0) bdraw1 <- c(bdraw1,bdraw[sum(x$N[1:6])+1:x$N[7]])
                                    Ymarg <- rep(0,maxmes*x$ng)
                                    out <- .Fortran(C_predictcont,
                                                    as.double(newdata1),
                                                    as.integer(idprob1),
                                                    as.integer(x$idea),
                                                    as.integer(x$idg),
                                                    as.integer(x$idcor),
                                                    as.integer(x$ng),
                                                    as.integer(x$N[7]),
                                                    as.integer(nv),
                                                    as.integer(maxmes),
                                                    as.integer(x$idiag),
                                                    as.integer(x$N[6]),
                                                    as.integer(npm1),
                                                    as.double(bdraw1),
                                                    as.double(x$epsY),
                                                    as.integer(x$linktype),
                                                    as.integer(nbzitr),
                                                    as.double(x$linknodes),
                                                    as.integer(nsim),
                                                    as.integer(methInteg),
                                                    as.integer(Ydiscret),
                                                    y=as.double(Ymarg))
                                    
                                    out$y[which(out$y==9999)] <- NA
                                    
                                    return(out$y)
                                }

                            ydraws <- replicate(ndraws,doOneDraw())
                        

                            f <- function(x)
                                {
                                    quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
                                }
                            
                            ydistr <- apply(ydraws,1,FUN=f)
                            Ypred_50 <- matrix(ydistr[2,],ncol=x$ng,byrow=FALSE)
                            Ypred_2.5 <- matrix(ydistr[1,],ncol=x$ng,byrow=FALSE)
                            Ypred_97.5 <- matrix(ydistr[3,],ncol=x$ng,byrow=FALSE)

                            res <- cbind(Ypred_50,Ypred_2.5,Ypred_97.5)


                            if (x$ng==1)
                                {
                                    colnames(res) <- c("Ypred_50","Ypred_2.5","Ypred_97.5")
                                }
                            if (x$ng>1)
                                {
                                    colnames(res) <- c(paste("Ypred_50_class",1:x$ng,sep=""),paste("Ypred_2.5_class",1:x$ng,sep=""),paste("Ypred_97.5_class",1:x$ng,sep=""))
                                }
                            
                        }

                   
                }

            res.list <- NULL
            res.list$pred <- res
            res.list$times <- times        

        }
    else
        {
            cat("Predictions can not be computed since the program stopped abnormally. \n")
            res.list <- list(pred=NA,times=NA)
        }
    

    class(res.list) <- "predictY"
    return(res.list)
}

