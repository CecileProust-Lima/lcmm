#' @rdname predictY
#' @export
#'
predictY.lcmm <- function(x,newdata,var.time,methInteg=0,nsim=20,draws=FALSE,ndraws=2000,na.action=1,...){


    if(missing(newdata)) stop("The argument newdata should be specified")
    if(missing(x)) stop("The argument x should be specified")
    if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")
                                        # ad 2/04/2012 Xnames2
    if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) {                               
        cat("newdata should at least include the following covariates: ", "\n")
        cat(x$Xnames2[-1], "\n")}
    if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) stop("see above")
    if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
    if (!(methInteg %in% c(0,1))) stop("The integration method must be either 0 for Gauss-Hermite or 1 for Monte-Carlo")
    if ((methInteg==0)&(!(nsim %in% c(5,7,9,15,20,30,40,50)))) stop("For Gauss-Hermite integration method, 'nsim' should be either 5,7,9,15,20,30,40 or 50")
                                        #if(missing(var.time)) stop("missing argument 'var.time'")
                                        #if(!(var.time %in% colnames(newdata))) stop("'var.time' should be included in newdata")


    call_fixed <- x$call$fixed[3]
    if(is.null(x$call$random)) {call_random <- -1} else call_random <- x$call$random[2]
    if(is.null(x$call$classmb)) {call_classmb <- -1} else call_classmb <- x$call$classmb[2]
    if(is.null(x$call$mixture)) {call_mixture <- -1} else call_mixture <- x$call$mixture[2]


    if(x$conv==1|x$conv==2|x$conv==3) {

                                        #------------> changement Cecile 10/04/2012
        ## add 12/04/2012
        if(x$Xnames2[1]!="intercept"){
            newdata1 <- newdata[,x$Xnames2]
            colnames(newdata1) <- x$Xnames
            newdata1 <- data.frame(newdata1)
        }else{
            newdata1 <- cbind(rep(1,length=length(newdata[,1])),newdata[,x$Xnames2[-1]])
            colnames(newdata1) <- c("intercept",x$Xnames2[-1])
            newdata1 <- data.frame(newdata1)
        }


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



        X1 <- NULL                                                              
        X2 <- NULL
        b1 <- NULL
        b2 <- NULL


        if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

        if(na.action==1){
            na.action=na.omit
        }else{
            na.action=na.fail
        }

### pour les facteurs
        ##donnees de l estimation
        if(!is.null(x$data))
        {
            olddata <- x$data
        }
        else
        {
            olddata <- eval(x$call$data)
        }

                                        #cas ou une variable du dataset est un facteur
        for(v in x$Xnames2[-1])
        {
            if (is.factor(olddata[,v]))
            {
                mod <- levels(olddata[,v])
                if (!(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
                newdata1[,v] <- factor(newdata1[,v], levels=mod)
            }
        }
        
                                        #cas ou on a factor() dans l'appel
        z <- all.names(call_fixed)
        ind_factor <- which(z=="factor")
        if(length(ind_factor))
        {
            nom.factor <- z[ind_factor+1]  
            for (v in nom.factor)
            {
                mod <- levels(as.factor(olddata[,v]))
                if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
                newdata1[,v] <- factor(newdata1[,v], levels=mod)
            }
        }
        call_fixed <- gsub("factor","",call_fixed)

        z <- all.names(call_random)
        ind_factor <- which(z=="factor")
        if(length(ind_factor))
        {
            nom.factor <- z[ind_factor+1]
            for (v in nom.factor)
            {
                mod <- levels(as.factor(olddata[,v]))
                if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
                newdata1[,v] <- factor(newdata1[,v], levels=mod)
            }
        }
        call_random <- gsub("factor","",call_random)
        
        z <- all.names(call_classmb)
        ind_factor <- which(z=="factor")
        if(length(ind_factor))
        {
            nom.factor <- z[ind_factor+1]
            for (v in nom.factor)
            {
                mod <- levels(as.factor(olddata[,v]))
                if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
                newdata1[,v] <- factor(newdata1[,v], levels=mod)
            }
        }
        call_classmb <- gsub("factor","",call_classmb)
        
        z <- all.names(call_mixture)
        ind_factor <- which(z=="factor")
        if(length(ind_factor))
        {
            nom.factor <- z[ind_factor+1]
            for (v in nom.factor)
            {
                mod <- levels(as.factor(olddata[,v]))
                if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
                newdata1[,v] <- factor(newdata1[,v], levels=mod)
            }
        }
        call_mixture <- gsub("factor","",call_mixture)   

        call_mixture <- formula(paste("~",call_mixture,sep=""))
        call_random <- formula(paste("~",call_random,sep=""))
        call_classmb <- formula(paste("~",call_classmb,sep=""))

### Traitement des donnees manquantes

                                        # permet de conserver que data=... dans lcmm ; mcall= objet de type call
                                        #mcall <- x$call[c(1,match(c("data"),names(x$call),0))]
        mcall <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
        mcall$na.action <- na.action
        mcall$data <- newdata1

                                        # fixed
        m <- mcall
        m$formula <- formula(paste("~",call_fixed,sep=""))
        m[[1]] <- as.name("model.frame")	
        m <- eval(m, sys.parent()) 
        na.fixed <- attr(m,"na.action")

                                        # mixture
        if((length(attr(terms(call_mixture),"term.labels"))+attr(terms(call_mixture),"intercept"))>0){
            id.X_mixture <- 1
            m <- mcall
            m$formula <- call_mixture
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.mixture <- attr(m,"na.action")
        }else{
            id.X_mixture <- 0
            na.mixture <- NULL
        }

                                        # random
        if((length(attr(terms(call_random),"term.labels"))+attr(terms(call_random),"intercept"))>0){
            id.X_random <- 1
            m <- mcall
            m$formula <- call_random
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.random <- attr(m,"na.action")
        }else{
            id.X_random <- 0
            na.random <- NULL
        }
                                        # classmb
        if((length(attr(terms(call_classmb),"term.labels"))+attr(terms(call_classmb),"intercept"))>0){ 
            id.X_classmb <- 1
            m <- mcall	
            m$formula <- call_classmb
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.classmb <- attr(m,"na.action")
        }else{
            id.X_classmb <- 0
            na.classmb <- NULL
        }
                                        # cor
        na.cor <- NULL
        if(length(x$N)>5)
        {
            if(x$N[6]>0)
            {
                z <- which(x$idcor0==1)
                var.cor <- newdata1[,x$Xnames[z]]
                na.cor <- which(is.na(var.cor))
            }
        }

                                        #var.time
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

        ## Table sans donnees manquante: newdata
        na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.cor))
        if(length(na.action)){
            newdata1 <- newdata1[-na.action,]
            times <- times[-na.action,,drop=FALSE]
        }






        ## Construction de nouvelles var eplicatives sur la nouvelle table
        ## fixed
	
	X_fixed <- model.matrix(formula(paste("~",call_fixed,sep="")),data=newdata1)
	if(colnames(X_fixed)[1]=="(Intercept)"){
            colnames(X_fixed)[1] <- "intercept"
            int.fixed <- 1
	}
        ## mixture
	if(id.X_mixture == 1){
            X_mixture <- model.matrix(call_mixture,data=newdata1)	
            if(colnames(X_mixture)[1]=="(Intercept)"){
                colnames(X_mixture)[1] <- "intercept"
                int.mixture <- 1
            }
	}	
        ## random
	if(id.X_random == 1){
            X_random <- model.matrix(call_random,data=newdata1)	
            if(colnames(X_random)[1]=="(Intercept)"){
                colnames(X_random)[1] <- "intercept"
                int.random <- 1
            }
	}	
        ## classmb
	if(id.X_classmb == 1){ 
            X_classmb <- model.matrix(call_classmb,data=newdata1)
            colnames(X_classmb)[1] <- "intercept"
	}
        ##cor	
        if(length(x$N)>5)
        {
            if(x$N[6]>0)  #on reprend la variable de temps de cor
            {
                z <- which(x$idcor0==1)
                var.cor <- newdata1[,x$Xnames[z]]
            }
        }

        ## Construction des var expli
        newdata1 <- X_fixed
        colX <- colnames(X_fixed)

        if(id.X_mixture == 1){
            for(i in 1:length(colnames(X_mixture))){
		if((colnames(X_mixture)[i] %in% colnames(newdata1))==FALSE){
                    newdata1 <- cbind(newdata1,X_mixture[,i])
                    colnames(newdata1) <- c(colX,colnames(X_mixture)[i])
                    colX <- colnames(newdata1)
		}
            }
        }
        if(id.X_random == 1){
            for(i in 1:length(colnames(X_random))){
		if((colnames(X_random)[i] %in% colnames(newdata1))==FALSE){
                    newdata1 <- cbind(newdata1,X_random[,i])
                    colnames(newdata1) <- c(colX,colnames(X_random)[i])
                    colX <- colnames(newdata1)
		}	 
            }
        }
        if(id.X_classmb == 1){
            for(i in 1:length(colnames(X_classmb))){
		if((colnames(X_classmb)[i] %in% colnames(newdata1))==FALSE){
                    newdata1 <- cbind(newdata1,X_classmb[,i])
                    colnames(newdata1) <- c(colX,colnames(X_classmb)[i])
                    colX <- colnames(newdata1)
		}	
            }
        }
        if(length(x$N)>5)
        {
            if(x$N[6]>0)
            { 
                if( x$idg0[z]==0 & x$idea0[z]==0 & x$idprob0[z]==0)
                {
                    newdata1 <- cbind(newdata1,var.cor)
                    colnames(newdata1) <- c(colX,x$Xnames[z])
                    colX <- colnames(newdata1)
                }
            }
        }




        nv <- length(x$idg0)
        maxmes <- length(newdata1[,1])
        npm <- length(x$best)
        best <- x$best
        if(x$idiag==0 & x$N[3]>0) best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])] <- x$cholesky
        if(x$idiag==1 & x$N[3]>0) best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])] <- sqrt(best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])])
        nwg <- x$N[4]
        ncor <- 0
        if (x$linktype!=3)
        { 
            if(length(x$N)>5) {ncor <- x$N[6]}
        }

        

### for linear trajectory 
        


        if (x$linktype==0){

            if (!draws) {

                                        # prediction
                X1 <- NULL
                X2 <- NULL
                b1 <- NULL
                b2 <- NULL

                kk<-0
                for(k in 1:length(x$idg0)){
                    if(x$idg0[k]==1){
                        X1 <- cbind(X1,newdata1[,k])
                        if (k==1) b1 <- c(b1,0)
                        if (k>1) {
                            place <- x$N[1]+kk
                            b1 <- c(b1,x$best[place+1])
                            kk <- kk+1
                        }
                    }

                    if(x$idg0[k]==2){
                        X2 <- cbind(X2,newdata1[,k])
                        if (k==1){
                            place1 <- x$N[1]+kk+1
                            place2 <- x$N[1]+kk+x$ng-1
                            b2 <- rbind(b2,c(0,x$best[place1:place2]))
                            kk <- kk+x$ng-1
                        }
                        if (k>1){
                            place1 <- x$N[1]+kk+1
                            place2 <- x$N[1]+kk+x$ng
                            b2 <- rbind(b2,x$best[place1:place2])
                            kk <- kk+x$ng
                        }
                    }
                }

                Ypred<-matrix(0,length(newdata1[,1]),x$ng)
                colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
                if (x$ng==1) colnames(Ypred) <- "Ypred"

                for(g in 1:x$ng){
                    if(length(b1) != 0){
                        Ypred[,g]<- X1 %*% b1 
                    }
                    if(length(b2) != 0){
                        Ypred[,g]<- Ypred[,g] + X2 %*% b2[,g]
                    }

                    Ypred[,g] <- Ypred[,g]*abs(x$best[(npm-ncor)])+x$best[(npm-1-ncor)]
                }
            }

            if (draws){


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



                for (j in 1:ndraws) {

                    bdraw <- rnorm(npm)
                    bdraw <- best + Chol %*% bdraw


                                        # prediction
                    X1 <- NULL
                    X2 <- NULL
                    b1 <- NULL
                    b2 <- NULL


                    kk<-0
                    for(k in 1:length(x$idg0)){
                        if(x$idg0[k]==1){
                            X1 <- cbind(X1,newdata1[,k])
                            if (k==1) b1 <- c(b1,0)
                            if (k>1) {
                                place <- x$N[1]+kk
                                b1 <- c(b1,bdraw[place+1])
                                kk <- kk+1
                            }
                        }

                        if(x$idg0[k]==2){
                            X2 <- cbind(X2,newdata1[,k])
                            if (k==1){
                                place1 <- x$N[1]+kk+1
                                place2 <- x$N[1]+kk+x$ng-1
                                b2 <- rbind(b2,c(0,bdraw[place1:place2]))
                                kk <- kk+x$ng-1
                            }
                            if (k>1){
                                place1 <- x$N[1]+kk+1
                                place2 <- x$N[1]+kk+x$ng
                                b2 <- rbind(b2,bdraw[place1:place2])
                                kk <- kk+x$ng
                            }
                        }
                    }

                    Ypred<-matrix(0,length(newdata1[,1]),x$ng)
                    colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
                    if (x$ng==1) colnames(Ypred) <- "Ypred"

                    for(g in 1:x$ng){
                        if(length(b1) != 0){
                            Ypred[,g]<- X1 %*% b1 
                        }
                        if(length(b2) != 0){
                            Ypred[,g]<- Ypred[,g] + X2 %*% b2[,g]
                        }
                        Ypred[,g] <- Ypred[,g]*abs(bdraw[(npm-ncor)])+bdraw[(npm-1-ncor)]
                    }
                    pred <- as.vector(Ypred)
                    ydraws <- cbind(ydraws,pred)
                }

                f <- function(x) {
                    quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
                }
                ydistr <- apply(ydraws,1,FUN=f)
                Ypred_50 <- matrix(ydistr[2,],ncol=x$ng,byrow=F)
                Ypred_2.5 <- matrix(ydistr[1,],ncol=x$ng,byrow=F)
                Ypred_97.5 <- matrix(ydistr[3,],ncol=x$ng,byrow=F)

                Ypred <- cbind(Ypred_50,Ypred_2.5,Ypred_97.5)


                if (x$ng==1){
                    colnames(Ypred) <- c("Ypred_50","Ypred_2.5","Ypred_97.5")
                }
                if (x$ng>1){
                    colnames(Ypred) <- c(paste("Ypred_50_class",1:x$ng,sep=""),paste("Ypred_2.5_class",1:x$ng,sep=""),paste("Ypred_97.5_class",1:x$ng,sep=""))
                }




            }
        }



### for threshold trajectory

        if (x$linktype==3){


            if(!draws) {

                                        # debut prediction
                X1 <- NULL
                X2 <- NULL
                X3 <- NULL
                b1 <- NULL
                b2 <- NULL

                kk<-0
                for(k in 1:length(x$idg0)){
                    if(x$idg0[k]==1){
                        X1 <- cbind(X1,newdata1[,k])
                        if (k==1) b1 <- c(b1,0)
                        if (k>1) {
                            place <- x$N[1]+kk
                            b1 <- c(b1,x$best[place+1])
                            kk <- kk+1
                        }
                    }

                    if(x$idg0[k]==2){
                        X2 <- cbind(X2,newdata1[,k])
                        if (k==1){
                            place1 <- x$N[1]+kk+1
                            place2 <- x$N[1]+kk+x$ng-1
                            b2 <- rbind(b2,c(0,x$best[place1:place2]))
                            kk <- kk+x$ng-1
                        }
                        if (k>1){
                            place1 <- x$N[1]+kk+1
                            place2 <- x$N[1]+kk+x$ng
                            b2 <- rbind(b2,x$best[place1:place2])
                            kk <- kk+x$ng
                        }
                    }

                    if(x$idea0[k]==1){
                        X3 <- cbind(X3,newdata1[,k])
                    }
                }

                Ypred<-matrix(0,length(newdata1[,1]),x$ng)
                colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
                if (x$ng==1) colnames(Ypred) <- "Ypred"

                for(g in 1:x$ng){
                    if(length(b1) != 0){
                        Ypred[,g]<- X1 %*% b1 
                    }
                    if(length(b2) != 0){
                        Ypred[,g]<- Ypred[,g] + X2 %*% b2[,g]
                    }
                }


                varpred <- 0
                nea <- sum(x$idea0)
                if(nea!=0){
                    if(nea==x$N[3]){
                        varpred <- X3 %*%  best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])]^2
                    }
                    if(nea!=x$N[3]){
                        U <- matrix(0,nrow=nea,ncol=nea)
                        U[upper.tri(U,diag=TRUE)] <- best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])]
                        varpred <-  X3 %*% t(U)
                        varpred <- varpred %*% t(varpred)
                    }
                    if(nea>1) varpred <- diag(varpred)
                }

                wg <- rep(1,x$ng)
                if(x$N[4]!=0&x$ng>1){
                    wg[1:(x$ng-1)] <- best[(x$N[1]+x$N[2]+x$N[3]+1):(x$N[1]+x$N[2]+x$N[3]+x$N[4])]
                    wg <- wg^2
                }

                ntrtot0 <- sum(x$ide==1) 
                seuils <- x$ide
                Nseuils <- length(x$ide)
                seuils[x$ide==1] <- as.vector(best[(npm-ntrtot0+1):npm])   #ncor=0 donc ok
                seuils[x$ide==0] <- 0
                if (Nseuils>=2){
                    cumseuils <- cumsum(seuils[2:Nseuils]*seuils[2:Nseuils])
                    seuils[2:Nseuils] <- rep(seuils[1],(Nseuils-1))+cumseuils
                }

                pred <- Ypred
                for(g in 1:x$ng){
                    Ypred[,g] <- rep(x$linknodes[2],maxmes)
                    for(i in 1:Nseuils)
                        Ypred[,g] <- Ypred[,g] - pnorm((seuils[i]-pred[,g])/sqrt(wg[g]*varpred+1))
                }

                if (x$ng>1) colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
                if (x$ng==1) colnames(Ypred) <- "Ypred"
            }

            if (draws) {

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


                for (j in 1:ndraws) {

                    bdraw <- rnorm(npm)
                    bdraw <- best + Chol %*% bdraw

                                        # debut prediction
                    X1 <- NULL
                    X2 <- NULL
                    X3 <- NULL
                    b1 <- NULL
                    b2 <- NULL

                    kk<-0
                    for(k in 1:length(x$idg0)){
                        if(x$idg0[k]==1){
                            X1 <- cbind(X1,newdata1[,k])
                            if (k==1) b1 <- c(b1,0)
                            if (k>1) {
                                place <- x$N[1]+kk
                                b1 <- c(b1,bdraw[place+1])
                                kk <- kk+1
                            }
                        }

                        if(x$idg0[k]==2){
                            X2 <- cbind(X2,newdata1[,k])
                            if (k==1){
                                place1 <- x$N[1]+kk+1
                                place2 <- x$N[1]+kk+x$ng-1
                                b2 <- rbind(b2,c(0,bdraw[place1:place2]))
                                kk <- kk+x$ng-1
                            }
                            if (k>1){
                                place1 <- x$N[1]+kk+1
                                place2 <- x$N[1]+kk+x$ng
                                b2 <- rbind(b2,bdraw[place1:place2])
                                kk <- kk+x$ng}
                        }

                        if(x$idea0[k]==1){
                            X3 <- cbind(X3,newdata1[,k])
                        }
                    }

                    Ypred<-matrix(0,length(newdata1[,1]),x$ng)
                    colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
                    if (x$ng==1) colnames(Ypred) <- "Ypred"

                    for(g in 1:x$ng){
                        if(length(b1) != 0){
                            Ypred[,g]<- X1 %*% b1 
                        }
                        if(length(b2) != 0){
                            Ypred[,g]<- Ypred[,g] + X2 %*% b2[,g]
                        }
                    }


                    varpred <- 0
                    nea <- sum(x$idea0)
                    if(nea!=0){
                        if(nea==x$N[3]){
                            varpred <- X3 %*%  bdraw[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])]^2 
                        }
                        if(nea!=x$N[3]){
                            U <- matrix(0,nrow=nea,ncol=nea)
                            U[upper.tri(U,diag=TRUE)] <- bdraw[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])]
                            varpred <- X3 %*% t(U) 
                            varpred <- varpred %*% t(varpred)
                        }
                        if(nea>1) varpred <- diag(varpred)
                    }

                    wg <- rep(1,x$ng)
                    if(x$N[4]!=0&x$ng>1){
                        wg[1:(x$ng-1)] <- bdraw[(x$N[1]+x$N[2]+x$N[3]+1):(x$N[1]+x$N[2]+x$N[3]+x$N[4])]
                        wg <- wg^2
                    }


                    ntrtot0 <- sum(x$ide==1) 
                    seuils <- x$ide
                    Nseuils <- length(x$ide)
                    seuils[x$ide==1] <- as.vector(bdraw[(npm-ntrtot0+1):npm])
                    seuils[x$ide==0] <- 0
                    if (Nseuils>=2){
                        cumseuils <- cumsum(seuils[2:Nseuils]*seuils[2:Nseuils])
                        seuils[2:Nseuils] <- rep(seuils[1],(Nseuils-1))+cumseuils
                    }
                    pred <- Ypred
                    for(g in 1:x$ng){
                        Ypred[,g] <- rep(x$linknodes[2],maxmes)
                        for(i in 1:Nseuils)
                            Ypred[,g] <- Ypred[,g] - pnorm((seuils[i]-pred[,g])/sqrt(wg[g]*varpred+1))
                    }
                    pred <- as.vector(Ypred)
                    ydraws <- cbind(ydraws,pred)
                }

                f <- function(x) {
                    quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
                }
                ydistr <- apply(ydraws,1,FUN=f)
                Ypred_50 <- matrix(ydistr[2,],ncol=x$ng,byrow=F)
                Ypred_2.5 <- matrix(ydistr[1,],ncol=x$ng,byrow=F)
                Ypred_97.5 <- matrix(ydistr[3,],ncol=x$ng,byrow=F)

                Ypred <- cbind(Ypred_50,Ypred_2.5,Ypred_97.5)

                if (x$ng==1){
                    colnames(Ypred) <- c("Ypred_50","Ypred_2.5","Ypred_97.5")
                }
                if (x$ng>1){
                    colnames(Ypred) <- c(paste("Ypred_50_class",1:x$ng,sep=""),paste("Ypred_2.5_class",1:x$ng,sep=""),paste("Ypred_97.5_class",1:x$ng,sep=""))
                }
            }
        }



### for splines or beta trajectory


        if (x$linktype %in% c(1,2)){


            nbzitr <- length(x$linknodes)
            epsY <- x$epsY
            Ymarg <- rep(0,maxmes*x$ng)

                                        #cat(c(nv,x$ng,nbzitr,epsY,nwg,nsim,methInteg,x$Ydiscrete),"\n")
                                        #cat(epsY,"\n")

            if (!draws){      
                out <- .Fortran(C_predictcont,
                                as.double(newdata1),
                                as.integer(x$idprob0),
                                as.integer(x$idea0),
                                as.integer(x$idg0),
                                as.integer(x$idcor0),
                                as.integer(x$ng),
                                as.integer(ncor),
                                as.integer(nv),
                                as.integer(maxmes),
                                as.integer(x$idiag),
                                as.integer(nwg),
                                as.integer(npm),
                                as.double(best),
                                as.double(epsY),
                                as.integer(x$linktype),
                                as.integer(nbzitr),
                                as.double(x$linknodes),
                                as.integer(nsim),
                                as.integer(methInteg),
                                as.integer(x$Ydiscrete),
                                Ymarg=as.double(Ymarg))
                
                out$Ymarg[out$Ymarg==9999] <- NA

                                        #cat(out$Ymarg)
                Ypred <- matrix(out$Ymarg,ncol=x$ng,byrow=FALSE)

                if (x$ng==1)colnames(Ypred) <- "Ypred"
                if (x$ng>1)colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
            }

########### ajout ndraws ###############################


            if (draws) {

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

                for (j in 1:ndraws) {

                    bdraw <- rnorm(npm)
                    bdraw <- best + Chol %*% bdraw
                    
                    out <- .Fortran(C_predictcont,
                                    as.double(newdata1),
                                    as.integer(x$idprob0),
                                    as.integer(x$idea0),
                                    as.integer(x$idg0),
                                    as.integer(x$idcor0),
                                    as.integer(x$ng),
                                    as.integer(ncor),
                                    as.integer(nv),
                                    as.integer(maxmes),
                                    as.integer(x$idiag),
                                    as.integer(nwg),
                                    as.integer(npm),
                                    as.double(bdraw),
                                    as.double(epsY),
                                    as.integer(x$linktype),
                                    as.integer(nbzitr),
                                    as.double(x$linknodes),
                                    as.integer(nsim),
                                    as.integer(methInteg),
                                    as.integer(x$Ydiscrete),
                                    Ymarg=as.double(Ymarg))

                    out$Ymarg[out$Ymarg==9999] <- NA
                    ydraws <- cbind(ydraws,out$Ymarg)
                }

                f <- function(x) {
                    quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
                }
                ydistr <- apply(ydraws,1,FUN=f)
                Ypred_50 <- matrix(ydistr[2,],ncol=x$ng,byrow=FALSE)
                Ypred_2.5 <- matrix(ydistr[1,],ncol=x$ng,byrow=FALSE)
                Ypred_97.5 <- matrix(ydistr[3,],ncol=x$ng,byrow=FALSE)

                Ypred <- cbind(Ypred_50,Ypred_2.5,Ypred_97.5)


                if (x$ng==1){
                    colnames(Ypred) <- c("Ypred_50","Ypred_2.5","Ypred_97.5")
                }
                if (x$ng>1){
                    colnames(Ypred) <- c(paste("Ypred_50_class",1:x$ng,sep=""),paste("Ypred_2.5_class",1:x$ng,sep=""),paste("Ypred_97.5_class",1:x$ng,sep=""))
                }
            }
        }


        res.list <- NULL
        res.list$pred <- Ypred
        res.list$times <- times
    }

    else  #cas xconv != 1 ou 2
    { 
        cat("Predictions can not be computed since the program stopped abnormally. \n")
        res.list <- list(pred=NA,times=NA)
    }

    class(res.list) <- "predictY"
    return(res.list)
}




#' Predictions (marginal and possibly subject-specific in some cases) of a \code{hlme},
#' \code{lcmm}, \code{multlcmm} or \code{Jointlcmm} object in the natural scale
#' of the longitudinal outcome(s) computed from a profile of covariates (marginal) or
#' individual data (subject specific in case of \code{hlme}).
#' 
#' For \code{hlme} and \code{Jointlcmm} objects, the function computes the
#' predicted values of the longitudinal marker (in each latent class of ng>1) for a
#' specified profile of covariates.  For \code{lcmm} and \code{multlcmm}
#' objects, the function computes predicted values in the natural scale of the
#' outcomes for a specified profile of covariates. For linear and threshold
#' links, the predicted values are computed analytically. For splines and Beta
#' links, a Gauss-Hermite or Monte-Carlo integration are used to numerically
#' compute the predictions. In addition, for any type of link function,
#' confidence bands (and median) can be computed by a Monte Carlo approximation
#' of the posterior distribution of the predicted values.
#' 
#' 
#' @param x an object inheriting from class \code{lcmm}, \code{hlme},
#' \code{Jointlcmm} or \code{multlcmm} representing a general latent class
#' mixed model.
#' @param newdata data frame containing the data from which predictions are to be
#' computed. The data frame should include at least all the covariates listed
#' in x$Xnames2. Names in the data frame should be exactly x$Xnames2 that are
#' the names of covariates specified in \code{lcmm}, \code{hlme},
#' \code{Jointlcmm} or \code{multlcmm} calls. For \code{hlme} object and marg=FALSE,
#' the grouping structure and values for the outcome should also be specified.
#' @param var.time A character string containing the name of the variable that
#' corresponds to time in the data frame (x axis in the plot).
#' @param methInteg optional integer specifying the type of numerical
#' integration required only for predictions with splines or Beta link
#' functions. Value 0 (by default) specifies a Gauss-Hermite integration which
#' is very rapid but neglects the correlation between the predicted values (in
#' presence of random-effects). Value 1 refers to a Monte-Carlo integration
#' which is slower but correctly account for the correlation between the
#' predicted values.
#' @param nsim For a \code{lcmm}, \code{multlcmm} or \code{Jointlcmm} object
#' only; optional number of points used in the numerical integration with
#' splines or Beta link functions. For methInteg=0, nsim should be chosen among
#' the following values: 5, 7, 9, 15, 20, 30, 40 or 50 (nsim=20 by default). If
#' methInteg=1, nsim should be relatively important (more than 200).
#' @param draws optional boolean specifying whether median and confidence bands
#' of the predicted values should be computed (TRUE) - whatever the type of
#' link function. For a \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}
#' object, a Monte Carlo approximation of the posterior distribution of the
#' predicted values is computed and the median, 2.5\% and 97.5\% percentiles
#' are given. Otherwise, the predicted values are computed at the point
#' estimate. By default, draws=FALSE.
#' @param ndraws For a \code{lcmm}, \code{multlcmm} or \code{Jointlcmm} object
#' only; if draws=TRUE, ndraws specifies the number of draws that should be
#' generated to approximate the posterior distribution of the predicted values.
#' By default, ndraws=2000.
#' @param marg Optional boolean specifying whether the      
#' predictions are marginal (the default) or subject-specific (marg=FALSE). marge=FALSE 
#' only works with \code{hlme} objects.
#' @param subject For a \code{hlme} object with marg=FALSE only, character specifying
#' the name of the grouping strucuture. If NULL (the default), the same as in the model
#' (argument x) will be used.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @param \dots further arguments to be passed to or from other methods. 
#' Only the argument 'median' will be used, other are ignored. 'median' should
#' be a logical indicating whether the median should be computed. By 
#' default, the mean value is computed.
#' @return An object of class \code{predictY} with values :
#' 
#' - \code{pred} : a matrix with the same rows (number and order) as in
#' newdata.
#' 
#' For \code{hlme} objects and \code{lcmm} or \code{Jointlcmm} with
#' \code{draws=FALSE}, returns a matrix with ng columns corresponding to the ng
#' class-specific vectors of predicted values computed at the point estimate
#' 
#' For objects of class \code{lcmm} or \code{Jointlcmm} with \code{draws=TRUE},
#' returns a matrix with ng*3 columns representing the ng class-specific 50\%,
#' 2.5\% and 97.5\% percentiles of the approximated posterior distribution of
#' the class-specific predicted values.
#' 
#' For objects of class \code{multlcmm} with \code{draws=FALSE}, returns a
#' matrix with ng+1 columns: the first column indicates the name of the outcome
#' which is predicted and the ng subsequent columns correspond to the ng
#' class-specific vectors of predicted values computed at the point estimate
#' 
#' For objects of class \code{multlcmm} with \code{draws=TRUE}, returns a
#' matrix with ng*3+1 columns: the first column indicates the name of the
#' outcome which is predicted and the ng*3 subsequent columns correspond to the
#' ng class-specific 50\%, 2.5\% and 97.5\% percentiles of the approximated
#' posterior distribution of the class-specific predicted values.
#' 
#' For objects of class \code{hlme} with \code{marg=FALSE}, returns a matrix
#' with 2+ng columns : the grouping structure, subject-specific predictions (pred_ss) averaged 
#' over classes and the class-specific subject-specific predictions (with the
#' number of the latent class: pred_ss_1,pred_ss_2,...)
#' 
#' - \code{times} : the \code{var.time} variable from \code{newdata}
#' @author Cecile Proust-Lima, Viviane Philipps, Sasha Cuau
#' @seealso \code{\link{lcmm}}, \code{\link{multlcmm}}, \code{\link{hlme}},
#' \code{\link{Jointlcmm}}
#' @examples
#' 
#' 
#' #### Prediction from a 2-class model with a Splines link function
#' \dontrun{
#' ## fitted model
#' m<-lcmm(Ydep2~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
#' subject='ID',ng=2,data=data_lcmm,link="splines",B=c(
#' -0.175,      -0.191,       0.654,      -0.443, 
#' -0.345,      -1.780,       0.913,       0.016, 
#'  0.389,       0.028,       0.083,      -7.349, 
#'  0.722,       0.770,       1.376,       1.653, 
#'  1.640,       1.285))
#' summary(m)
#' ## predictions for times from 0 to 5 for X1=0
#' newdata<-data.frame(Time=seq(0,5,length=100),
#' X1=rep(0,100),X2=rep(0,100),X3=rep(0,100))
#' pred0 <- predictY(m,newdata,var.time="Time")
#' head(pred0)
#' ## Option draws=TRUE to compute a MonteCarlo 
#' # approximation of the predicted value distribution 
#' # (quite long with ndraws=2000 by default)
#' \dontrun{
#' pred0MC <- predictY(m,newdata,draws=TRUE,var.time="Time")
#' }
#' ## predictions for times from 0 to 5 for X1=1
#' newdata$X1 <- 1
#' pred1 <- predictY(m,newdata,var.time="Time")
#' ## Option draws=TRUE to compute a MonteCarlo 
#' # approximation of the predicted value distribution 
#' # (quite long with ndraws=2000 by default)
#' \dontrun{
#' pred1MC <- predictY(m,newdata,draws=TRUE,var.time="Time")
#' }
#' }
#' 
#' @export
#' 
predictY <- function(x,newdata,var.time,...){
  dots <- list(...)
  median <- FALSE
  if(length(dots$median)) median <- as.logical(eval(dots$median))
  if(median==TRUE) 
  {
    .predYmedian(m=x,newdata=newdata,var.time=var.time,...)
  }
  else
    {
       UseMethod("predictY")
    }
}
