#' @rdname predictY
#' @export
#'
predictY.hlme <- function(x, newdata, var.time, draws=FALSE, na.action=1, marg=TRUE, subject=NULL, ...){

    if(missing(newdata)) stop("The argument newdata should be specified")
    if(missing(x)) stop("The argument x should be specified")
    if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
                                        # ad 2/04/2012 Xnames2
    if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) {
        cat("newdata should at least include the following covariates: ", "\n")
        cat(x$Xnames2[-1], "\n")}
    if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) stop("see above")
    if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
                                        #if(missing(var.time)) stop("missing argument 'var.time'")
                                        #if(!(var.time %in% colnames(newdata))) stop("'var.time' should be included in newdata")
    if(!marg)
    {
        if(draws==TRUE) stop("No confidence intervals are provided for subject-specific prediction")
        
        Yname <- as.character(x$call$fixed[[2]])
        subj <- subject
        if(is.null(subj)) subj <- as.character(x$call$subject)

        if(!(Yname %in% colnames(newdata))) stop(paste("newdata should include the outcome", Yname))

        if(!(subj %in% colnames(newdata))) stop(paste("newdata should include the grouping structure", subj))
    }

    call_fixed <- x$call$fixed[3]
    if(is.null(x$call$random)) {call_random <- -1} else call_random <- x$call$random[2]
    if(is.null(x$call$classmb)) {call_classmb <- -1} else call_classmb <- x$call$classmb[2]
    if(is.null(x$call$mixture)) {call_mixture <- -1} else call_mixture <- x$call$mixture[2]


    if(x$conv==1|x$conv==2) {
        res.list <- NULL

        ##------------> changement Cecile 10/04/2012
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


        X1 <- NULL                                                              
        X2 <- NULL
        b1 <- NULL
        b2 <- NULL
        Z <- NULL


        if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

        if(na.action==1){
            na.action=na.omit
        }else{
            na.action=na.fail
        }

### pour les facteurs
        ##donnees de l estimation
        ##    if(!is.null(x$data))
        ##    {
        ##        olddata <- x$data
        ##    }
        ##    else
        ##    {
        ##        olddata <- eval(x$call$data)
        ##    }

        ## #cas ou une variable du dataset est un facteur
        ##  for(v in x$Xnames2[-1])
        ## {
        ##  if (is.factor(olddata[,v]))
        ##  {
        ##   mod <- levels(olddata[,v])
        ##   if (!(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
        ##   newdata1[,v] <- factor(newdata1[,v], levels=mod)
        ##  }
        ## }

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
        
                                        #cas ou on a factor() dans l'appel
        ## z <- all.names(call_fixed)
        ## ind_factor <- which(z=="factor")
        ## if(length(ind_factor))
        ## {
        ##  nom.factor <- z[ind_factor+1]  
        ##  for (v in nom.factor)
        ##  {
        ##   mod <- levels(as.factor(olddata[,v]))
        ##   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
        ##   newdata1[,v] <- factor(newdata1[,v], levels=mod)
        ##  }
        ## }
        call_fixed <- gsub("factor","",call_fixed)

        ## z <- all.names(call_random)
        ## ind_factor <- which(z=="factor")
        ## if(length(ind_factor))
        ## {
        ##  nom.factor <- z[ind_factor+1]
        ##  for (v in nom.factor)
        ##  {
        ##   mod <- levels(as.factor(olddata[,v]))
        ##   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
        ##   newdata1[,v] <- factor(newdata1[,v], levels=mod)
        ##  }
        ## }
        call_random <- gsub("factor","",call_random)
        
        ## z <- all.names(call_classmb)
        ## ind_factor <- which(z=="factor")
        ## if(length(ind_factor))
        ## {
        ##  nom.factor <- z[ind_factor+1]
        ##  for (v in nom.factor)
        ##  {
        ##   mod <- levels(as.factor(olddata[,v]))
        ##   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
        ##   newdata1[,v] <- factor(newdata1[,v], levels=mod)
        ##  }
        ## }
        call_classmb <- gsub("factor","",call_classmb)
        
        ## z <- all.names(call_mixture)
        ## ind_factor <- which(z=="factor")
        ## if(length(ind_factor))
        ## {
        ##  nom.factor <- z[ind_factor+1]
        ##  for (v in nom.factor)
        ##  {
        ##   mod <- levels(as.factor(olddata[,v]))
        ##   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
        ##   newdata1[,v] <- factor(newdata1[,v], levels=mod)
        ##  }
        ## }
        call_mixture <- gsub("factor","",call_mixture)   
        
        call_mixture <- formula(paste("~",call_mixture,sep=""))
        call_random <- formula(paste("~",call_random,sep=""))
        call_classmb <- formula(paste("~",call_classmb,sep=""))

### Traitement des donnees manquantes

        ## permet de conserver que data=... dans lcmm ; mcall= objet de type call
                                        #mcall <- x$call[c(1,match(c("data"),names(x$call),0))]
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

        ## random
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
        ## classmb
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
        ## cor
        na.cor <- NULL
        if(length(x$N)>4)
        {
            if(x$N[5]>0)
            {
                z <- which(x$idcor0==1)
                var.cor <- newdata1[,x$Xnames[z]]
                na.cor <- which(is.na(var.cor))
            }
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

        ## Table sans donnees manquante: newdata
        na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.cor))
        
        if(length(na.action)){
            newdata1 <- newdata1[-na.action,]
            times <- times[-na.action,,drop=FALSE]
        }

        if(marg)
        {
            
            ## create one data frame for each formula (useful with factors)
            newdata1fixed <- newdata1
            for(v in colnames(newdata1fixed))
            {0
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
            

            ## Construction de nouvelles var eplicatives sur la nouvelle table
            ## fixed
            
            X_fixed <- model.matrix(formula(paste("~",call_fixed,sep="")),data=newdata1fixed)
            if(colnames(X_fixed)[1]=="(Intercept)"){
		colnames(X_fixed)[1] <- "intercept"
		int.fixed <- 1
            }
            ## mixture
            if(id.X_mixture==1){
		X_mixture <- model.matrix(call_mixture,data=newdata1mixture)	
		if(colnames(X_mixture)[1]=="(Intercept)"){
                    colnames(X_mixture)[1] <- "intercept"
                    int.mixture <- 1
                }
            }
            
            ## random
            if(id.X_random == 1){
		X_random <- model.matrix(call_random,data=newdata1random)	
		if(colnames(X_random)[1]=="(Intercept)"){
                    colnames(X_random)[1] <- "intercept"
                    int.random <- 1
                }
            }	
            ## classmb
            if(id.X_classmb == 1){ 
                X_classmb <- model.matrix(call_classmb,data=newdata1classmb)
                colnames(X_classmb)[1] <- "intercept"
            }

            ##cor
            if(length(x$N)>4)
            {	
                if(x$N[5]>0)  #on reprend la variable de temps de cor (sans NA)
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
            if(length(x$N)>4)
            {
                if(x$N[5]>0)
                { 
                    if( x$idg0[z]==0 & x$idea0[z]==0 & x$idprob0[z]==0)
                    {
                        newdata1 <- cbind(newdata1,var.cor)
                        colnames(newdata1) <- c(colX,x$Xnames[z])
                        colX <- colnames(newdata1)
                    }
                }
            }

            placeV <- list() #places pour les variances
            placeV$commun <- NA
            for(i in 1:x$ng)
            {
                placeV[paste("class",i,sep="")] <- NA
            }

            kk<-0
            for(k in 1:length(x$idg0)) 
            {
                if(x$idg0[k]==1)
                {
                    X1 <- cbind(X1,newdata1[,k])  
                    place <- x$N[1]+kk
                    b1 <- c(b1,x$best[place+1])
                    placeV$commun <- c(placeV$commun,place+1)
                    kk <- kk+1
                }
                
                if(x$idg0[k]==2) 
                {
                    X2 <- cbind(X2,newdata1[,k])
                    place1 <- x$N[1]+kk+1
                    place2 <- x$N[1]+kk+x$ng
                    b2 <- rbind(b2,x$best[place1:place2])
                    for(i in 1:x$ng)
                    {
                        placeV[[paste("class",i,sep="")]] <- c(placeV[[paste("class",i,sep="")]],x$N[1]+kk+i)
                    }
                    kk <- kk+x$ng
                }
            }
            


            Y <- matrix(0,length(newdata1[,1]),x$ng) 
            for(g in 1:x$ng){
                if(length(b1) != 0){
                    Y[,g] <- X1 %*% b1 
                }
                if(length(b2) != 0){
                    Y[,g] <- Y[,g] + X2 %*% b2[,g]
                }
            }
            



            if(draws==TRUE)
            {
                ##extraction de Var(beta)
                Vbeta <- matrix(0,x$N[2],x$N[2])
                npm <- length(x$best)
                indice <- 1:npm * (1:npm+1) /2
                indtmp <- indice[(x$N[1]+1):(x$N[1]+x$N[2])]
                indtmp <- cbind(indtmp-0:(length(indtmp)-1),indtmp)
                
                indV <- NULL
                for(i in 1:nrow(indtmp))
                {
                    indV <- c(indV,seq(indtmp[i,1],indtmp[i,2]))
                }
                
                Vbeta[upper.tri(Vbeta, diag=TRUE)] <- x$V[indV]
                Vbeta <- t(Vbeta)
                Vbeta[upper.tri(Vbeta,diag=TRUE)] <- x$V[indV] 
                
                
                ##IC pour les predictions 
                lower <- matrix(0,nrow(Y),ncol(Y))  
                upper <- matrix(0,nrow(Y),ncol(Y))
                colnames(lower) <- paste("lower.class",1:x$ng,sep="")
                colnames(upper) <- paste("upper.class",1:x$ng,sep="") 
                
                if(x$ng==1)
                { 
                    varpred <- apply(X1,1,function(x) matrix(x,nrow=1) %*% Vbeta %*% matrix(x,ncol=1)) 
                    lower[,1] <- Y[,1] -1.96 * sqrt(varpred)
                    upper[,1] <- Y[,1] +1.96 * sqrt(varpred)    
                }
                else
                {
                    for(g in 1:x$ng)
                    {
                        ind <- na.omit(c(placeV[["commun"]],placeV[[paste("class",g,sep="")]]))
                        X12 <- cbind(X1,X2)
                        X12 <- X12[,order(ind)]
                        
                        varclass <- Vbeta[sort(ind)-x$N[1],sort(ind)-x$N[1]]
                        varpred <- diag(X12 %*% varclass %*% t(X12))
                        
                        lower[,g] <- Y[,g] -1.96 * sqrt(varpred)
                        upper[,g] <- Y[,g] +1.96 * sqrt(varpred)    
                    }
                }
                
            }

            if(draws==TRUE)
            {
                res <- cbind(Y,lower,upper)
                if(x$ng==1) colnames(res) <- c("Ypred","lower.Ypred","upper.Ypred")
                if(x$ng>1) colnames(res) <- c(paste("Ypred_class",1:x$ng,sep=""),paste("lower.Ypred_class",1:x$ng,sep=""),paste("upper.Ypred_class",1:x$ng,sep=""))

                res.list$pred <- res
                res.list$times <- times
            }
            if(draws==FALSE)
            {
                if (x$ng==1){
                    colnames(Y) <- c("Ypred")
                }
                if (x$ng>1){
                    colnames(Y) <- c(paste("Ypred_class",1:x$ng,sep=""))
                }
                
                res.list$pred <- Y
                res.list$times <- times
            }
        }
        else
        {
            ##if(draws==FALSE)
            ##{
                arguments <- as.list(x$call)
                argfunction <- as.character(arguments[[1]]) 
                arguments[[1]] <- NULL
                arguments[["data"]] <- newdata
                arguments[["B"]] <- x$best 
                arguments[["maxiter"]] <- 0
                arguments[["verbose"]] <- FALSE
                if(!is.null(subject)){
                    arguments[['subject']] <- subject
                }
                newmodel <- do.call(argfunction, c(arguments))
                
                res.list$pred <- newmodel$pred[,c(1,4,6+x$ng+1:x$ng)]
                res.list$times <- NA #times
            ##}
            ##else
            ##{
            ##     doone <- function(bdraw)
            ##     {

            ##     }

            ##     npm <- length(x$best)
            ##     bdraw <- replicate(ndraws, rnorm(npm))

            ##     res <- apply(bdraw, 2, doone)

            ##     res.list$pred <- apply(res,1, quantile, prob=c(0.5, 0.025, 0.975))
            ##     res.list$times <- times
            ##}

        }

    }
    else
    {
        cat("Predictions can not be computed since the program stopped abnormally. \n")
        res.list <- list(pred=NA,times=NA)
    }

    class(res.list) <- "predictY"
    return(res.list)
}       


