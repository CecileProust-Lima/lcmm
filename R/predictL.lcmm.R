
predictL.lcmm <- function(x,newdata,var.time,na.action=1,confint=FALSE,...)
{
if(missing(newdata)) stop("The argument newdata should be specified")
if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")
# ad 2/04/2012 Xnames2
if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) {
stop(paste(c("newdata should at least include the following covariates: ","\n",x$Xnames2[-1]),collapse=" "))}
if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
#if(missing(var.time)) stop("missing argument 'var.time'")
#if(!(var.time %in% colnames(newdata))) stop("'var.time' should be included in newdata")

#if(is.null(x$call$random)) x$call$random <- ~-1
#if(is.null(x$call$classmb)) x$call$classmb <- ~-1
#if(is.null(x$call$mixture)) x$call$mixture <- ~-1
#
call_fixed <- x$call$fixed[3]
if(is.null(x$call$random)) {call_random <- ~-1} else call_random <- x$call$random
if(is.null(x$call$classmb)) {call_classmb <- ~-1} else call_classmb <- x$call$classmb
if(is.null(x$call$mixture)) {call_mixture <- ~-1} else call_mixture <- x$call$mixture


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

 #cas ou une variable du dataset est un facteur
 olddata <- eval(x$call$data)
  for(v in x$Xnames2[-1])
 {
  if (is.factor(olddata[,v]) & !(is.factor(newdata[,v])))
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
if(!is.null(x$call$mixture)){
	m <- mcall
	m$formula <- formula(paste("~",call_mixture,sep=""))
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
	na.mixture <- attr(m,"na.action")
}else{
	na.mixture <- NULL
}

# random
if(!is.null(x$call$random)){
	m <- mcall
	m$formula <- formula(paste("~",call_random,sep=""))
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.random <- attr(m,"na.action")
}else{
	na.random <- NULL
}
# classmb
if(!is.null(x$call$classmb)){ 
	m <- mcall	
	m$formula <- formula(paste("~",call_classmb,sep=""))
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.classmb <- attr(m,"na.action")
}else{
	na.classmb <- NULL
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
na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb))
if(!is.null(na.action)){
	newdata1 <- newdata1[-na.action,]
        times <- times[-na.action]
}




## Construction de nouvelles var eplicatives sur la nouvelle table
## fixed
	
	X_fixed <- model.matrix(formula(paste("~",call_fixed,sep="")),data=newdata1)
	if(colnames(X_fixed)[1]=="(Intercept)"){
		colnames(X_fixed)[1] <- "intercept"
		int.fixed <- 1
	}
## mixture
	if(!is.null(x$call$mixture)){
		X_mixture <- model.matrix(formula(paste("~",call_mixture,sep="")),data=newdata1)	
		if(colnames(X_mixture)[1]=="(Intercept)"){
			colnames(X_mixture)[1] <- "intercept"
			int.mixture <- 1
		}
		id.X_mixture <- 1
	}else{
		id.X_mixture <- 0
	}	
## random
	if(!is.null(x$call$random)){
		X_random <- model.matrix(formula(paste("~",call_random,sep="")),data=newdata1)	
		if(colnames(X_random)[1]=="(Intercept)"){
			colnames(X_random)[1] <- "intercept"
			int.random <- 1
		}
		id.X_random <- 1
	}else{
		id.X_random <- 0
	}	
## classmb
	if(!is.null(x$call$classmb)){ 
		X_classmb <- model.matrix(formula(paste("~",call_classmb,sep="")),data=newdata1)
		colnames(X_classmb)[1] <- "intercept"
		id.X_classmb <- 1
	}else{
		id.X_classmb <- 0
	}
##cor	
if(x$N[6]>0)  #on ajoute la variable de temps de cor
{
 z <- which(x$idcor0==1)
 var.cor <- newdata1[,x$Xnames[z]]
}


## Construction des var expli
newdata1 <- X_fixed
colX <- colnames(X_fixed)

if(id.X_mixture == 1){
	for(i in 1:length(colnames(X_mixture))){
		if((colnames(X_mixture)[i] %in% colnames(newdata1))==F){
			newdata1 <- cbind(newdata1,X_mixture[,i])
                        colnames(newdata1) <- c(colX,colnames(X_mixture)[i])
                        colX <- colnames(newdata1)
		}
	}
}
if(id.X_random == 1){
	for(i in 1:length(colnames(X_random))){
		if((colnames(X_random)[i] %in% colnames(newdata1))==F){
			newdata1 <- cbind(newdata1,X_random[,i])
                        colnames(newdata1) <- c(colX,colnames(X_random)[i])
                        colX <- colnames(newdata1)
		}	 
	}
}
if(id.X_classmb == 1){
	for(i in 1:length(colnames(X_classmb))){
		if((colnames(X_classmb)[i] %in% colnames(newdata1))==F){
			newdata1 <- cbind(newdata1,X_classmb[,i])
                        colnames(newdata1) <- c(colX,colnames(X_classmb)[i])
                        colX <- colnames(newdata1)
		}	
	}
}
if(x$N[6]>0)
{ 
 if( x$idg0[z]==0 & x$idea0[z]==0 & x$idprob0[z]==0)
     {
         newdata1 <- cbind(newdata1,var.cor)
         colnames(newdata1) <- c(colX,x$Xnames[z])
         colX <- colnames(newdata1)
     }
}


### calcul des predicitons

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
    if (k==1) b1 <- c(b1,0)
    if (k>1) 
    {
     place <- x$N[1]+kk
     b1 <- c(b1,x$best[place+1])
     placeV$commun <- c(placeV$commun,place+1)     
     kk <- kk+1
    }
   }
   
   if(x$idg0[k]==2)
   {
    X2 <- cbind(X2,newdata1[,k])
    if (k==1)
    {
     place1 <- x$N[1]+kk+1
     place2 <- x$N[1]+kk+x$ng-1
     b2 <- rbind(b2,c(0,x$best[place1:place2]))
     for(i in 2:x$ng)
     {
      placeV[[paste("class",i,sep="")]] <- c(placeV[[paste("class",i,sep="")]],x$N[1]+kk+i-1)
     } 
     kk <- kk+x$ng-1
    }
    if (k>1)
    {
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
  }
  
  
  Y<-matrix(0,length(newdata1[,1]),x$ng)
  colnames(Y) <- paste("class",1:x$ng,sep="") 
  for(g in 1:x$ng){
  if(length(b1) != 0){
  Y[,g]<- X1 %*% b1 
  }
  if(length(b2) != 0){
  Y[,g]<- Y[,g] + X2 %*% b2[,g]
  }
  }

  #extraction de Var(beta)
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
  Vbeta[upper.tri(Vbeta, diag=TRUE)] <- x$V[indV] 
  
  
  #IC pour les predictions 
  lower <- matrix(0,nrow(Y),ncol(Y))  
  upper <- matrix(0,nrow(Y),ncol(Y))
  colnames(lower) <- paste("lower.class",1:x$ng,sep="")
  colnames(upper) <- paste("upper.class",1:x$ng,sep="") 
   
  if(x$ng==1)
  {
   varpred <- apply(X1[,-1,drop=FALSE],1,function(x) matrix(x,nrow=1) %*% Vbeta %*% matrix(x,ncol=1)) 
   lower[,1] <- Y[,1] -1.96 * sqrt(varpred)
   upper[,1] <- Y[,1] +1.96 * sqrt(varpred)
  }
  else
  {
   for(g in 1:x$ng)
   {
    ind <- na.omit(c(placeV[["commun"]],placeV[[paste("class",g,sep="")]]))

    if(g==1)
    {
     if(x$idg0[1]==1)
     {
      X12 <- X12 <- cbind(X1[,-1,drop=FALSE],X2)
     }
     
     if(x$idg0[1]==2)
     {
      X12 <- X12 <- cbind(X1,X2[,-1,drop=FALSE])
     }
    }
    else
    {
     X12 <- cbind(X1,X2)    
    }
    
    X12 <- X12[,order(ind),drop=FALSE]

    varclass <- Vbeta[sort(ind)-x$N[1],sort(ind)-x$N[1]]
    varpred <- apply(X12,1,function(x) matrix(x,nrow=1) %*% varclass %*% matrix(x,ncol=1)) 
    
    lower[,g] <- Y[,g] -1.96 * sqrt(varpred)
    upper[,g] <- Y[,g] +1.96 * sqrt(varpred)    
   }


}  


  if(confint==TRUE)
      {
          res <- cbind(Y,lower,upper)
          if(x$ng==1) colnames(res) <- c("pred","lower.pred","upper.pred")
          if(x$ng>1) colnames(res) <- c(paste("pred_class",1:x$ng,sep=""),paste("lower.pred_class",1:x$ng,sep=""),paste("upper.pred_class",1:x$ng,sep=""))

          res.list <- NULL
          res.list$pred <- res
          res.list$times <- times
      }
  if(confint==FALSE)
      {
          if(x$ng==1) colnames(Y) <- "pred"
          if(x$ng>1) colnames(Y) <- paste("pred_class",1:x$ng,sep="")

          res.list <- NULL
          res.list$pred <- Y
          res.list$times <- times
      }

}
else{
cat("Output can not be produced since the program stopped abnormally.")
 res.list <- list(pred=NA,times=NA)
}

  class(res.list) <- "predictL"
  return(res.list)

}         


predictL <- function(x,newdata,var.time,na.action=1,confint=FALSE,...) UseMethod("predictL")
