#' @export
#'
predictL.multlcmm <- function(x,newdata,var.time,na.action=1,confint=FALSE,...)
{
 if(missing(newdata)) stop("The argument newdata should be specified")
 if(missing(x)) stop("The argument x should be specified")
 if (!inherits(x, "multlcmm")) stop("use only with \"multlcmm\" objects")
 # ad 2/04/2012 Xnames2
 if (!all(x$Xnames2 %in% colnames(newdata))) {
 stop(paste(c("newdata should at least include the following covariates: ","\n",x$Xnames2),collapse=" "))}
 if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
# if(missing(var.time)) stop("missing argument 'var.time'")
# if(!(var.time %in% colnames(newdata))) stop("'var.time' should be included in newdata")

 
 if(x$conv==1 | x$conv==2 | x$conv==3)
 { 
  if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")

    call_fixed <- x$call$fixed[3]
    if(is.null(x$call$random)) {call_random <- -1} else call_random <- x$call$random[2]
    if(is.null(x$call$classmb)) {call_classmb <- -1} else call_classmb <- x$call$classmb[2]
  if(is.null(x$call$mixture)) {call_mixture <- -1} else call_mixture <- x$call$mixture[2]
    
  if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 
  
  if(na.action==1){
      na.action=na.omit
  }else{
      na.action=na.fail
  }

  ## transform to factor is the variable appears in levels$levelsdata
  for(v in colnames(newdata))
  {
      if(v %in% names(x$levels$levelsdata))
      {
          if(!is.null(x$levels$levelsdata[[v]]))
          {
              newdata[,v] <- factor(newdata[,v], levels=x$levels$levelsdata[[v]])
          }
      }
  }
  
  call_fixed <- gsub("factor","",call_fixed)
  call_fixed <- gsub("contrast","",call_fixed)
  call_random <- gsub("factor","",call_random)
  call_classmb <- gsub("factor","",call_classmb)
  call_mixture <- gsub("factor","",call_mixture) 
  
  call_mixture <- formula(paste("~",call_mixture,sep=""))
  call_random <- formula(paste("~",call_random,sep=""))
  call_classmb <- formula(paste("~",call_classmb,sep=""))

  ## Traitement des donnees manquantes

  mcall <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
  mcall$na.action <- na.action
  mcall$data <- newdata

  ## fixed
  m <- mcall
  m$formula <- formula(paste("~",call_fixed,sep=""))
  m[[1]] <- as.name("model.frame")	
  m <- eval(m, sys.parent()) 
  na.fixed <- attr(m,"na.action")
  ## mixture
  na.mixture <- NULL
  if(call_mixture != -1)
  {
      m <- mcall
      m$formula <- call_mixture
      m[[1]] <- as.name("model.frame")	
      m <- eval(m, sys.parent()) 
      na.mixture <- attr(m,"na.action")
  }
  ## random
  na.random <- NULL
  if(call_random != -1)
  {
      m <- mcall
      m$formula <- call_random
      m[[1]] <- as.name("model.frame")	
      m <- eval(m, sys.parent()) 
      na.random <- attr(m,"na.action")
  }
  ## classmb
  na.classmb <- NULL
  if(call_classmb != -1)
  {
      m <- mcall
      m$formula <- call_classmb
      m[[1]] <- as.name("model.frame")	
      m <- eval(m, sys.parent()) 
      na.classmb <- attr(m,"na.action")
  }
  ## cor
  na.cor <- NULL
  if(x$N[7]>0)
  {
      z <- which(x$idcor0==1)
      var.cor <- newdata[,x$Xnames[z]]
      na.cor <- which(is.na(var.cor))
  }

  ##var.time
  if(!missing( var.time))
  {
      if(!(var.time %in% colnames(newdata))) stop("'var.time' should be included in newdata")
      times <- newdata[,var.time,drop=FALSE]
      
  }
  else
    {
        times <- newdata[,1,drop=FALSE]
    }
  
  ## Table sans donnees manquante: newdata
  na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.cor))
  
  if(length(na.action)){
      newdata1 <- newdata[-na.action,]
      times <- times[-na.action,,drop=FALSE]
  } else {
      newdata1 <- newdata
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
  

  ## Construction de nouvelles var eplicatives sur la nouvelle table
  X_fixed <- X_mixture <- X_random <- X_classmb <- NULL
  ## fixed
  X_fixed <- model.matrix(formula(paste("~",call_fixed,sep="")),data=newdata1fixed)
  if(colnames(X_fixed)[1]=="(Intercept)"){
      colnames(X_fixed)[1] <- "intercept"
  }
  ## mixture
  if(call_mixture != ~-1){
      X_mixture <- model.matrix(call_mixture,data=newdata1mixture)	
      if(colnames(X_mixture)[1]=="(Intercept)"){
          colnames(X_mixture)[1] <- "intercept"
      }
  }            
  ## random
  if(call_random != ~-1){
      X_random <- model.matrix(call_random,data=newdata1random)	
      if(colnames(X_random)[1]=="(Intercept)"){
          colnames(X_random)[1] <- "intercept"
      }
  }	
  ## classmb
  if(call_classmb != ~-1){ 
      X_classmb <- model.matrix(call_classmb,data=newdata1classmb)
      colnames(X_classmb)[1] <- "intercept"
  }
  
  ##cor
  if(x$N[7]>0)  #on reprend la variable de temps de cor (sans NA)
  {
      z <- which(x$idcor0==1)
      var.cor <- newdata1[,x$Xnames[z]]
  }
  
  ## pour mettre les var dans le bon ordre
  newdata1 <- X_fixed
  colX <- colnames(X_fixed)
  
  if(!is.null(X_mixture)){
      for(i in 1:length(colnames(X_mixture))){
          if((colnames(X_mixture)[i] %in% colnames(newdata1))==FALSE){
              newdata1 <- cbind(newdata1,X_mixture[,i])
              colnames(newdata1) <- c(colX,colnames(X_mixture)[i])
              colX <- colnames(newdata1)
          }
      }
  }
  if(!is.null(X_random)){
      for(i in 1:length(colnames(X_random))){
          if((colnames(X_random)[i] %in% colnames(newdata1))==FALSE){
              newdata1 <- cbind(newdata1,X_random[,i])
              colnames(newdata1) <- c(colX,colnames(X_random)[i])
              colX <- colnames(newdata1)
          }	 
      }
  }
  if(!is.null(X_classmb)){
      for(i in 1:length(colnames(X_classmb))){
          if((colnames(X_classmb)[i] %in% colnames(newdata1))==FALSE){
              newdata1 <- cbind(newdata1,X_classmb[,i])
              colnames(newdata1) <- c(colX,colnames(X_classmb)[i])
              colX <- colnames(newdata1)
          }	
      }
  }
  if(x$N[7]>0)
  { 
      if( x$idg0[z]==0 & x$idea0[z]==0 & x$idprob0[z]==0)
      {
          newdata1 <- cbind(newdata1,var.cor)
          colnames(newdata1) <- c(colX,x$Xnames[z])
          colX <- colnames(newdata1)
      }
  }


  
  ###calcul des predictions
     
  X1 <- NULL
  X2 <- NULL
  b1 <- NULL
  b2 <- NULL
  
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
  for(g in 1:x$ng)
  {
   if(length(b1) != 0)
   {
    Y[,g]<- X1 %*% b1
   }
   if(length(b2) != 0)
   {
    Y[,g]<- Y[,g] + X2 %*% b2[,g]
   }
  }

  
  ny <- length(x$Ynames)

  #extraction de Var(beta) 
  Vbeta <- matrix(0,x$N[3]-x$N[2]-x$N[1],x$N[3]-x$N[2]-x$N[1])
  npm <- length(x$best)
  indice <- 1:npm * (1:npm+1) /2
  indtmp <- indice[(x$N[1]+1):(x$N[3]-x$N[2])]
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
   #varpred <- diag(X1[,-1,drop=FALSE] %*% Vbeta %*% t(X1[,-1,drop=FALSE]))
   varpred <- apply(X1[,-1,drop=FALSE],1,function(x) matrix(x,nrow=1) %*% Vbeta %*% matrix(x,ncol=1)) 
     # browser()
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
 else
 {
  cat("Output can not be produced since the program stopped abnormally.")
  res.list <- list(pred=NA,times=NA)
 }

  class(res.list) <- "predictL"
  return(res.list)
}

