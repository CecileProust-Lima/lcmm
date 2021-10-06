#' @rdname predictY
#' @export
#'
predictY.multlcmm <- function(x,newdata,var.time,methInteg=0,nsim=20,draws=FALSE,ndraws=2000,na.action=1,...)
{
if(missing(newdata)) stop("The argument newdata should be specified")
if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "multlcmm")) stop("use only with \"lcmm\" or \"multlcmm\" objects")
if (!all(x$Xnames2 %in% colnames(newdata))) stop(paste(c("newdata should at least include the following covariates: ","\n",x$Xnames2),collapse=" "))
if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
#if(missing(var.time)) stop("missing argument 'var.time'")
#if(!(var.time %in% colnames(newdata))) stop("'var.time' should be included in newdata")

if(any(x$linktype==3))
{
     if(methInteg==0) stop("predictions for ordinal outcomes are only available with MC method. Please use methInteg=1 and set nsim argument.")
#     if(any(x$linktype!=3)) stop("predictions for mixed outcomes (ordinal and continuous) are not available yet.")
#     if(any(x$idcor>0)) stop("predictions with BM or AR correlations are not available yet.")
}

if(x$conv==1 | x$conv==2 | x$conv==3) 
{
  if(x$conv>1 & draws==TRUE)
  {
   cat("No confidence interval will be provided since the program did not converge properly \n")
   draws <- FALSE
  }
  
  if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")

    call_fixed <- x$call$fixed[3]
    if(is.null(x$call$random)) {call_random <- -1} else call_random <- x$call$random[2]
    if(is.null(x$call$classmb)) {call_classmb <- -1} else call_classmb <- x$call$classmb[2]
    if(is.null(x$call$mixture)) {call_mixture <- -1} else call_mixture <- x$call$mixture[2]

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

  
  #arguments pour predictMult
   maxmes <- dim(newdata1)[1]
   X0 <- as.vector(newdata1)
   ncor <- x$N[7]
   nalea <- x$N[6]
   nv <- dim(newdata1)[2]                          
   ny <- x$N[8]
   best <- x$best
   nef <- x$N[3]
   nvc <- x$N[4]
   if(nvc>0) 
   {
    if(x$idiag==0) best[nef+1:nvc] <- x$cholesky[-1]
    else best[nef+1:nvc] <- x$cholesky[c((1:(nvc+1)*2:(nvc+2))/2)[-1]]
   }
   npm <- length(best)
   nbzitr <- rep(2,ny)
   nbzitr[which(x$linktype==2)] <- x$nbnodes
   Ymarg <- matrix(0,maxmes*ny,x$ng)
   
   #if(verbose==TRUE) print(head(newdata1))
  
  if(!draws)
  { 
   #if(verbose==TRUE) cat("ny=",ny,"nvc=",nvc,"ncontr=",x$N[2],"nv=",nv,"\n",x$Xnames,"\n",head(newdata1),"\n","idg",x$idg0)
   
      out <- .Fortran(C_predictmult,
                      as.double(X0),
                      as.integer(x$idprob0),
                      as.integer(x$idea0),
                      as.integer(x$idg0),
                      as.integer(x$idcor0),
                      as.integer(x$idcontr0),
                      as.integer(x$ng),
                      as.integer(ncor),
                      as.integer(nalea),
                      as.integer(nv),
                      as.integer(ny),
                      as.integer(maxmes),
                      as.integer(x$idiag),
                      as.integer(x$N[5]),
                      as.integer(npm),
                      as.double(best),
                      as.double(x$epsY),
                      as.integer(x$linktype),
                      as.integer(nbzitr),
                      as.double(x$linknodes),
                      as.integer(unlist(x$modalites)),
                      as.integer(x$nbmod),
                      as.integer(nsim),
                      as.integer(methInteg),
                      Ymarg=as.double(Ymarg))
  
   out$Ymarg[out$Ymarg==9999] <- NA
   #Ypred <- matrix(out$Ymarg,ncol=x$ng,byrow=FALSE)
   Ypred <- data.frame(rep(x$Ynames,each=maxmes),matrix(out$Ymarg,ncol=x$ng,byrow=FALSE))
  
   if (x$ng==1) colnames(Ypred) <- c("Yname","Ypred")
   if (x$ng>1) colnames(Ypred) <- c("Yname",paste("Ypred_class",1:x$ng,sep=""))
  
  }
  else  #draws
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
   
  

    
   for (j in 1:ndraws)
   {   #cat("boucle sur ndraws j=",j,"\n")
    bdraw <- rnorm(npm)
    bdraw <- best + Chol %*% bdraw
   
    out <- .Fortran(C_predictmult,
                    as.double(X0),
                    as.integer(x$idprob0),
                    as.integer(x$idea0),
                    as.integer(x$idg0),
                    as.integer(x$idcor0),
                    as.integer(x$idcontr0),
                    as.integer(x$ng),
                    as.integer(ncor),
                    as.integer(nalea),
                    as.integer(nv),
                    as.integer(ny),
                    as.integer(maxmes),
                    as.integer(x$idiag),
                    as.integer(x$N[5]),
                    as.integer(npm),
                    as.double(bdraw),
                    as.double(x$epsY),
                    as.integer(x$linktype),
                    as.integer(nbzitr),
                    as.double(x$linknodes),
                    as.integer(unlist(x$modalites)),
                    as.integer(x$nbmod),
                    as.integer(nsim),
                    as.integer(methInteg),
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
  
   #Ypred <- cbind(Ypred_50,Ypred_2.5,Ypred_97.5)
   Ypred <- data.frame(rep(x$Ynames,each=maxmes),Ypred_50,Ypred_2.5,Ypred_97.5)
   
   if (x$ng==1) colnames(Ypred) <- c("Yname","Ypred_50","Ypred_2.5","Ypred_97.5")
   if (x$ng>1) colnames(Ypred) <- c("Yname",c(paste("Ypred_50_class",1:x$ng,sep=""),paste("Ypred_2.5_class",1:x$ng,sep=""),paste("Ypred_97.5_class",1:x$ng,sep="")))

  }

  res.list <- NULL
  res.list$pred <- Ypred
  res.list$times <- times
}

else
{
 cat(" The program stopped abnormally. No prediction can be computed.\n")

 res.list <- list(pred=NA,times=NA)
}


 class(res.list) <- "predictY"
 return(res.list)
}
