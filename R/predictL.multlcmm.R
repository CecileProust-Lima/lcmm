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

  
  ## var.time
  if(!missing( var.time))
      {
          if(!(var.time %in% colnames(newdata))) stop("'var.time' should be included in newdata")
          times <- newdata[,var.time,drop=FALSE]         
      }
  else
      {
          times <- newdata[,1,drop=FALSE]
      }
 
  ### Traitement des donnees manquantes
  newdata <- newdata[,x$Xnames2]
  newdata <- as.data.frame(newdata)
  colnames(newdata) <- x$Xnames2
  
  linesNA <- apply(newdata,2,function(v) which(is.na(v)))
  linesNA <- unique(unlist(linesNA))
  
  if(length(linesNA) & na.action==2) stop("newdata contains missing values")
  if(length(linesNA) & na.action==1) 
  {
   newdata <- as.data.frame(newdata[-linesNA,])
   colnames(newdata) <- x$Xnames2
   times <- times[-linesNA]
  }
  

  fixed <- as.formula(x$call$fixed)
  random <- as.formula(x$call$random)
  mixture <- as.formula(x$call$mixture)
  classmb <- as.formula(x$call$classmb)
  
  if(x$ng==1) classmb <- ~-1
  if(x$ng>1 & is.null(x$call$classmb)) classmb <- ~1
  if(x$ng==1) mixture <- ~-1
  
  ##pour les facteurs
  ##donnees de l estimation
  if(!is.null(x$data))
  {
    olddata <- x$data
  }
  else
  {
    olddata <- eval(x$call$data)
  }
  
  #cas ou une variable dans le dataset du modele est un facteur
  for(v in x$Xnames2)
  {
    if (is.factor(olddata[,v]))
    {
      mod <- levels(olddata[,v])
      if (!(levels(as.factor(newdata[,v])) %in% mod)) stop(paste("invalid level in factor", v))
      newdata[,v] <- factor(newdata[,v], levels=mod)
      
    }
  }
  
  #cas ou on a factor() dans l'appel de la fonction
  dans_appel <- c(all.names(x$call$fixed),all.names(x$call$random),all.names(x$call$mixture),all.names(x$call$classmb))
  ind_factor <- which(dans_appel=="factor")
  if(length(ind_factor))
  {
    nom.factor <- dans_appel[ind_factor+1]
    for (v in nom.factor)
    {
      mod <- levels(as.factor(olddata[,v]))
      if (!all(levels(as.factor(newdata[,v])) %in% mod)) stop(paste("invalid level in factor", v))
      newdata[,v] <- factor(newdata[,v], levels=mod)
    }
  }
  fixed <- gsub("factor","",fixed)
  
  z <- all.names(random)
  ind_factor <- which(z=="factor")
  if(length(ind_factor))
  {
    nom.factor <- z[ind_factor+1]
    for (v in nom.factor)
    {
      mod <- levels(as.factor(olddata[,v]))
      if (!all(levels(as.factor(newdata[,v])) %in% mod)) stop(paste("invalid level in factor", v))
      newdata[,v] <- factor(newdata[,v], levels=mod)
    }
  }
  random <- reformulate(gsub("factor","",random))
  
  z <- all.names(classmb)
  ind_factor <- which(z=="factor")
  if(length(ind_factor))
  {
    nom.factor <- z[ind_factor+1]
    for (v in nom.factor)
    {
      mod <- levels(as.factor(olddata[,v]))
      if (!all(levels(as.factor(newdata[,v])) %in% mod)) stop(paste("invalid level in factor", v))
      newdata[,v] <- factor(newdata[,v], levels=mod)
    }
  }
  classmb <- reformulate(gsub("factor","",classmb))
  
  z <- all.names(mixture)
  ind_factor <- which(z=="factor")
  if(length(ind_factor))
  {
    nom.factor <- z[ind_factor+1]
    for (v in nom.factor)
    {
      mod <- levels(as.factor(olddata[,v]))
      if (!all(levels(as.factor(newdata[,v])) %in% mod)) stop(paste("invalid level in factor", v))
      newdata[,v] <- factor(newdata[,v], levels=mod)
    }
  }
  mixture <- reformulate(gsub("factor","",mixture))
  
  
  
  afixed <- terms(reformulate(fixed), specials=c("factor","contrast"))
  
  
  ##fixed sans contrast
  fixed2 <- gsub("contrast","",fixed)
  fixed2 <- formula(paste(fixed2[2],fixed2[3],sep="~")) 
  
  ##contrast
  contr <- ~-1
  if(!is.null(attr(afixed,"specials")$contrast))
  {
    vcontr <- attr(afixed,"term.labels")[setdiff(attr(afixed,"specials")$contrast-1,untangle.specials(afixed,"contrast",2)$terms)]
    vcontr <- gsub("contrast","",vcontr)
    contr <- as.formula(paste("~-1+",paste(vcontr,collapse="+")))
  }
  
  ###creation de X0 ( var + interactions)
  
  Xfixed <- model.matrix(fixed2[-2], data=newdata)
  Xmixture <- model.matrix(mixture, data=newdata)
  Xrandom <- model.matrix(random, data=newdata)
  Xclassmb <- model.matrix(classmb, data=newdata)
  Xcontr <- model.matrix(contr,data=newdata)
  
  z.fixed <- strsplit(colnames(Xfixed),split=":",fixed=TRUE)
  z.fixed <- lapply(z.fixed,sort)
  
  z.random <- strsplit(colnames(Xrandom),split=":",fixed=TRUE)
  z.random <- lapply(z.random,sort)
  
  #if(mixture != ~-1)
  if(ncol(Xmixture)>0)
  {
    z.mixture <- strsplit(colnames(Xmixture),split=":",fixed=TRUE)
    z.mixture <- lapply(z.mixture,sort)
  }
  else
  {
    z.mixture <- list()
  }
  
  #if(classmb != ~-1)
  if(ncol(Xclassmb)>0)
  {
    z.classmb <- strsplit(colnames(Xclassmb),split=":",fixed=TRUE)
    z.classmb <- lapply(z.classmb,sort)
  }
  else
  {
    z.classmb <- list()
  }
  
  if(contr != ~-1)
  {
    z.contr <- strsplit(colnames(Xcontr),split=":",fixed=TRUE)
    z.contr <- lapply(z.contr,sort)
  }
  else
  {
    z.contr <- list()
  }
  
  X0 <- cbind(Xfixed, Xrandom, Xclassmb)        
  nom.unique <- unique(colnames(X0))
  X0 <- X0[,nom.unique,drop=FALSE]
  
  if (x$N[7]>0)
  {
    cor.var.time <- x$Xnames[which(x$idcor0==1)]
    if(!(cor.var.time %in% colnames(X0)))
    {
      X0 <- cbind(X0, newdata[,cor.var.time])
      colnames(X0) <- c(nom.unique, cor.var.time)
      nom.unique <- c(nom.unique,cor.var.time)
    }
  }
  
  newdata1 <- as.matrix(X0)
  
  
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

