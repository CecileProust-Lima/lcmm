
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

    ##pour les facteurs
  olddata <- eval(x$call$data)
  termes <- x$Xnames[-1]
  
   #cas ou une variable dans le dataset du modele est un facteur
   for(v in x$Xnames2[-1])
   {                                                          
    if (is.factor(olddata[,v]))
    {
     mod <- levels(olddata[,v])
     if (!(levels(as.factor(newdata[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     newdata[,v] <- factor(newdata[,v], levels=mod)
  
     for(m in mod)
     {
      termes <- gsub(paste(v,m,sep=""),v,termes)
     }
     
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
     
     factorv <- paste("factor(",v,")",sep="")
     for(m in mod)
     {
      termes <- gsub(paste(factorv,m,sep=""),v,termes,fixed=TRUE)
     }
    }
   }



  ## pour poly()
  if(any(grep("poly",termes)))
      {
          split_termes <- strsplit(termes,split="")
          for(j in 1:length(split_termes))
              {
                  splitj <- split_termes[[j]]
                  
                  if(paste(splitj[c(1:5,length(splitj)-1)],collapse="")=="poly()")
                      {
                          termes[j] <- paste(splitj[1:(length(splitj)-1)],collapse="")
                      }

              }
              
      }


  
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
  

  
  
  ###matrice avec toutes les var et toutes les interactions
  newdata1 <- model.matrix(as.formula(paste("~",paste(termes,collapse="+"))),data=newdata)
  #remettre les termes dans le bon ordre
   Xnames <- x$Xnames[-1]
   z <- grep("factor\\(",Xnames) 
   if (length(z))
   {
    Xnames <- gsub("factor\\(","",Xnames)
    Xnames[z] <- gsub("\\)","",Xnames[z])
   }
   newdata1 <- newdata1[,c("(Intercept)",Xnames),drop=FALSE]
  
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

