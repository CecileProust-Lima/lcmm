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
  
  if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")
  

 
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
   times <- times[-linesNA,,drop=FALSE]
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
  
  
  #newdata1 <- model.matrix(as.formula(paste("~",paste(x$Xnames[-1],collapse="+"))),data=newdata)
  
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
