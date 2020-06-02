.predYmedian <- function(m,newdata,var.time,nsim=1000,ndraws=0,cl=NULL,seed=NULL,...)
{
  ## resultat au format predictY
  res.list <- predictY(m,newdata=newdata,var.time=var.time,methInteg=1,nsim=1,draws=as.logical(ndraws>0),ndraws=1)
  
  ## avec les IC ##
  if(ndraws>0)
  {
    ##graines
    if(!is.null(seed))
    {
      if(length(seed)!=ndraws) stop("seed should be of length ndraws")
    }
    else
    {
      seed <- c(1:ndraws)
    }
    
    ## fonction qui donne la mediane pour 1 draw 
    doone <- function(m,newdata,nsim,seed){
      f <- function(s){
        set.seed(s)
        as.vector(predictY(m,newdata=newdata,var.time=var.time,methInteg=1,nsim=1,draws=TRUE,ndraws=1)$pred[,1:m$ng])
      }
      predb <- replicate(nsim,f(seed)) 
      medb <- apply(predb,1,median,na.rm=TRUE)
      return(medb)
    }
    
    ## en parallele
    if(!is.null(cl))
    {
      if(!inherits(cl,"cluster"))
      {
        if(!is.numeric(cl)) stop("argument cl should be either a cluster or a numeric value indicating the number of cores")
        
        ncl <- cl
        cl <- makeCluster(ncl)
      }
      
      ## export other arguments
      clusterExport(cl, list("m", "newdata", "nsim"), envir = environment())
      
      ## get and export loaded packages
      pck <- .packages()
      dir0 <- find.package()
      dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
      clusterExport(cl,list("pck","dir"),envir=environment())
      clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))
      
      ## faire les ndraws repliques
      allpred <- parSapply(cl,seed,doone,m=m,newdata=newdata,nsim=nsim)
      
      if(!is.null(ncl)) stopCluster(cl)
    }
    else ## en sequentiel
    {
      allpred <- sapply(seed,doone,m=m,newdata=newdata,nsim=nsim)
    }
    
    ## prendre les quantiles pour avoir l'IC
    med_ic <- apply(allpred,1,quantile,probs=c(0.5,0.025,0.975),na.rm=TRUE)
    
    ## remettre en matrice avec ng colonnes
    pred <- matrix(as.vector(t(med_ic)),nrow=nrow(newdata),ncol=3*m$ng)
  }
  else
  {
    ### sans IC ###
    
    ## fonction pour faire 1 simu 
    onesim <- function(m,newdata){
      as.vector(predictY(m,newdata=newdata,var.time=var.time,methInteg=1,nsim=1,draws=FALSE)$pred[,1:m$ng])
    }
    
    ## faire les nsim simus
    pred <- replicate(nsim,onesim(m=m,newdata=newdata)) 
    
    ## mediane, en supprimant les NA ou NaN
    med <- apply(pred,1,median,na.rm=TRUE)
    
    ## remettre en matrice avec ng colonnes
    pred <- matrix(med,nrow=nrow(newdata),ncol=m$ng)
  }
  
  colnames(pred) <- colnames(res.list$pred)  
  res.list$pred <- pred
  
  return(res.list)
}


