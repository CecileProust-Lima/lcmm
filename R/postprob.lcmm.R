postprob.lcmm <- function(x,threshold=c(0.7,0.8,0.9),...)
{
 if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")
 if(x$ng==1)
 { 
  cat("Postprob function can only be used when ng > 1 \n")
  res <- list()
 } 
 else 
 {
  classif<-NULL
  cl.table<-NULL
  thr.table <- NULL 
  for (g in 1:x$ng) 
  {
   temp<- subset(x$pprob,x$pprob[,2]==g)
   temp1<-apply(temp[,3:(x$ng+2)],2,mean)
   cl.table<-rbind(cl.table,temp1)
   classif<-cbind(classif,c(as.integer(length(temp[,2])),as.double((length(temp[,2])/x$ns*100))))
   if(!is.null(threshold)) thr.table <- cbind(thr.table,sapply(threshold,function(x) length(which(temp[,2+g]>x)))/length(temp[,1]))   
  }
  
  rownames(cl.table)<-paste("class",1:x$ng,sep="")
  colnames(cl.table)<-paste("prob",1:x$ng,sep="")
  colnames(classif)<-paste("class",1:x$ng,sep="")
  rownames(classif)<-c("N","%")
  
  res <- list(round(classif,2),round(cl.table,4))

  if(!is.null(thr.table))
  {
   thr.table <- 100*thr.table  
    
   rownames(thr.table) <- paste("prob>",threshold,sep="") 
   colnames(thr.table) <- paste("class",1:x$ng,sep="") 
    
   res <- list(round(classif,2),round(cl.table,4),round(thr.table,2))  
  }   
   
  cat(" \n")
  cat("Posterior classification:", "\n")
  print(round(classif,2))
  cat(" \n")
  
  cat("Posterior classification table:", "\n")
  cat("     --> mean of posterior probabilities in each class", "\n")
  print(round(cl.table,4))
  cat(" \n")
   
  if(!is.null(thr.table))
  {
   cat("Posterior probabilities above a threshold (%):", "\n")
   print(round(thr.table,2))
   cat(" \n")    
  }   
 }
 
 return(invisible(res))
}


postprob <- function(x,threshold=c(0.7,0.8,0.9),...) UseMethod("postprob")
