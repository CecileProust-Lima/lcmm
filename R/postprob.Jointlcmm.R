postprob.Jointlcmm <- function(x,threshold=c(0.7,0.8,0.9),...)
{
 if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
 
 res <- list()
 
 if(x$ng==1) 
 { 
  cat("Postprob function can only be used when ng > 1 \n")
 } 
 else 
 {
  if(x$conv==1|x$conv==2|x$conv==3) 
  {  
    classif <- NULL
    classifY <- NULL
    cl.table <- NULL
    thr.table <- NULL
    pprob <- x$pprob[,-1]
    pprobY <- x$pprobY[,-1]
    for (g in 1:x$ng) 
    {
    temp<- subset(pprob,pprob[,1]==g)
    tempY<- subset(pprobY,pprobY[,1]==g)
    temp1<-apply(temp[,1+1:x$ng],2,mean)
    temp1Y<-apply(tempY[,1+1:x$ng],2,mean)
    cl.table<-rbind(cl.table,temp1)
    classif<-cbind(classif,length(temp[,1]))
    classifY<-cbind(classifY,length(tempY[,1]))
    if(!is.null(threshold)) thr.table <- cbind(thr.table,sapply(threshold,function(x) length(which(temp[,1+g]>x)))/length(temp[,1]))  
   }
   
   classif <- rbind(classif,100*classif/x$ns)
   classifY <- rbind(classifY,100*classifY/x$ns)
   
   rownames(cl.table)<-paste("class",1:x$ng,sep="")
   colnames(cl.table)<-paste("prob",1:x$ng,sep="")
   colnames(classif)<-paste("class",1:x$ng,sep="")
   rownames(classif)<-c("N","%")
   colnames(classifY)<-paste("class",1:x$ng,sep="")
   rownames(classifY)<-c("N","%")
    
   if(!is.null(thr.table))
   {
    thr.table <- 100*thr.table  
    
    rownames(thr.table) <- paste("prob>",threshold,sep="") 
    colnames(thr.table) <- paste("class",1:x$ng,sep="") 
   }
    
   if(sum(is.na(pprob[1+1:x$ng]))==0)
   {
    cat(" \n")
    cat("Posterior classification based on longitudinal and time-to-event data:", "\n")
    print(round(classif,2))
    cat(" \n")
    
    cat("Posterior classification table:", "\n")
    cat("     --> mean of posterior probabilities in each class", "\n")
    print(round(cl.table,4))
    #print(cl.table)
    cat(" \n")
    
    if(!is.null(thr.table))
    {
     cat("Posterior probabilities above a threshold (%):", "\n")
     print(round(thr.table,2))
     cat(" \n")    
    }

    res <- c(res,list(round(classif,2),round(cl.table,4)))
    if(!is.null(thr.table)) res <- c(res,list(round(thr.table,2)))
   }
   else
   {
    cat("Error in the computation of posterior class-membership probabilities given all the information")
   }
   
   
   if(sum(is.na(pprobY[,1+1:x$ng]))==0)
   {
    cat(" \n")
    cat("Posterior classification based only on longitudinal data:", "\n")
    print(round(classifY,2))
    cat(" \n")
    
    res <- c(res,list(round(classifY,2)))
   }
   else
   {
    cat("Error in the computation of posterior class-membership probabilities given the repeated measures of the marker")
   }
  }
  else
  {
   cat("Output can not be produced since the program stopped abnormally.")
  }
 }
 
 return(invisible(res))
}


postprob <- function(x,threshold=c(0.7,0.8,0.9),...) UseMethod("postprob")
