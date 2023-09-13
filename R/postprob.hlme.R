#' @export
postprob.hlme <- function(x,threshold=c(0.7,0.8,0.9),...)
{
 if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
 if(x$ng==1)  
 {
  cat("Postprob function can only be used when ng > 1 \n")
  res <- list()
 } 
 else 
 {
  classif <- NULL
  cl.table <- NULL
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




#' Posterior classification stemmed from a \code{hlme}, \code{lcmm},
#' \code{multlcmm} or \code{Jointlcmm} estimation
#' 
#' This function provides informations about the posterior classification
#' stemmed from a \code{hlme}, \code{lcmm}, \code{multlcmm}, \code{Jointlcmm},
#' \code{mpjlcmm}, \code{externSurv} or \code{externX} object.
#' 
#' This function provides the number of subjects classified a posteriori in
#' each latent class, the percentage of subjects classified with a posterior
#' probability above a certain threshold, and the classification table that
#' contains the mean of the posterior probability of belonging to each latent
#' class over the subjects classified in each of the latent classes. This table
#' aims at evaluating the quality of the posterior classification.  For
#' \code{hlme}, \code{lcmm} objects, the posterior classification and the
#' classification table are derived from the posterior class-membership
#' probabilities given the vector of repeated measures that are contained in
#' pprob output matrix. For a \code{Jointlcmm} object, the first posterior
#' classification and the classification table are derived from the posterior
#' class-membership probabilities given the vector of repeated measures and the
#' time-to-event information (that are contained in columns probYT1, probYT2,
#' etc in pprob output matrix). The second posterior classification is derived
#' from the posterior class-membership probabilities given only the vector of
#' repeated measures (that are contained in columns probY1, probY2, etc in
#' pprob output matrix).
#' 
#' @aliases postprob postprob.hlme postprob.lcmm postprob.Jointlcmm
#' postprob.multlcmm postprob.mpjlcmm postprob.externSurv postprob.externX
#' @param x an object inheriting from classes \code{hlme}, \code{lcmm},
#' \code{Jointlcmm} or \code{multlcmm} representing respectively a fitted
#' latent class linear mixed-effects model, a more general latent class mixed
#' model, a joint latent class model or a multivariate general latent class
#' mixed model.
#' @param threshold optional vector of thresholds for the posterior
#' probabilities
#' @param \dots further arguments to be passed to or from other methods.  They
#' are ignored in this function.
#' @return A list containing the posterior classification, the posterior
#' classification table and the percentage of subjects classified with a
#' posterior probability above the given thresholds.
#' @note This function can only be used with latent class mixed models and
#' joint latent class mixed models that include at least 2 latent classes
#' @author Cecile Proust-Lima, Benoit Liquet and Viviane Philipps
#' @seealso \code{\link{Jointlcmm}}, \code{\link{lcmm}}, 
#' \code{\link{hlme}}, \code{\link{plot.lcmm}}
#' @examples
#' 
#' 
#' m<-lcmm(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
#' subject='ID',ng=2,data=data_hlme,B=c(0.41,0.55,-0.18,-0.41,
#' -14.26,-0.34,1.33,13.51,24.65,2.98,1.18,26.26,0.97))
#' postprob(m)
#' 
#' 
#' @export
#' 
postprob <- function(x,threshold=c(0.7,0.8,0.9),...) UseMethod("postprob")
