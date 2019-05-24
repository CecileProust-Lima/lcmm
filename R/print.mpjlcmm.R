#' @export
#'
print.mpjlcmm <- function(x,...)
{
 cat("Multivariate joint latent class model for quantitative outcome and competing risks", "\n")
 cat("     fitted by maximum likelihood method", "\n")

 cl <- x$call
 cl$B <- NULL
 if(is.data.frame(cl$data))
     {
         cl$data <- NULL
         x$call$data <- NULL    
     }
 cat(" \n")
 dput(cl)
 cat(" \n")

 posfix <- eval(cl$posfix)
 nbevt <- x$nbevt
 K <- x$K
 
 cat("Statistical Model:", "\n")
 cat(paste("     Dataset:", x$call$data),"\n")
 cat(paste("     Number of subjects:", x$ns),"\n")
 cat(paste("     Number of longitudinal models:", x$K),"\n")
 cat(paste("     Number of observations:", paste(x$N[11+1:K],collapse=" ")),"\n")
 cat(paste("     Number of latent classes:", x$ng), "\n")
 cat(paste("     Number of parameters:", length(x$best))," \n")
 if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")

 if(nbevt>0)
 {
     for(ke in 1:nbevt)
     {

         cat(paste("     Event",ke,": \n"))
         cat(paste("        Number of events: ", x$N[11+K+ke],"\n",sep=""))
         if(x$ng>1)
         {
             if (x$hazardtype[ke]=="Specific") cat("        Class-specific hazards and \n")
             if (x$hazardtype[ke]=="PH") cat("        Proportional hazards over latent classes and \n")
             if (x$hazardtype[ke]=="Common") cat("        Common hazards over classes and \n")
         }
         
         if (x$typrisq[ke]==2)
         {
             cat("        Weibull baseline risk function \n")
         }
         if (x$typrisq[ke]==1)
         {
             cat("        Piecewise constant baseline risk function with nodes \n")
             cat("        ",x$hazardnodes[1:x$nz[ke],ke]," \n")
         }
         if (x$typrisq[ke]==3)
         {
             cat("        M-splines constant baseline risk function with nodes \n")
             cat("        ",x$hazardnodes[1:x$nz[ke],ke]," \n")
         }
         
         
     }
 }

 Ynames <- x$Names$Yname
 if(any(x$linktype!=-1))
 {
     cat("     Link functions: ")
     for (yk in 1:sum(x$ny))
     {
         if (x$linktype[yk]==0)
         {
             if (yk>1) cat("                     ")
             cat("Linear for",Ynames[yk]," \n")
         }
         if (x$linktype[yk]==1)
         {
             if (yk>1) cat("                     ")
             cat("Standardised Beta CdF for",Ynames[yk]," \n")
         }
         if (x$linktype[yk]==2)
         {
             if (yk>1) cat("                     ")
             cat("Quadratic I-splines with nodes", x$linknodes[1:x$nbzitr[yk],yk]," for ",Ynames[yk], "\n")
         }
     }
 }
 
 cat(" \n")
 cat("Iteration process:", "\n")

 if(x$conv==1) cat("     Convergence criteria satisfied")
 if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
 if(x$conv==3) cat("     Convergence with restrained Hessian matrix")
 if(x$conv==4|x$conv==12)
 {
  cat("     The program stopped abnormally. No results can be displayed.\n")
 }
 else
 {
  cat(" \n")
  cat("     Number of iterations: ", x$niter, "\n")
  cat("     Convergence criteria: parameters=", signif(x$gconv[1],2), "\n")
  cat("                         : likelihood=", signif(x$gconv[2],2), "\n")
  cat("                         : second derivatives=", signif(x$gconv[3],2), "\n")
  cat(" \n")
  cat("Goodness-of-fit statistics:", "\n")
  cat(paste("     maximum log-likelihood:", round(x$loglik,2))," \n")
  cat(paste("     AIC:", round(x$AIC,2))," \n")
  cat(paste("     BIC:", round(x$BIC,2))," \n")
  cat(" \n")
 }

}
