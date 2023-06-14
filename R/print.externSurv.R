#' @export
#'
print.externSurv <- function(x,...){
  if (!inherits(x, "externSurv")) stop("use only with \"externSurv\" objects")
  
  cat("Secondary survival model", "\n")
  cat("     fitted by maximum likelihood method", "\n")
  if(x$varest == "none") cat("     primary model variance not accounted for", "\n")
  if(x$varest == "Hessian") cat("     primary model variance accounted for through the hessian of the joint likelihood", "\n")
  if(x$varest == "paramBoot") cat("     primary model variance accounted for through parametric boostrap", "\n")
  
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
  
  cat("Statistical Model:", "\n")
  cat(paste("     Dataset:", x$call$data),"\n")
  cat(paste("     Number of subjects:", x$ns),"\n")
  cat(paste("     Number of longitudinal models:", x$K),"\n")
  cat(paste("     Number of latent classes:", x$ng), "\n")
  cat(paste("     Number of parameters:", length(x$best))," \n")
  if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")
  
  nbevt <- x$nbevt
  nprisq <- rep(NA,nbevt)
  nrisq <- rep(NA,nbevt)
  
  for(ke in 1:nbevt)
  {
    if(x$typrisq[ke]==1) nprisq[ke] <- x$nz[ke]-1
    if(x$typrisq[ke]==2) nprisq[ke] <- 2
    if(x$typrisq[ke]==3) nprisq[ke] <- x$nz[ke]+2
    
    nrisq[ke] <- x$Nprm[1+ke]
    
    cat(paste("     Event",ke,": \n"))
    cat(paste("        Number of events: ", x$N[2+ke],"\n",sep=""))
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
      if(x$varest == "paramBoot") {
        cat("     Proportion of convergence on bootstrap iterations (%)=", x$Mconv, "\n")
      } else {
        cat("     Number of iterations: ", x$niter, "\n")
        cat("     Convergence criteria: parameters=", signif(x$gconv[1],2), "\n")
        cat("                         : likelihood=", signif(x$gconv[2],2), "\n")
        cat("                         : second derivatives=", signif(x$gconv[3],2), "\n")
      }
      cat(" \n")
      cat("Goodness-of-fit statistics:", "\n")
      cat(paste("     maximum log-likelihood:", round(x$loglik,2))," \n")
      cat(paste("     AIC:", round(x$AIC,2))," \n")
      cat(paste("     BIC:", round(x$BIC,2))," \n")
      cat(" \n")
      
    }
  }
}