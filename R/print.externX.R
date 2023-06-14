#' @export
#'
print.externX = function(x, ...){
  if(!inherits(x, "externX")) stop("use only with\"externX\" objects")
  
  cat("Secondary multinomial model for external class predictor", "\n")
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
  cat(paste("     Number of latent classes:", x$ng), "\n")
  cat(paste("     Number of parameters:", length(x$best))," \n")
  
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