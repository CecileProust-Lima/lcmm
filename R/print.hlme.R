print.hlme <- function(x,...){
if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")

#cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")


cat("Heterogenous linear mixed model", "\n")
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

cat("Statistical Model:", "\n")
cat(paste("     Dataset:", x$dataset),"\n")
cat(paste("     Number of subjects:", x$ns),"\n")
cat(paste("     Number of observations:", length(x$pred[,1])),"\n")
if(length(x$na.action))cat(paste("     Number of observations deleted:",length(x$na.action)),"\n")
cat(paste("     Number of latent classes:", x$ng), "\n")
cat(paste("     Number of parameters:", length(x$best))," \n")
 if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")
cat(" \n")
cat("Iteration process:", "\n")

if(x$conv==1) cat("     Convergence criteria satisfied")
if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
if(x$conv==4|x$conv==5)  {
cat("     The program stopped abnormally. No results can be displayed.\n")
}else{



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
cat(" \n")

}
}
