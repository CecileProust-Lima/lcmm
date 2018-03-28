print.multlcmm <- function(x,...){
if (!inherits(x, "multlcmm")) stop("use only with \"multlcmm\" objects")

cat("General latent class mixed model", "\n")
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
cat(paste("     Dataset:", x$call$data),"\n")
cat(paste("     Number of subjects:", x$ns),"\n")

cat(paste("     Number of observations:", x$N[9]),"\n")
#if(length(x$linesNA))cat(paste("     Number of observations deleted:",length(x$linesNA)),"\n")
cat(paste("     Number of latent classes:", x$ng), "\n")
cat(paste("     Number of parameters:", length(x$best))," \n")
 if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")

ntrtot <- rep(NA,x$N[8])
numSPL <- 0
cat("     Link functions: ")
for (yk in 1:x$N[8])
{
 if (x$linktype[yk]==0) {
 ntrtot[yk] <- 2
 if (yk>1) cat("                     ")
 cat("linear for",x$Ynames[yk]," \n")
 }
 if (x$linktype[yk]==1)
 {
 ntrtot[yk] <- 4
 if (yk>1) cat("                     ")
 cat("Standardised Beta CdF for",x$Ynames[yk]," \n")
 }
 if (x$linktype[yk]==2) {
 numSPL <- numSPL+1
 ntrtot[yk] <- x$nbnodes[numSPL]+2
 if (yk>1) cat("                     ")
 cat("Quadratic I-splines with nodes", x$linknodes[1:x$nbnodes[numSPL],yk]," for ",x$Ynames[yk], "\n")
 }
}

cat(" \n")
cat("Iteration process:", "\n")

if(x$conv==1) cat("     Convergence criteria satisfied")
if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
if(x$conv==3) cat("     Convergence with restrained Hessian matrix")
if(x$conv==4|x$conv==12) {
cat("     The program stopped abnormally. No results can be displayed.\n")
}
else{

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
