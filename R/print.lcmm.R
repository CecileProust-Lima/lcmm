#' Brief summary of a \code{hlme}, \code{lcmm},
#' \code{Jointlcmm},\code{multlcmm}, \code{epoce} or \code{Diffepoce} objects
#' 
#' The function provides a brief summary of \code{hlme},
#' \code{lcmm},\code{multlcmm} or \code{Jointlcmm} estimations, and
#' \code{epoce} or \code{Diffepoce} computations.
#' 
#' 
#' @aliases print.hlme print.lcmm print.Jointlcmm print.multlcmm print.epoce
#' print.Diffepoce print.externSurv print.externX print.mpjlcmm
#' @param x an object inheriting from classes \code{hlme}, \code{lcmm}, \code{multlcmm}
#' for fitted latent class mixed-effects, or class \code{Jointlcmm}, code{mpjclmm} for
#' a Joint latent class mixed model or \code{epoce} for predictive accuracy
#' computations or \code{externSurv}, \code{externX} for secondary regression models.
#' @param \dots further arguments to be passed to or from other methods.  They
#' are ignored in this function.
#' @author Cecile Proust-Lima, Viviane Philipps, Amadou Diakite and Benoit
#' Liquet
#' @seealso \code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{Jointlcmm}},
#' \code{epoce}, \code{Diffepoce}
#' @keywords print
#' 
#' @export
#'
print.lcmm <-
function(x,...){
if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")

if(inherits(x, "externVar")){
  cat("Secondary linear mixed model", "\n")
} else {
  cat("General latent class mixed model", "\n")
}
cat("     fitted by maximum likelihood method", "\n")
if(inherits(x, "externVar")){
  if(x$varest == "none") cat("     primary model variance not accounted for", "\n")
  if(x$varest == "Hessian") cat("     primary model variance accounted for through the hessian of the joint likelihood", "\n")
  if(x$varest == "paramBoot") cat("     primary model variance accounted for through parametric boostrap", "\n")
}

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
cat(paste("     Number of observations:", x$N[5]),"\n")
if(length(x$na.action))cat(paste("     Number of observations deleted:",length(x$na.action)),"\n")
cat(paste("     Number of latent classes:", x$ng), "\n")
cat(paste("     Number of parameters:", length(x$best))," \n")
 if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")
if (x$linktype==0) {
ntrtot <- 1
cat("     Link function: linear"," \n")
}
if (x$linktype==1)
{
ntrtot <- 3
cat("     Link function: Standardised Beta CdF"," \n")
}
if (x$linktype==2) {
ntrtot <- length(x$linknodes)+1
cat("     Link function: Quadratic I-splines with nodes"," \n")
cat(     x$linknodes," \n")
}

cat(" \n")
cat("Iteration process:", "\n")

if(x$conv==1) cat("     Convergence criteria satisfied")
if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
if(x$conv==3) cat("     Convergence with restrained Hessian matrix")
if(x$conv==4|x$conv==12) {
cat("     The program stopped abnormally. No results can be displayed.\n")
}else{

cat(" \n")
cat("     Number of iterations: ", x$niter, "\n")
if(inherits(x, "externVar")) {
  if(x$varest == "paramBoot"){
    cat("     Proportion of convergence on bootstrap iterations (%)=", x$Mconv, "\n")
  } else {
    cat("     Number of iterations: ", x$niter, "\n")
    cat("     Convergence criteria: parameters=", signif(x$gconv[1],2), "\n")
    cat("                         : likelihood=", signif(x$gconv[2],2), "\n")
    cat("                         : second derivatives=", signif(x$gconv[3],2), "\n")
  }
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
if (x$Ydiscrete==1&x$N[6]==0){
cat(paste("     Discrete posterior log-likelihood:", round(x$discrete_loglik,2))," \n")
cat(paste("     Discrete AIC:", round(-2*x$discrete_loglik+2*(length(x$best)-length(posfix)),2))," \n")
cat(" \n")

cat(paste("     Mean discrete AIC per subject:",round((-x$discrete_loglik+length(x$best)-length(posfix))/as.double(x$ns),4))," \n")
cat(paste("     Mean UACV per subject:",round(x$UACV,4))," \n")
cat(paste("     Mean discrete LL per subject:",round(x$discrete_loglik/as.double(x$ns),4))," \n")
}
cat(" \n")
}}
