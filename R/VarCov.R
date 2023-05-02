#' Variance-covariance of the estimates
#' 
#' This function provides the variance-covariance matrix of the estimates. vcov
#' is an alias for it.
#' 
#' 
#' @param x an object of class \code{hlme}, \code{lcmm}, \code{multlcmm},
#' \code{Jointlcmm} or \code{mpjlcmm}
#' @return a matrix containing the variance-covariance of the estimates. For
#' the parameters of the matrix of variance-covariance of the random effects,
#' the Cholesky transformed parameters are considered so that VarCov provides
#' the covariance matrix of function \code{estimates} with cholesky=TRUE.
#' @author Cecile Proust-Lima, Viviane Philipps
#' @seealso \code{\link{estimates}}
#' 
#' @export
#' 
VarCov <- function(x)
{
 if(missing(x)) stop("The argument x should be specified")
 if (!inherits(x, c("hlme","lcmm","multlcmm","Jointlcmm","mpjlcmm","externX","externSurv"))) stop("use only with hlme, lcmm, multlcmm, Jointlcmm, mpjlcmm, externX or externSurv objects")

 if(x$conv %in% c(1,2,3))
 {
  res <- matrix(0,length(x$best),length(x$best))
  res[upper.tri(res,diag=TRUE)] <- x$V
  res <- t(res)
  res[upper.tri(res,diag=TRUE)] <- x$V

  noms <- sub("varcov","cholesky",names(x$best))
  colnames(res) <- noms
  rownames(res) <- noms
 }
 else
 {
  res <- NA
  cat("Output can not be produced since the program stopped abnormally. \n")
 }

 return(res)
}
