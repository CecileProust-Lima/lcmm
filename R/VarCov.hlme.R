#' @export
#'
VarCov.hlme <- function(x)
{
 if(missing(x)) stop("The argument x should be specified")
 if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")

 if(x$conv==1 | x$conv==2)
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



#' Variance-covariance of the estimates
#' 
#' This function provides the variance-covariance matrix of the estimates. vcov
#' is an alias for it.
#' 
#' 
#' @aliases VarCov VarCov.hlme VarCov.lcmm VarCov.Jointlcmm VarCov.multlcmm
#' @param x an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or
#' \code{Jointlcmm}
#' @return a matrix containing the variance-covariance of the estimates. For
#' the parameters of the matrix of variance-covariance of the random effects,
#' the Cholesky transformed parameters are considered so that VarCov provides
#' the covariance matrix of function \code{estimates} with cholesky=TRUE.
#' @author Cecile Proust-Lima, Viviane Philipps
#' @seealso \code{\link{estimates}}, \code{\link{hlme}}, \code{\link{lcmm}},
#' \code{\link{multlcmm}}, \code{\link{Jointlcmm}}
#' 
#' @export
#' 
VarCov <- function(x) UseMethod("VarCov")
