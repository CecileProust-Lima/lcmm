#' @export
estimates.hlme <- function(x,cholesky=TRUE)
{
 if(missing(x)) stop("The argument x should be specified")
 if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
 if(is.na(as.logical(cholesky))) stop("cholesky should be TRUE or FALSE")

 if(x$conv==1 | x$conv==2)
 {
  res <- x$best
  cholesky <- as.logical(cholesky)

  if(isTRUE(cholesky) & x$N[3]>0)
  {
   if(!isTRUE(x$call$idiag) | x$N[3]==1)
   {
    res[x$N[1]+x$N[2]+1:x$N[3]] <- x$cholesky
   }
   else
   {
    nea <- sum(x$idea0)
    res[x$N[1]+x$N[2]+1:x$N[3]] <- x$cholesky[-setdiff(1:(nea*(nea+1)/2),1:nea*(1:nea+1)/2)]
   }
    
   names(res) <- sub("varcov","cholesky",names(res))
  }
 }
 else
 {
  res <- NA
  cat("Output can not be produced since the program stopped abnormally. \n")
 }
 
 return(res)
}



#' Maximum likelihood estimates
#' 
#' This function provides the vector of maximum likelihood estimates of a model
#' estimated with \code{hlme}, \code{lcmm}, \code{multlcmm},
#' \code{Jointlcmm}, \code{mpjlcmm}, \code{externSurv}, or \code{externX}.
#' 
#' 
#' @aliases estimates estimates.hlme estimates.lcmm estimates.Jointlcmm
#' estimates.multlcmm estimates.externSurv estimates.externX estimates.mpjlcmm
#' @param x an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or
#' \code{Jointlcmm}
#' @param cholesky optional logical indicating if the parameters of
#' variance-covariance of the random effets should be displayed instead of
#' their cholesky transformations used in the estimation process.
#' @return a vector with all estimates of the model.
#' @author Cecile Proust-Lima, Viviane Philipps
#' @seealso \code{\link{VarCov}}, \code{\link{hlme}}, \code{\link{lcmm}},
#' \code{\link{multlcmm}}, \code{\link{Jointlcmm}}
#' 
#' @export
#' 
estimates <- function(x,cholesky=TRUE) UseMethod("estimates")
