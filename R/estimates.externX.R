#' @export
estimates.externX <- function(x, cholesky)
{
  if(missing(x)) stop("The argument x should be specified")
  if (!inherits(x, "externX")) stop("use only with \"externX\" objects")
  if(!missing(cholesky)) warning("No cholesky is defined in 'externSurv' objects")
  
  if(x$conv==1 | x$conv==2 | x$conv==3)
  {
    res <- x$best
  }
  else
  {
    res <- NA
    cat("Output can not be produced since the program stopped abnormally. \n")
  }
  
  return(res)
}
  
