estimates.Jointlcmm <- function(x,cholesky=TRUE)
{
 if(missing(x)) stop("The argument x should be specified")
 if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")
 if(is.na(as.logical(cholesky))) stop("cholesky should be TRUE or FALSE")

 if(x$conv==1 | x$conv==2 | x$conv==3)
 {
  res <- x$best
  cholesky <- as.logical(cholesky)

  if(isTRUE(cholesky) & x$N[5]>0)
  {
   res[sum(x$N[1:4])+1:x$N[5]] <- x$cholesky

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

estimates <- function(x,cholesky=TRUE) UseMethod("estimates")
