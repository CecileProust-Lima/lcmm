
estimates.multlcmm <- function(x,cholesky=TRUE)
{
 if(missing(x)) stop("The argument x should be specified")
 if (!inherits(x, "multlcmm")) stop("use only with \"multlcmm\" objects")
 if(is.na(as.logical(cholesky))) stop("cholesky should be TRUE or FALSE")

 if(x$conv==1 | x$conv==2 | x$conv==3)
 {
  res <- x$best
  cholesky <- as.logical(cholesky)

  if(isTRUE(cholesky) & x$N[4]>0)
  {
   cholesky <- x$cholesky[-1] # le premier n'est pas estime
   if(!isTRUE(x$call$idiag))
   {
    res[x$N[3]+1:x$N[4]] <- cholesky
   }
   else #cas diagonal : on ne met pas les 0 de cholesky
   {
    nea <- sum(x$idea0)
    res[x$N[3]+1:x$N[4]] <- cholesky[-(setdiff(1:(nea*(nea+1)/2),1:nea*(1:nea+1)/2)-1)]
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

estimates <- function(x,cholesky=TRUE) UseMethod("estimates")
