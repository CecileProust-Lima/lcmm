VarCov.Jointlcmm <- function(x)
{
 if(missing(x)) stop("The argument x should be specified")
 if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")

 if(x$conv==1 | x$conv==2 | x$conv==3)
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

VarCov <- function(x) UseMethod("VarCov")
