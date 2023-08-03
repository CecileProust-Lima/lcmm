hessienne <- function(x, method="numDeriv", nproc=1)
{
  if(!inherits(x, "mpjlcmm")) stop("Use only with mpjlcmm objects")
  
  z <- x$call
  z$posfix <- NULL
  z$B <- x$best
  z[[1]] <- as.name("argsmpj")
  argsloglik <- eval(z)
  
  if(method=="numDeriv")
  {
    names(argsloglik)[which(names(argsloglik)=="b")] <- "x"
    res <- do.call(numDeriv::hessian, c(list(func=loglikmpjlcmm), argsloglik))
  }
  else
  {
      if(nproc > 1)
      {
          clustpar <- parallel::makeCluster(nproc)
          doParallel::registerDoParallel(clustpar)
      }
      derivees <- do.call(marqLevAlg::deriva, c(list(funcpa=loglikmpjlcmm, nproc=nproc), argsloglik))
      if(nproc > 1)  parallel::stopCluster(clustpar)
      
    npm <- length(argsloglik$b)
    d2 <- -derivees$v[1:(npm*(npm+1)/2)]
    res <- matrix(0,npm,npm)
    res[upper.tri(res, diag=TRUE)] <- d2
    res <- t(res)
    res[upper.tri(res, diag=TRUE)] <- d2
  }
  
  return(res)
}
