
#' @export
loglikF <- function(b0,npm0,Y0,X0,prior0,idprob0,idea0,idg0,idcor0,iderr0,ns0,ng0,nv0,nobs0,nea0,nmes0,idiag0,nwg0,ncor0,nerr0,nfix0,bfix0,fix0)
{
  res<-0
  # enlever 
  #PACKAGE="lcmmMLA" --> "lcmm"
  loglikhlme <- 0
  loglikhlme <- .Fortran(C_loglikhlmevarhetero,
                         as.double(Y0),
                         as.double(X0),
                         as.integer(prior0),
                         as.integer(idprob0),
                         as.integer(idea0),
                         as.integer(idg0),
                         as.integer(idcor0),
                         as.integer(iderr0),
                         as.integer(ns0),
                         as.integer(ng0),
                         as.integer(nv0),
                         as.integer(nobs0),
                         as.integer(nea0),
                         as.integer(nmes0),
                         as.integer(idiag0),
                         as.integer(nwg0),
                         as.integer(ncor0),
                         as.integer(nerr0),
                         as.integer(npm0),
                         as.double(b0),
                         as.integer(nfix0),
                         as.double(bfix0),
                         as.integer(fix0),
                         loglik=as.double(res))$loglik
  return(loglikhlme)
}