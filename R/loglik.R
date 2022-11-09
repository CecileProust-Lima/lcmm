#' Wrapper to the Fortran subroutines computing the log-likelihood
#'
#' Log-likelihood of hlme, lcmm, multlcmm, Jointlcmm and mpjlcmm models.
#' The argument's specification is not straightforward, so these functions are usually not directly used.
#'
#' @name loglik
#' @aliases loglik loglikhlme logliklcmm loglikmultlcmm loglikJointlcmm loglikmpjlcmm
#'
#' @param b the vector of estimated parameters (length npm0)
#' @param Y0 the observed values of the outcome(s) (length nobs0)
#' @param X0 the observed values of all covariates included in the model (dim nob0 * nv0)
#' @param prior0 the prior latent class (length ns0)
#' @param pprior0 the prior probabilty of each latent class (dim ns0 * ng0)
#' @param idprob0 indicator of presence in the class membership submodel (length nv0)
#' @param idea0 indicator of presence in the random part of the longitudinal submodel (length nv0)
#' @param idg0 indicator of presence in the fixed part of the longitudinal submodel (length nv0)
#' @param idcor0 indicator of presence in the correlation part of the longitudinal submodel (length nv0)
#' @param ns0 number of subjects
#' @param ng0 number of latent classes
#' @param nv0 number of covariates
#' @param nobs0 number of observations
#' @param nea0 number of random effects
#' @param nmes0 number of mesures for each subject (length ns0 or dom ns0*ny0)
#' @param idiag0 indicator of diagonal variance matrix of the random effects
#' @param nwg0 number of parameters for proportional random effects over latent classes
#' @param ncor0 number of parameters for the correlation
#' @param npm0 total number of parameters
#' @param fix0 indicator of non estimated parameter (length npm0+nfix0)
#' @param nfix0 number of non estimated parameters
#' @param bfix0 vector of non estimated parameters
#' @param epsY0 epsY values for Beta transformations
#' @param idlink0 type of transformation
#' @param nbzitr0 number of nodes for the transformations
#' @param zitr0 nodes for the transformations
#' @param minY0 minimum value for the longitudinal outcome
#' @param maxY0 maximum value for the longitudinal outcome
#' @param ide0 indicator of observed values for ordinal outcomes
#' @param idcontr0 indicator of presence as contrast in the fixed part of the longitudinal submodel (length nv0)
#' @param ny0 number of longitudinal outcomes
#' @param nalea0 number of parameters f the outcome specific random effect
#' @param uniqueY0 unique values of the longitudinal outcomes
#' @param indiceY0 correspondance between Y0 and uniqueY0
#' @param nvalSPLORD0 number of unique values for outcomes modeled with splines transformations or as ordinal outcome
#' @param methInteg0 type of integration
#' @param nMC0 number of nodes for Monte Carlo integration
#' @param dimMC0 dimension of the integration
#' @param seqMC0 sequence of integration nodes
#' @param chol0 indicator of Cholesky parameterization
#' @param tentr0 entry time for the survival submodel
#' @param tevt0 event time for the survival submodel
#' @param devt0 indicator of event for the survival submodel
#' @param ind_survint0 indicator of risk change
#' @param idcom0 indicator of presence in the survival submodel with common effect
#' @param idspecif0 indicator of presence in the survival submodel with cause-specific or class specific effect
#' @param idtdv0 indicator of 'TimeDepVar' covariate
#' @param nvalSPL0 number of unique values for outcomes modeled with splines transformations
#' @param typrisq0 type of baseline risk
#' @param risqcom0 specification of baseline risk across latent classes
#' @param nz0 number of nodes for the baseline
#' @param zi0 nodes for the baseline
#' @param nbevt0 number of events
#' @param idtrunc0 indicator of left truncation
#' @param logspecif0 indicator of logarithm parameterization
#' @param K0 number of latent processes
#' @param Xns0 the observed values of the covariates included in the survival submodel (dim ns0*nv20)
#' @param nv20 number of covariates in Xns0
#' @param Tentr0 entry time for the survival submodel (length ns0)
#' @param Tevt0 event time for the survival submodel (length ns0)
#' @param Devt0 indicator of event for the survival submodel (length ns0)
#' @param idnv0 indicator of presence in each subpart of the longitudinal models (length 4*sum(nv0))
#' @param idnv20 indicator of presence in each subpart of the survival models (length 3*nv20)
#' @param nw0 number of parameters for proportional random effects over latent classes
#' @param contrainte0 type of identifiability constraints
#'
#' @return the log-likelihood
#' @author Cecile Proust-Lima, Viviane Philipps

#'@export
loglikhlme <- function(b,Y0,X0,prior0,pprior0,idprob0,idea0,idg0,idcor0,
                       ns0,ng0,nv0,nobs0,nea0,nmes0,idiag0,nwg0,ncor0,
                       npm0,fix0,nfix0,bfix0)
{
    res <- 0
    ppi0 <- rep(0,ns0*ng0)
    resid_m <- rep(0,nobs0)
    resid_ss <- rep(0,nobs0)
    pred_m_g <- rep(0,nobs0*ng0)
    pred_ss_g <- rep(0,nobs0*ng0)
    predRE <- rep(0,ns0*nea0)
    varRE <- rep(0,ns0*nea0*(nea0+1)/2)
    estim0 <- 1
    .Fortran(C_loglikhlme,as.double(Y0),as.double(X0),as.integer(prior0),as.double(pprior0),as.integer(idprob0),as.integer(idea0),as.integer(idg0),as.integer(idcor0),as.integer(ns0),as.integer(ng0),as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg0),as.integer(ncor0),as.integer(npm0),as.double(b),as.double(ppi0),as.double(resid_m),as.double(resid_ss),as.double(pred_m_g),as.double(pred_ss_g),as.double(predRE),as.double(varRE),as.integer(fix0),as.integer(nfix0),as.double(bfix0),as.integer(estim0),loglik=as.double(res))$loglik
}


#'@rdname loglik
#'@export
logliklcmm <- function(b,Y0,X0,prior0,idprob0,idea0,idg0,idcor0,ns0,ng0,nv0,nobs0,
                       nea0,nmes0,idiag0,nwg0,ncor0,npm0,epsY0,idlink0,nbzitr0,zitr0,
                       minY0,maxY0,ide0,fix0,nfix0,bfix0)
{
    res <- 0
    ppi0 <- rep(0,ns0*ng0)
    resid_m <- rep(0,nobs0)
    resid_ss <- rep(0,nobs0)
    pred_m_g <- rep(0,nobs0*ng0)
    pred_ss_g <- rep(0,nobs0*ng0)
    predRE <- rep(0,ns0*nea0)
    Yobs <- rep(0,nobs0)
    Ydiscret <- 0
    vraisdiscret <- 0
    UACV <- 0
    rlindiv <- rep(0,ns0)
    v <- rep(0,(npm0*(npm0+1))/2)
    estim0 <- 1

    if(idlink0==3)
    {
        marker <- rep(0,2*(maxY0-minY0+1))
        transfY <- rep(0,2*(maxY0-minY0+1))
        ll <- .Fortran(C_logliklcmmord,as.double(Y0),as.double(X0),as.integer(prior0),as.integer(idprob0),as.integer(idea0),as.integer(idg0),as.integer(ns0),as.integer(ng0),as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg0),as.integer(npm0),as.double(b),as.double(ppi0),as.double(resid_m),as.double(resid_ss),as.double(pred_m_g),as.double(pred_ss_g),as.double(predRE),as.integer(minY0),as.integer(maxY0),as.integer(ide0),as.double(marker),as.double(transfY),as.double(UACV),as.double(rlindiv),as.double(v),as.integer(fix0),as.integer(nfix0),as.double(bfix0),as.integer(estim0),loglik=as.double(res))$loglik
    }
    else
    {
        nsim0 <- 0
        marker <- rep(0,nsim0)
        transfY <- rep(0,nsim0)
        
        ll <- .Fortran(C_logliklcmmcont,as.double(Y0),as.double(X0),as.integer(prior0),as.integer(idprob0),as.integer(idea0),as.integer(idg0),as.integer(idcor0),as.integer(ns0),as.integer(ng0),as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg0),as.integer(ncor0),as.integer(npm0),as.double(b),as.double(ppi0),as.double(resid_m),as.double(resid_ss),as.double(pred_m_g),as.double(pred_ss_g),as.double(predRE),as.double(epsY0),as.integer(idlink0),as.integer(nbzitr0),as.double(zitr0),as.double(marker),as.double(transfY),as.integer(nsim0),as.double(Yobs),as.integer(Ydiscret),as.double(vraisdiscret),as.double(UACV),as.double(rlindiv),as.double(v),as.integer(fix0),as.integer(nfix0),as.double(bfix0),as.integer(estim0),loglik=as.double(res))$loglik
    }

    return(ll)
}




#'@rdname loglik
#'@export
loglikmultlcmm <- function(b,Y0,X0,prior0,idprob0,idea0,idg0,idcor0,idcontr0,ny0,ns0,ng0,
                           nv0,nobs0,nea0,nmes0,idiag0,nwg0,ncor0,nalea0,npm0,
                           #ppi0,resid_m,resid_ss,pred_m_g,pred_ss_g,pred_RE,pred_RE_Y,
                           epsY0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0,nvalSPLORD0,
                           #marker,transfY,nsim0,Yobs,Ydiscret,vraisdiscret,UACV,rlindiv,
                           fix0,nfix0,bfix0,methInteg0,nMC0,dimMC0,seqMC0,chol0)#,estim0)
{
    res <- 0
    ppi0 <- rep(0,ns0*ng0)
    resid_m <- rep(0,nobs0)
    resid_ss <- rep(0,nobs0)
    pred_m_g <- rep(0,nobs0*ng0)
    pred_ss_g <- rep(0,nobs0*ng0)
    predRE <- rep(0,ns0*nea0)
    predRE_Y <- rep(0,ns0*nalea0)
    nsim0 <- 0
    marker <- rep(0,nsim0*ny0)
    transfY <- rep(0,nsim0*ny0)
    Yobs <- rep(0,nobs0)
    Ydiscret <- 0
    vraisdiscret <- 0
    UACV <- 0
    rlindiv <- rep(0,ns0)
    estim0 <- 1
    .Fortran(C_loglikmultlcmm,as.double(Y0),as.double(X0),as.integer(prior0),as.integer(idprob0),as.integer(idea0),as.integer(idg0),as.integer(idcor0),as.integer(idcontr0),as.integer(ny0),as.integer(ns0),as.integer(ng0),as.integer(nv0),as.integer(nobs0),as.integer(nea0),as.integer(nmes0),as.integer(idiag0),as.integer(nwg0),as.integer(ncor0),as.integer(nalea0),as.integer(npm0),as.double(b),as.double(ppi0),as.double(resid_m),as.double(resid_ss),as.double(pred_m_g),as.double(pred_ss_g),as.double(predRE),as.double(predRE_Y),as.double(epsY0),as.integer(idlink0),as.integer(nbzitr0),as.double(zitr0),as.double(uniqueY0),as.integer(indiceY0),as.integer(nvalSPLORD0),as.double(marker),as.double(transfY),as.integer(nsim0),as.double(Yobs),as.integer(Ydiscret),as.double(vraisdiscret),as.double(UACV),as.double(rlindiv),as.integer(fix0),as.integer(nfix0),as.double(bfix0),as.integer(methInteg0),as.integer(nMC0),as.integer(dimMC0),as.double(seqMC0),as.integer(chol0),as.integer(estim0),loglik=as.double(res))$loglik
}



#' @rdname loglik
#' @export
loglikJointlcmm <- function(b,Y0,X0,prior0,tentr0,tevt0,devt0,ind_survint0,idprob0,idea0,idg0,idcor0,idcom0,idspecif0,idtdv0,idlink0,epsY0,nbzitr0,zitr0,uniqueY0,nvalSPL0,indiceY0,typrisq0,risqcom0,nz0,zi0,ns0,ng0,nv0,nobs0,nmes0,nbevt0,nea0,nwg0,ncor0,idiag0,idtrunc0,logspecif0,npm0,fix0,nfix0,bfix0)
{
    res <- 0
    ppi0 <- rep(0,ns0*ng0)
    ppitest0 <- rep(0,ns0*ng0)
    resid_m <- rep(0,nobs0)
    resid_ss <- rep(0,nobs0)
    pred_m_g <- rep(0,nobs0*ng0)
    pred_ss_g <- rep(0,nobs0*ng0)
    predRE <- rep(0,ns0*nea0)
    nsim0 <- 0
    marker <- rep(0,nsim0)
    transfY <- rep(0,nsim0)
    Yobs <- rep(0,nobs0)
    time <- rep(0,nsim0)
    risq_est <- rep(0,nsim0*ng0*nbevt0)
    risqcum_est <- rep(0,nsim0*ng0*nbevt0)
    statglob <- 0
    statevt <- rep(0,nbevt0)
    estim0 <- 1
    .Fortran(C_loglikjointlcmm,
             as.double(Y0),
             as.double(X0),
             as.integer(prior0),
             as.double(tentr0),
             as.double(tevt0),
             as.integer(devt0),
             as.integer(ind_survint0),
             as.integer(idprob0),
             as.integer(idea0),
             as.integer(idg0),
             as.integer(idcor0),
             as.integer(idcom0),
             as.integer(idspecif0),
             as.integer(idtdv0),
             as.integer(idlink0),
             as.double(epsY0),
             as.integer(nbzitr0),
             as.double(zitr0),
             as.double(uniqueY0),
             as.integer(nvalSPL0),
             as.integer(indiceY0),
             as.integer(typrisq0),
             as.integer(risqcom0),
             as.integer(nz0),
             as.double(zi0),
             as.integer(ns0),
             as.integer(ng0),
             as.integer(nv0),
             as.integer(nobs0),
             as.integer(nmes0),
             as.integer(nbevt0),
             as.integer(nea0),
             as.integer(nwg0),
             as.integer(ncor0),
             as.integer(idiag0),
             as.integer(idtrunc0),
             as.integer(logspecif0),
             as.integer(npm0),
             best=as.double(b),
             ppi=as.double(ppi0),
             ppitest=as.double(ppitest0),
             resid_m=as.double(resid_m),
             resid_ss=as.double(resid_ss),
             pred_m_g=as.double(pred_m_g),
             pred_ss_g=as.double(pred_ss_g),
             predRE=as.double(predRE),
             time=as.double(time),
             risq_est=as.double(risq_est),
             risqcum_est=as.double(risqcum_est),
             marker=as.double(marker),
             transfY=as.double(transfY),
             as.integer(nsim0),
             Yobs=as.double(Yobs),
             statglob=as.double(statglob),
             statevt=as.double(statevt),
             as.integer(fix0),
             as.integer(nfix0),
             as.double(bfix0),
             as.integer(estim0),
             loglik=as.double(res))$loglik
}



#'@rdname loglik
#'@export
loglikmpjlcmm <- function(b,K0,ny0,nbevt0,ng0,ns0,Y0,nobs0,X0,nv0,Xns0,nv20,prior0,
                          Tentr0,Tevt0,Devt0,ind_survint0,idnv0,idnv20,idspecif0,idlink0,
                          epsY0,nbzitr0,zitr0,uniqueY0,nvalSPL0,indiceY0,typrisq0,
                          risqcom0,nz0,zi0,nmes0,nea0,nw0,ncor0,nalea0,idiag0,idtrunc0,
                          logspecif0,npm0,fix0,contrainte0,nfix0,bfix0)
{
    res <- 0
    estim0 <- 1
    ppi0 <- rep(0,ns0*ng0)
    ppitest0 <- rep(0,ns0*ng0)
    resid_m <- rep(0,nobs0)
    resid_ss <- rep(0,nobs0)
    pred_m_g <- rep(0,nobs0*ng0)
    pred_ss_g <- rep(0,nobs0*ng0)
    predRE <- rep(0,ns0*sum(nea0))
    predRE_Y <- rep(0,ns0*sum(nalea0))
    nsim0 <- 0
    time <- rep(0,nsim0)
    risq_est <- rep(0,nsim0*ng0*nbevt0)
    risqcum_est <- rep(0,nsim0*ng0*nbevt0)
    marker <- rep(0,nsim0*sum(ny0))
    transfY <- rep(0,nsim0*sum(ny0))
    Yobs <- rep(0,nobs0)
    statscoretest <- rep(0,1+nbevt0)

    .Fortran(C_loglikmpjlcmm,
             as.integer(K0),
             as.integer(ny0),
             as.integer(nbevt0),
             as.integer(ng0),
             as.integer(ns0),
             as.double(Y0),
             as.integer(nobs0),
             as.double(X0),
             as.integer(nv0),
             as.double(Xns0),
             as.integer(nv20),
             as.integer(prior0),
             as.double(Tentr0),
             as.double(Tevt0),
             as.integer(Devt0),
             as.integer(ind_survint0),
             as.integer(idnv0),
             as.integer(idnv20),
             as.integer(idspecif0),
             as.integer(idlink0),
             as.double(epsY0),
             as.integer(nbzitr0),
             as.double(zitr0),
             as.double(uniqueY0),
             as.integer(nvalSPL0),
             as.integer(indiceY0),
             as.integer(typrisq0),
             as.integer(risqcom0),
             as.integer(nz0),
             as.double(zi0),
             as.integer(nmes0),
             as.integer(nea0),
             as.integer(nw0),
             as.integer(ncor0),
             as.integer(nalea0),
             as.integer(idiag0),
             as.integer(idtrunc0),
             as.integer(logspecif0),
             as.integer(npm0),
             best=as.double(b),
             ppi=as.double(ppi0),
             ppitest=as.double(ppitest0),
             resid_m=as.double(resid_m),
             resid_ss=as.double(resid_ss),
             pred_m_g=as.double(pred_m_g),
             pred_ss_g=as.double(pred_ss_g),
             predRE=as.double(predRE),
             predRE_Y=as.double(predRE_Y),
             time=as.double(time),
             risq_est=as.double(risq_est),
             risqcum_est=as.double(risqcum_est),
             marker=as.double(marker),
             transfY=as.double(transfY),
             as.integer(nsim0),
             Yobs=as.double(Yobs),
             statscoretest=as.double(statscoretest),
             as.integer(fix0),
             as.integer(contrainte0),
             as.integer(nfix0),
             as.double(bfix0),
             as.integer(estim0),
             loglik=as.double(res),
             NAOK=TRUE)$loglik
}
