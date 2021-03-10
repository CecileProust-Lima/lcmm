## probaIRT <- function(Y,X,Tentr,Tevt,Devt,ind_survint
##      idea,idg,idcor,idcontr,idsurv,idtdv,
##      typrisq,nz,zi,nbevt,idtrunc,logspecif,
##      ny,nv,nobs,nea,nmes,idiag,ncor,nalea,
##      epsY,idlink,nbzitr,zitr,uniqueY,indiceY,
##      nvalSPLORD,fix,methInteg,nMC,dimMC,seqMC,npmtot,btot)
## {
##     proba <- 0
##     .Fortran("proba_irtsre",as.double(Y),as.double(X),as.double(Tentr),as.double(Tevt),as.integer(Devt),as.integer(ind_survint)
##      as.integer(idea),as.integer(idg),as.integer(idcor),as.integer(idcontr),as.integer(idsurv),as.integer(idtdv),
##      as.integer(typrisq),as.integer(nz),as.double(zi),as.integer(nbevt),as.integer(idtrunc),as.integer(logspecif),
##      as.integer(ny),as.integer(nv),as.integer(nobs),as.integer(nea),as.integer(nmes),as.integer(idiag),as.integer(ncor),as.integer(nalea),
##      as.double(epsY),as.integer(idlink),as.integer(nbzitr),as.double(zitr),as.integer(uniqueY),as.integer(indiceY),
##      as.integer(nvalSPLORD),as.integer(fix),as.integer(methInteg),as.integer(nMC),as.integer(dimMC),as.double(seqMC),as.integer(npmtot),as.double(btot),res=as.double(proba))$res
## }


#' @export
expectancy <- function(x, event, cond=NULL, newdata, var.time, start=0, nMC=1000, subdivisions=100L, rel.tol=.Machine$double.eps^0.25)
{
    if(missing(x)) stop("the model (argument x) is missing")
    if(!inherits(x,"irt")) stop("use only with irt model")

    if(missing(event)) stop("argument event is missing")
    if(!missing(cond) & (start==0)) stop("argument cond should only be used with start > 0")
    if(missing(var.time)) stop("argument var.time is missing")
    if(!(var.time %in% x$Names$Xvar)) stop("var.time does not appear in the model")
    if(missing(newdata) & length(setdiff(x$Names$Xvar,var.time))) stop("argument newdata is missing")
    if(!all(setdiff(x$Names$Xvar,var.time) %in% colnames(newdata))) stop(paste("newdata should include variables",paste(setdiff(x$Names$Xvar, var.time),collapse=" ")))
    
    if(any(x$idtdv==1))
    {
        if(is.na(newdata[,x$Names$TimeDepVar])) newdata[,x$Names$TimeDepVar] <- Inf
    }
    
    newdata1 <- na.omit(newdata[1,setdiff(x$Names$Xvar,c(var.time)),drop=FALSE])


    ny <- x$N[12]
    idlink <- x$linktype

    nmes <- rep(0,ny)
    nmescond <- rep(0,ny)
    Yevent <- rep(NA,ny)
    indiceYevent <- rep(NA,ny)
    Ycond <- rep(NA,ny)
    indiceYcond <- rep(NA,ny)
    nvalSPLORD <- rep(1,ny)
    uniqueY <- NULL
    for(k in 1:ny)
    {
        if(x$Names$Ynames[k] %in% names(cond))
        {
            if(idlink[k] != 3) stop("cond should only include ordinal outcomes")
            Ycond[k] <- cond[[x$Names$Ynames[k]]]
            nmes[k] <- nmes[k] + 1
            nmescond[k] <- nmescond[k] + 1
            nvalSPLORD[k] <- length(x$mod[[k]])
            indiceYcond[k] <- which(x$mod[[k]] == Ycond[k])
        }
        
        if(x$Names$Ynames[k] %in% names(event))
        {
            if(idlink[k] != 3) stop("event should only include ordinal outcomes")
            Yevent[k] <- event[[x$Names$Ynames[k]]]
            nmes[k] <- nmes[k] + 1
            nvalSPLORD[k] <- length(x$mod[[k]])
            indiceYevent[k] <- which(x$mod[[k]] == Yevent[k])
        }

        if(idlink[k]==3) uniqueY <- c(uniqueY,x$mod[[k]])
        if(idlink[k]==2) uniqueY <- c(uniqueY,x$zitr[1,k])
        
    }
    
    Y <- na.omit(as.vector(rbind(Ycond,Yevent)))
    indiceY <- na.omit(as.vector(rbind(indiceYcond,indiceYevent)))
    
    
    fixed <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$fixed[2])
    random <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$random[2])
    contr <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$contr[2])
    surv <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$form.commun[2])
    survcause <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$form.cause[2])
    cor <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",x$form$form.cor[2])

    fixed <- formula(paste("~",fixed))
    random <- formula(paste("~",random))
    contr <- formula(paste("~",contr))
    surv <- formula(paste("~",surv))
    survcause <- formula(paste("~",survcause))
    cor <- formula(paste("~",cor))

    Tentr <- 0
    Devt <- 0
    idea <- x$idea
    idg <- x$idg
    idcor <- x$idcor
    idcontr <- x$idcontr
    idsurv <- x$idsurv
    idtdv <- x$idtdv
    typrisq <- x$typrisq
    nz <- x$nz
    zi <- x$hazardnodes
    nbevt <- length(x$nevent)
    idtrunc <- 0
    logspecif <- x$logspecif
    nv <- length(x$Names$Xnames)
    nobs <- length(Y)
    nea <- sum(idea)
    idiag <- x$idiag
    ncor <- x$N[9]
    nalea <- x$N[11]
    epsY <- x$epsY
    nbzitr <- rep(2,ny)
    nbzitr[which(idlink==2)] <- x$nbnodes
    zitr <- x$linknodes
    fix <- rep(0,length(x$best))
    posfix <- eval(x$call$posfix)
    if(length(posfix)) fix[posfix] <- 1
    methInteg <- 3
    dimMC <- nea + nalea + sum(nmes)*as.numeric(ncor>0)
    seqMC <- randtoolbox::sobol(n=nMC,dim=dimMC,normal=TRUE,scrambling=1)
    npmtot <- length(x$best)
    btot <- x$best
    
    
    
    fctprob <- function(t, s, x, newdata, Y, fixed, random, contr, surv, survcause, cor,
                        Tentr,Devt,
                        idea,idg,idcor,idcontr,idsurv,idtdv,
                        typrisq,nz,zi,nbevt,idtrunc,logspecif,
                        ny,nv,nobs,nea,nmes,idiag,ncor,nalea,
                        epsY,idlink,nbzitr,zitr,uniqueY,indiceY,
                        nvalSPLORD,fix,methInteg,nMC,dimMC,seqMC,npmtot,btot)
    {
        if(s>0)
        {
            time <- unlist(lapply(as.vector(nmes),function(x) rep(c(s,t),length.out=x))) # suppose que si on a une seule mesure d'un Y, c'est au temps s
        }
        else
        {
            time <- rep(t, sum(nmes))
        }
        
        newdata <- data.frame(newdata, t=time, row.names=NULL)
        
        Xfixed <- model.matrix(fixed, data=newdata)
        Xrandom <- model.matrix(random, data=newdata)
        Xcontr <- model.matrix(contr,data=newdata)
        Xsurv <- model.matrix(surv,data=newdata)
        Xsurvcause <- model.matrix(survcause,data=newdata)
        Xcor <- model.matrix(cor,data=newdata)   
        
        X0 <- cbind(Xfixed, Xrandom, Xsurv, Xsurvcause, Xcor)

        Xnames <- x$Names$Xnames
        Xnames[1] <- "(Intercept)"
        Xnames <- gsub(paste("\\b",var.time,"\\b",sep=""),"t",Xnames)
        X0 <- X0[,Xnames,drop=FALSE]

        Tevt <- t
        ind_survint <- 0
        if(any(idtdv==1)) ind_survint <- as.numeric(newdata[1,x$Names$TimeDepVar] < t)

        proba <- 0
        
        res <- .Fortran(C_proba_irtsre,
                        as.double(Y),
                        as.double(X0),
                        as.double(Tentr),
                        as.double(Tevt),
                        as.integer(Devt),
                        as.integer(ind_survint),
                        as.integer(idea),
                        as.integer(idg),
                        as.integer(idcor),
                        as.integer(idcontr),
                        as.integer(idsurv),
                        as.integer(idtdv),
                        as.integer(typrisq),
                        as.integer(nz),
                        as.double(zi),
                        as.integer(nbevt),
                        as.integer(idtrunc),
                        as.integer(logspecif),
                        as.integer(ny),
                        as.integer(nv),
                        as.integer(nobs),
                        as.integer(nea),
                        as.integer(nmes),
                        as.integer(idiag),
                        as.integer(ncor),
                        as.integer(nalea),
                        as.double(epsY),
                        as.integer(idlink),
                        as.integer(nbzitr),
                        as.double(zitr),
                        as.double(uniqueY),
                        as.integer(indiceY),
                        as.integer(nvalSPLORD),
                        as.integer(fix),
                        as.integer(methInteg),
                        as.integer(nMC),
                        as.integer(dimMC),
                        as.double(seqMC),
                        as.integer(npmtot),
                        as.double(btot),
                        res=as.double(proba))$res

        if(!is.na(res)){ if(res == -1E-9) res <- NA }

        return(exp(res))
    }
    
    fctprobVect <- Vectorize(fctprob,"t")
    
    res1 <- integrate(f=fctprobVect, lower=start, upper=Inf, s=start, x=x, newdata=newdata1, Y=Y,
                      fixed=fixed, random=random, contr=contr, surv=surv, survcause=survcause, cor=cor,
                      Tentr=Tentr,Devt=Devt,
                      idea=idea,idg=idg,idcor=idcor,idcontr=idcontr,idsurv=idsurv,idtdv=idtdv,
                      typrisq=typrisq,nz=nz,zi=zi,nbevt=nbevt,idtrunc=idtrunc,logspecif=logspecif,
                      ny=ny,nv=nv,nobs=nobs,nea=nea,nmes=nmes,idiag=idiag,ncor=ncor,nalea=nalea,
                      epsY=epsY,idlink=idlink,nbzitr=nbzitr,zitr=zitr,uniqueY=uniqueY,indiceY=indiceY,
                      nvalSPLORD=nvalSPLORD,fix=fix,methInteg=methInteg,nMC=nMC,dimMC=dimMC,seqMC=seqMC,npmtot=npmtot,btot=btot, rel.tol=rel.tol, subdivisions=subdivisions)

    result <- res1$value
    
    if(start > 0)
    {
        nobscond <- length(na.omit(Ycond))
        Ycond <- na.omit(Ycond)
        indiceYcond <- na.omit(indiceYcond)

        res2 <- fctprob(t=start, s=0, x=x, newdata=newdata1, Y=Ycond,
            fixed=fixed, random=random, contr=contr, surv=surv, survcause=survcause, cor=cor,
            Tentr=Tentr,Devt=Devt,
            idea=idea,idg=idg,idcor=idcor,idcontr=idcontr,idsurv=idsurv,idtdv=idtdv,
            typrisq=typrisq,nz=nz,zi=zi,nbevt=nbevt,idtrunc=idtrunc,logspecif=logspecif,
            ny=ny,nv=nv,nobs=nobscond,nea=nea,nmes=nmescond,idiag=idiag,ncor=ncor,nalea=nalea,
            epsY=epsY,idlink=idlink,nbzitr=nbzitr,zitr=zitr,uniqueY=uniqueY,indiceY=indiceYcond,
            nvalSPLORD=nvalSPLORD,fix=fix,methInteg=methInteg,nMC=nMC,dimMC=dimMC,seqMC=seqMC,
            npmtot=npmtot,btot=btot)

        result <- result/res2
    }

    return(result)
}



