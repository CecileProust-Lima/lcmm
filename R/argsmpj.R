argsmpj <- function(longitudinal,subject,classmb,ng,survival,
                    hazard="Weibull",hazardtype="Specific",hazardnodes=NULL,TimeDepVar=NULL,
                    data,B,prior,logscale=FALSE,subset=NULL,na.action=1,posfix=NULL,
                    partialH=FALSE, ...)
{
  
  ptm <- proc.time()
  
  cl <- match.call()
  
  if(!is.list(longitudinal)) stop("longitudinal should be a list of estimated models")
  longclass <- sapply(longitudinal,class)
  if(any(!(longclass %in% c("hlme","lcmm","multlcmm")))) stop("longitudinal should only contain hlme, lcmm or multlcmm objects")
  #if(length(longclass)!=1) stop("longitudinal should only contain objects of the same class")
  if(!missing(classmb) & ng==1) stop("No classmb can be specified with ng=1")
  #if(missing(classmb) & ng==1) classmb <- ~-1
  #if(missing(classmb) & ng>1) classmb <- ~1
  if(missing(classmb)) classmb <- ~1
  if(missing(survival)) survival <- NULL
  if(ng==1) hazardtype <- "Specific"
  if(any(!(hazardtype %in% c("Common","Specific","PH")))) stop("'hazardtype' should be either 'Common' or 'Specific' or 'PH'")
  
  if(!inherits(classmb,"formula")) stop("The argument classmb must be a formula")
  if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
  if(nrow(data)==0) stop("Data should not be empty")
  if(missing(subject)){ stop("The argument subject must be specified in any model")}
  
  nom.subject <- as.character(subject)
  if(!isTRUE(nom.subject %in% colnames(data))) stop(paste("Data should contain variable",nom.subject))
  
  nom.prior <- NULL
  if(!missing(prior))
  {
    nom.prior <- as.character(prior)
    if(!isTRUE(nom.prior %in% colnames(data))) stop(paste("Data should contain variable",nom.prior))
  }
  
  nom.timedepvar <- NULL
  if(!missing(TimeDepVar))
  {
    if(!is.null(TimeDepVar))
    {  
      nom.timedepvar <- as.character(TimeDepVar)
      if(!isTRUE(nom.timedepvar %in% colnames(data))) stop(paste("Data should contain variable",nom.timedepvar))
    }
  }
  
  if(!(na.action %in% c(1,2))) stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")
  
  ##        if(length(posfix) & missing(B)) stop("A set of initial parameters must be specified if some parameters are not estimated")
  
  if(!isTRUE(all.equal(as.character(cl$subset),character(0))))
  {
    cc <- cl
    cc <- cc[c(1,which(names(cl)=="subset"))]
    cc[[1]] <- as.name("model.frame")
    cc$formula <- formula(paste("~",paste(colnames(data),collapse="+")))
    cc$data <- data
    cc$na.action <- na.pass
    data <- eval(cc)
  }
  
  attributes(data)$terms <- NULL
  
  ## contrainte selon le type de modele longitudinal
  contrainte <- as.numeric(sapply(longclass, function(x){switch(x,"hlme"=0,"lcmm"=1,"multlcmm"=2)}))
  
  ## ### donnees Y ### ##
  K <- length(longitudinal)
  dataY <- NULL
  ny <- rep(NA,K)
  Ynames <- vector("list",K)
  Xnames <- vector("list",K)
  nomsX <- unique(unlist(sapply(longitudinal,function(x) setdiff(x$Xnames2,"intercept"))))
  longicall <- vector("list",K)
  timeobs <- NULL
  
  for(k in 1:K)
  {
    ## modele k
    modk <- longitudinal[[k]]
    z <- modk$call
    z$data <- data
    z$maxiter <- 0
    z$B <- modk$best
    z$verbose <- FALSE
    if(any(modk$linktype == 2))
    {
      if(inherits(modk, "lcmm"))
      {
        ## lcmm
        nb <- length(modk$linknodes)
        z$link <- paste(nb,"-manual-splines",sep="")
        z$range <- modk$linknodes[c(1,nb)]
        z$intnodes <- modk$linknodes[c(2:(nb-1))]
      }
      else
      {
        ## multlcmm
        nspl <- 0
        link <- rep(NA,length(modk$linktype))
        intnodes <- NULL
        range <- NULL
        for (yk in 1:length(modk$linktype))
        {
          if(modk$linktype[yk] == 2)
          {
            nspl <- nspl+1
            nb <- modk$nbnodes[nspl]
            link[yk] <- paste(nb,"-manual-splines",sep="")
            intnodes <- c(intnodes, modk$linknodes[2:(nb-1),yk])
            range <- c(range, modk$linknodes[c(1,nb),yk])
          }
          else
          {
            if(modk$linktype[yk]==0) link[yk] <- "linear"
            if(modk$linktype[yk]==1)
            {
              link[yk] <- "beta"
              range <- c(range, modk$linknodes[c(1,2), yk])
            }
          }
        }
        z$link <- link
        z$intnodes <- intnodes
        z$range <- range
      }
    }
    mod <- eval(z)
    longicall[[k]] <- mod$call
    longicall[[k]]$data <- cl$data
    
    ##mod <- longitudinal[[k]]
    assign(paste("mod",k,sep=""),mod)
    
    subject <- mod$call$subject
    if(k>1){if(subject != colnames(dataY)[1]) stop("Subject variable should be the same for all longitudinal models")}
    
    ## pas de classmb
    if(mod$N[1]>(ng-1)) stop("No classmb should be specified in the longitudinal models")
    if(mod$ng!=ng) stop(paste("The longitudinal model (number ",k,") does not define the correct number of latent classes",sep=""))
    
    
    if(longclass[k]=="multlcmm")
    {
      Ynames[[k]] <- mod$Ynames
      Xnames[[k]] <- mod$Xnames
      ny[k] <- length(mod$Ynames)
    }
    else
    {
      Ynames[[k]] <- as.character(mod$call$fixed[2])
      Xnames[[k]] <- mod$Xnames
      ny[k] <- 1
    }
    
    ## donnees km
    for(m in 1:ny[k])
    {
      ## data frame de l'outcome m
      colx <- c(subject,nomsX,Ynames[[k]][m])
      if(longclass[k]=="multlcmm")
      {
        if(length(mod$na.action[[m]]))
        {
          datam <- data[-mod$na.action[[m]],colx,drop=FALSE]
        }
        else
        {
          datam <- data[,colx,drop=FALSE]
        }
      }
      else
      {
        if(length(mod$na.action))
        {
          datam <- data[-mod$na.action,colx,drop=FALSE]
        }
        else
        {
          datam <- data[,colx,drop=FALSE]
        }
      }
      
      datam <- datam[order(datam[,1]),,drop=FALSE]
      old <- colnames(datam)
      datam$processK <- k
      datam$outcomeM <- m
      datam <- datam[,c(old[1],"processK","outcomeM",old[-1])]
      
      if((k==1) & (m==1))
      {
        dataY <- datam
        colnames(dataY)[which(colnames(dataY)==Ynames[[k]][m])] <- "measureY"
      }
      else
      {
        colnames(datam)[which(colnames(datam)==Ynames[[k]][m])] <- "measureY"
        Xplus <- setdiff(colnames(datam),colnames(dataY))
        if(length(Xplus))
        {
          for(l in 1:length(Xplus))
          {
            old <- colnames(dataY)
            dataY <- cbind(dataY,NA)
            colnames(dataY) <- c(old,Xplus[l])
          }
        }
        
        Xmqt <- setdiff(colnames(dataY),colnames(datam))
        if(length(Xmqt))
        {
          for(l in 1:length(Xmqt))
          {
            old <- colnames(datam)
            datam <- cbind(datam,NA)
            colnames(datam) <- c(old,Xmqt[l])
          }
          datam <- datam[,colnames(dataY),drop=FALSE]
        }
        
        dataY <- rbind(dataY,datam)
      }
      
      ## save var.time
      if(!is.null(mod$var.time))
      {
        if(longclass[k]=="multlcmm")
        {
          timeobs <- c(timeobs,mod$pred[which(mod$pred[,2]==Ynames[[k]][m]),ncol(mod$pred)])
        }
        else
        {
          timeobs <- c(timeobs,mod$pred[,ncol(mod$pred)])
        }
      }
      else
      {
        timeobs <- c(timeobs, rep(NA, nrow(datam)))
      }
    }
    
  }
  
  dataY <- data.frame(dataY, var.time.timeobs=as.numeric(timeobs))
  
  if(is.null(survival))
  {
    nbevt <- 0
    idtrunc <- 0
    nprisq <- 0
    nrisq <- 0
    nrisqtot <- 0
    nvarxevt <- 0
    nvarxevt2 <- 0
    typrisq <- 0
    risqcom <- 0
    nom.Tentry <- NULL
    nom.Tevent <- NULL
    nom.Event <- NULL
    form.commun <- ~-1
    form.mixture <- ~-1
    form.cause <- ~-1
    survival <- ~-1
    nz <- 0
    zi <- 0
    minT <- 0
    maxT <- 0
  }
  else
  {
    ## objet Surv
    surv <- cl$survival[[2]]
    
    if(length(surv)==3) #censure droite sans troncature gauche
    {
      idtrunc <- 0 
      
      nom.Tevent <- as.character(surv[2])
      nom.Event <- as.character(surv[3])
      nom.Tentry <- NULL #si pas de troncature, Tentry=0
      
      noms.surv <-  c(nom.Tevent,nom.Event) 
    }
    
    if(length(surv)==4) #censure droite et troncature
    {
      idtrunc <- 1 
      
      nom.Tentry <- as.character(surv[2])
      nom.Tevent <- as.character(surv[3])
      nom.Event <- as.character(surv[4])
      
      noms.surv <-  c(nom.Tentry,nom.Tevent,nom.Event)
    }  
    
    ## nombre d'evenement concurrents
    Tevent <- getElement(object=data,name=nom.Tevent)
    Event <- getElement(object=data,name=nom.Event)  
    nbevt <- length(attr(do.call("Surv",list(time=Tevent,event=Event,type="mstate")),"states"))
    ##nbevt <- length(which(names(table(data[,nom.Event]))>0))    #length(unique(Event))-1   
    if(nbevt<1) nbevt <- 1 #stop("No observed event in the data")
    
    
    ## pour la formule pour survivial, creer 3 formules : 
    ## une pour les covariables en mixture, une pour les covariables avec effet specifique a la cause, et une pour les effets communs.  
    form.surv <- cl$survival[3]
    
    noms.form.surv <- all.vars(attr(terms(formula(paste("~",form.surv))),"variables"))
    if(length(noms.form.surv)==0)
    {
      form.cause <- ~-1
      form.causek <- vector("list",nbevt)
      for(k in 1:nbevt) form.causek[[k]] <- ~-1
      form.mixture <- ~-1
      form.commun <- ~-1
      asurv <- terms(~-1)
    }
    else
    {
      ##creer la formula pour cause
      form1 <- gsub("mixture","",form.surv)
      form1 <- formula(paste("~",form1))
      asurv1 <- terms(form1,specials="cause")  
      ind.cause <- attr(asurv1,"specials")$cause
      if(length(ind.cause))
      {
        form.cause <- paste(labels(asurv1)[ind.cause],collapse="+")
        form.cause <- gsub("cause","",form.cause)
        form.cause <- formula(paste("~",form.cause))
      }
      else
      {
        form.cause <- ~-1 
      }
      
      ## formules pour causek
      form.causek <- vector("list",nbevt)
      for(k in 1:nbevt)
      {
        formk <- gsub("mixture","",form.surv)
        for(kk in 1:nbevt)
        {
          if(kk != k) formk <- gsub(paste("cause",kk,sep=""),"",formk)
        }
        
        asurvk <- terms(formula(paste("~",formk)),specials=paste("cause",k,sep=""))
        ind.causek <- attr(asurvk,"specials")$cause
        
        if(length(ind.causek))
        {
          formcausek <- paste(labels(asurvk)[ind.causek],collapse="+")
          formcausek <- gsub(paste("cause",k,sep=""),"",formcausek)
          formcausek <- formula(paste("~",formcausek))
          form.causek[[k]] <- formcausek
        }
        else
        {
          form.causek[[k]] <- ~-1
        }
      }
      
      
      
      ##creer la formule pour mixture
      form2 <- form.surv
      for( k in 1:nbevt)
      {
        form2 <- gsub(paste("cause",k,sep=""),"",form2)
      }
      form2 <- gsub("cause","",form2)
      form2 <- formula(paste("~",form2))         
      asurv2 <- terms(form2,specials="mixture") 
      ind.mixture <- attr(asurv2,"specials")$mixture
      if(length(ind.mixture))
      {
        form.mixture <- paste(labels(asurv2)[ind.mixture],collapse="+")
        form.mixture <- gsub("mixture","",form.mixture)
        form.mixture <- formula(paste("~",form.mixture))
      }
      else
      {
        form.mixture <- ~-1 
      }  
      
      ## creer la formule pour ni cause ni mixture
      asurv <- terms(formula(paste("~",form.surv)),specials=c("cause","mixture",paste("cause",1:nbevt,sep="")))
      ind.commun <- setdiff(1:length(labels(asurv)),unlist(attr(asurv,"specials")))
      if(length(ind.commun))
      {
        form.commun <- paste(labels(asurv)[ind.commun],collapse="+")
        form.commun <- gsub("mixture","",form.commun) #si X1*mixture(X2), alors X1:mixture(X2) dans form.commun
        form.commun <- gsub("cause","",form.commun)   # si X1:cause(X2)
        form.commun <- formula(paste("~",form.commun))  
        ##NB: si mixture(X1)*cause(X2), X1:X2 en commun
      }
      else
      {
        form.commun <- ~-1 
      }
    }
  }
  
  ## attributs classmb
  aclassmb <- terms(classmb)
  
  ##verifier si toutes les variables sont dans data
  varSurvClas <- unique(c(all.vars(terms(survival)),all.vars(aclassmb)))
  if(!is.null(nom.timedepvar)){if(!(nom.timedepvar %in% all.vars(terms(survival)))) stop("Variable in 'TimeDepVar' should also appear as a covariate in the 'survival' argument")}  
  if(!all(varSurvClas %in% colnames(data))) stop(paste("Data should contain the variables",paste(varSurvClas,collapse=" ")))
  
  
  ##subset de data avec les variables utilisees pr S et C
  newdata <- data[,unique(c(nom.subject,varSurvClas,nom.prior)),drop=FALSE]
  
  ## remplacer les NA de prior par 0  
  if(!is.null(nom.prior))
  {
    prior <- newdata[,nom.prior]
    newdata[which(is.na(prior)),nom.prior] <- 0
  }
  
  if(nbevt>0)
  {
    ## remplacer les NA de TimeDepVar par Tevent
    Tint <- newdata[,nom.Tevent]
    Tevent <- newdata[,nom.Tevent]
    if(is.null(nom.Tentry)) Tentry <- rep(0,length(Tevent))
    else Tentry <- newdata[,nom.Tentry]
    nvdepsurv <- 0  
    if(!is.null(nom.timedepvar))
    {
      Tint <- newdata[,nom.timedepvar]
      Tint[(is.na(Tint))] <- Tevent[(is.na(Tint))]
      Tint[Tint>Tevent] <- Tevent[Tint>Tevent]
      Tint[Tint<Tentry] <- Tentry[Tint<Tentry]
      nvdepsurv <- 1
      if (length(Tint[Tint<Tevent])==0)
      {
        stop("TimeDepVar is always greater than Time of Event. \n")
        nvdepsurv <- 0
      }
      if (length(Tint[Tint>Tentry])==0)
      {
        Tint <- Tevent
        stop("TimeDepVar is always lower than Time of Entry (0 by default). \n")
        nvdepsurv  <- 0
      }
      
      newdata[,nom.timedepvar] <- Tint
    }
  }
  
  ##enlever les NA de survie et classmb
  linesNA <- apply(newdata[,varSurvClas,drop=FALSE],2,function(v) which(is.na(v)))
  linesNA <- unique(unlist(linesNA))  
  
  if(length(linesNA))
  {
    if(na.action==2) stop("Data contain missing values.")
    if(na.action==1)
    {
      newdata <- data.frame(newdata[-linesNA,nom.subject],
                            newdata[-linesNA,c(varSurvClas,nom.prior),drop=FALSE])
      colnames(newdata) <- unique(c(nom.subject,varSurvClas,nom.prior))
    }
  }
  
  
  ## prendre les sujets dans dataY et dans newdata
  selectid <- intersect(unique(dataY[,subject]),unique(newdata[,nom.subject]))
  ns <- length(selectid)
  dataY <- dataY[which(dataY[,nom.subject] %in% selectid),,drop=FALSE]
  newdata <- newdata[which(newdata[,nom.subject] %in% selectid),,drop=FALSE]
  nsdata <- unique(newdata[,c(nom.subject,varSurvClas,nom.prior),drop=FALSE])
  if(nrow(nsdata)!=ns) stop("No time-dependant variable should appear in survival nor in classmb")
  nsdata <- nsdata[order(nsdata[,1]),,drop=FALSE] # tri
  prior <- nsdata[,nom.prior]
  if(is.null(nom.prior)) prior <- rep(0,ns)
  if(nbevt>0)
  {
    if(length(nom.Tentry)) Tentry <- nsdata[,nom.Tentry]
    else Tentry <- rep(0,ns)
    Tevent <- nsdata[,nom.Tevent]
    Event <- nsdata[,nom.Event]
    if(!is.null(nom.timedepvar)) Tint <- nsdata[,nom.timedepvar]
    else Tint <- Tevent
    ind_survint <- (Tint<Tevent) + 0
  }
  else
  {
    Tevent <- 0
    Event <- 0
    Tentry <- 0
    ind_survint <- 0
  }
  
  ## Y0
  Y0 <- dataY$measureY
  
  ## var.time
  timeobs <- dataY$var.time.timeobs[order(dataY[,nom.subject])]
  
  ## X0 pour longitudinal
  nomxk <- vector("list",K)
  nv <- rep(NA,K)
  idlink <- rep(NA,sum(ny))
  nobs <- rep(NA,K)
  idiag <- rep(NA,K)
  npmtot <- rep(NA,K)
  ncor <- rep(NA,K)
  nvc <- rep(NA,K)
  ncontr <- rep(0,K)
  nalea <- rep(0,K)
  p1 <- rep(NA,K)
  p2 <- rep(NA,K)
  q <- rep(NA,K)
  ctr <- rep(NA,K)
  nalea <- rep(0,K)
  epsY <- rep(0,sum(ny))
  nw <- rep(NA,K)
  namesmod <- NULL
  idg <- NULL
  idea <- NULL
  idcontr <- NULL
  idcor <- NULL
  nodes <- NULL
  nbzitr <- rep(0,sum(ny))
  ntr <- rep(0,sum(ny))
  zitr <- vector("list", sum(ny))
  levelsFRM <- vector("list",K)
  for(k in 1:K)
  {   
    mod <- get(paste("mod",k,sep=""))
    namesmod <- c(namesmod,names(mod$best))
    
    ## levels fixed random and mixture
    levelsFRM[[k]] <- mod$levels
    
    ## formule k
    formf <- gsub("contrast","",mod$call$fixed[3])
    formk <- paste("processK+outcomeM",paste(formf,collapse="+"), sep="+")
    xfk <- model.matrix(formula(paste("~",formk)), data=dataY[which(dataY$processK==k),,drop=FALSE])
    xrk <- model.matrix(formula(paste("~",mod$call$random[2])), data=dataY[which(dataY$processK==k),,drop=FALSE])
    ## X0 dans le bon ordre
    xk <- cbind(xfk, xrk)[,union(colnames(xfk),colnames(xrk)),drop=FALSE]
    if( mod$N[5+contrainte[k]]>0)
    {
      namescor <- as.character(mod$call$cor)[2]
      xck <- model.matrix(formula(paste("~-1+",namescor)), data=dataY[which(dataY$processK==k),,drop=FALSE])
      if(!(namescor %in% colnames(xk)))
      {
        xk <- cbind(xk,xck)
      }
    }
    
    nomxk[[k]] <- colnames(xk)
    
    ## X0 merge (range par K)
    if(k>1)
    {
      #X0 <- merge(X0,xk,all=TRUE,sort=FALSE)[,union(colnames(X0),colnames(xk))]
      xh <- cbind(X0,matrix(NA,nrow=nrow(X0),ncol=ncol(xk)))
      xb <- cbind(matrix(NA,nrow=nrow(xk),ncol=ncol(X0)),xk)
      colKM <- c(colKM,ncol(X0)+which(colnames(xk) %in% c("processK","outcomeM")))
      old <- colnames(X0)
      X0 <- rbind(xh,xb)
      colnames(X0) <- c(old,colnames(xk))
    }
    else
    {
      X0 <- xk
      colKM <- which(colnames(xk) %in% c("processK","outcomeM"))
    }
    
    
    ## idx
    nv[k] <- ncol(xk)-2 # enlever process et outcome
    idg <- c(idg,mod$idg)
    idea <- c(idea,mod$idea)
    if(is.null(mod$idcontr)) mod$idcontr <- rep(0, length(mod$idg))
    idcontr <- c(idcontr,mod$idcontr)
    idcor <- c(idcor,mod$idcor)
    
    ## N
    nobs[k] <- nrow(xk)
    idiag[k] <- mod$idiag
    npmtot[k] <- length(mod$best)-mod$N[1]
    p1[k] <- sum(mod$idg==1)
    p2[k] <- sum(mod$idg==2)
    ctr[k] <- sum(mod$idcontr)
    q[k] <- sum(mod$idea)
    nvc[k] <- ifelse(idiag[k]==1,q[k],q[k]*(q[k]+1)/2)
    ncor[k] <- mod$N[5+contrainte[k]]
    nw[k] <- ifelse(contrainte[k]==2,mod$N[5],mod$N[4])
    
    ## link
    if(contrainte[k]==0) idlink[sum(ny[1:k])-ny[k]+1:ny[k]] <- -1
    else idlink[sum(ny[1:k])-ny[k]+1:ny[k]] <- mod$linktype
    
    if(contrainte[k]==2)
    {
      ncontr[k] <- mod$N[2]
      nalea[k] <- mod$N[6]
      
      nbtmp <- rep(2,ny[k])
      nbtmp[which(mod$linktype==2)] <- mod$nbnodes
      nbzitr[sum(ny[1:k])-ny[k]+1:ny[k]] <- nbtmp
      nodes <- c(nodes,as.vector(mod$linknodes))
      epsY[sum(ny[1:k])-ny[k]+1:ny[k]] <- mod$epsY
      mspl <- 0
      for (m in 1:ny[k])
      {
        if(mod$linktype[m]==0) ntr[sum(ny[1:k])-ny[k]+m] <- 2
        if(mod$linktype[m]==1) ntr[sum(ny[1:k])-ny[k]+m] <- 4
        if(mod$linktype[m]==2)
        {
          mspl <- mspl +1
          ntr[sum(ny[1:k])-ny[k]+m] <- mod$nbnodes[mspl]+2
          zitr[[sum(ny[1:k])-ny[k]+m]] <- mod$linknodes[1:mod$nbnodes[mspl],m]
        }
        else
        {
          zitr[[sum(ny[1:k])-ny[k]+m]] <- mod$linknodes[1:2,m]
        }
      }
    }
    if(contrainte[k]==1)
    {
      nbzitr[sum(ny[1:k])-ny[k]+1] <- length(mod$linknodes)
      nodes <- c(nodes,as.vector(mod$linknodes))
      zitr[[sum(ny[1:k])-ny[k]+1]] <- mod$linknodes
      epsY[sum(ny[1:k])-ny[k]+1] <- mod$epsY
      ntr[sum(ny[1:k])-ny[k]+1] <- ifelse(idlink[k]==0,2,nbzitr[k]+2)
    }
  }
  
  ## zitr en matrice
  zzitr <- matrix(0, max(sapply(zitr, length)), sum(ny))
  for(k in 1:K)
  {
    for(m in 1:ny[k])
    {
      km <- sum(ny[1:k])-ny[k]+m
      if(nbzitr[km]>0) zzitr[1:nbzitr[km],km] <- zitr[[km]]
    }
  }
  zitr <- zzitr
  
  ## refaire les idg etc pour tenir compte de toutes les var
  colnames(X0)[which(colnames(X0)=="(Intercept)")] <- "intercept"
  #X0 <- X0[,setdiff(colnames(X0),c("processK","outcomeM")),drop=FALSE]
  X0 <- X0[,-colKM,drop=FALSE]
  ## idg0 <- matrix(0,nrow=K,ncol=ncol(X0))
  ## idcontr0 <- matrix(0,nrow=K,ncol=ncol(X0))
  ## idea0 <- matrix(0,nrow=K,ncol=ncol(X0))
  ## idcor0 <- matrix(0,nrow=K,ncol=ncol(X0))
  ## for (k in 1:K)
  ## {
  ##     idg0[k,match(Xnames[[k]],colnames(X0))] <- idg[sum(nv[1:k])-nv[k]+1:nv[k]]
  ##     if(ctr[k]>0) idcontr0[k,match(Xnames[[k]],colnames(X0))] <- idcontr[sum(nv[1:k])-nv[k]+1:nv[k]]
  ##     if(q[k]>0) idea0[k,match(Xnames[[k]],colnames(X0))] <- idea[sum(nv[1:k])-nv[k]+1:nv[k]]
  ##     if(ncor[k]>0) idcor0[k,match(Xnames[[k]],colnames(X0))] <- idcor[sum(nv[1:k])-nv[k]+1:nv[k]]
  ## }
  idg0 <- idg
  idea0 <- idea
  idcontr0 <- idcontr
  if(is.null(idcontr)) idcontr0 <- rep(0,ncol(X0))
  idcor0 <- idcor
  nv0 <- ncol(X0)
  
  if(any(idlink==3))  stop("The link function thresholds is not available yet")
  
  ## X0 pour survie et classmb
  Xclassmb <- model.matrix(classmb, data=nsdata)
  Xsurv <- model.matrix(form.commun,data=nsdata)
  Xsurvmix <- model.matrix(form.mixture,data=nsdata)
  Xsurvcause <- model.matrix(form.cause,data=nsdata)
  if(nbevt>0)
  {
    for (k in 1:nbevt)
    {
      assign(paste("Xsurvcause",k,sep=""),model.matrix(form.causek[[k]],data=nsdata))
    }
  }
  
  Xns0 <- cbind(Xclassmb,Xsurv,Xsurvmix,Xsurvcause)
  if(nbevt>0)
  {
    for(k in 1:nbevt)
    {
      Xns0 <- cbind(Xns0,get(paste("Xsurvcause",k,sep="")))
    }
  }
  
  nom.unique <- unique(colnames(Xns0))
  Xns0 <- Xns0[,nom.unique,drop=FALSE]  
  Xns0 <- as.matrix(Xns0)
  
  
  if(classmb != ~-1)
  {
    z.classmb <- strsplit(colnames(Xclassmb),split=":",fixed=TRUE)
    z.classmb <- lapply(z.classmb,sort)
  }
  else
  {
    z.classmb <- list()
  }
  
  if(form.commun != ~-1)
  {
    z.surv <- strsplit(colnames(Xsurv),split=":",fixed=TRUE)
    z.surv <- lapply(z.surv,sort)
  }
  else
  {
    z.surv <- list() 
  }
  
  if(form.mixture != ~-1)
  {
    z.survmix <- strsplit(colnames(Xsurvmix),split=":",fixed=TRUE)
    z.survmix <- lapply(z.survmix,sort)
  }
  else
  {
    z.survmix <- list() 
  }  
  
  if(form.cause != ~-1)
  {
    z.survcause <- strsplit(colnames(Xsurvcause),split=":",fixed=TRUE)
    z.survcause <- lapply(z.survcause,sort)
  }
  else
  {
    z.survcause <- list() 
  }
  
  if(nbevt>0)
  {
    for(k in 1:nbevt)
    {
      if(form.causek[[k]] != ~-1)
      {
        assign(paste("z.survcause",k,sep=""),strsplit(colnames(get(paste("Xsurvcause",k,sep=""))),split=":",fixed=TRUE))
        assign(paste("z.survcause",k,sep=""),lapply(get(paste("z.survcause",k,sep="")),sort))
      }
      else
      {
        assign(paste("z.survcause",k,sep=""),list())
      }
    }
  }
  
  
  ###uniqueY0 et indiceY0
  uniqueY0 <- NULL
  indiceY0 <- NULL
  nvalSPL0 <- NULL
  nb <- 0
  for (k in 1:K)
  {
    for (m in 1:ny[k])
    {
      if(idlink[sum(ny[1:k])-ny[k]+m]!=2)
      {
        indiceY0 <- c(indiceY0, rep(0,length(dataY$measureY[which((dataY$processK==k) & (dataY$outcomeM==m))])))
        next
      }
      
      ym <- dataY$measureY[which((dataY$processK==k) & (dataY$outcomeM==m))]
      uniqueTemp <- sort(unique(ym))
      permut <- order(order(ym))  # sort(y)[order(order(y))] = y
      if(length(as.vector(table(ym)))==length(uniqueTemp))
      {
        indice <- rep(1:length(uniqueTemp), as.vector(table(ym)))
        indiceTemp <- nb + indice[permut]
        
        nb <- nb + length(uniqueTemp)
        uniqueY0 <- c(uniqueY0, uniqueTemp)
        indiceY0 <- c(indiceY0, indiceTemp)
        nvalSPL0 <- c(nvalSPL0,length(uniqueTemp))
      }
      else
      {
        uniqueY0 <- c(uniqueY0, ym)
        indiceY0 <- c(indiceY0, c(1:length(ym)))
        nb <- nb + length(ym)
        nvalSPL0 <- c(nvalSPL0, length(ym))
      }
    }
  }
  if(is.null(nvalSPL0)) nvalSPL0 <- 0
  
  
  ###ordonner les mesures par individu
  IND <- dataY[,nom.subject] 
  if(!length(indiceY0)) indiceY0 <- rep(0,length(Y0))  
  matYX <- data.frame(IND,processK=dataY$processK,outcomeM=dataY$outcomeM,
                      Y0,indiceY0,X0)
  matYXord <- matYX[order(IND),]
  Y0 <- as.numeric(matYXord[,4])
  old <- colnames(X0)
  X0 <- apply(matYXord[,-c(1:5),drop=FALSE],2,as.numeric)
  colnames(X0) <- old
  IND <- matYXord[,1]
  indiceY0 <- as.numeric(matYXord[,5])
  
  nmes <- as.vector(table(matYXord[,1]))
  #nmesM <- matrix(NA,ns,sum(ny))
  nmesM <- data.frame(id=unique(IND))
  for  (k in 1:K)
  {
    for (m in 1:ny[k])
    {
      ##nmesM[,sum(ny[1:k])-ny[k]+m] <- table(matYXord[which(matYXord$processK==k & matYXord$outcomeM==m),1]) # va pas si nb sujets differents par outcome
      temp <- rle(matYXord[which(matYXord$processK==k & matYXord$outcomeM==m),1])
      tempnmesm <- data.frame(id=temp[[2]], nm=temp[[1]])
      if(nrow(tempnmesm)!=nrow(nmesM)) warning(paste("Some subjects have no measure for outcome",m,"in process",k))
      old <- colnames(nmesM)
      nmesM <- merge(nmesM, tempnmesm, by="id", all=TRUE)
      colnames(nmesM) <- c(old, paste("k", k, "m", m, sep=""))
    }
  }
  nmesM <- as.matrix(nmesM[,-1])
  nmesM[which(is.na(nmesM))] <- 0
  
  
  if(nbevt>0)
  {
    ## test de hazard
    arghaz <- hazard
    hazard <- rep(hazard,length.out=nbevt)
    if(any(hazard %in% c("splines","Splines")))
    {
      hazard[which(hazard %in% c("splines","Splines"))] <- "5-quant-splines" 
    }
    if(any(hazard %in% c("piecewise","Piecewise")))
    {
      hazard[which(hazard %in% c("piecewise","Piecewise"))] <- "5-quant-piecewise" 
    }
    
    haz13 <- strsplit(hazard[which(!(hazard=="Weibull"))],"-")
    if(any(sapply(haz13,length)!=3)) stop("Invalid argument hazard")
    
    ## si plusieurs splines, noeuds doivent etre identiques
    nbspl <- length(which(sapply(haz13,getElement,3) %in% c("splines","Splines")))
    if(nbspl>1)
    {
      haz3 <- haz13[which(sapply(haz13,getElement,3) %in% c("splines","Splines"))]
      if(length(unique(sapply(haz3,getElement,1)))>1) stop("The nodes location of all splines hazard functions must be identical")
      
      if(length(unique(sapply(haz3,getElement,2)))>1) stop("The nodes location of all splines hazard functions must be identical")
    }
    
    nz <- rep(2,nbevt) 
    locnodes <- NULL  
    typrisq <- rep(2,nbevt)   
    nprisq <- rep(2,nbevt) 
    
    nbnodes <- 0 #longueur de hazardnodes
    ii <- 0
    dejaspl <- 0
    if(any(hazard!="Weibull"))
    {
      
      for (i in 1:nbevt)
      {
        if(hazard[i]=="Weibull") next;
        
        ii <- ii+1
        
        nz[i] <- as.numeric(haz13[[ii]][1])
        if(nz[i]<3) stop("At least 3 nodes are required")  
        typrisq[i] <- ifelse(haz13[[ii]][3] %in% c("splines","Splines"),3,1)
        nprisq[i] <- ifelse(haz13[[ii]][3] %in% c("splines","Splines"),nz[i]+2,nz[i]-1)  
        locnodes <- c(locnodes, haz13[[ii]][2])
        if(!(haz13[[ii]][3] %in% c("splines","Splines","piecewise","Piecewise"))) stop("Invalid argument hazard")
        
        if((haz13[[ii]][2]=="manual"))
        {
          if(typrisq[i]==1 | dejaspl==0)
          {
            if(length(arghaz)>1 | i==1 )
            {
              nbnodes <- nbnodes + nz[i]-2
              
            }
          }
          if(typrisq[i]==3) dejaspl <- 1
        }
        
        if(!all(locnodes %in% c("equi","quant","manual"))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")      
      }
      
      if(!is.null(hazardnodes))
      {
        if(!any(locnodes == "manual"))  stop("hazardnodes should be NULL if the nodes are not chosen manually")
        
        if(length(hazardnodes) != nbnodes) stop(paste("Vector hazardnodes should be of length",nbnodes)) 
      }  
    }
    else
    {
      if(!is.null(hazardnodes)) stop("hazardnodes should be NULL if Weibull baseline risk functions are chosen")
    }
    
    
    if(nbevt>1 & length(arghaz)==1 & nbnodes>0)
    {
      hazardnodes <- rep(hazardnodes,length.out=nbnodes*nbevt)
    }
    
    zi <- matrix(0,nrow=max(nz),ncol=nbevt)
    nb <- 0   
    
    minT1 <- 0
    maxT1 <- max(Tevent)
    
    if(idtrunc==1)
    {
      minT1 <- min(Tevent,Tentry)
      maxT1 <- max(Tevent,Tentry)
    }
    
    ## arrondir
    minT2 <- round(minT1,3)
    if(minT1<minT2) minT2 <- minT2-0.001
    minT <- minT2
    
    maxT2 <- round(maxT1,3)
    if(maxT1>maxT2) maxT2 <- maxT2+0.001
    maxT <- maxT2
    
    
    ii <- 0  
    for(i in 1:nbevt)
    {
      if(typrisq[i]==2)
      {
        zi[1:2,i] <- c(minT,maxT)
      }
      else
      {
        ii <- ii+1
        
        if(locnodes[ii]=="manual") 
        {
          zi[1:nz[i],i] <- c(minT,hazardnodes[nb+1:(nz[i]-2)],maxT)
          nb <- nb + nz[i]-2 
        } 
        if(locnodes[ii]=="equi")
        {
          zi[1:nz[i],i] <- seq(minT,maxT,length.out=nz[i]) 
        }
        if(locnodes[ii]=="quant")
        {
          #pi <- seq(0,1,length.out=nz[i])
          #pi <- pi[-length(pi)]
          #pi <- pi[-1]
          pi <- c(1:(nz[i]-2))/(nz[i]-1)
          qi <- quantile(Tevent,prob=pi)
          zi[1,i] <- minT
          zi[2:(nz[i]-1),i] <- qi
          zi[nz[i],i] <- maxT
        }
      }   
    }
    
    hazardtype <- rep(hazardtype,length.out=nbevt)  
    risqcom <- (hazardtype=="Common") + (hazardtype=="PH")*2   
    nrisq <- (risqcom==1)*nprisq + (risqcom==0)*nprisq*ng + (risqcom==2)*(nprisq+ng-1)  
    nrisqtot <- sum(nrisq)
    
    ## pour Brandom
    wBrandom <- NULL
    b0Brandom <- NULL
    nn <- 0
    for(ke in 1:nbevt)
    {
      if(hazardtype[ke]=="Common")
      {
        wBrandom <- c(wBrandom,nn+1:nprisq[ke])
        nn <- nn + nprisq[ke]
      }
      if(hazardtype[ke]=="PH")
      {
        wBrandom <- c(wBrandom,nn+1:nprisq[ke],rep(0,ng-1))
        b0Brandom <- c(b0Brandom,rep(0,ng-1))
        nn <- nn + nprisq[ke]
      }
      if(hazardtype[ke]=="Specific")
      {
        wBrandom <- c(wBrandom,rep(nn+1:nprisq[ke],ng))
      }
      nn <- nn + nprisq[ke]
    }
  }
  
  
  
  ## idprob a partir de Xns
  nv2 <- ncol(Xns0)
  z.Xns0 <- strsplit(colnames(Xns0),split=":",fixed=TRUE)
  z.Xns0 <- lapply(z.Xns0,sort)
  
  idprob <- z.Xns0 %in% z.classmb + 0
  
  
  ## indicateurs pr variables survie
  idtdv <- z.Xns0 %in% nom.timedepvar + 0
  idcause <- rep(0,nv2)
  idcom <- rep(0,nv2)
  idspecif <- matrix(0,ncol=nv2,nrow=nbevt)
  nn <- sum(nprisq)
  if(nbevt>0)
  {
    for(k in 1:nbevt)
    {
      idcause <- idcause + (z.Xns0 %in% get(paste("z.survcause",k,sep="")) )
    }
    
    for(j in 1:nv2)
    {
      if((z.Xns0[j] %in% z.surv) &  !(z.Xns0[j] %in% z.survcause) & (idcause[j]==0))
      {
        idcom[j] <- 1
        idspecif[,j] <- 1
        if(j>1)
        {
          wBrandom <- c(wBrandom,nn+1)
          nn <- nn + 1
        }
        #print("X")
      }
      
      if((z.Xns0[j] %in% z.survmix) & ( !(z.Xns0[j] %in% z.survcause) & (idcause[j]==0)))
      {
        idcom[j] <- 1
        idspecif[,j] <- 2
        if(j>1)
        {
          wBrandom <- c(wBrandom,rep(nn+1,ng))
          nn <- nn + 1
        }
        #print("mixture(X)")
      }    
      if((z.Xns0[j] %in% z.survmix) & (z.Xns0[j] %in% z.survcause))
      {
        idcom[j] <- 0
        idspecif[,j] <- 2
        if(j>1)
        {
          wBrandom <- c(wBrandom,rep(nn+1:nbevt,each=ng))
          nn <- nn + nbevt
        }
        #print("cause(mixture(X))")
      }
      
      if((z.Xns0[j] %in% z.survcause) & (!(z.Xns0[j] %in% z.survmix)))
      {
        idcom[j] <- 0
        idspecif[,j] <- 1
        if(j>1)
        {
          wBrandom <- c(wBrandom,nn+1:nbevt)
          nn <- nn + nbevt
        }
        #print("cause(X)")
      }              
      
      if(idcause[j]!=0)
      {
        if(z.Xns0[j] %in% z.survmix)
        {
          for(k in 1:nbevt)
          {
            if(z.Xns0[j] %in% get(paste("z.survcause",k,sep="")))
            {
              idcom[j] <- 0
              idspecif[k,j] <- 2
              if(j>1)
              {
                wBrandom <- c(wBrandom,rep(nn+1,ng))
                nn <- nn + 1
              }
              #cat("causek(mixture(X)) ,k=",k,"\n")
            }
          }
        }
        else
        {
          for(k in 1:nbevt)
          {
            if(z.Xns0[j] %in% get(paste("z.survcause",k,sep="")))
            {
              idcom[j] <- 0
              idspecif[k,j] <- 1
              if(j>1)
              {
                wBrandom <- c(wBrandom,nn+1)
                nn <- nn + 1
              }
              #cat("causek(X) ,k=",k,"\n")
            }                                        
          }
        }
      }
      
    }
    
    
    ## mettre des 0 pour l'intercept
    idcom[1] <- 0
    idspecif[,1] <- 0
    
    ## nb coef pr survie
    nvarxevt <- 0
    nvarxevt2 <- 0 # pr valeurs initiales calculees avec Fortran
    for(j in 1:nv2)
    {
      if(idcom[j]==1)
      {
        if(all(idspecif[,j]==1))
        {
          nvarxevt <- nvarxevt + 1
          nvarxevt2 <- nvarxevt2 + 1
        }
        if(all(idspecif[,j]==2))
        {
          nvarxevt <- nvarxevt + ng
          nvarxevt2 <- nvarxevt2 + 1
        }
      }
      
      if(idcom[j]==0)
      {
        if(all(idspecif[,j]==0)) next
        for(k in 1:nbevt)
        {
          if(idspecif[k,j]==1)
          {
            nvarxevt <- nvarxevt + 1
            nvarxevt2 <- nvarxevt2 + 1
          }
          if(idspecif[k,j]==2)
          {
            nvarxevt <- nvarxevt + ng
            nvarxevt2 <- nvarxevt2 + 1
          }
        }
      }
    }
    
  }
  
  
  
  #nea <- apply(idea0,1,sum)
  nea <- q
  predRE <- rep(0,ns*sum(nea))
  predRE_Y <- rep(0,ns*sum(nalea))
  
  nef <- rep(NA,K)
  nerr <- rep(NA,K)
  for(k in 1:K)
  {
    if(contrainte[k]==0) nef[k] <- p1[k]+ng*p2[k] else nef[k] <- p1[k]+ng*p2[k]-1
    if(contrainte[k]==2) nvc[k] <- nvc[k]-1 else nvc[k] <- nvc[k]
    if(contrainte[k]==1) nerr[k] <- 0 else nerr[k] <- ny[k]
  }
  
  
  ## nb prm
  nprob <- (ng-1)*sum(idprob)
  neftot <- sum(nef)
  ncontrtot <- sum(ncontr)
  nvctot <- sum(nvc)
  nwtot <- sum(nw)
  ncortot <- sum(ncor)
  nerrtot <- sum(nerr)
  naleatot <- sum(nalea)
  ntrtot <- sum(ntr)
  
  
  
  ##nombre total de parametres
  NPM <- nprob + nrisqtot + nvarxevt +
    neftot + ncontrtot + nvctot + nwtot + ncortot +
    naleatot + ntrtot + nerrtot
  
  V <- rep(0, NPM*(NPM+1)/2)  #pr variance des parametres
  
  ## prm fixes
  fix <- rep(0,NPM)
  if(length(posfix))
  {
    if(any(!(posfix %in% 1:NPM))) stop("Indexes in posfix are not correct")
    
    fix[posfix] <- 1
  }        
  if(length(posfix)==NPM) stop("No parameters to estimate")
  
  
  ## pour H restreint
  Hr <- as.numeric(partialH)
  pbH <- rep(0,NPM)
  if(is.logical(partialH))
  {
    if(partialH)
    {
      if(any(typrisq %in% c(1,3)))
      {
        for(k in 1:nbevt)
        {
          if(typrisq[k] %in% c(1,3))
          {
            pbH[nprob+sum(nrisq[1:k])-nrisq[k]+1:nrisq[k]] <- 1
            if(risqcom[k]==2)
            {
              pbH[nprob+sum(nrisq[1:k])] <- 0
            }
          }
        }
      }
      if(any(idlink==2))
      {
        sumnpm <- 0
        for(k in 1:K)
        {
          summ <- nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]
          for(m in 1:ny[k])
          {
            ym <- sum(ny[1:k])-ny[k]+m
            if(idlink[ym]==2)
            {
              pbH[nprob+nrisqtot+nvarxevt+sumnpm+summ+1:ntr[ym]] <- 1
            }
            summ <- summ+ntr[ym]
          }
          sumnpm <- sumnpm + npmtot[k]
        }
      }
    }
    if(sum(pbH)==0 & Hr==1) stop("No partial Hessian matrix can be defined")
  }
  else
  {
    if(!all(Hr %in% 1:NPM)) stop("Indexes in partialH are not correct")
    pbH[Hr] <- 1
  }
  indexHr <- NULL
  if(sum(pbH)>0)
  {
    if(length(posfix)) pbH1 <- pbH[-posfix] else pbH1 <- pbH
    indexHr <- which(pbH1==1)
  }
  
  ## gestion de B=random(mod)
  
  Brandom <- FALSE
  if(!missing(B))
  {
    B <- try(eval(B),silent=TRUE)
    if(inherits(B,"try-error"))
    {
      if(length(cl$B)==1) stop(B)
      if(!inherits(eval(cl$B[[2]], parent.env(environment())),"mpjlcmm")) stop("The model specified in B should be of class mpjlcmm")
      if(as.character(cl$B[1])!="random") stop("Please use random() to specify random initial values")
      
      Brandom <- TRUE
      B <- eval(cl$B[[2]], parent.env(environment()))
      if(B$conv != 1) stop("Model in argument B did not converge properly")
      #if(length(posfix)) stop("Argument posfix is not compatible with random intial values")
    }
  }
  ## valeurs initiales :
  ## faire le modele longitudinal avec ng=1 : m0
  ## faire le modele longitudinal avec ng=ng, B=random(m0), maxiter=0 : m
  ## donner m en entree -> valeurs initiales deja pretes
  ## -> faire random slt pour survie ? ou refaire tout ?
  
  logspecif <- as.numeric(logscale)
  nobs0 <- length(Y0)
  
  
  ##valeurs initiales
  if(!(missing(B)))
  {
    if(is.vector(B))
    {
      if (length(B)==NPM) b <- B
      else stop(paste("Vector B should be of length",NPM))
      
      ## remplacer varcov par cholesky
      sumnpm <- 0
      for(k in 1:K)
      {
        if(nvc[k]>0)
        {
          if(idiag[k]==1)
          {
            b[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]] <- sqrt(b[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]])
          }
          else
          {
            mvc <- matrix(0,nea[k],nea[k])
            if(contrainte[k]==2)
            {
              mvc[upper.tri(mvc,diag=TRUE)] <- c(1,b[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]])
              mvc <- t(mvc)
              mvc[upper.tri(mvc,diag=TRUE)] <- c(1,b[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]])        
            }
            else
            {
              mvc[upper.tri(mvc,diag=TRUE)] <- b[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]]
              mvc <- t(mvc)
              mvc[upper.tri(mvc,diag=TRUE)] <- b[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]]
            }
            
            ch <- chol(mvc)
            
            if(contrainte[k]==2)
            {
              b[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]] <- ch[upper.tri(ch,diag=TRUE)][-1]
            }
            else
            {
              b[nprob+nrisqtot+nvarxevt+sumnpm+nef[k]+ncontr[k]+1:nvc[k]] <- ch[upper.tri(ch,diag=TRUE)]
            }
          }
        }
        sumnpm <- sumnpm + npmtot[k]
      }
    }
    else
    {
      if(!inherits(B,"mpjlcmm")) stop("B should be either a vector or an object of class mpjlcmm")
      nef2 <- p1+p2
      ##if(contrainte!=0) nef2 <- p1+p2-1
      nef2 <- sapply(1:K, function(k){ifelse(contrainte[k]!=0,p1[k]+p2[k]-1,p1[k]+p2[k])})
      NPM2 <- sum(nprisq)+nvarxevt2+sum(nef2)+sum(ncontr)+sum(nvc)+
        sum(ncor)+sum(nerr)+sum(nalea)+sum(ntr)
      if(length(B$best)!=NPM2) stop(paste("B is not correct. The number of parameters should be",NPM2))
      
      if(Brandom==FALSE) # B deterministe
      {
        b <- rep(0,nprob) # pour nprob
        
        if(nbevt>0) #survie
        {
          bSurv <- Brandom(theta0=B$best[1:(sum(nprisq)+nvarxevt2)],v0=matrix(0,nrow=sum(nprisq)+nvarxevt2,ncol=sum(nprisq)+nvarxevt2),b0=b0Brandom,w=wBrandom)
          b <- c(b,bSurv)
        }
        
        for (k in 1:K) #m1,..,mK
        {
          mi1 <- get(as.character(B$call$longitudinal[[1+k]]),envir=sys.frame())
          z <- update(longitudinal[[k]],B=mi1,verbose=FALSE,evaluate=FALSE)
          b1 <- eval(z)$best[-c(1:(ng-1))]
          
          b <- c(b,b1)
        }
      }
      else # B random
      {
        ## variance de minit
        var0 <- matrix(0,nrow=length(B$best),ncol=length(B$best))
        var0[upper.tri(var0,diag=TRUE)] <- B$V
        var0 <- t(var0)
        var0[upper.tri(var0,diag=TRUE)] <- B$V
        
        theta0 <- B$best
        
        fix1 <- which(diag(var0)==0)
        
        ##b <- rep(0,nprob)
        w <- rep(0,nprob)
        b0 <- rep(0,nprob)
        
        if(nbevt>0) #survie
        {
          w <- c(w,wBrandom)
          b0 <- c(b0, b0Brandom)
        }
        
        sumnpm <- 0
        sumch <- 0
        sumnpmG <- nprob+nrisqtot+nvarxevt
        cholRandom <- vector("list",K)
        for (k in 1:K) #m1,..,mK
        {
          mk <- get(paste("mod",k,sep=""))
          if(inherits(mk,"multlcmm"))
          {
            multRandom <- TRUE
            if(mk$N[4]>0) cholRandom[[k]] <- sumnpmG+1:mk$N[4]
            if(idiag[k]==1)
            {
              names(cholRandom)[k] <- "diag"
            }
            else
            {
              names(cholRandom)[k] <- "full"
            }
          }
          else
          {
            multRandom <- FALSE
            if(mk$N[3]>0) cholRandom[[k]] <- sumnpmG+mk$N[2]+1:mk$N[3]
            if(idiag[k]==1)
            {
              names(cholRandom)[k] <- "diag"
            }
            else
            {
              names(cholRandom)[k] <- "full"
            }
          }
          
          ## remplacer varcov par cholesky
          avt <- 3
          if(nbevt>1) avt <- 2+nbevt
          if(B$Nprm[avt+2*K+k]>0) theta0[sum(nprisq)+nvarxevt2+sumnpm+B$Nprm[avt+k]+B$Nprm[avt+K+k]+1:B$Nprm[avt+2*K+k]] <- B$cholesky[sumch+1:B$Nprm[avt+2*K+k]]
          
          mkw <- max(w)+mk$wRandom
          mkw[which(mk$wRandom==0)] <- 0
          w <- c(w,mkw[-c(1:(ng-1))])
          
          b0 <- c(b0,mk$b0Random[-c(1:(ng-1))])
          
          sumnpm <- sumnpm + B$npmK[k]
          sumch <- sumch + B$Nprm[3+2*K+k]
          sumnpmG <- sumnpmG + npmtot[k]
        }
        
        
        
        ## gerer les posfix de B
        ww <- w
        for(j in fix1)
        {
          ww[which(w>j)] <- ww[which(w>j)]-1
          b0 <- c(b0,rep(theta0[j],length(which(w==j))))
        }
        ww[which(w %in% fix1)] <- 0
        theta1 <- theta0[setdiff(1:length(theta0),fix1)]
        var1 <- var0[setdiff(1:length(B$best),fix1),setdiff(1:length(B$best),fix1)]
        b <- Brandom(theta0=theta1,v0=var1,w=ww,b0=b0,chol=NULL, mult=NULL)
        
      } # fin random
      
    } # fin B de classe mpjlcmm
    ## cas ng>1 et B$ng=1, B deterministe ou random 
  }
  else # B missing
  {
    b <- rep(0,NPM)
    
    ## pr classmb : prm=0 par defaut
    
    if(nbevt>0)
    {
      ## prm initiaux par defaut pr risques :
      for(i in 1:nbevt)
      {
        if (typrisq[i]==2)
        {
          if(logspecif==1)
          {  
            b[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- c(rep(c(log(sum(Event==i)/sum(Tevent[Event==i])),0),ifelse(risqcom[i]==0,ng,1)),rep(1,(ng-1)*(risqcom[i]==2)))  
          }
          else
          {
            b[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- c(rep(c(sqrt(sum(Event==i)/sum(Tevent[Event==i])),1),ifelse(risqcom[i]==0,ng,1)),rep(1,(ng-1)*(risqcom[i]==2)))
          }   
        }
        else
        {
          
          if(logspecif==1)
          {  
            b[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- c(rep(log(1/nprisq[i]),ifelse(risqcom[i]==0,ng*nprisq[i],nprisq[i])),rep(1,(ng-1)*(risqcom[i]==2)))
          }
          else
          {
            b[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- c(rep(sqrt(1/nprisq[i]),ifelse(risqcom[i]==0,ng*nprisq[i],nprisq[i])),rep(1,(ng-1)*(risqcom[i]==2)))
          }   
          
          
        }
      }
      ## pr covariable survie : prm=0 par defaut
    }
    
    ## prm initiaux pr MM :
    tmp <- 0
    for (k in 1:K)
    {
      mod <- get(paste("mod",k,sep=""))
      if(ng>1)
      {
        b[nprob+nrisqtot+nvarxevt+tmp+1:(length(mod$best)-(ng-1))] <- estimates(mod)[-(1:(ng-1))] # enlever les intercept de classmb
        tmp <- tmp + length(mod$best)-(ng-1)
      }
      else
      {
        b[nprob+nrisqtot+nvarxevt+tmp+1:length(mod$best)] <- estimates(mod)
        tmp <- tmp + length(mod$best)
      }
    }
    
  }
  # B par defaut ok!
  
  ### nom au vecteur best
  
  nom.X0 <- colnames(X0)
  nom.X0[which(nom.X0=="(Intercept)")] <- "intercept"
  nom.Xns0 <- colnames(Xns0)
  nom.Xns0[which(nom.Xns0=="(Intercept)")] <- "intercept"
  
  ##prm classmb   
  if(ng>=2)
  {
    nom <- rep(nom.Xns0[idprob==1],each=ng-1)
    nom1 <- paste(nom," class",c(1:(ng-1)),sep="")
    names(b)[1:nprob] <- nom1
  }
  
  ##prm fct de risque
  if(nbevt>0)
  {
    if(isTRUE(logscale))
    {
      for(i in 1:nbevt)
      {
        nom1 <- rep(paste("event",i,sep=""),nrisq[i])
        if(typrisq[i]==2)
        {
          nom2 <- paste(nom1[1:2]," log(Weibull",1:2,")",sep="")  
          nom1[1:(2*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
          if(risqcom[i]==2) nom1[2+1:(ng-1)] <- paste(nom1[2+1:(ng-1)],"SurvPH") 
          names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1 
        }
        if(typrisq[i]==1)  
        {
          nom2 <- paste(nom1[1:(nz[i]-1)]," log(piecewise",1:(nz[i]-1),")",sep="")  
          nom1[1:((nz[i]-1)*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
          if(risqcom[i]==2) nom1[nz[i]-1+1:(ng-1)] <- paste(nom1[nz[i]-1+1:(ng-1)],"SurvPH")  
          names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1  
        }
        if(typrisq[i]==3)  
        {
          nom2 <- paste(nom1[1:(nz[i]-1)]," log(splines",1:(nz[i]+2),")",sep="")  
          nom1[1:((nz[i]+2)*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
          if(risqcom[i]==2) nom1[nz[i]+2+1:(ng-1)] <- paste(nom1[nz[i]+2+1:(ng-1)],"SurvPH")
          names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1  
        }  
      }
    }
    else
    {
      for(i in 1:nbevt)
      {
        nom1 <- rep(paste("event",i,sep=""),nrisq[i])
        if(typrisq[i]==2)
        {
          nom2 <- paste(nom1[1:2]," +/-sqrt(Weibull",1:2,")",sep="")  
          nom1[1:(2*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
          if(risqcom[i]==2) nom1[2+1:(ng-1)] <- paste(nom1[2+1:(ng-1)],"SurvPH")
          names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1  
        }
        if(typrisq[i]==1)  
        {
          nom2 <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(piecewise",1:(nz[i]-1),")",sep="")  
          nom1[1:((nz[i]-1)*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
          if(risqcom[i]==2) nom1[nz[i]-1+1:(ng-1)] <- paste(nom1[nz[i]-1+1:(ng-1)],"SurvPH")
          names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1  
        }
        if(typrisq[i]==3)  
        {
          nom2 <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(splines",1:(nz[i]+2),")",sep="")  
          nom1[1:((nz[i]+2)*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
          if(risqcom[i]==2) nom1[nz[i]+2+1:(ng-1)] <- paste(nom1[nz[i]+2+1:(ng-1)],"SurvPH") 
          names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1  
        }  
      }   
    }
    if(ng>1)
    {
      for(i in 1:nbevt)
      {
        if(risqcom[i]==1) next;
        if(risqcom[i]==0)
        {
          names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- paste(names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]],paste("class",rep(1:ng,each=nprisq[i]))) 
        }
        if(risqcom[i]==2)
        {
          names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- paste(names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]],c(rep("",nprisq[i]),paste(" class",1:(ng-1),sep="")),sep="") 
        }    
      }
    }
    
    ##prm covariables survival
    nom1 <- NULL  
    for(j in 1:nv2)
    {
      if(idcom[j]==0 & all(idspecif[,j]==0)) next
      
      if(idcom[j]==1 & all(idspecif[,j]==1)) #X
      {
        if(idtdv[j]==1)
        {
          nom1 <- c(nom1,paste("I(t>",nom.timedepvar,")",sep=""))
        }
        else
        {
          nom1 <- c(nom1,nom.Xns0[j])
        }
        next
      }
      
      if(idcom[j]==1 & all(idspecif[,j]==2)) #mixture(X)
      {
        if(idtdv[j]==1)
        {
          nom1 <- c(nom1,paste("I(t>",nom.timedepvar,") class",1:ng,sep=""))
        }
        else
        {                  
          nom1 <- c(nom1,paste(nom.Xns0[j],paste("class",1:ng,sep="")))
        }
        next
      }
      
      if(idcom[j]==0 & all(idspecif[,j]==1)) #cause(X)
      {
        if(idtdv[j]==1)
        {
          nom1 <- c(nom1,paste("I(t>",nom.timedepvar,") event",1:nbevt,sep=""))
        }
        else
        {                  
          nom1 <- c(nom1,paste(nom.Xns0[j],paste("event",1:nbevt,sep="")))
        }
        next
      }
      
      if(idcom[j]==0 & all(idspecif[,j]==2)) #cause(mixture(X))
      {
        if(idtdv[j]==1)
        {
          xevt <- paste("I(t>",nom.timedepvar,") event",1:nbevt,sep="")
          classg <- paste("class",1:ng,sep="")
          nom1 <- c(nom1,paste(rep(xevt,each=ng),rep(classg,nbevt)))
        }
        else
        {
          xevt <- paste(nom.Xns0[j],paste("event",1:nbevt,sep=""))
          classg <- paste("class",1:ng,sep="")
          nom1 <- c(nom1,paste(rep(xevt,each=ng),rep(classg,nbevt)))
        }
        next
      }
      
      
      
      if(idcom[j]==0 & idcause[j]!=0) #causek
      {
        for(k in 1:nbevt)
        {
          if(idspecif[k,j]==0) next
          
          if(idspecif[k,j]==1)
          {
            if(idtdv[j]==1)
            {
              nom1 <- c(nom1,paste("I(t>",nom.timedepvar,") event",k,sep=""))
            }
            else
            {
              nom1 <- c(nom1,paste(nom.Xns0[j],paste("event",k,sep="")))
            }
            next
          }
          
          if(idspecif[k,j]==2)
          {
            if(idtdv[j]==1)
            {
              xevtk <- paste("I(t>",nom.timedepvar,") event",k,sep="")
              classg <- paste("class",1:ng,sep="")
              nom1 <- c(nom1,paste(xevtk,classg))
            }
            else
            {
              xevtk <- paste(nom.Xns0[j],paste("event",k,sep=""))
              classg <- paste("class",1:ng,sep="")
              nom1 <- c(nom1,paste(xevtk,classg))
            }
            next
          }
        }
        
      }                
    }
    
    if(nvarxevt>0) names(b)[nprob+nrisqtot+1:nvarxevt] <- nom1 
    ##NB : pour chaque variable, coef ranges par evenement puis par classe  
  }
  
  ## noms MM
  tmp <- 0
  for (k in 1:K)
  {
    mod <- get(paste("mod",k,sep=""))
    if(ng>1)
    {
      names(b)[nprob+nrisqtot+nvarxevt+tmp+1:npmtot[k]] <- names(mod$best[-(1:(ng-1))])
    }
    else
    {
      names(b)[nprob+nrisqtot+nvarxevt+tmp+1:npmtot[k]] <- names(mod$best)
    }
    tmp <- tmp + npmtot[k]
  }
  names_best <- names(b)
  
  
  #print(b)
  #### estimation
  #idnv0 <- c(as.vector(t(idg0)),as.vector(t(idcontr0)),as.vector(t(idea0)),as.vector(t(idcor0)))
  idnv0 <- c(idg0,idcontr0,idea0,idcor0)
  idnv2 <- c(idprob,idcom,idtdv)
  idspecif <- as.vector(t(idspecif))
  #return(b)
  
  
  nfix <- sum(fix)
  bfix <- 0
  if(nfix>0)
  {
    bfix <- b[which(fix==1)]
    b <- b[which(fix==0)]
    NPM <- NPM-nfix
  }
  
  res <- list(b=b,K0=K,ny0=ny,nbevt0=nbevt,ng0=ng,ns0=ns,Y0=Y0,nobs0=nobs0,X0=X0,nv0=nv,
              Xns0=Xns0,nv20=nv2,prior0=prior,Tentr0=Tentry,Tevt0=Tevent,Devt0=Event,
              ind_survint0=ind_survint,idnv0=idnv0,idnv20=idnv2,idspecif0=idspecif,
              idlink0=idlink,epsY0=epsY,nbzitr0=nbzitr,zitr0=zitr,uniqueY0=uniqueY0,
              nvalSPL0=nvalSPL0,indiceY0=indiceY0,typrisq0=typrisq,risqcom0=risqcom,
              nz0=nz,zi0=zi,nmes0=nmesM,nea0=nea,nw0=nw,ncor0=ncor,nalea0=nalea,
              idiag0=idiag,idtrunc0=idtrunc,logspecif0=logspecif,npm0=NPM,fix0=fix,
              contrainte0=contrainte,nfix0=nfix,bfix0=bfix)
  
  return(res)
}