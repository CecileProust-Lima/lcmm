#' Data simulation according to models from lcmm package
#'
#' This function simulates a sample according to a model estimated with \code{hlme},
#' \code{lcmm}, \code{multlcmm} or \code{Jointlcmm} functions.
#' 
#' @param object an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or
#' \code{Jointlcmm}
#' @param nsim not used (for compatibility with stats::simulate). The function simulates only one sample
#' @param seed the random seed
#' @param  times  either a data frame with 2 columns containing IDs and measurement times, or a vector of length 4 specifying the minimal and maximum measurement times, the spacing between 2 consecutive visits and the margin around this spacing
#' @param  tname the name of the variable representing the measurement times in \code{object}.
#' Default to the second column's name of times if it is a data frame, and to object$var.time otherwise.
#' @param n number of subjects to simulate. Required only if times is not a data frame.
#' @param  Xbin an optional named list giving the probabilities of the binary
#' covariates to simulate. The list's names should match the binary covariate's names
#' used in \code{object}.
#' @param  Xcont  an optional named list giving the mean and standard deviation
#' of the Gaussian covariates to simulate. The list's names should match the
#' continuous covariate's names used in \code{object}.
#' @param entry expression to simulate a subject's entry time. Default to 0.
#' @param dropout expression to simulate a subject's time to dropout. Default to NULL,
#' no dropout is considered.
#' @param pMCAR optional numeric giving an observation's probability to be missing.
#' Default to 0, no missing data are introduced.
#' @param \dots additionnal options. None is used yet.
#' @return a data frame with one line per observation and one column per variable. Variables appears in the following order : subject id, measurement time, entry time, binary covariates, continuous covariates, longitudinal outcomes, latent class, entry time, survival time, event indicator.
#'
#' @author Viviane Philipps and Cecile Proust-Lima
#'
#' @examples
#' ## estimation of a 2 classes mixed model
#' m2 <- hlme(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,subject='ID',
#'          ng=2,data=data_hlme,B=c(0.11,-0.74,-0.07,20.71,
#'                                  29.39,-1,0.13,2.45,-0.29,4.5,0.36,0.79,0.97))
#'
#' ## simulate according to model m2 with same number of subjects and
#' ## same measurement times as in data_lcmm. Binary covariates X1 and X2 are simulated
#' ## according to a Bernoulli distribution with probability p=0.5, continuous covariate
#' ## X3 is simulated according to a Gaussian distribution with mean=1 and sd=1 :
#' dsim1 <- simulate(m2, times=data_hlme[,c("ID","Time")],
#'                   Xbin=list(X1=0.5, X2=0.5), Xcont=list(X3=c(1,1)))
#'
#' ## simulate a dataset of 300 subjects according to the same model
#' ## with new observation times, equally spaced and ranging from 0 to 3 :
#' dsim2 <- simulate(m2, times=c(0,3,0.5,0), n=300, tname="Time",
#'                   Xbin=list(X1=0.5, X2=0.5), Xcont=list(X3=c(1,1)))
#'
#' 
#'
#'@export
simulate.lcmm <- function(object, nsim, seed, times, tname=NULL, n,
                          Xbin=NULL, Xcont=NULL, entry=0,
                          dropout=NULL, pMCAR=0, ...)
{
    
####   simulation des donnees    ####
    if(missing(object)) stop("The argument object should be specified")
    if(missing(times)) stop("The argument times should be specified")
    if(!missing(nsim))
    {
        if(nsim != 1) stop("The function is only available with nsim=1")
    }
    nsim <- 1 # not used
    if(!missing(seed)) set.seed(seed)
    options(warn=-1)
    on.exit(options(warn=0))

    ## si model est un modele estime, prendre les formules
    modele <- list(fixed=formula(paste("~",object$call$fixed[3])),
                   random=formula(paste("~",c(as.character(object$call$random[2]),"-1")[1])),
                   mixture=formula(paste("~",c(as.character(object$call$mixture[2]),"-1")[1])),
                   classmb=formula(paste("~",c(as.character(object$call$classmb[2]),"-1")[1])))

    nRE <- sum(object$idea0)
    ng <- object$ng
    if(any(object$linktype==1)) stop("Beta link functions are not implemented yet")
    
    if(class(object) == "multlcmm")
    {
        nprob <- object$N[1]
        ncontr <- object$N[2]
        nef <- object$N[3]-ncontr-nprob
        nvc <- object$N[4]
        nw <- object$N[5]
        nalea <- object$N[6]
        ncor <- object$N[7]
        ny <- object$N[8]
        nrisqtot <- 0
        nvarxevt <- 0
        
        ntr <- rep(0, ny)
        ntr[which(object$linktype==0)] <- 2
        ntr[which(object$linktype==1)] <- 4
        ntr[which(object$linktype==2)] <- object$nbnodes + 2
        ntr[which(object$linktype==3)] <- object$nbmod[which(object$linktype==3)] - 1
        
        bfixed <- c(0, object$best[nprob+1:nef])
        B.random <- c(1, object$best[nprob+nef+ncontr+1:nvc])
        modalites <- object$modalites
        names(modalites) <- sapply(as.character(object$linktype), function(x){switch(x, "0"="linear", "2"="spline", "3"="thres")})
        seuils <- vector("list", length=ny)
        sumntr <- 0
        for(k in 1:ny)
        {
            seuils[[k]] <- object$best[nprob+nef+ncontr+nvc+nw+ncor+nalea+ny+sumntr+1:ntr[k]]
            if(object$linktype[k]==3)
            {
                seuils[[k]] <- c(seuils[[k]][1], seuils[[k]][1]+cumsum(seuils[[k]][-1]^2))
            }
            if(object$linktype[k]==2)
            {
                modalites[[k]] <- object$linknodes[1:(ntr[k]-2),k]
            }
            sumntr <- sumntr + ntr[k]
        }
        
        sigma <- object$best[nprob+nef+ncontr+nvc+nw+ncor+1:ny]
        Ynames <- object$Ynames
        nbevt <- 0
        if(ncor>0) corr <- object$best[nprob+nef+ncontr+nvc+nw+1:ncor]
    }

    if(class(object)=="hlme")
    {
        nprob <- object$N[1]
        nef <- object$N[2]
        nvc <- object$N[3]
        nw <- object$N[4]
        ncor <- object$N[5]
        ny <- 1
        ntr <- 0
        nalea <- 0 
        nrisqtot <- 0
        nvarxevt <- 0
        ncontr <- 0
        
        bfixed <- object$best[nprob+1:nef]
        if(nRE>0) B.random <- object$best[nprob+nef+1:nvc]
        modalites <- list(none=0)
        seuils <- list(none=0)
        sigma <- object$best[nprob+nef+nvc+nw+ncor+1]
        Ynames <- as.character(object$call$fixed[2])
        nbevt <- 0
        if(ncor>0) corr <- object$best[nprob+nef+nvc+nw+1:ncor]
    }

    if(class(object)=="lcmm")
    {
        nprob <- object$N[1]
        nrisqtot <- 0
        nvarxevt <- 0
        nef <- object$N[2]
        nvc <- object$N[3]
        nw <- object$N[4]
        ncor <- object$N[6]
        ntr <- 2
        modalites <- list(linear=0)
        seuils <- list(none=object$best[nprob+nef+nvc+nw+1:2])
        if(object$linktype==2)
        {
            ntr <- length(object$linknodes)+2
            modalites <- list(spline=object$linknodes)
            seuils <- list(none=object$best[nprob+nef+nvc+nw+1:ntr])
        }
        if(object$linktype==3)
        {
            ntr <- sum(object$ide)
            modalites <- list(thres=seq(object$linknodes[1], object$linknodes[2])[as.logical(object$ide)])
            seuils <- list(none=object$best[nprob+nef+nvc+nw+1:ntr])
            seuils[[1]] <- c(seuils[[1]][1], seuils[[1]][1]+cumsum(seuils[[1]][-1]^2))
        }
        ny <- 1
        nalea <- 0 
        ncontr <- 0
        nbevt <- 0
        
        bfixed <- c(0,object$best[nprob+1:nef])
        if(nRE>0) B.random <- object$best[nprob+nef+1:nvc]
        sigma <- 1
        Ynames <- as.character(object$call$fixed[2])
        nbevt <- 0       
        if(ncor>0) corr <- object$best[nprob+nef+nvc+nw+ntr+1:ncor]
    }
    
    if(class(object)=="Jointlcmm")
    {
        nprob <- object$N[1]
        nrisqtot <- object$N[2]
        nvarxevt <- object$N[3]
        nef <- object$N[4]
        nvc <- object$N[5]
        nw <- object$N[6]
        ncor <- object$N[7]
        ntr <- object$N[8]
        ny <- 1
        nalea <- 0 
        ncontr <- 0
        nbevt <- length(object$N)-9

        if(ntr==1)# pas de link function
        {
            bfixed <- object$best[nprob+nrisqtot+nvarxevt+1:nef]
            sigma <- object$best[length(object$best)]
        }
        else
        {
            bfixed <- c(0,object$best[nprob+nrisqtot+nvarxevt+1:nef])
            sigma <- 1
        }
        if(nRE>0) B.random <- object$best[nprob+nrisqtot+nvarxevt+nef+1:nvc]
        modalites <- list(none=0)
        seuils <- list(none=0)
        if(object$linktype==2) modalites <- list(spline=object$linknodes)
        if(ntr>1) seuils <- list(object$best[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+1:ntr])
        Ynames <- as.character(object$call$fixed[2])
        if(ncor>0) corr <- object$best[nprob+nrisqtot+nvarxevt+nef+nvc+nw+1:ncor]

        
        noms.surv <- object$Names$Tnames #t0,t,d ou t,d
        if(length(noms.surv)==3 & !length(match.call()$entry)) stop("entry should be specified")
        
        typerisq <- object$hazard[[1]] # 1=piecewise 2=weibull 3=splines
        risqcom <- object$hazard[[2]] #0=specific 1=common 2=PH
        nz <- object$hazard[[4]]
        #if(any(typerisq==3)) stop("splines risk functions are not implemented yet")

        nprisq <- rep(2,nbevt)
        nprisq[which(typerisq==1)] <- nz[which(typerisq==1)]-1
        nprisq[which(typerisq==3)] <- nz[which(typerisq==3)]+2
        nrisq <- (risqcom=="Common")*nprisq + (risqcom=="Specific")*nprisq*ng + (risqcom=="PH")*(nprisq+ng-1)
        prisq <- vector("list", length=nbevt)
        ph <- matrix(NA, nbevt, ng)
        sumke <- 0
        for(ke in 1:nbevt)
        {
            r <- matrix(NA, nprisq[ke], ng)
            sumg <- 0
            for(g in 1:ng)
            {
                if(risqcom[ke]!="Specific")
                {
                    r[,g] <- object$best[nprob+sumke+sumg+1:nprisq[ke]]
                    if(risqcom[ke]=="Common")
                    {
                        ph[ke,] <- rep(0,ng)
                    }
                    else
                    {
                        ph[ke,] <- c(object$best[nprob+sumke+nprisq[ke]+1:(ng-1)],0)
                    }
                }
                else
                {
                    r[,g] <- object$best[nprob+sumke+sumg+1:nprisq[ke]]
                    sumg <- sumg + nprisq[ke]
                    ph[ke,] <- rep(0, ng)
                }

            }

            if(object$logspecif==1) r <- exp(r) else r <- r^2
            prisq[[ke]] <- r

            sumke <- sumke + nrisq[ke]
        }

        
        fsurv <- gsub("mixture","", object$call$survival[3])
        fsurv <- gsub("cause","", fsurv)
        for(ke in 1:nbevt) fsurv <- gsub(paste("cause",ke,sep=""),"", fsurv)
        fsurv <- formula(paste("~",fsurv))
        fsurv <- update(fsurv, ~.-1)
        
        idspecif <- matrix(object$idspecif, nrow=nbevt, byrow=TRUE)
        vsurv <- as.numeric(unlist(apply(idspecif,1,function(x){which(x>0)})))
        bsurv <- array(0, c(length(unique(vsurv)),ng,nbevt))
        sumnv <- 0
        jj <- 1
        for(k in 1:length(object$idcom))
        {
            if(object$idtdv[k]==1) stop("No time dependent variable can be included in the survival model")
            
            if(object$idcom[k]==1)
            {
                if(all(idspecif[,k]==1))
                {
                    bsurv[jj,,] <- object$best[nprob+nrisqtot+sumnv+1]
                    sumnv <- sumnv+1
                    jj <- jj+1
                }
                if(all(idspecif[,k]==2))
                {
                    bsurv[jj,1:ng,] <- object$best[nprob+nrisqtot+sumnv+1:ng]
                    sumnv <- sumnv+ng
                    jj <- jj+1
                }               
            }
            else
            {
                if(all(idspecif[,k]==0)) next
                for(ke in 1:nbevt)
                {
                    if(idspecif[ke,k]==1)
                    {
                        bsurv[jj,,ke] <- object$best[nprob+nrisqtot+sumnv+1]
                        sumnv <- sumnv+1  
                        jj <- jj+1                      
                    }
                    if(idspecif[ke,k]==2)
                    {
                        bsurv[jj,1:ng,ke] <- object$best[nprob+nrisqtot+sumnv+1:ng]
                        sumnv <- sumnv+ng
                        jj <- jj+1
                    }
                }
            }
        }
    }
    
        

    ## passer des parametres aux seuils
    ##seuils <- lapply(seuils, function(x){c(x[1], x[1]+cumsum(x[-1]^2))})
    names(seuils) <- names(modalites)
    
    
    if(ncontr>0)
    {
        fixed <- gsub("contrast","", modele$fixed)
        fixed <- formula(paste("~",fixed[2]))
        
        bcontr <- matrix(object$best[nprob+nef+1:ncontr], nrow=ncontr/(ny-1), ncol=ny-1)
        bcontr <- cbind(bcontr, apply(bcontr,1, function(x) {-sum(x)}))
    }
    
    if(ng>1)
    {
        bclassmb <- matrix(object$best[1:nprob],ncol=ng-1,byrow=TRUE)
        bclassmb <- cbind(bclassmb,0)
        
        beta <- matrix(NA, sum(object$idg != 0), ng)
        jj <- 0
        for(j in 1:length(object$idg))
        {
            if(object$idg[j]==1)
            {
                beta[j,] <- rep(bfixed[jj+1], ng)
                jj <- jj+1
            }
            if(object$idg[j]==2)
            {
                beta[j,] <- bfixed[jj+1:ng]
                jj <- jj+ng
            }
        }
    }
    else
    {
        beta <- matrix(bfixed, ncol=1)
    }
    
    wg <- rep(1, ng)
    if(nw>0)
    {
        wg <- c(object$best[nprob+nrisqtot+nvarxevt+nef+ncontr+nvc+1:(ng-1)],1)
    }
   

   
    if(is.data.frame(times))
    {
        ## times : data frame avec numero et temps de mesure
        times <- times[order(times[,1], times[,2]),]
        num <- unique(times[,1])
        n <- length(num)
        nmes <- as.numeric(table(times[,1]))
        tname <- c(tname,colnames(times)[2],object$var.time,"t")[1]
        idname <- c(colnames(times)[1],as.character(object$call$subject), "ID")[1]
    }
    else
    {
        ## times : tmin, tmax, ecart, marge
        if(length(times) != 4) stop("times should be either a data frame or a vector of length 4")

        t.min <- times[1]
        t.max <- times[2]
        t.ecart <- times[3]
        t.marge <- times[4]
        tname <- c(tname,object$var.time,"t")[1]
        idname <- c(as.character(object$call$subject),"ID")[1]

        if(missing(n)) stop("the number of subjects should be specified")
    }
    
    ## simuler les donnees pour n sujets
    res <- NULL
    i <- 1
    sumnmes <- 0
    while(i<=n)
    {

        if(is.data.frame(times))
        {
            t <- times[sumnmes + 1:nmes[i], 2]

            tcensure <- Inf
            if(length(match.call()$dropout))  tcensure <- eval(match.call()$dropout)

            age0 <- 0
            if(length(match.call()$entry)) age0 <- eval(match.call()$entry)

        }
        else
        {
            ## simuler les temps
            t <- seq(t.min,t.max,t.ecart)
            t <- t + runif(length(t), -abs(t.marge), abs(t.marge))
            if(t.min==0) t[1] <- t.min

            ## simuler le temps d'entree 
            age0 <- 0
            if(length(match.call()$entry))
            {
                age0 <- eval(match.call()$entry)
            }
            
            ## censure
            tcensure <- t.max
            if(length(match.call()$dropout))
            {
                tcensure <- eval(match.call()$dropout)
            }
            if(tcensure>t.min)
            {
                if(any(t>tcensure)) t <- t[-which(t>tcensure)] 
            }
            else
            {
                if(any(t<tcensure)) t <- t[-which(t<tcensure)]
            }
            
            ## ajouter temps d'entree aux temps simules
            t <- t + age0
        }
        
       
        ##initialisation du data frame
        donnees <- data.frame(1,t=t, age0=age0)
        colnames(donnees) <- c("intercept", tname, "entry")
        
        
        ## simuler la variable expl
        if(length(Xbin))
        {
            for(p in 1:length(Xbin))
            {
                X <- rbinom(1,size=1,prob=Xbin[[p]])
            
                ## ajouter X dans donnees
                old <- colnames(donnees)
                donnees <- cbind(donnees,X)
                colnames(donnees) <- c(old, names(Xbin)[p])
            }
        }

        Xc <- NULL
        if(length(Xcont))
        {
            for(p in 1:length(Xcont))
            {
                Xc <- rnorm(1,mean=Xcont[[p]][1],sd=Xcont[[p]][2])
                
                ## ajouter X dans donnees
                old <- colnames(donnees)
                donnees <- cbind(donnees,Xc)
                colnames(donnees) <- c(old, names(Xcont)[p])
            }
        }
        

        ## simuler la classe
        class <- 1
        if(ng>1)
        {
            Xclassmb <- model.matrix(modele$classmb, data=donnees[1,])
            proba <- as.vector(exp(t(Xclassmb %*% bclassmb)))
            proba <- proba/sum(proba)
            tmp <- runif(1,0,1)
            class <- as.numeric(cut(tmp, breaks=c(0,cumsum(proba))))
        }

        ## simuler les effets aleatoires
        if(nRE>0)
        {
            varRE <- matrix(0,nrow=nRE,ncol=nRE)
            if(length(B.random)==nRE)
            {
                diag(varRE) <- B.random
            }
            else
            {
                varRE[upper.tri(varRE,diag=TRUE)] <- B.random
                varRE <- t(varRE)
                varRE[upper.tri(varRE,diag=TRUE)] <- B.random
            }
            
            brandom <- as.numeric(rmvnorm(1, mean=rep(0,nRE), sigma=wg[class]^2*varRE)) 
        }
        
        ## simuler BM ou AR
        omega <- 0
        if(ncor>0)
        { 
            if(ncor==1)
            {
                Vcorr <- outer(t,t,Vectorize(function(t1,t2){corr[1]^2*min(t1,t2)}))
            }
            if(ncor==2)
            {
                Vcorr <- outer(t,t,Vectorize(function(t1,t2){corr[2]^2*exp(-corr[1]*abs(t1-t2))}))
            }

            omega <- as.numeric(rmvnorm(1, mean=rep(0,length(t)), sigma=Vcorr))
        }
        
        ## calculer la valeur du processus latent
        if(ncontr>0)
        {
            Xfixed <-  model.matrix(fixed, data=donnees)
            Xcontr <- Xfixed[,which(object$idcontr==1), drop=FALSE]
            Xbeta <- Xfixed %*% beta[,class]
        }
        else
        {
            Xbeta <- model.matrix(modele$fixed, data=donnees) %*% beta[,class]
        }
        

        if(nRE==0)
        {
            Zu <- rep(0, nrow(Xbeta))
        }
        else
        {
            Zu <- model.matrix(modele$random, data=donnees) %*% brandom
        }

        platent <- Xbeta + Zu + omega
        
        ## simuler les erreurs de mesures
        erreurs <- matrix(sapply(sigma,rnorm,n=length(t),mean=0),nrow=length(t),ncol=ny)

        ## calculer Lambda+erreurs
        ytransf <- matrix(sweep(erreurs,STATS=platent,FUN="+",MARGIN=1),nrow=length(t),ncol=ny)
        
        ## contrasts
        if(ncontr>0)
        {
            for(k in 1:ny)
            {
                ytransf[,k] <- ytransf[,k] + Xcontr %*% bcontr[,1,drop=FALSE]
            }
        }

        ## randomY
        if(nalea>0)
        {
            alpha <- rnorm(ny,0,object$best[nprob+nef+ncontr+nvc+nw+ncor+ny+1:ny])
            for(k in 1:ny)
            {
                ytransf[,k] <- ytransf[,k] + alpha[k]
            }            
        }

        ##fonction pour passer de H(Y) a Y
        transfinv <- function(ytilde,k)
        {
            nam <- names(seuils)
            nam <- tolower(nam)
            namk <- substr(nam[k],1,6)

            linktype <- pmatch(namk,c("linear","spline","thresh","none"))

            y <- NA
            
            if(linktype==3) ## thresholds
            {
                v <- c(ytilde,seuils[[k]])
                indic <- c(1,rep(0,length(seuils[[k]])))
                indic <- indic[order(v)]
                pos <- which(indic==1)
                y <- modalites[[k]][pos]
            }
            
            if(linktype==1) ##linear
            {
                y <- seuils[[k]][2]*ytilde+seuils[[k]][1]
            }

            if(linktype==2) ##splines
            {
                z <- modalites[[k]] #spline nodes
                ff <- function(x,hy,z,b){transfo_spl(x,z,b)-hy}
                y <- uniroot(f=ff,lower=z[1],upper=z[length(z)],
                             hy=ytilde,z=z,b=seuils[[k]])$root
            }

            if(linktype==4) #no transfo
            {
                y <- ytilde
            }

            
            return(y)
        }
        
        ## prendre l'inverse pour avoir des simulations de Y
        y <- NULL
        ok <- TRUE
        for(k in 1:ny)
        {
            a <- try(yk <- sapply(ytransf[,k],transfinv,k=k), silent=TRUE)
            if(class(a)=="try-error")
            {
                ok <- FALSE
                break
            }
            else
            {
                y <- c(y,yk)
            }
        }

        if(ok==FALSE) next
        
        ## garder t,X,y,age0
        if(is.data.frame(times))
        {
            d <- cbind(times[sumnmes + 1:nmes[i], ], donnees[,-c(1:2)], matrix(y,nrow=length(t),ncol=ny), class)
            sumnmes <- sumnmes + nmes[i]
        }
        else
        {
            d <- cbind(rep(i, length(t)), donnees[,-1], matrix(y,nrow=length(t),ncol=ny), class)
        }
        colnames(d) <- c(idname, tname, colnames(donnees)[-c(1:2)], Ynames, "class")


        ## survie
        if(nbevt>0)
        {
            ## variables explicatives dans a survie
            Xsurv <- model.matrix(fsurv, data=donnees[1,])
            xbsurv <- apply(bsurv,3,function(b,c){Xsurv%*%b[,c]},c=class)

            ## simuler un temps pour chaque evenement (on est en competitif)
            tevt <- rep(-Inf, nbevt)
            for(ke in 1:nbevt)
            {
                while(tevt[ke]<age0)
            {
                if(typerisq[ke]==2) # Weibull
                {
                    weib <- prisq[[ke]][,class]
                    
                    if(object$logspecif==0) #en carre, (w1*t)^w2
                    {
                        tevt[ke] <- rweibull(1, shape=weib[2], scale=1/(weib[1]*exp(xbsurv+ph[ke,class])^(1/weib[2])))
                    }
                    else #en expo, w1 * t^w2
                    {
                        tevt[ke] <- rweibull(1, shape=weib[2], scale=1/((weib[1]*exp(xbsurv+ph[ke,class]))^(1/weib[2])))
                    }
                }
                else
                {
                    if(typerisq[ke]==1) # piecewise constant
                    {
                        bl <- prisq[[ke]][,class]
                        zl <- object$hazard[[3]][1:nz[ke],ke]
                        
                        surv <- function(t,zl,bl,expXb,p)
                        {
                            j <- which.max(zl[which(zl<=t)])
                            if(j==1) som <- 0 
                            else som <- sum(bl[1:(j-1)]*(zl[2:j]-zl[1:(j-1)]))
                            
                            if(j<length(zl)) surv <- exp(-(som+bl[j]*(t-zl[j]))*expXb)
                            else surv <- exp(-som*expXb)
                            
                            return(surv-p)
                        }
                        
                        unif <- runif(1)
                        zero <- try(uniroot(surv,interval=c(zl[1],zl[length(zl)]),
                                            zl=zl,bl=bl,
                                            expXb=exp(xbsurv+ph[ke,class]),p=unif),
                                    silent=TRUE)
                        if(class(zero)!="try-error") tevt[ke] <- zero$root
                    }
                    else
                    {
                        if(typerisq[ke]==3)# splines risk
                        {
                            bl <- prisq[[ke]][,class]
                            zl <- object$hazard[[3]][1:nz[ke],ke]

                            surv <- function(t,zl,bl,expXb,p)
                            {
                                temp <- risqcum_spl(t,zl,bl)
                                surv <- exp(-temp*expXb)
                                
                                return(surv-p)
                            }
                            
                            unif <- runif(1)
                            zero <- try(uniroot(surv,interval=c(zl[1],zl[length(zl)]),
                                                zl=zl,bl=bl,
                                                expXb=exp(xbsurv+ph[ke,class]),p=unif),
                                        silent=TRUE)
                            if(class(zero)!="try-error") tevt[ke] <- zero$root
                        }
                    }
                }
            }
            }

            ## le premier evt qui arrive
            Devent <- which.min(tevt)
            Tevent <- min(tevt)

            ## censure
            if(Tevent > tcensure)
            {
                Devent <- 0
                Tevent <- tcensure
            }

            
            oldnames <- colnames(d)
            if(length(noms.surv)==3)
            {
                d <- cbind(d, rep(age0,length(t)), rep(Tevent,length(t)), rep(Devent,length(t)))
            }
            else
            {
                d <- cbind(d, rep(Tevent,length(t)), rep(Devent,length(t)))
            }
            colnames(d) <- c(oldnames, noms.surv)
                        
            ## enlever les observations de Y apres Tevent
            if(any(t > Tevent))
            {
                ## ?? on supp qu'on a tjrs la meme echelle de temps pour Y et T
            }
        }
           
        i <- i+1
        
        res <- rbind(res,d)
        
    }  
    

    
    
    ## introduire les NA en MCAR
    if(pMCAR>0)
    {
        nbmestot <- length(unlist(res[,Ynames]))
        isna <- rbinom(n=nbmestot,size=1,prob=pMCAR)
        idMCAR <- which(isna==1)
        y <- unlist(res[,Ynames])
        y[idMCAR] <- NA
        y <- matrix(y,ncol=length(Ynames))
        data <- res
        data[,Ynames] <- y
        
    }
    else
    {
        data <- res
    }
    
    
    return(data)
}


#'@export
simulate.hlme <- simulate.lcmm
#'@export
simulate.multlcmm <- simulate.lcmm
#'@export
simulate.Jointlcmm <- simulate.lcmm
