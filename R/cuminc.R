#' Predicted cumulative incidence of event according to a profile of covariates
#' 
#' This function computes the predicted cumulative incidence of each cause of
#' event according to a profile of covariates from a joint latent class model.
#' Confidence bands can be computed by a Monte-Carlo method.
#' 
#' 
#' @param x an object inheriting from class \code{Jointlcmm} or \code{mpjlcmm}
#' @param time a vector of times at which the cumulative incidence is
#' calculated
#' @param draws optional boolean specifying whether a Monte Carlo approximation
#' of the posterior distribution of the cumulative incidence is computed and
#' the median, 2.5\% and 97.5\% percentiles are given. Otherwise, the predicted
#' cumulative incidence is computed at the point estimate. By default,
#' draws=FALSE.
#' @param ndraws if draws=TRUE, ndraws specifies the number of draws that
#' should be generated to approximate the posterior distribution of the
#' predicted cumulative incidence. By default, ndraws=2000.
#' @param integrateOptions optional list specifying the subdivisions, rel.tol
#' and stop.on.error options (see ?integrate).
#' @param \dots further arguments, in particular values of the covariates
#' specified in the survival part of the joint model.
#' @return An object of class \code{cuminc} containing as many matrices as
#' profiles defined by the covariates values. Each of these matrices contains
#' the event-specific cumulative incidences in each latent class at the
#' different times specified.
#' @author Viviane Philipps and Cecile Proust-Lima
#' @examples
#' m2 <- Jointlcmm(fixed= Ydep1~Time*X1,mixture=~Time,random=~Time,
#' classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~X1+mixture(X2),
#' hazard="3-quant-splines",hazardtype="PH",ng=2,data=data_lcmm,
#' B=c(0.64,-0.62,0,0,0.52,0.81,0.41,0.78,0.1,0.77,-0.05,10.43,11.3,-2.6,
#' -0.52,1.41,-0.05,0.91,0.05,0.21,1.5))
#'
#' par(mfrow=c(1,2))
#' plot(cuminc(m2,time=seq(0,20),X1=0,X2=0), ylim=c(0,1))
#' plot(cuminc(m2,time=seq(0,20),X1=0,X2=1), ylim=c(0,1))
#' 
#' @seealso
#' \code{\link{Jointlcmm}}, \code{\link{plot.Jointlcmm}}, \code{\link{plot.cuminc}}
#' @export
cuminc <- function(x,time,draws=FALSE,ndraws=2000,integrateOptions=NULL,...)
    {   
        if(!inherits(x,c("Jointlcmm","mpjlcmm"))) stop("The argument 'x' must be a'Jointlcmm' or 'mpjlcmm' object")
        if(isTRUE(draws) & x$conv!=1) stop("No confidence interval can be provided since the model did not converge properly") 
       
        ## infos sur le modele
        if(inherits(x,"Jointlcmm"))
        {
            nbevt <- length(x$hazard[[1]])
            typrisq <- x$hazard[[1]]
            risqcom <- x$hazard[[2]]
            zi <- x$hazard[[3]]
            nz <- x$hazard[[4]]
        }
        else
        {
            nbevt <- x$nbevt
            if(nbevt==0) stop("No survival part in the model.")
            typrisq <- x$typrisq
            risqcom <- x$hazardtype
            zi <- x$hazardnodes
            nz <- x$nz
        }
        
        carre <- 0
        if(any(typrisq==2)) carre <- 1-x$logspecif
        ng <- x$ng

        nprisq <- sapply(1:nbevt, function(k) 2*(typrisq[k]==2)+(nz[k]-1)*(typrisq[k]==1)+(nz[k]+2)*(typrisq[k]==3))

        nrisq <- (risqcom %in% c(1,"Common"))*nprisq + (risqcom %in% c(0,"Specific"))*nprisq*ng + (risqcom %in% c(2,"PH"))*(nprisq+ng-1) 

        idspecif <- matrix(x$idspecif,nbevt,length(x$Names$Xnames),byrow=TRUE)
        idsurv <- (x$idcom!=0)+(apply(idspecif,2,function(x) any(x!=0)))
        idcom <- x$idcom
        idtdv <- x$idtdv 

        nvdepsurv <- length(x$Names$TimeDepVar.name)

        ## chercher les covariables dans ...
        Xnames <- x$Names$Xnames[which(idsurv!=0)]
        dots <- list(...)
        Xdots <- dots[Xnames] #une liste
        Xdots <- Xdots[!sapply(Xdots,is.null)]
        namesXdots <- names(Xdots)
        
        if(length(namesXdots)!=length(Xnames)) #stop(paste("Please indicate values for each of the following covariates:",paste(Xnames,collapse=" ")))
            {
                Xdefaut <- setdiff(Xnames,namesXdots)
                Xdots <- c(Xdots,as.list(rep(0,length(Xdefaut))))
                names(Xdots) <- c(namesXdots,Xdefaut)
            }
        
        ## si timedepvar est NA
        if(nvdepsurv>0)
            {
                itimedepvar <- which(names(Xdots)==x$Names$TimeDepVar.name)
                Xtimedepvar <- Xdots[[itimedepvar]]
                Xdots <- Xdots[-itimedepvar]
                if(any(is.na(Xtimedepvar))) Xtimedepvar[which(is.na(Xtimedepvar))] <- max(time)+1
                Xdots <- c(Xdots,list(Xtimedepvar))
                names(Xdots)[length(Xdots)] <- x$Names$TimeDepVar.name
            }
        
        ## matrice avec 1 profil par ligne
        Xprofil <- as.matrix(do.call("expand.grid",Xdots))
        if(nrow(Xprofil)==0) Xprofil <- matrix(0,1,1)
        
        ## fonction d'integration
        io <- list(subdivisions=100L, rel.tol=.Machine$double.eps^0.25, stop.on.error=TRUE)
        io[names(integrateOptions)] <- integrateOptions
        integrate2 <- function(...) return(integrate(..., subdivisions=io$subdivisions, rel.tol=io$rel.tol, stop.on.error=io$stop.on.error)$value)

        calculincid <- function(idraw)
            { 
                ## sous-vecteur de best avec que les prm des risques/covariables
                if(idraw>0)
                    {
                        b <- estimates(x)
                        V <- VarCov(x)
                        bdraw <- rmvnorm(1,mean=b,sigma=V)

                        brisqtot <- as.vector(bdraw[x$N[1]+1:x$N[2]])
                        if(x$N[3]>0) bvarxevt <- as.vector(bdraw[sum(x$N[1:2])+1:x$N[3]]) 
                    }
                else
                    {
                        brisqtot <- as.vector(x$best[x$N[1]+1:x$N[2]])
                        if(x$N[3]>0) bvarxevt <- as.vector(x$best[sum(x$N[1:2])+1:x$N[3]])
                    }
                

                ## declaration du resultat
                if(idraw==0) list.res <- vector(mode="list",length=nrow(Xprofil))
                if(idraw>0) array.res <- array(NA,c(length(time)*nbevt,ng+2,nrow(Xprofil)))

                
                ##  calcul selon profil puis classe :
                for (i in 1:nrow(Xprofil))
                    {
                        if(length(x$Names$TimeDepVar.name))
                            {
                                Xevt <- as.vector(Xprofil[i,])
                                tint <- Xevt[length(Xevt)]
                                Xevt <- Xevt[-length(Xevt)]
                            }
                        else
                            {
                                Xevt <- as.vector(Xprofil[i,])
                                tint <- max(time)+1
                            }
                        mat.res <- matrix(NA,nrow=length(time)*nbevt,ncol=2+ng)
                        mat.res[,1] <- rep(1:nbevt,each=length(time))
                        mat.res[,2] <- rep(time,nbevt)
                        colnames(mat.res) <- c("event","time",paste("class",1:ng,sep=""))
                        
                        for (g in 1:ng)
                            {
                                ## coef pour la classe g (1 evt par colonne)
                                bevt <- matrix(0,length(which(idsurv!=0)),nbevt)
                                bevtint <- rep(0,nbevt)
                                l <- 0
                                kcurr <- 0
                                for(k in 1:length(idcom))
                                    {
                                        if(idsurv[k]==0) next

                                        if(idcom[k]==1 & idspecif[1,k]==1)
                                            {
                                                if(idtdv[k]==1)
                                                    {
                                                        bevtint[1:nbevt] <- bvarxevt[kcurr+1]
                                                    }
                                                else
                                                    {
                                                        l <- l+1
                                                        bevt[l,1:nbevt] <- bvarxevt[kcurr+1]
                                                    }
                                                kcurr <- kcurr+1
                                            }
                                        
                                        if(idcom[k]==1 & idspecif[1,k]==2)
                                            {
                                                if(idtdv[k]==1)
                                                    {
                                                        bevtint[1:nbevt] <- bvarxevt[kcurr+g]
                                                    }
                                                else
                                                    {
                                                        l <- l+1
                                                        bevt[l,1:nbevt] <- bvarxevt[kcurr+g]
                                                    }
                                                kcurr <- kcurr+ng
                                            }

                                        if(idcom[k]==0)
                                            {
                                                if(idtdv[k]==1)
                                                    {
                                                        kecurr <- 0
                                                        for(ke in 1:nbevt)
                                                            {
                                                                if(idspecif[ke,k]==0)
                                                                    {
                                                                        bevtint[ke] <- 0
                                                                    }
                                                                if(idspecif[ke,k]==1)
                                                                    {
                                                                        bevtint[ke] <- bvarxevt[kcurr+kecurr+1]
                                                                        kecurr <- kecurr+1
                                                                    }
                                                                if(idspecif[ke,k]==2)
                                                                    {
                                                                        bevtint[ke] <- bvarxevt[kcurr+kecurr+g]
                                                                        kecurr <- kecurr+ng
                                                                    }
                                                            }
                                                        kcurr <- kcurr+kecurr
                                                    }
                                                else
                                                    {
                                                        l <- l+1
                                                        
                                                        kecurr <- 0
                                                        for(ke in 1:nbevt)
                                                            {
                                                                if(idspecif[ke,k]==0)
                                                                    {
                                                                        bevt[l,ke] <- 0
                                                                    }
                                                                if(idspecif[ke,k]==1)
                                                                    {
                                                                        bevt[l,ke] <- bvarxevt[kcurr+kecurr+1]
                                                                        kecurr <- kecurr+1
                                                                    }
                                                                if(idspecif[ke,k]==2)
                                                                    {
                                                                        bevt[l,ke] <- bvarxevt[kcurr+kecurr+g]
                                                                        kecurr <- kecurr+ng
                                                                    }
                                                            }
                                                        kcurr <- kcurr+kecurr
                                                    }
                                            }
                                    }
                                if(nrow(bevt)==0) #pas de variables dans survie
                                    {
                                        Xevt <- 0
                                        bevt <- matrix(0,1,nbevt)
                                    }                                
                                
                                
                                ##brisq de la classe g                   
                                brisqtmp <- sapply(1:nbevt,function(k) brisqtot[sum(nrisq[1:k])-nrisq[k] + (g-1)*ifelse(risqcom[k] %in% c(0,"Specific"),nprisq[k],0) + 1:nprisq[k]])
                                


                                if(is.matrix(brisqtmp))
                                    {
                                        if(x$logspecif==1) brisq <- exp(brisqtmp)
                                        if(x$logspecif==0) brisq <- brisqtmp**2
                                    }
                                if(is.list(brisqtmp))
                                    {
                                        lmax <- max(sapply(brisqtmp,length))

                                        if(x$logspecif==1)
                                            {                
                                                brisq <- lapply(brisqtmp, function(l) c(exp(l),rep(0,lmax-length(l))))
                                                brisq <- matrix(unlist(brisq),nrow=lmax,ncol=nbevt)
                                            }

                                        if(x$logspecif==0)
                                            {                
                                                brisq <- lapply(brisqtmp, function(l) c(l**2,rep(0,lmax-length(l))))
                                                brisq <- matrix(unlist(brisq),nrow=lmax,ncol=nbevt)
                                            }                              
                                        
                                    }

                                
                                ## coef PH de la classe g
                                bPH <- rep(1,nbevt)
                                if(any(risqcom %in% c(2,"PH")))
                                    {
                                        if(g<ng)
                                            {
                                                bPH <- sapply(1:nbevt,function(k) ifelse(risqcom[k] %in% c(2,"PH"),brisqtot[sum(nrisq[1:k])-(ng-1)+g],1))

                                                bPH <- exp(bPH)
                                            }                                
                                    }

                                
########## fonction de risque ##########
                                
                                risq <- function(t,evt,typrisq,brisq,zi,nz,bPH,carre)
                                    {
                                        if(typrisq[evt]==2)
                                            {
                                                if(carre==0)
                                                    {
                                                        res <- brisq[2,evt]*brisq[1,evt]*t**(brisq[2,evt]-1)
                                                    }
                                                if(carre==1)
                                                    {
                                                        res <- brisq[1,evt]*brisq[2,evt]*(brisq[1,evt]*t)**(brisq[2,evt]-1)
                                                    }
                                            }

                                        if(typrisq[evt]==1)
                                            {
                                                res <- NULL
                                                for(i in 1:length(t))
                                                    {  
                                                        zz <- zi[1:nz[evt],evt]
                                                        ordtzz <- order(c(t[i],zz))
                                                        itzz <- c(1,rep(0,length(zz)))
                                                        j <- which(itzz[ordtzz]==1)-1
                                                        if(t[i]==zz[1]) j <- 1

                                                        res <- c(res,brisq[j,evt])
                                                    }
                                            }

                                        if(typrisq[evt]==3)
                                            {
                                                res <- risq_spl(t=t,z=zi[1:nz[evt],evt],b=brisq[,evt]) 
                                            }

                                        return(res*bPH[evt])
                                    }
########## fin fonction de risque


########### fonction de risque cumule ##########
                                risqcum <- function(t,evt,typrisq,brisq,zi,nz,bPH,carre)
                                    {
                                        if(typrisq[evt]==2)
                                            {
                                                if(carre==0)
                                                    {
                                                        res <-  brisq[1,evt]*t**brisq[2,evt]
                                                    }
                                                if(carre==1)
                                                    {
                                                        res <- (brisq[1,evt]*t)**brisq[2,evt]
                                                    }
                                            }

                                        if(typrisq[evt]==1)
                                            {
                                                res <- NULL
                                                for(i in 1:length(t))
                                                    {
                                                        zz <- zi[1:nz[evt],evt]
                                                        ordtzz <- order(c(t[i],zz))
                                                        itzz <- c(1,rep(0,length(zz)))

                                                        if(t[i]==zz[1])
                                                            {
                                                                j <- 1
                                                                som <- 0
                                                            }
                                                        else
                                                            {
                                                                j <- which(itzz[ordtzz]==1)-1
                                                                if(j==1) som <- 0
                                                                else som <- sum(brisq[1:(j-1),evt]*diff(zz)[1:(j-1)])
                                                            }
                                                        
                                                        res <- c(res,som+brisq[j,evt]*(t[i]-zz[j]))
                                                    }
                                            }

                                        if(typrisq[evt]==3)
                                            {
                                                res <- risqcum_spl(t=t,z=zi[1:nz[evt],evt],b=brisq[,evt])   
                                            }

                                        return(res*bPH[evt])
                                        
                                    }
########### fin fonction de risque cumule

                                

#### fonction a integrer pour avoir l'incidence cumulee ####
                                fct_incid <- function(t,evt,tint,typrisq,brisq,zi,nz,bPH,Xevt,bevt,nbevt,bevtint,carre)
                                    {
                                        risq_tevt <- risq(t,evt,typrisq,brisq,zi,nz,bPH,carre)*exp(sum(Xevt*bevt[,evt]))*(exp(bevtint[evt])**(t>tint))

                                        somme <- 0
                                        for(ke in 1:nbevt)
                                            {
                                                somme <- somme + risqcum(t,ke,typrisq,brisq,zi,nz,bPH,carre)*exp(sum(Xevt*bevt[,ke]))*(exp(bevtint[ke])**(t>tint))
                                                
                                            }

                                        return(risq_tevt * exp(-somme) )               
                                    }


                                
                                for(ke in 1:nbevt)
                                    {
                                        mat.res[length(time)*(ke-1)+1:length(time),2+g] <- sapply(time,integrate2,f=fct_incid,lower=0,evt=ke,tint=tint,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH,Xevt=Xevt,bevt=bevt,nbevt=nbevt,bevtint=bevtint,carre=carre)
                                        
                                    }


                                
                                
                            } #fin boucle g

                        if(idraw==0)
                            {
                                list.res[[i]] <- mat.res
                                if(nvdepsurv==0)
                                    {
                                        names(list.res)[i] <- paste(colnames(Xprofil),"=",Xevt,collapse=",")
                                    }
                                else
                                    {
                                        if(tint>max(time)) tint <- NA
                                        names(list.res)[i] <- paste(colnames(Xprofil),"=",c(Xevt,tint),collapse=",")
                                    }
                            }
                        else
                            {
                                array.res[,,i] <- mat.res
                            }
                        
                    } # fin boucle i (profil)

                if(idraw==0) return(list.res)
                if(idraw>0) return(array.res)
            }

        if(!isTRUE(draws))
            {
                res <- calculincid(0)
                class(res) <- "cuminc"
                return(res)
            }
        else
            {
                res.draws <- replicate(ndraws,calculincid(1)) ## array a 4 dim
                res.draws2 <- res.draws[,-c(1,2),,,drop=FALSE]
                res.med <- apply(res.draws2,c(1,2,3),median)
                res.q025 <- apply(res.draws2,c(1,2,3),quantile,prob=0.025)
                res.q975 <- apply(res.draws2,c(1,2,3),quantile,prob=0.975)

                res <- vector("list",nrow(Xprofil))
                for(i in 1:nrow(Xprofil))
                    {
                        res[[i]] <- cbind(res.draws[,c(1,2),i,1],res.med[,,i],res.q025[,,i],res.q975[,,i])
                        colnames(res[[i]]) <- c("event","time",paste("50_class",1:ng,sep=""),paste("2.5_class",1:ng,sep=""),paste("97.5_class",1:ng,sep=""))

                        if(nvdepsurv>0 & Xprofil[i,ncol(Xprofil)]>max(time))
                            {
                                Xprofil[i,ncol(Xprofil)] <- NA
                            }
                        names(res)[i] <- paste(colnames(Xprofil),"=",Xprofil[i,],collapse=",")
                    }

                class(res) <- "cuminc"
                return(res)
            }
    }
