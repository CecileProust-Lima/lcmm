#' Individual dynamic predictions from a joint latent class model
#' 
#' This function computes individual dynamic predictions and 95\% confidence
#' bands. Given a joint latent class model, a landmark time s, a horizon time t
#' and measurements until time s, the predicted probability of event in the
#' window [s,s+t] is calculated. Confidence bands can be provided using a Monte
#' Carlo method.
#' 
#' 
#' @param model an object inheriting from class \code{Jointlcmm}.
#' @param newdata a data frame containing the data from which predictions are
#' computed. This data frame must contain all the model's covariates, the
#' observations of the longitudinal and survival outcomes, the subject
#' identifier and if necessary the variables specified in prior and TimeDepVar
#' argumentsfrom Jointlcmm.
#' @param event integer giving the event for which the prediction is to be
#' calculated
#' @param landmark a numeric vector containing the landmark times.
#' @param horizon a numeric vector containing the horizon times.
#' @param var.time a character indicating the time variable in \code{newdata}
#' @param fun.time an optional function. This is only required if the time
#' scales in the longitudinal part of the model and the survival part are
#' different. In that case, \code{fun.time} is the function that translates the
#' times from the longitudinal part into the time scale of the survival part.
#' The default is the identity function which means that the two time scales
#' are the same.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @param draws optional boolean specifying whether median and confidence bands
#' of the predicted values should be computed (TRUE). IF TRUE, a Monte Carlo
#' approximation of the posterior distribution of the predicted values is
#' computed and the median, 2.5\% and 97.5\% percentiles are given. Otherwise,
#' the predicted values are computed at the point estimate. By default,
#' draws=FALSE.
#' @param ndraws if draws=TRUE, ndraws specifies the number of draws that
#' should be generated to approximate the posterior distribution of the
#' predicted values. By default, ndraws=2000.
#' @return A list containing : \item{pred}{a matrix with 4 columns if
#' draws=FALSE and 6 columns if draws=TRUE, containing the subjects identifier,
#' the landmark times, the horizon times, the predicted probability (if
#' draws=FALSE) or the median, 2.5\% and 97.5 \% percentiles of the 'ndraws'
#' probabilities calculated (if draws=TRUE). If a subject has no measurement
#' before time s or if the event has already occured at time s, his probability
#' is NA.} \item{newdata}{a data frame obtained from argument newdata
#' containing time measurements and longitudinal observations used to compute
#' the predictions}
#' @author Cecile Proust-Lima, Viviane Philipps
#' @seealso
#' \code{\link{plot.dynpred}}, \code{\link{Jointlcmm}}, \code{\link{predictY}}, \code{\link{plot.predict}}
#' @references Proust-Lima, Sene, Taylor and Jacqmin-Gadda (2014). Joint latent
#' class models of longitudinal and time-to-event data: a review. Statistical
#' Methods in Medical Research 23, 74-90.
#' @examples
#' 
#' 
#' ## Joint latent class model with 2 classes :
#' m32 <- Jointlcmm(Ydep1~Time*X1,mixture=~Time,random=~Time,subject="ID",
#' classmb=~X3,ng=2,survival=Surv(Tevent,Event)~X1+mixture(X2),
#' hazard="3-quant-splines",hazardtype="PH",data=data_lcmm,
#' B = c(0.641, -0.6217, 0, 0, 0.5045, 0.8115, -0.4316, 0.7798, 0.1027, 
#' 0.7704, -0.0479, 10.4257, 11.2972, -2.5955, -0.5234, 1.4147, 
#' -0.05, 0.9124, 0.0501, 0.2138, 1.5027))
#' 
#' ## Predictions at landmark 10 and 12 for horizon 3, 5 and 10 for two subjects :
#' 
#' dynpred(m32,landmark=c(10,12),horizon=c(3,5,10),var.time="Time",
#' fun.time=function(x){10*x},newdata=data_lcmm[1:8,])
#' \dontrun{
#' dynpred(m32,landmark=c(10,12),horizon=c(3,5,10),var.time="Time",
#' fun.time=function(x){10*x},newdata=data_lcmm[1:8,],draws=TRUE,ndraws=2000)
#' }
#' 
#' @export
#' 
#' 
dynpred <- function(model,newdata,event=1,landmark,horizon,var.time,
                    fun.time=identity,na.action=1,draws=FALSE,ndraws=2000)    
{                                                       
    if(missing(model)) stop("The argument model must be specified")
    if(class(model)!="Jointlcmm") stop("The argument model must be a 'Jointlcmm' object")
    if(missing(newdata)) stop("The argument newdata should be specified")
    if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
    if(length(event)>1) stop("Please specify only one event")
    if(!(event %in% c(1:(length(model$N)-9)))) stop("Argument 'event' is not correct")
    if(missing(landmark)) stop("Please specify at least one landmark time")
    if(missing(horizon)) stop("Please specify at least one horizon time")
    if(any(horizon<=0)) stop("horizon must be positive times")
    if(missing(var.time)) stop("The argument 'var.time' should be specified")
    if(!is.character(var.time)) stop("'var.time' should be a character")
    if(!(var.time %in% colnames(newdata))) stop("'var.time' should be a variable included in 'newdata'") 
    if(!is.function(fun.time)) stop("'fun.time' should be a function")



    
    if(model$conv==1 | model$conv==2 | model$conv==3)
        {
            if(model$conv==2 & draws==TRUE)
                {
                    cat("No confidence interval will be provided since the program did not converge properly \n")
                    draws <- FALSE
                }
            if(model$conv==3 & draws==TRUE)
                {
                    cat("No confidence interval will be provided since the program did not converge properly \n")
                    draws <- FALSE
                }
            
            
            nbland <- length(landmark)
            nbhoriz <- length(horizon) 

            nbevt <- length(model$N)-9
            idprob <- model$idprob
            idea <- model$idea
            idg <- model$idg
            idcor <- model$idcor
            idcom <- model$idcom
            idspecif <- matrix(model$idspecif,nbevt,length(idg),byrow=TRUE)
            idspecif2 <- model$idspecif
            idsurv <- (idcom!=0)+(apply(idspecif,2,function(x) any(x!=0)))
            idtdv <- model$idtdv
            idiag <- model$idiag
            nv <- length(idprob)
            nwg <- model$N[6]
            ng <- model$ng
            ncor <- model$N[7]
            nz <- model$hazard[[4]]
            zi <- model$hazard[[3]]
            typrisq <-  model$hazard[[1]]
            risqcom <- model$hazard[[2]]
            logspecif <- model$logspecif
            nvdepsurv <- length(model$Names$TimeDepVar.name)
            nvarxevt <- model$N[3]
            best <- model$best
            npm <- length(best) 
            nprisq <- sapply(1:nbevt, function(k) 2*(typrisq[k]==2)+(nz[k]-1)*(typrisq[k]==1)+(nz[k]+2)*(typrisq[k]==3))
            nrisq <- (risqcom %in% c(1,"Common"))*nprisq + (risqcom %in% c(0,"Specific"))*nprisq*ng + (risqcom %in% c(2,"PH"))*(nprisq+ng-1) 
            zitr <- model$linknodes
            nbzitr <- length(zitr)
            idlink <- model$linktype
            epsY <- model$epsY
            for(ke in 1:nbevt)
                {
                    if(risqcom[ke]=="Specific") risqcom[ke] <- 0
                    if(risqcom[ke]=="Common") risqcom[ke] <- 1
                    if(risqcom[ke]=="PH") risqcom[ke] <- 2
                }
            
            
            
            ##mettre cholesky a la place de varcov des effets aleatoires
            if(model$N[5]>0)
                {   
                    best[sum(model$N[1:4])+1:model$N[5]] <- na.omit(model$cholesky)
                }
            
            call_fixed <- model$call$fixed[3]
            if(is.null(model$call$random)) {call_random <- ~-1} else call_random <- model$call$random
            if(is.null(model$call$classmb)) {call_classmb <- ~-1} else call_classmb <- model$call$classmb
            if(is.null(model$call$survival)) {call_survival <- ~-1} else call_survival <- model$call$survival[3]

            call_survival <- gsub("mixture","",call_survival)
            for(ke in 1:nbevt)
                {
                    call_survival <- gsub(paste("cause",ke,sep=""),"",call_survival)
                }
            call_survival <- gsub("cause","",call_survival)
            call_survival <- call(call_survival) 


            if(!(na.action%in%c(1,2))) stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")

            if(na.action==1)
                {
                    na.action=na.omit
                }
            else
                {
                    na.action=na.fail
                }

            modelNames <- c(model$Names$Xnames2, # variables expl (MM et survie)
                            model$Names$Yname,         # nom de l'outcome du MM
                            model$Names$ID,            # identifiant sujets
                            model$Names$prior.name)    # nom prior
                            #model$Names$TimeDepVar.name)     # nom TimeDepVar
            
            if(length(model$Names$Tnames)>2)
                {
                    modelNames <- c(modelNames,model$Names$Tnames[1]) #Tentry
                }                


            if(model$Names$Xnames[1]=="intercept")
                {
                    newdata1 <- cbind(intercept=rep(1,length=length(newdata[,1])),newdata[,setdiff(colnames(newdata),model$Names$Xnames[1]),drop=FALSE])
                    colnames(newdata1) <- c("intercept",setdiff(colnames(newdata),model$Names$Xnames[1]))
                    newdata1 <- data.frame(newdata1)
                } 
            
            ##verif des variables dans newdata1 #mettre des valeurs par defaut dans prior et TimedepVar ?
            if(!all(modelNames %in% colnames(newdata1))) stop(paste(c("newdata should at least include the following covariates: ","\n",modelNames),collapse=" "))

            
            ## ordonner selon numero et temps
            newdata1 <- newdata1[order(newdata1[,model$Names$ID],newdata1[,var.time]),,drop=FALSE]  #est ce que fun.time peut modifier l'ordre????

### pour les facteurs

            Xnames2 <- model$Names$Xnames2

            ##donnees de l estimation
            if(!is.null(model$data))
            {
                olddata <- model$data
            }
            else
            {
                olddata <- eval(model$call$data)
            }
            
            ##cas ou une variable du dataset est un facteur
            for(v in Xnames2[-1])
                {
                    if (is.factor(olddata[,v]) & !(is.factor(newdata[,v])))
                        {
                            mod <- levels(olddata[,v])
                            if (!(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
                            newdata1[,v] <- factor(newdata1[,v], levels=mod)
                        }
                }

            ##cas ou on a factor() dans l'appel
            z <- all.names(as.formula(paste("~",call_fixed)))
            ind_factor <- which(z=="factor")
            if(length(ind_factor))
                {
                    nom.factor <- z[ind_factor+1]
                    for (v in nom.factor)
                        {
                            mod <- levels(as.factor(olddata[,v]))
                            if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
                            newdata1[,v] <- factor(newdata1[,v], levels=mod)

                            fv <- paste("factor\\(",v,"\\)",sep="")
                            if(length(grep(fv,model$Names$Xnames)))
                                {
                                    model$Names$Xnames <- gsub(fv,v,model$Names$Xnames)
                                }
                        }
                }
            call_fixed <- gsub("factor","",call_fixed)

            z <- all.names(as.formula(paste("~",call_random)))
            ind_factor <- which(z=="factor")
            if(length(ind_factor))
                {
                    nom.factor <- z[ind_factor+1]
                    for (v in nom.factor)
                        {
                            mod <- levels(as.factor(olddata[,v]))
                            if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
                            newdata1[,v] <- factor(newdata1[,v], levels=mod)

                            fv <- paste("factor\\(",v,"\\)",sep="")
                            if(length(grep(fv,model$Names$Xnames)))
                                {
                                    model$Names$Xnames <- gsub(fv,v,model$Names$Xnames)
                                }
                        }
                }
            call_random <- gsub("factor","",call_random)

            z <- all.names(as.formula(paste("~",call_classmb)))
            ind_factor <- which(z=="factor")
            if(length(ind_factor))
                {
                    nom.factor <- z[ind_factor+1]
                    for (v in nom.factor)
                        {
                            mod <- levels(as.factor(olddata[,v]))
                            if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
                            newdata1[,v] <- factor(newdata1[,v], levels=mod)

                            fv <- paste("factor\\(",v,"\\)",sep="")
                            if(length(grep(fv,model$Names$Xnames)))
                                {
                                    model$Names$Xnames <- gsub(fv,v,model$Names$Xnames)
                                }
                        }
                }
            call_classmb <- gsub("factor","",call_classmb)

            z <- all.names(as.formula(paste("~",call_survival)))
            ind_factor <- which(z=="factor")
            if(length(ind_factor))
                {
                    nom.factor <- z[ind_factor+1]
                    for (v in nom.factor)
                        {
                            mod <- levels(as.factor(olddata[,v]))
                            if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
                            newdata1[,v] <- factor(newdata1[,v], levels=mod)

                            fv <- paste("factor\\(",v,"\\)",sep="")
                            if(length(grep(fv,model$Names$Xnames)))
                                {
                                    model$Names$Xnames <- gsub(fv,v,model$Names$Xnames)
                                }
                        }
                }
            call_survival <- gsub("factor","",call_survival)



### Traitement des donnees manquantes
            mcall <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
            mcall$na.action <- na.action
            mcall$data <- newdata1

            ## fixed
            m <- mcall
            m$formula <- formula(paste("~",call_fixed,sep=""))
            m[[1]] <- as.name("model.frame")
            m <- eval(m, sys.parent())
            na.fixed <- attr(m,"na.action")

            ## random
            if(!is.null(model$call$random))
                {
                    m <- mcall
                    m$formula <- formula(paste("~",call_random,sep=""))
                    m[[1]] <- as.name("model.frame")
                    m <- eval(m, sys.parent())
                    na.random <- attr(m,"na.action")
                }
            else
                {
                    na.random <- NULL
                }

            ## classmb
            if(!is.null(model$call$classmb))
                {
                    m <- mcall
                    m$formula <- formula(paste("~",call_classmb,sep=""))
                    m[[1]] <- as.name("model.frame")
                    m <- eval(m, sys.parent())
                    na.classmb <- attr(m,"na.action")
                }
            else
                {
                    na.classmb <- NULL
                }

            ##survival
            if(!is.null(model$call$survival))
                {
                    m <- mcall
                    m$formula <- formula(paste("~",call_survival,sep=""))
                    m[[1]] <- as.name("model.frame")
                    m <- eval(m, sys.parent())
                    na.survival <- attr(m,"na.action")
                }
            else
                {
                    na.survival <- NULL
                }
            
            ##cor
            na.cor <- NULL
            if(model$N[7]>0)
                {
                    z <- which(model$idcor==1)
                    var.cor <- newdata1[,model$Names$Xnames[z]]
                    na.cor <- which(is.na(var.cor))
                }
       
            
            ##var.time
            na.time <- NULL
            if(!(var.time %in% model$Names[[1]]))
                {
                    na.time <- which(is.na(newdata1[,var.time]))
                }
            
            ## toutes les variables pas var expl
            dtmp <- newdata1[,setdiff(modelNames,Xnames2),drop=FALSE]

            m <- model.frame(formula=formula(paste(model$Names$Yname,"~",paste(setdiff(modelNames,Xnames2),collapse="+"),sep="")),data=dtmp,na.action=na.action)
            na.other <- attr(m,"na.action")
            
            ## fun.time 
            if(length(match.call()$fun.time))
                {   
                    if(!isTRUE(as.character(match.call()[[which(names(match.call())=="fun.time")]])=="identity"))
                        {
                            z <- match.call()$fun.time 
                            
                            if(isTRUE(as.character(z) %in% ls(.GlobalEnv)))
                                {
                                    bodyfun <- as.expression(body(get(as.character(z))))
                                }
                            else
                                {
                                    bodyfun <- as.expression(z[[3]])
                                }
                            vars <- intersect(colnames(newdata1),all.names(bodyfun))
                            getvars <- paste("getElement(newdata1,'",vars,"')",sep="")
                            
                            if(length(vars))
                                {
                                    for(i in 1:length(vars))
                                        {
                                            bodyfun <- sub(vars[i],getvars[i],bodyfun)
                                            bodyfun <- sub(as.character(match.call()$newdata),"",bodyfun)
                                            tmp <- strsplit(bodyfun,split="")[[1]]
                                            dollard <- which(tmp=="$")
                                            if(length(dollard))
                                                {
                                                    tmp <- tmp[-dollard]
                                                }
                                            bodyfun <- paste(tmp,collapse="")
                                        }
                                }

                            headfun <- paste("function(",names(formals(fun.time)),")",collapse="")
                            ff1 <- paste(headfun,bodyfun)
                            ff2 <- parse(text=ff1)
                            ff3 <- eval(ff2)

                            timemes <- do.call("ff3",list(newdata1[,var.time]))
                        }
                    else
                        {
                            timemes <- newdata1[,var.time]
                        } 
                }
            else
                {
                    timemes <- newdata1[,var.time]
                }
            na.fun <- which(is.na(timemes))
            
            ## Table sans donnees manquante: newdata1
            na.action <- unique(c(na.other,na.fixed,na.random,na.classmb,na.survival,na.cor,na.time,na.fun))
            if(length(na.action))
                {
                    newdata1 <- newdata1[-na.action,]
                }
            
            ## si Tevent et Devent dans newdata
            Tevent <- NULL
            Devent <- NULL
            
            if(length(model$Names$Tnames)==2)  # ordre : Tevent,Event 
                {
                    if(all(model$Names$Tnames[1:2] %in% colnames(newdata1)))
                        {
                            Tevent <- newdata1[,model$Names$Tnames[1]]
                            Devent <- newdata1[,model$Names$Tnames[2]]
                            
                            if(length(na.action))
                                {
                                    Tevent <- Tevent[-na.action]
                                    Devent <- Devent[-na.action]
                                }
                            
                            indNA <- which(is.na(Tevent) | is.na(Devent))
                            if(length(indNA)) 
                                {
                                    Tevent[indNA] <- 0
                                    Devent[indNA] <- 0
                                } 
                        }
                }

            if(length(model$Names$Tnames)==3)  # ordre : Tentry,Tevent,Event 
                {
                    if(all(model$Names$Tnames[2:3] %in% colnames(newdata1)))
                        {
                            Tevent <- newdata1[,model$Names$Tnames[2]]
                            Devent <- newdata1[,model$Names$Tnames[3]]
                            
                            if(!is.null(na.action))
                                {
                                    Tevent <- Tevent[-na.action]
                                    Devent <- Devent[-na.action]
                                }  
                            
                            indNA <- which(is.na(Tevent) | is.na(Devent))
                            if(length(indNA)) 
                                {
                                    Tevent[indNA] <- 0
                                    Devent[indNA] <- 0
                                }  
                        }
                }
            
            ## nb de sujets
            id.subject <- newdata1[,model$Names$ID]
            if(class(id.subject)=="factor") id.subject <- factor(id.subject,levels=unique(id.subject))
            ns <- length(unique(id.subject))
            
            ## vecteur Y et indiceY,uniqueY,nvalSPL
            Y <- newdata1[,model$Names$Yname]


            uniqueY <- 0
            indiceY <- 0
            nvalSPL <- 0

            if(idlink == 2)
                {
                    uniqueY <- sort(unique(Y))
                    permut <- order(order(Y))  # sort(y)[order(order(y))] = y
                    indice <- rep(1:length(uniqueY), as.vector(table(Y)))
                    indiceY <- indice[permut]
                    nvalSPL <- length(uniqueY)
                }
            
            ## nb d'observations
            nobs <- length(Y)
            
            ## nb de mesures par sujet
            nmes <- as.vector(table(id.subject))
            
            ## prior
            prior <- rep(0,ns)
            if(length(model$Names$prior.name))
                {
                    prior <- newdata1[cumsum(nmes),model$Names$prior.name]
                }

            
            ## age entree si type = counting
            tsurv0 <- rep(0,ns)
            idtrunc <- 0
            if(length(model$Names$Tnames)>2)
                {
                    idtrunc <- 1
                    tsurv0 <- newdata1[cumsum(nmes),model$Names$Tnames[1]] 
                }

                           
            ## reduire Tevent et Devent au nb de sujets
            if(!is.null(Tevent))
                {
                    Tevent <- Tevent[cumsum(nmes)]
                    Devent <- Devent[cumsum(nmes)]   
                } 


            ## Construction de nouvelles var explicatives sur la nouvelle table :

            X_intercept <- model.matrix(~1,data=newdata1)
            
            ## fixed

            X_fixed <- model.matrix(formula(paste("~",call_fixed,sep="")),data=newdata1)
            if(colnames(X_fixed)[1]=="(Intercept)")
                {
                    colnames(X_fixed)[1] <- "intercept"
                }

            ## mixture pas besoin car les variables sont dans fixed
            
            ## random
            if(!is.null(model$call$random))
                {
                    X_random <- model.matrix(formula(paste("~",call_random,sep="")),data=newdata1)
                    if(colnames(X_random)[1]=="(Intercept)")
                        {
                            colnames(X_random)[1] <- "intercept"
                        }
                }
            else
                {
                    X_random <- NULL
                }
            
            ## classmb
            if(!is.null(model$call$classmb))
                {
                    X_classmb <- model.matrix(formula(paste("~",call_classmb,sep="")),data=newdata1)
                    colnames(X_classmb)[1] <- "intercept"
                }
            else
                {
                    X_classmb <- NULL
                }

            ## survival
            if(!is.null(model$call$survival))
                {
                    X_survival <- model.matrix(formula(paste("~",call_survival,sep="")),data=newdata1)
                    if(ncol(X_survival)) colnames(X_survival)[1] <- "intercept"
                }
            else
                {
                   X_survival <- NULL
                }

            ##cor
            if(model$N[7]>0)  #on reprend la variable de temps de cor
                {
                    z <- which(model$idcor==1)
                    X_cor <- newdata1[,model$Names$Xnames[z]]
                }
            else
                {
                    X_cor <- NULL
                }
            
            
            ## var.time
            if(length(na.action)) timemes <- timemes[-na.action]
                                        #  temps deja dans la bonne echelle
            

            ## Construction de newdata1 dans le bon ordre
            X <- cbind(X_intercept,X_fixed,X_random,X_classmb,X_survival,X_cor)
            colX <- strsplit(colnames(X),split=":",fixed=TRUE)
            colX <- lapply(colX,sort)
            colX <- lapply(colX,paste,collapse=":")
            colnames(X) <- unlist(colX)

            Xnames <- strsplit(model$Names$Xnames,split=":",fixed=TRUE)
            Xnames <- lapply(Xnames,sort)
            Xnames <- lapply(Xnames,paste,collapse=":")
            Xnames <- unlist(Xnames)
            
            X <- X[,Xnames,drop=FALSE]


            
########## fonction de risque ##########
            
            risq0 <- function(t,evt,typrisq,brisq,zi,nz,bPH,logspecif)
                {
                    if(typrisq[evt]==2)
                        {
                            if(logspecif==1)
                                {
                                    res <- brisq[2,evt]*brisq[1,evt]*t**(brisq[2,evt]-1)
                                }
                            if(logspecif==0)
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
            risqcum0 <- function(t,evt,typrisq,brisq,zi,nz,bPH,logspecif)
                {
                    if(typrisq[evt]==2)
                        {
                            if(logspecif==1)
                                {
                                    res <-  brisq[1,evt]*t**brisq[2,evt]
                                }
                            if(logspecif==0)
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
            fct_incid <- function(t,evt,typrisq,brisq,zi,nz,bPH,Xevt,bevt,nbevt,bevtint,logspecif)
                {
                    risq_tevt <- risq0(t,evt,typrisq,brisq,zi,nz,bPH,logspecif)*exp(sum(Xevt*bevt[,evt]))*exp(bevtint[evt])

                    somme <- 0
                    for(ke in 1:nbevt)
                        {
                            somme <- somme + risqcum0(t,ke,typrisq,brisq,zi,nz,bPH,logspecif)*exp(sum(Xevt*bevt[,ke]))*exp(bevtint[ke])
                        }

                    return(risq_tevt * exp(-somme) )               
                }
######## fin fct incid

            

#### matrice Xevt a ns lignes

            Xevt_ns <- X[cumsum(nmes),which(idsurv!=0),drop=FALSE]

### ok differents profils reperes

            
            ## fonction d'integration
            integrate2 <- function(...) return(integrate(...)$value)


            
            if(!isTRUE(draws))
                {

                    ## calcul des proba a posteriori pour tous les temps s   
                    ppi <- rep(0,ns*nbland*ng)
                    
                    out <- .Fortran(C_postprob2,
                                    as.double(Y),
                                    as.double(X),
                                    as.integer(ns),
                                    as.integer(nmes),
                                    as.integer(nobs),
                                    as.integer(ng),
                                    as.integer(nv),
                                    as.integer(idiag),
                                    as.integer(nwg),
                                    as.integer(ncor),
                                    as.integer(logspecif),
                                    as.double(zi),
                                    as.integer(idea),
                                    as.integer(idg),
                                    as.integer(idprob),
                                    as.integer(idcor),
                                    as.integer(idcom),
                                    as.integer(idspecif2),
                                    as.integer(idtdv),
                                    as.integer(risqcom),
                                    as.integer(nvdepsurv),
                                    as.integer(typrisq),
                                    as.integer(nz),
                                    as.integer(nprisq),
                                    as.integer(nbzitr),
                                    as.double(zitr),
                                    as.integer(idlink),
                                    as.integer(indiceY),
                                    as.double(uniqueY),
                                    as.integer(nvalSPL),
                                    as.double(epsY),
                                    as.integer(nbevt),
                                    as.double(tsurv0),
                                    as.integer(idtrunc),
                                    as.double(best),
                                    as.integer(npm),
                                    as.double(timemes),
                                    as.double(landmark),
                                    as.integer(nbland),
                                    ppi=as.double(ppi))
                    
                    out$ppi[which(out$ppi==-5)] <- NA
                    ppi <- matrix(out$ppi,ns*nbland,ng)



                    ## incid et surv a remplir
                    incid <- array(NA,c(ns*nbland,ng,length(horizon)))
                    #surv <- array(1,c(ns*nbland,ng))
        
                    ## calcul de incidcum et survie 
                    for(ks in 1:nbland)
                        {
                            for(i in 1:nrow(Xevt_ns))
                                {
                                    ## si sujet a 0 mesure, suivant
                                    if(nmes[i]==0) next
                                    
                                    if(nvdepsurv>0)
                                        {
                                            itimedepvar <- which(colnames(Xevt_ns)==model$Names$TimeDepVar.name)
                                            tint <- Xevt[i,itimedepvar]
                                            Xevt <- as.vector(Xevt_ns[i,-itimedepvar])
                                        }
                                    else
                                        {
                                            Xevt <- as.vector(Xevt_ns[i,])
                                            tint <- max(landmark)+1
                                        }
                                    
                                    brisqtot <- as.vector(best[model$N[1]+1:model$N[2]])
                                    bvarxevt <- as.vector(best[sum(model$N[1:2])+1:model$N[3]])
                       
                                    for (g in 1:ng)
                                        {
                                            ## coef pour la classe g
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
                                                    if(logspecif==1) brisq <- exp(brisqtmp)
                                                    if(logspecif==0) brisq <- brisqtmp**2
                                                }

                                            if(is.list(brisqtmp))
                                                {
                                                    lmax <- max(sapply(brisqtmp,length))
                                                    
                                                    if(model$logspecif==1)
                                                        {                
                                                            brisq <- lapply(brisqtmp, function(l) c(exp(l),rep(0,lmax-length(l))))
                                                            brisq <- matrix(unlist(brisq),nrow=lmax,ncol=nbevt)
                                                        }

                                                    if(model$logspecif==0)
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


                                           
                                            ## incidcum entre s et s+t (ou difference des survies)
                                            if(nbevt==1)
                                                {
                                                    bevtint2 <- bevtint
                                                    if(nvdepsurv!=0 & tint>landmark[ks]) bevtint2 <- rep(0,nbevt)

                                                    surv_s <- exp(-risqcum0(t=landmark[ks],evt=1,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH,logspecif=logspecif)*exp(sum(Xevt*bevt[,1])))
                                                    surv_tint <- exp(-risqcum0(t=tint,evt=1,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH,logspecif=logspecif)*exp(sum(Xevt*bevt[,1])))                                     

                                                    if(nvdepsurv!=0 & tint<landmark[ks])
                                                        {
                                                            incid[ns*(ks-1)+i,g,] <- sapply(landmark[ks]+horizon,function(h) (surv_tint * (surv_s-surv_tint)*exp(bevtint[1])) - surv_tint*(exp(-risqcum0(t=h,evt=1,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH,logspecif=logspecif)*exp(sum(Xevt*bevt[,1]))*exp(bevtint[1]))-surv_tint*exp(bevtint[1])))
                                                        }
                                                    else
                                                        {
                                                            incid[ns*(ks-1)+i,g,] <- sapply(landmark[ks]+horizon,function(h) surv_s-exp(-risqcum0(t=h,evt=1,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH,logspecif=logspecif)*exp(sum(Xevt*bevt[,1]))))
                                                        }
                                                }
                                            else
                                                {
                                                    for(ke in 1:nbevt)
                                                        {
                                                            bevtint2 <- bevtint
                                                            if(nvdepsurv!=0 & tint>landmark[ks]) bevtint2 <- rep(0,nbevt)
                                                            if(ke==event)
                                                                {
                                                                    incid[ns*(ks-1)+i,g,] <- sapply(landmark[ks]+horizon,integrate2,f=fct_incid,lower=landmark[ks],evt=ke,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH,Xevt=Xevt,bevt=bevt,nbevt=nbevt,bevtint=bevtint2,logspecif=logspecif)
                                                                }
                                                            
                                        #if(nvdepsurv!=0 & tint<landmark[ks])
                                        #    {
                                        #        surv[ns*(ks-1)+i,g] <- surv[ns*(ks-1)+i,g]*exp(-risqcum0(t=tint,evt=ke,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH)*exp(sum(Xevt*bevt[,ke])))*(exp(-risqcum0(t=landmark[ks],evt=ke,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH)*exp(sum(Xevt*bevt[,ke]))*exp(bevtint[ke]))-exp(-risqcum0(t=tint,evt=ke,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH)*exp(sum(Xevt*bevt[,ke]))*exp(bevtint[ke])))
                                        #   }
                                                    #else
                                        #    {
                                        #        surv[ns*(ks-1)+i,g] <- surv[ns*(ks-1)+i,g]*exp(-risqcum0(t=landmark[ks],evt=ke,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH)*exp(sum(Xevt*bevt[,ke])))
                                        #    }
                                                        }
                                                }
                                            
                                        } # fin boucle g
                                   
                                }# fin boucle i
                                   
                        } # fin boucle ks


                    res <- data.frame(rep(sort(unique(id.subject)),each=nbland*nbhoriz))
                    res <- cbind(res,rep(landmark,each=nbhoriz))
                    res <- cbind(res,horizon)
                    res <- cbind(res,NA)

                    colnames(res) <- c(model$Names$ID,"landmark","horizon","pred")

                    pred <- NULL
                    #pred_g <- incid*array(ppi/surv,c(ns*nbland,ng,nbhoriz))
                    pred_g <- incid*array(ppi,c(ns*nbland,ng,nbhoriz))
                    pred_evt <- apply(pred_g,c(1,3),sum) # res de ie horiz dans ie colonne
                    pred_evt <- pred_evt[order(rep(1:ns,nbland)),]
                    pred_evt <- as.vector(t(pred_evt))
                    pred <- c(pred,pred_evt)
                            
                    
                    res[,4] <- pred

                    ## mettre NA si Devent!=0 et Tevent<s
                    if(!is.null(Tevent))
                        {
                            matevt <- cbind(rep(Tevent,each=nbland*nbhoriz),rep(Devent,each=nbland*nbhoriz))
                            indevt <- which( (matevt[,2]!=0) & (matevt[,1]<res[,"landmark"]) )
                            if(length(indevt)) res[indevt,"pred"] <- NA
                        }

                }
            else   # ie avec draws
                {
                    posfix <- eval(model$call$posfix)

                    Mat <- matrix(0,ncol=npm,nrow=npm)
                    Mat[upper.tri(Mat,diag=TRUE)]<- model$V
                    if(length(posfix))
                        {
                            Mat2 <- Mat[-posfix,-posfix]
                            Chol2 <- chol(Mat2)
                            Chol <- matrix(0,npm,npm)
                            Chol[setdiff(1:npm,posfix),setdiff(1:npm,posfix)] <- Chol2
                            Chol <- t(Chol)
                        }
                    else
                        {
                            Chol <- chol(Mat)
                            Chol <- t(Chol)
                        }
                	
                    
                    doOneDraw <- function()
                        {
                            bdraw <- rnorm(npm)
                            bdraw <- best + Chol %*% bdraw
                            
                            ## calcul des proba a posteriori pour tous les temps s   
                            ppi <- rep(0,ns*nbland*ng)
                            
                            out <- .Fortran(C_postprob2,
                                            as.double(Y),
                                            as.double(X),
                                            as.integer(ns),
                                            as.integer(nmes),
                                            as.integer(nobs),
                                            as.integer(ng),
                                            as.integer(nv),
                                            as.integer(idiag),
                                            as.integer(nwg),
                                            as.integer(ncor),
                                            as.integer(logspecif),
                                            as.double(zi),
                                            as.integer(idea),
                                            as.integer(idg),
                                            as.integer(idprob),
                                            as.integer(idcor),
                                            as.integer(idcom),
                                            as.integer(idspecif2),
                                            as.integer(idtdv),
                                            as.integer(risqcom),
                                            as.integer(nvdepsurv),
                                            as.integer(typrisq),
                                            as.integer(nz),
                                            as.integer(nprisq),
                                            as.integer(nbzitr),
                                            as.double(zitr),
                                            as.integer(idlink),
                                            as.integer(indiceY),
                                            as.double(uniqueY),
                                            as.integer(nvalSPL),
                                            as.double(epsY),
                                            as.integer(nbevt),
                                            as.double(tsurv0),
                                            as.integer(idtrunc),
                                            as.double(bdraw),
                                            as.integer(npm),
                                            as.double(timemes),
                                            as.double(landmark),
                                            as.integer(nbland),
                                            ppi=as.double(ppi)) 
                            
                            out$ppi[which(out$ppi==-5)] <- NA
                            if(all(is.na(out$ppi))) return(rep(NA,ns*nbland*nbhoriz*nbevt))
                            ppi <- matrix(out$ppi,ns*nbland,ng)


                            ## incid et surv a remplir
                            incid <- array(NA,c(ns*nbland,ng,length(horizon)))
                            #surv <- array(1,c(ns*nbland,ng))
                            
                            ## calcul de incidcum et survie 
                            for(ks in 1:nbland)
                                {
                                    for(i in 1:nrow(Xevt_ns))
                                        {
                                            ## si sujet a 0 mesure,suivant
                                            if(nmes[i]==0) next

                                            if(nvdepsurv>0)
                                                {
                                                    itimedepvar <- which(colnames(Xevt_ns)==model$Names$TimeDepVar.name)
                                                    tint <- Xevt[i,itimedepvar]
                                                    Xevt <- as.vector(Xevt_ns[i,-itimedepvar])
                                                }
                                            else
                                                {
                                                    Xevt <- as.vector(Xevt_ns[i,])
                                                    tint <- max(landmark)+1
                                                }
                                            
                                            brisqtot <- as.vector(bdraw[model$N[1]+1:model$N[2]])
                                            bvarxevt <- as.vector(bdraw[sum(model$N[1:2])+1:model$N[3]])


                                            for (g in 1:ng)
                                                {
                                                    ## coef pour la classe g
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
                                                            if(logspecif==1) brisq <- exp(brisqtmp)
                                                            if(logspecif==0) brisq <- brisqtmp**2                                                   
                                                        }
                                                    
                                                    if(is.list(brisqtmp))
                                                        {
                                                            lmax <- max(sapply(brisqtmp,length))

                                                            if(model$logspecif==1)
                                                                {                
                                                                    brisq <- lapply(brisqtmp, function(l) c(exp(l),rep(0,lmax-length(l))))
                                                                    brisq <- matrix(unlist(brisq),nrow=lmax,ncol=nbevt)
                                                                }

                                                            if(model$logspecif==0)
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


                                                    
                                                    ## incidcum entre s et s+t et survie en s
                                                    if(nbevt==1)
                                                        {
                                                            bevtint2 <- bevtint
                                                            if(nvdepsurv!=0 & tint>landmark[ks]) bevtint2 <- rep(0,nbevt)
                                                            
                                                            surv_s <- exp(-risqcum0(t=landmark[ks],evt=1,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH,logspecif=logspecif)*exp(sum(Xevt*bevt[,1])))
                                                            surv_tint <- exp(-risqcum0(t=tint,evt=1,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH,logspecif=logspecif)*exp(sum(Xevt*bevt[,1])))
               
                                                    
                                                            if(nvdepsurv!=0 & tint<landmark[ks])
                                                                {
                                                                    incid[ns*(ks-1)+i,g,] <- sapply(landmark[ks]+horizon,function(h) (surv_tint * (surv_s-surv_tint)*exp(bevtint[1])) - surv_tint*(exp(-risqcum0(t=h,evt=1,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH,logspecif=logspecif)*exp(sum(Xevt*bevt[,1]))*exp(bevtint[1]))-surv_tint*exp(bevtint[1])))
                                                                }
                                                            else
                                                                {
                                                                    incid[ns*(ks-1)+i,g,] <- sapply(landmark[ks]+horizon,function(h) surv_s-exp(-risqcum0(t=h,evt=1,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH,logspecif=logspecif)*exp(sum(Xevt*bevt[,1]))))
                                                                }
                                                        }
                                                    else
                                                        {
                                                            for(ke in 1:nbevt)
                                                                {
                                                                    bevtint2 <- bevtint
                                                                    if(nvdepsurv!=0 & tint>landmark[ks]) bevtint2 <- rep(0,nbevt)
                                                                    if(ke==event)
                                                                        {
                                                                            incid[ns*(ks-1)+i,g,] <- sapply(landmark[ks]+horizon,integrate2,f=fct_incid,lower=landmark[ks],evt=ke,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH,Xevt=Xevt,bevt=bevt,nbevt=nbevt,bevtint=bevtint2,logspecif=logspecif)
                                                                        }
                                                            
                                                            #if(nvdepsurv!=0 & tint<landmark[ks])
                                                            #    {
                                                            #        surv[ns*(ks-1)+i,g] <- surv[ns*(ks-1)+i,g]*exp(-risqcum0(t=tint,evt=ke,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH)*exp(sum(Xevt*bevt[,ke])))*(exp(-risqcum0(t=landmark[ks],evt=ke,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH)*exp(sum(Xevt*bevt[,ke]))*exp(bevtint[ke]))-exp(-risqcum0(t=tint,evt=ke,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH)*exp(sum(Xevt*bevt[,ke]))*exp(bevtint[ke])))
                                                            #    }
                                                            #else
                                                            #    {
                                                            #        surv[ns*(ks-1)+i,g] <- surv[ns*(ks-1)+i,g]*exp(-risqcum0(t=landmark[ks],evt=ke,typrisq=typrisq,brisq=brisq,zi=zi,nz=nz,bPH=bPH)*exp(sum(Xevt*bevt[,ke])))
                                                             #   }
                                                                }
                                                        }
                                                    
                                                } # fin boucle g
                                            
                                        }# fin boucle i
                                    
                                } # fin boucle ks
                            



                            
                            res <- data.frame(rep(sort(unique(id.subject)),each=nbland*nbhoriz))
                            res <- cbind(res,rep(landmark,each=nbhoriz))
                            res <- cbind(res,horizon)
                            res <- cbind(res,NA)

                            colnames(res) <- c(model$Names$ID,"landmark","horizon","pred")

                            pred <- NULL
                            #pred_g <- incid*array(ppi/surv,c(ns*nbland,ng,nbhoriz))
                            pred_g <- incid*array(ppi,c(ns*nbland,ng,nbhoriz))
                            pred_evt <- apply(pred_g,c(1,3),sum) # res de ie horiz dans ie colonne
                            pred_evt <- pred_evt[order(rep(1:ns,nbland)),]
                            pred_evt <- as.vector(t(pred_evt))
                            pred <- c(pred,pred_evt)
                                    

                            res[,4] <- pred

                            ## mettre NA si Devent!=0 et Tevent<s
                            if(!is.null(Tevent))
                                {
                                    matevt <- cbind(rep(Tevent,each=nbland*nbhoriz),rep(Devent,each=nbland*nbhoriz))
                                    indevt <- which( (matevt[,2]!=0) & (matevt[,1]<res[,"landmark"]) )
                                    if(length(indevt)) res[indevt,"pred"] <- NA
                                }

                            
                            return(res[,"pred"])                          
                        }

                    
                    ndraws <- as.integer(ndraws)
                    
                    resdraws <- replicate(ndraws,doOneDraw())

                    if(is.vector(resdraws))
                        {
                            resdraws <- matrix(resdraws,nrow=1)
                        }
                    #probs <- resdraws[-1,,drop=FALSE]
                    #pb <- sum(resdraws[1,])
                    
                    #if(pb>0) warning("Infinite probabilities have been found. Confidence intervals are based on ", ndraws-pb ," simulations.")
                    
                    med <- apply(resdraws,1,median,na.rm=TRUE)
                    qmin <- apply(resdraws,1,quantile,prob=0.025,na.rm=TRUE)
                    qmax <- apply(resdraws,1,quantile,prob=0.975,na.rm=TRUE)


                    res <- data.frame(rep(sort(unique(id.subject)),each=nbland*nbhoriz))
                    res <- cbind(res,rep(landmark,each=nbhoriz))
                    res <- cbind(res,horizon)
                    res <- cbind(res,med,qmin,qmax)

                    colnames(res) <- c(model$Names$ID,"landmark","horizon","pred_50","pred_2.5","pred_97.5")                   
                }
        }
    else  # ie conv != 1 ou 2 ou 3
        {
            cat("Output can not be produced since the program stopped abnormally. \n")
            res <- NA
            id.subject <- NA
            nmes <- NA
            Y <- NA
            timemes <- NA 
        }
    
    prmatrix(res,quote=FALSE)  
    
    res.list <- list(pred=res,newdata=data.frame(id=rep(unique(id.subject),nmes)[which(timemes<=max(landmark))],y=Y[which(timemes<=max(landmark))],time=timemes[which(timemes<=max(landmark))]))
    class(res.list) <- "dynpred"
    
    return(invisible(res.list))
}
