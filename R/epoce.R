##  -> epoce implemente uniquement pr 1 evt sans transfo





#' Estimators of the Expected Prognostic Observed Cross-Entropy (EPOCE) for
#' evaluating predictive accuracy of joint latent class models estimated using
#' \code{Jointlcmm}
#' 
#' This function computes estimators of the Expected Prognostic Observed
#' Cross-Entropy (EPOCE) for evaluating the predictive accuracy of joint latent
#' class models estimated using \code{Jointlcmm}. On the same data as used for
#' estimation of the \code{Jointlcmm} object, this function computes both the
#' Mean Prognostic Observed Log-Likelihood (MPOL) and the Cross-Validated
#' Observed Log-Likelihood (CVPOL), two estimators of EPOCE. The latter
#' corrects the MPOL estimate for over-optimism by approximated
#' cross-validation. On external data, this function only computes the Mean
#' Prognostic Observed Log-Likelihood (MPOL).
#' 
#' This function does not apply for the moment with multiple causes of event
#' (competing risks).
#' 
#' EPOCE assesses the prognostic information of a joint latent class model. It
#' relies on information theory.
#' 
#' MPOL computed at time s equals minus the mean individual contribution to the
#' conditional log-likelihood of the time to event given the longitudinal data
#' up to the time of prediction s and given the subject is still at risk of
#' event in s.
#' 
#' CVPOL computed at time s equals MPOL at time s plus a penalty term that
#' corrects for over-optimism when computing predictive accuracy measures on
#' the same dataset as used for estimation. This penalty term is computed from
#' the inverse of the Hessian of the joint log-likelihood and the product of
#' the gradients of the contributions to respectively the joint log-likelihood
#' and the conditional log-likelihood.
#' 
#' The theory of EPOCE and its estimators MPOL and CVPOL is given in Commenges
#' et al. (2012), and further detailed and illustrated for joint models in
#' Proust-Lima et al. (2013).
#' 
#' @param model an object inheriting from class \code{Jointlcmm}
#' @param pred.times Vector of times of prediction, from which predictive
#' accuracy is evaluated (only subjects still at risk at the time of prediction
#' are included in the computation, and only information before the time of
#' prediction is considered.
#' @param var.time Name of the variable indicating time in the dataset
#' @param fun.time an optional function. This is only required if the time
#' scales in the longitudinal part of the model and the survival part are
#' different. In that case, \code{fun.time} is the function that translates the
#' times from the longitudinal part into the time scale of the survival part.
#' The default is the identity function which means that the two time scales
#' are the same.
#' @param newdata optional. When missing, the data used for estimating the
#' \code{Jointlcmm} object are used, and CVPOL and MPOL are computed (internal
#' validation). When newdata is specified, only MPOL is computed on this
#' newdataset (external validation).
#' @param subset a specification of the rows to be used: defaults to all rows.
#' This can be any valid indexing vector for the rows of data or if that is not
#' supplied, a data frame made up of the variable used in formula.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @return \item{call.Jointlcmm}{the \code{Jointlcmm} call}
#' \item{call.epoce}{the matched call} \item{EPOCE}{Dataframe containing, for
#' each prediction time s, the number of subjects still at risk at s (and with
#' at least one measure before s), the number of events after time s, the MPOL,
#' and the CVPOL when computation is done on the dataset used for
#' \code{Jointlcmm} estimation} \item{IndivContrib}{Individual contributions to
#' the prognostic observed log-likelihood at each time of prediction. Used for
#' computing tracking intervals of EPOCE differences between models.}
#' \item{new.data}{a boolean for internal use only, which is FALSE if
#' computation is done on the same data as for \code{Jointlcmm} estimation, and
#' TRUE otherwise.}
#' @author Cecile Proust-Lima and Amadou Diakite
#' @seealso
#' \code{\link{Jointlcmm}}, \code{\link{print.epoce}}, \code{\link{summary.epoce}}, \code{\link{plot.epoce}}
#' @references Commenges, Liquet and Proust-Lima (2012). Choice of prognostic
#' estimators in joint models by estimating differences of expected conditional
#' Kullback-Leibler risks. Biometrics 68(2), 380-7.
#' 
#' Proust-Lima, Sene, Taylor and Jacqmin-Gadda (2014). Joint latent class
#' models of longitudinal and time-to-event data: a review. Statistical Methods
#' in Medical Research 23, 74-90.
#' @examples
#' 
#' \dontrun{
#' ## estimation of a joint latent class model with 2 latent classes (ng=2)
#' # (see the example section of Jointlcmm for details about
#' #  the model specification)
#' 
#' m <- Jointlcmm(fixed= Ydep1~Time*X1,random=~Time,mixture=~Time,subject='ID'
#' ,survival = Surv(Tevent,Event)~ X1+X2 ,hazard="Weibull"
#' ,hazardtype="PH",ng=2,data=data_lcmm,logscale=TRUE,
#' B=c(0.7608, -9.4974 , 1.0242,  1.4331 , 0.1063 , 0.6714, 10.4679, 11.3178,
#'  -2.5671, -0.5386,  1.4616, -0.0605,  0.9489,  0.1020 , 0.2079,  1.5045))
#' summary(m)
#' 
#' ## Computation of the EPOCE on the same dataset as used for
#' # estimation of m with times at predictions from 1 to 15 
#' VecTime <- c(1,3,5,7,9,11,13,15)
#' cvpl <- epoce(m,var.time="Time",pred.times=VecTime)
#' summary(cvpl)
#' plot(cvpl,bty="l",ylim=c(0,2))
#' }
#' 
#' 
#' @export
#' 
#' 
#' 
epoce <- function(model,pred.times,var.time,fun.time=identity,newdata=NULL,subset=NULL,na.action=1)
{
    cl <- match.call()
    if(missing(var.time)) stop("The argument var.time should be specified")
    if (!inherits(var.time, "character")) stop("the class of var.time should be character")
    if(missing(model)) stop("The argument model must be specified")
    if(!inherits(model,"Jointlcmm")) stop("The argument model must be a class 'Jointlcmm'")
    if(!is.function(fun.time)) stop("'fun.time' should be a function")

    if(!missing(newdata))
        {
            if(!is.data.frame(newdata)) stop("The argument newdata must be a 'data.frame'")
            nomsvar <- unique(c(model$Names$Xnames2,model$Names$Ynames,model$Names$ID,model$Names$Tnames,model$Names$prior.name,model$Names$TimeDepVar.name))
            if(!all(nomsvar %in% colnames(newdata))) stop(paste("newdata should include the following covariates: \n",paste(nomsvar,collapse=" ")))
            if(!(var.time %in% names(newdata))) stop("The new dataset should contain the var.time argument")	
            
            new.data <- TRUE
        } 


    if(!(na.action%in%c(1,2))) stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 


    ## pas implemente pr competitif ni pr link non lineaire
    if(length(model$N)>10) stop("EPOCE is not implemented for competing risks yet")
    if(model$linktype!=-1) stop("EPOCE is not implemented for models with a link function yet")



    call_fixed <- model$call$fixed[3]
    if(is.null(model$call$random)) {call_random <- -1} else call_random <- model$call$random[2]
    if(is.null(model$call$classmb)) {call_classmb <- -1} else call_classmb <- model$call$classmb[2]
    if(is.null(model$call$mixture)) {call_mixture <- -1} else call_mixture <- model$call$mixture[2]
    if(is.null(model$call$survival)) {call_survival <- -1} else call_survival <- model$call$survival[3]

    nbevt <- length(model$hazard[[1]])
    call_survival <- gsub("mixture","",call_survival)
    for(k in 1:nbevt)
        {
            call_survival <- gsub(paste("cause",k,sep=""),"",call_survival)
        }
    call_survival <- gsub("cause","",call_survival)
    call_survival <- call(call_survival) 




######### RECUP SPECIFICATION ##########################
    nT <- length(pred.times)
    rl_cond <- rep(0,nT)
    epoir <- rep(0,nT)
    vopt <- as.double(model$V)
    best <- as.double(model$best)
    NPM <- length(best)
    ns_vect <- rep(0,nT)
    nevt_vect <- rep(0,nT)

### specification recup de args
    idprob0 <- model$idprob
    idprob0[1] <- 0
    idea0 <- model$idea
    idg0 <- model$idg
    idcor0 <- model$idcor
    idiag0 <- model$idiag
    nv0 <- length(idprob0)
    nwg0 <- model$N[6]
    ng0 <- model$ng
    ncor0 <- model$N[7]
    nz0 <- model$hazard[[4]]
    zi0 <- as.vector(model$hazard[[3]])
    typrisq0 <-  model$hazard[[1]]
    risqcom0 <- switch(model$hazard[[2]],"Specific"=0,"PH"=2,"Common"=1)
    logspecif <- model$logspecif
    idcom <- model$idcom
    idspecif <- model$idspecif
    idtdv <- model$idtdv

                                        #cat(paste("  risqcom ",risqcom0)," \n")
                                        #cat(paste("  idprob0 ",idprob0)," \n")
                                        #cat(paste("  idea0 ",idea0)," \n")
                                        #cat(paste("  idg0 ",idg0)," \n")
                                        #cat(paste("  idpxevt ",idxevt)," \n")

##############   RECUP DATA   ##########################


### new data or data from estimation

    if(missing(newdata))
    {
        if(!is.null(model$data))
        {
            data <- model$data
        }
        else
        {
            data <- eval(model$call$data)
        }
        if(!(var.time %in% names(data))) stop("The Jointlcmm data must contain the var.time variable")
            new.data <- FALSE
            ##  var.time est dedans
        }
    else
        {
            data <- newdata
        }

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


    
### pour les facteurs
    ## transform to factor is the variable appears in levels$levelsdata
    for(v in colnames(data))
    {
        if(v %in% names(model$levels$levelsdata))
        {
            if(!is.null(model$levels$levelsdata[[v]]))
            {
                data[,v] <- factor(data[,v], levels=model$levels$levelsdata[[v]])
            }
        }
    }
    
    ## ##cas ou on a factor() dans l'appel
    ## z <- all.names(call_fixed)
    ## ind_factor <- which(z=="factor")
    ## if(length(ind_factor))
    ##     {
    ##         nom.factor <- z[ind_factor+1]  
    ##         for (v in nom.factor)
    ##             {
    ##                 mod <- levels(as.factor(olddata[,v]))
    ##                 if (!all(levels(as.factor(data[,v])) %in% mod)) stop(paste("invalid level in factor", v))
    ##                 data[,v] <- factor(data[,v], levels=mod)
    ##             }
    ##     }
    call_fixed <- gsub("factor","",call_fixed)

    ## z <- all.names(call_random)
    ## ind_factor <- which(z=="factor")
    ## if(length(ind_factor))
    ##     {
    ##         nom.factor <- z[ind_factor+1]
    ##         for (v in nom.factor)
    ##             {
    ##                 mod <- levels(as.factor(olddata[,v]))
    ##                 if (!all(levels(as.factor(data[,v])) %in% mod)) stop(paste("invalid level in factor", v))
    ##                 data[,v] <- factor(data[,v], levels=mod)
    ##             }
    ##     }
    call_random <- gsub("factor","",call_random)
    
    ## z <- all.names(call_classmb)
    ## ind_factor <- which(z=="factor")
    ## if(length(ind_factor))
    ##     {
    ##         nom.factor <- z[ind_factor+1]
    ##         for (v in nom.factor)
    ##             {
    ##                 mod <- levels(as.factor(olddata[,v]))
    ##                 if (!all(levels(as.factor(data[,v])) %in% mod)) stop(paste("invalid level in factor", v))
    ##                 data[,v] <- factor(data[,v], levels=mod)
    ##             }
    ##     }
    call_classmb <- gsub("factor","",call_classmb)
    
    ## z <- all.names(call_survival)
    ## ind_factor <- which(z=="factor")
    ## if(length(ind_factor))
    ##     {
    ##         nom.factor <- z[ind_factor+1]
    ##         for (v in nom.factor)
    ##             {
    ##                 mod <- levels(as.factor(olddata[,v]))
    ##                 if (!all(levels(as.factor(data[,v])) %in% mod)) stop(paste("invalid level in factor", v))
    ##                 data[,v] <- factor(data[,v], levels=mod)
    ##             }
    ##     }
    
    call_survival <- gsub("factor","",call_survival)


    

    
    ## objet Surv
    surv <- model$call$survival[[2]]
    
    if(length(surv)==3) #censure droite sans troncature gauche
        {
            idtrunc <- 0 
            
            Tevent <- getElement(object=data,name=as.character(surv[2]))
            Event <- getElement(object=data,name=as.character(surv[3]))  
            Tentry <- rep(0,length(Tevent)) #si pas de troncature, Tentry=0
            
            noms.surv <-  c(as.character(surv[2]),as.character(surv[3]))
            
            surv <- do.call("Surv",list(time=Tevent,event=Event,type="mstate")) 
        }
    
    if(length(surv)==4) #censure droite et troncature
        {
            idtrunc <- 1 
            
            Tentry <- getElement(object=data,name=as.character(surv[2]))
            Tevent <- getElement(object=data,name=as.character(surv[3]))
            Event <- getElement(object=data,name=as.character(surv[4]))  
            
            noms.surv <-  c(as.character(surv[2]),as.character(surv[3]),as.character(surv[4]))   
            
            surv <- do.call("Surv",list(time=Tentry,time2=Tevent,event=Event,type="mstate"))   
        }   



    ##fun.time
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
                    vars <- intersect(colnames(data),all.names(bodyfun))
                    getvars <- paste("getElement(data,'",vars,"')",sep="")
                            
                    if(length(vars))
                        {
                            for(i in 1:length(vars))
                                {
                                    bodyfun <- sub(vars[i],getvars[i],bodyfun)
                                    if(!is.null(newdata)) bodyfun <- sub(as.character(match.call()$newdata),"",bodyfun)
                                    else bodyfun <- sub(as.character(model$call$data),"",bodyfun)
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
                    
                    Time <- do.call("ff3",list(data[,var.time]))
                }
            else
                {
                    Time <- data[,var.time]
                } 
        }
    else
        {
            Time <- data[,var.time]
        }


    
###subset de data avec les variables utilisees

            newdata1 <- data[,unique(c(model$Names$ID,model$Names$Yname,noms.surv,model$Names$Xnames2,model$Names$prior.name,model$Names$TimeDepVar.name,var.time)),drop=FALSE]


    ## remplacer les NA de prior par 0  
    if(!is.null(model$Names$prior.name))
        {
            prior <- newdata[,model$Names$prior.name]
            newdata1[which(is.na(prior)),model$Names$prior.name] <- 0
        }

    ## remplacer les NA de TimeDepVar par Tevent
    Tint <- Tevent
    nvdepsurv <- 0

    if(!is.null(model$Names$TimeDepVar.name))
        {
            Tint <- newdata1[,model$Names$TimeDepVar.name]
            Tint[(is.na(Tint))] <- Tevent[(is.na(Tint))]
            Tint[Tint>Tevent] <- Tevent[Tint>Tevent]
            Tint[Tint<Tentry] <- Tentry[Tint<Tentry]
            nvdepsurv <- 1
            if (length(Tint[Tint>Tentry])==0)
                {
                    Tint <- Tevent
                    cat("TimeDepVar will be ignored since it is always lower than Time of Entry (0 by default). \n")
                    nvdepsurv  <- 0
                }
            if (length(Tint[Tint<Tevent])==0)
                {
                    cat("TimeDepVar will be ignored since it is always greater than Time of Event. \n")
                    nvdepsurv <- 0
                }

            newdata1[,model$Names$TimeDepVar.name] <- Tint 
        }


    ##enlever les NA
    linesNA <- apply(newdata1,2,function(v) which(is.na(v)))
    linesNA <- unique(unlist(linesNA))  
    na.fun <- which(is.na(Time))
    if(length(na.fun)) linesNA <- unique(c(linesNA,na.fun))
    
    if(length(linesNA))
        {
            if(na.action==1) newdata1 <- newdata1[-linesNA,,drop=FALSE] 
            if(na.action==2) stop("Data contain missing values.")
            Tentry <- Tentry[-linesNA]  
            Tevent <- Tevent[-linesNA] 
            Event <- Event[-linesNA]
            Tint <- Tint[-linesNA]
            Time <- Time[-linesNA]
        }


    ## create one data frame for each formula (useful with factors)
    newdata1fixed <- newdata1
    for(v in colnames(newdata1fixed))
    {
        if(v %in% names(model$levels$levelsfixed))
        {
            if(!is.null(model$levels$levelsfixed[[v]]))
            {
                newdata1fixed[,v] <- factor(newdata1fixed[,v], levels=model$levels$levelsfixed[[v]])
                if(any(is.na(newdata1fixed[,v]))) stop(paste("Wrong factor level in variable",v))
            }
        }
    }
    newdata1random <- newdata1
    for(v in colnames(newdata1random))
    {
        if(v %in% names(model$levels$levelsrandom))
        {
            if(!is.null(model$levels$levelsrandom[[v]]))
            {
                newdata1random[,v] <- factor(newdata1random[,v], levels=model$levels$levelsrandom[[v]])
                if(any(is.na(newdata1random[,v]))) stop(paste("Wrong factor level in variable",v))
            }
        }
    }
    newdata1classmb <- newdata1
    for(v in colnames(newdata1classmb))
    {
        if(v %in% names(model$levels$levelsclassmb))
        {
            if(!is.null(model$levels$levelsclassmb[[v]]))
            {
                newdata1classmb[,v] <- factor(newdata1classmb[,v], levels=model$levels$levelsclassmb[[v]])
                if(any(is.na(newdata1classmb[,v]))) stop(paste("Wrong factor level in variable",v))
            }
        }
    }
    newdata1surv <- newdata1
    for(v in colnames(newdata1surv))
    {
        if(v %in% names(model$levels$levelssurv))
        {
            if(!is.null(model$levels$levelssurv[[v]]))
            {
                newdata1surv[,v] <- factor(newdata1surv[,v], levels=model$levels$levelssurv[[v]])
                if(any(is.na(newdata1surv[,v]))) stop(paste("Wrong factor level in variable",v))
            }
        }
    }
    
    
###creation de X0 (ttes les var + interactions)

    X_intercept <- model.matrix(~1,data=newdata1)
    colnames(X_intercept) <- "intercept"
    
    ## fixed

    X_fixed <- model.matrix(formula(paste("~",call_fixed,sep="")),data=newdata1fixed)
    if(colnames(X_fixed)[1]=="(Intercept)")
        {
            X_fixed <- X_fixed[,-1,drop=FALSE]        
        }

    ## random
    if(!is.null(model$call$random))
        {
            X_random <- model.matrix(formula(paste("~",call_random,sep="")),data=newdata1random)
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
            X_classmb <- model.matrix(formula(paste("~",call_classmb,sep="")),data=newdata1classmb)
            if(colnames(X_classmb)[1]=="(Intercept)")
                {
                    colnames(X_classmb)[1] <- "intercept"
                }
        }
    else
        {
            X_classmb <- NULL
        }

    ## survival
    if(!is.null(model$call$survival))
        {
            X_survival <- model.matrix(formula(paste("~",call_survival,sep="")),data=newdata1surv)
            if(colnames(X_survival)[1]=="(Intercept)")
                {
                    colnames(X_survival)[1] <- "intercept"
                }
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
    
    
    ## Construction de X0 dans le bon ordre
    X0 <- cbind(X_intercept,X_fixed)
    colX <- c("intercept",colnames(X_fixed))
    if(!is.null(X_random))
    {
        for(i in 1:length(colnames(X_random)))
        {
            if((colnames(X_random)[i] %in% colnames(X0))==FALSE)
            {
                X0 <- cbind(X0,X_random[,i])
                colnames(X0) <- c(colX,colnames(X_random)[i])
                colX <- colnames(X0)
            }	 
        }
    }
    if(!is.null(X_classmb))
    {
        for(i in 1:length(colnames(X_classmb)))
        {
            if((colnames(X_classmb)[i] %in% colnames(X0))==FALSE)
            {
                X0 <- cbind(X0,X_classmb[,i])
                colnames(X0) <- c(colX,colnames(X_classmb)[i])
                colX <- colnames(X0)
            }	
        }
    }
    if(!is.null(X_survival))
    {
        for(i in 1:length(colnames(X_survival)))
        {
            if((colnames(X_survival)[i] %in% colnames(X0))==FALSE)
            {
                X0 <- cbind(X0,X_survival[,i])
                colnames(X0) <- c(colX,colnames(X_survival)[i])
                colX <- colnames(X0)
            }
        }
    }
    
    if(model$N[7]>0)
    {
        idspecif <- matrix(model$idspecif,nbevt,length(model$idg),byrow=TRUE)
        if(model$idea[z]==0 & model$idprob[z]==0 & model$idg[z]==0 & model$idcom[z]==0 & all(idspecif[,z]==0))
        {
            X0 <- cbind(X0,X_cor)
            colnames(X0) <- c(colX,model$Names$Xnames[z])
            colX <- colnames(X0)
        }
        
    }
    
    ## X <- cbind(X_intercept,X_fixed,X_random,X_classmb,X_survival,X_cor)
    ## colX <- strsplit(colnames(X),split=":",fixed=TRUE)
    ## colX <- lapply(colX,sort)
    ## colX <- lapply(colX,paste,collapse=":")
    ## colnames(X) <- unlist(colX)
    
    ## Xnames <- strsplit(model$Names$Xnames,split=":",fixed=TRUE)
    ## Xnames <- lapply(Xnames,sort)
    ## Xnames <- lapply(Xnames,paste,collapse=":")
    ## Xnames <- unlist(Xnames)
    
    ## X0 <- X[,Xnames]
###X0 fini  


    ## DEP VAR
    Y0 <- newdata1[,model$Names$Yname]
    nobs0 <- length(Y0)


    if ((max(Tevent)>max(zi0))&(nz0!=2)) stop("The maximal time of event in the new dataset should not be greater than the maximal time of event in the dataset used in Jointlcmm.")


    ##identifiant sujets
    IND <- newdata1[,model$Names$ID]
    #IDnum <- as.numeric(IND)




    ## INCLUSION PRIOR 
    if(is.null(model$Names$prior.name))
        { 
            PRIOR <- rep(0,length(IND))
        } 
    else
        {
            PRIOR <- newdata1[,model$Names$prior.name]
            PRIOR[(is.na(PRIOR))] <- 0
        }


### DATA SORTING on IND variable
    matYX <- cbind(IND,PRIOR,Tentry,Tevent,Event,Tint,Y0,Time,X0)
    matYXord <- matYX[order(IND),]
    Y0 <- as.numeric(matYXord[,7])
    Time <- as.numeric(matYXord[,9])
    X0 <- apply(matYXord[,-c(1,2,3,4,5,6,7,8),drop=FALSE],2,as.numeric)
    #IDnum <- matYXord[,1]
    IND <- matYXord[,1]
    PRIOR <- as.numeric(matYXord[,2])
    PRIOR <- as.integer(as.vector(PRIOR))

### Tevent, Tentry et Event de dim ns  
    nmes0 <- as.vector(table(IND))
    data.surv <- data.frame(IND,apply(matYXord[,c(3,4,5,6)],2,as.numeric))
    data.surv <- data.surv[cumsum(nmes0),]
    
    #data.surv <- unique(data.surv)
    #if(nrow(data.surv) != length(unique(IND))) stop("Subjects cannot have several times to event.")

    tsurv0 <- data.surv[,2] 
    tsurv <- data.surv[,3]
    devt <- data.surv[,4]
    tsurvint <- data.surv[,5]
    ind_survint <- (tsurvint<tsurv) + 0 


    ns0 <- length(nmes0)

    INDuniq <- unique(IND)

    prior0 <- unique(data.frame(IND,PRIOR))[,2]


    contribt <- rep(0,length=ns0*nT)

    fix0 <- rep(0,NPM)
    posfix <- eval(model$call$posfix)
    if(length(posfix)) fix0[posfix] <- 1
    npmssfix <- NPM-length(posfix)

################ FORTRAN FUNCTION CALL #####################
    
    ptm<-proc.time()
    cat("Be patient, epoce function is running ... \n")


    out <- .Fortran(C_cvpl,
                    as.double(Y0),
                    as.double(X0),
                    as.integer(prior0),
                    as.integer(idprob0),
                    as.integer(idea0),
                    as.integer(idg0),
                    as.integer(idcor0),
                    as.integer(idcom),
                    as.integer(idspecif),
                    as.integer(idtdv),
                    as.integer(ns0),
                    as.integer(ng0),
                    as.integer(ncor0),
                    as.integer(nv0),
                    as.integer(nobs0),
                    as.integer(nmes0),
                    as.integer(idiag0),
                    as.integer(nwg0),
                    as.integer(NPM),
                    as.double(Time),
                    as.integer(typrisq0),
                    as.integer(idtrunc),
                    as.integer(risqcom0),
                    as.integer(nz0),
                    as.double(zi0),
                    as.double(tsurv0),
                    as.double(tsurv),
                    as.integer(devt),
                    as.integer(ind_survint),
                    as.double(vopt),
                    as.integer(nT),
                    as.double(pred.times),
                    best=as.double(best),
                    epoir=as.double(epoir),
                    rl_cond=as.double(rl_cond),
                    ns_vect=as.integer(ns_vect),
                    nevt_vect=as.integer(nevt_vect),
                    contribt=as.double(contribt),
                    as.integer(logspecif),
                    as.integer(fix0),
                    as.integer(npmssfix))


    ## construction de la matrice contribt
    contrib <- matrix(out$contribt,nrow=ns0,ncol=nT)
    namesContrib <- as.vector(apply(matrix(pred.times,nrow=1),MARGIN=2,FUN=function(x){paste("IndivContrib_time_",x,sep="")}))
    colnames(contrib) <- namesContrib
    contrib <- data.frame(INDuniq,contrib)
    contrib[contrib==0] <- NA
    ## remplacer les 0 (vrais 0 par NA) (le nombre de non nul pour un temps = ns_vect(de ce temps)
    ## colnames : permier colonne = le nom de la variable IND = colnames(model$pprob)[1]
    ##          : colonnes suivantes = IndivContrib_time_x, x etant la valeur du temps de prediction (ce qu'il y a dans pred.times)

    if (!is.null(newdata))
        {
            out$epoir <- rep(NA,length(pred.times))
        }
    out$epoir[out$epoir==1.e9] <- NA
    out$rl[out$rl==-1.e9] <- NA


    cvpl <- cbind(pred.times,out$ns_vect,out$nevt_vect,-out$rl,out$epoir)
    colnames(cvpl) <- c("pred. times"," N at risk"," N events","MPOL","CVPOL")
    rownames(cvpl) <- rep(" ",length(cvpl[,1]))

    ## sortie des resultats
    res <- list(call.Jointlcmm=model$call,call.epoce=cl,EPOCE=cvpl,IndivContrib=contrib,new.data=new.data)

    class(res) <-c("epoce")
    cost<-proc.time()-ptm
    cat("The program took", round(cost[3],2), "seconds \n")
    res
}



















