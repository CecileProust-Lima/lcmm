#' Estimation of secondary regression models after the estimation of a primary latent class model
#' 
#' This function fits regression models to relate a latent class structure (stemmed 
#' from a latent class model estimated within \code{lcmm} package) with either an external
#'  outcome or external class predictors. 
#'  Two inference techniques are implemented to account for the classification error: 
#'  
#'  - a 2-stage estimation of the joint likelihood of the primary latent class model 
#'  and the secondary/ external regression;
#'  
#'  - a regression between the posterior latent class assignment and the external variable 
#'  which internally corrects for the assignment misclassification. 
#'  
#' It returns an object from one of the \code{lcmm} package classes.
#' 
#' A. DATA STRUCTURE
#' 
#' The \code{data} argument must follow specific structure for individual variables,
#' i.e. variables with a unique constant value for each subject. For an individual variable
#' given as external outcome, data value must be present only once per subject,
#' independently of any time variable used in the primary latent class.
#' For an individual variable given as external class predictor,
#' data values must be given for every row of every individual (as usual)
#' 
#' B. VARIANCE ESTIMATION
#' 
#' Not taking into account first stage variance with specifing \code{"none"} may lead to
#' underestimation of the final variance. When possible, Method \code{"Hessian"} 
#' which relies on the combination of Hessians from the primary and secondary
#' model is recommended. However, it may become numerically intensive in the event 
#' of very high number of parameters in the primary latent class model. As an 
#' alternative, especially in situations with a complex primary model but rather 
#' parcimonious secondary model, method \code{"paramBoot"} which implements a 
#' parametric bootstrap can be used.
#' 
#' @param model an object inheriting from class \code{lcmm}, 
#' \code{Jointlcmm}, \code{multlcmm} or \code{mpjlcmm} giving the primary latent
#'  class model.
#' @param fixed optional two sided linear formula object for specifying the
#' fixed-effects in the secondary model with an external outcome variable.
#' The response outcome is on the left of \code{~} and the covariates are separated
#' by \code{+} on the right of the \code{~}. By default, an intercept is included.
#' @param mixture optional one-sided formula object for the class-specific fixed effects
#' in the model for the external outcome. Among the list of covariates included in fixed,
#' the covariates with class-specific regression parameters are entered in
#' mixture separated by \code{+}. By default, an intercept is included.
#' If no intercept, \code{-1} should be the first term included.
#' @param random optional one-sided linear formula object for specifying the
#' random-effects on external outcome in the secondary model, if appropriate. 
#' By default, no random effect is included.
#' @param subject name of the covariate representing the grouping structure.
#' Even in the absence of a hierarchical structure. By default, the function will try to retrieve it
#' from \code{model} argument.
#' @param classmb optional one-sided formula specifying the external predictors of 
#' latent class membership to be modelled in the secondary class-membership multinomial 
#' logistic model. Covariates are separated by \code{+} on the right of the \code{~}. 
#' By default, an intercept is included. 
#' @param varest optional character indicating the method to be used to compute the variance of 
#' the regression estimates. \code{"none"} does not account for the uncertainty in the primary latent 
#' class model, \code{"paramBoot"} computes the total variance using a parametric bootstrap technique, 
#' \code{"Hessian"} computes the total Hessian of the joint likelihood (implemented for \code{"twoStageJoint"}
#' method only). Default to \code{"Hessian"} for \code{"twoStageJoint"} method
#' @param M option integer indicating the number of draws for the parametric boostrap when \code{varest="paramBoot"}.
#' Default to 200.
#' @param idiag if appropriate, optional logical for the structure of the variance-covariance
#' matrix of the random-effects in the secondary model. 
#' If \code{FALSE}, a non structured matrix of
#' variance-covariance is considered (by default). If \code{TRUE} a diagonal
#' matrix of variance-covariance is considered.
#' @param nwg if appropriate, optional logical indicating if the variance-covariance of the
#' random-effects in the secondary model is class-specific. If \code{FALSE} the variance-covariance
#' matrix is common over latent classes (by default). If \code{TRUE} a
#' class-specific proportional parameter multiplies the variance-covariance
#' matrix in each class (the proportional parameter in the last latent class
#' equals 1 to ensure identifiability).
#' @param link optional family of parameterized link functions for the external outcome if appropriate.
#' Defaults to NULL, corresponding to continuous Gaussian distribution (hlme function).
#' @param intnodes optional vector of interior nodes. This argument is only
#' required for a I-splines link function with nodes entered manually.
#' @param epsY optional definite positive real used to rescale the marker in (0,1)
#' when the beta link function is used. By default, epsY=0.5.
#' @param data Data frame containing the variables named in
#' \code{fixed}, \code{mixture}, \code{random}, \code{classmb} and \code{subject},
#' for both the current function arguments and the primary model arguments
#' Check \code{details} to get information on the data structure, especially with
#' external outcomes.
#' @param method character indicating the inference technique to be used:
#' only "twoStageJoint" for 2-stage estimation from the joint likelihood 
#' is implemented at the moment.
#' @param B optional vector of initial parameter values for the secondary model. 
#' If external outcome, the vector has the same structure as a latent class model
#' estimated in the other functions of \code{lcmm} package for the same type of 
#' outcome. If external class predictors (of size p), the vector is of length 
#' (ng-1)*(1+p). If \code{B=NULL} (by default), internal initial values are selected. 
#' @param convB optional threshold for the convergence criterion based on the
#' parameter stability. By default, convB=0.0001.
#' @param convL optional threshold for the convergence criterion based on the
#' log-likelihood stability. By default, convL=0.0001.
#' @param convG optional threshold for the convergence criterion based on the
#' derivatives. By default, convG=0.0001.
#' @param longitudinal only with \code{mpjlcmm} primary models and "twoStageJoint" method:
#' mandatory list containing the longitudinal submodels used in the primary latent class model. 
#' @param maxiter optional maximum number of iterations for the secondary model estimation using 
#' Marquardt iterative algorithm. Defaults to 100
#' @param posfix optional vector specifying indices in parameter vector B the 
#' secondary model that should not be estimated. Default to NULL, all the 
#' parameters of the secondary regression are estimated. 
#' @param partialH optional logical for Piecewise and Splines baseline risk functions and
#' Splines link functions only. Indicates whether the parameters of the baseline risk or
#' link functions can be dropped from the Hessian matrix to define convergence criteria
#' (can solve non convergence due to estimates at the boundary of the parameter space - usually 0).
#' @param verbose logical indicating whether information about computation should be
#' reported. Default to FALSE.
#' @param nproc the number cores for parallel computation. Default to 1 (sequential mode).
#' 
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' 
#' ###### Estimation of the primary latent class model                   ######
#' 
#' set.seed(1234)
#' PrimMod <- hlme(Ydep1~Time,random=~Time,subject='ID',ng=1,data=data_lcmm)
#' PrimMod2 <- hlme(Ydep1~Time,mixture=~Time,random=~Time,subject='ID',
#'                  ng=2,data=data_lcmm,B=random(PrimMod))
#' 
#' ###### Example 1: Relationship between a latent class structure and         #
#' #                   external class predictors                          ######
#'                   
#' # estimation of the secondary multinomial logistic model with total variance
#' # computed with the Hessian
#' 
#' XextHess <- externVar(PrimMod2,
#'                       classmb = ~X1 + X2 + X3 + X4, 
#'                       subject = "ID",
#'                       data = data_lcmm,
#'                       method = "twoStageJoint") 
#' summary(XextHess)
#' 
#' # estimation of a secondary multinomial logistic model with total variance
#' # computed with parametric Bootstrap (much longer). When using the bootstrap 
#' # estimator, we recommend running first the analysis with option varest = "none" 
#' # which is faster but which underestimates the variance. And then use these values
#' # as initial values when running the model with varest = "paramBoot" to obtain 
#' # a valid variance of the parameters. 
#' 
#' XextNone <- externVar(PrimMod2,
#'                       classmb = ~X1 + X2 + X3 + X4, 
#'                       subject = "ID",
#'                       data = data_lcmm,
#'                       varest = "none",
#'                       method = "twoStageJoint") 
#' 
#' XextBoot <- externVar(PrimMod2,
#'                       classmb = ~X1 + X2 + X3 + X4, 
#'                       subject = "ID",
#'                       data = data_lcmm,
#'                       varest = "paramBoot",
#'                       method = "twoStageJoint",
#'                       B = XextNone$best) 
#' summary(XextBoot)
#' 
#'  
#' ###### Example 2: Relationship between a latent class structure and         #
#' #                external outcome (repeatedly measured over time)     ######
#'                   
#' # estimation of the secondary linear mixed model with total variance
#' # computed with the Hessian
#' 
#' YextHess = externVar(PrimMod2,   #primary model
#'                      fixed = Ydep2 ~ Time*X1,  #secondary model
#'                      random = ~Time, #secondary model
#'                      mixture = ~Time,  #secondary model
#'                      subject="ID",
#'                      data=data_lcmm,
#'                      method = "twoStageJoint")
#'                      
#' 
#' # estimation of a secondary linear mixed model with total variance
#' # computed with parametric Bootstrap (much longer). When using the bootstrap 
#' # estimator, we recommend running first the analysis with option varest = "none" 
#' # which is faster but which underestimates the variance. And then use these values
#' # as initial values when running the model with varest = "paramBoot" to obtain 
#' # a valid variance of the parameters. 
#' 
#' YextNone = externVar(PrimMod2,   #primary model
#'                      fixed = Ydep2 ~ Time*X1,  #secondary model
#'                      random = ~Time, #secondary model
#'                      mixture = ~Time,  #secondary model
#'                      subject="ID",
#'                      data=data_lcmm,
#'                      varest = "none",
#'                      method = "twoStageJoint")
#' 
#' YextBoot = externVar(PrimMod2,   #primary model
#'                      fixed = Ydep2 ~ Time*X1,  #secondary model
#'                      random = ~Time, #secondary model
#'                      mixture = ~Time,  #secondary model
#'                      subject="ID",
#'                      data=data_lcmm,
#'                      method = "twoStageJoint",
#'                      B = YextNone$best,
#'                      varest= "paramBoot")
#' 
#' summary(YextBoot) 
#' 
#' }
#'
#'
#'
#' 
#' @export
#' 
#' 
#' 


externVar = function(model,
                     fixed,
                     mixture,
                     random,
                     subject,
                     classmb,
                     survival,
                     hazard = "Weibull",
                     hazardtype = "Specific",
                     hazardnodes = NULL,
                     TimeDepVar = NULL,
                     nsim = 100,
                     logscale = FALSE,
                     idiag = FALSE,
                     nwg = FALSE,
                     link = NULL,
                     intnodes=NULL,
                     epsY = 0.5,
                     data,
                     longitudinal,
                     method,
                     varest = "Hessian",
                     M = 200,
                     B,
                     convB = 0.0001,
                     convL = 0.0001,
                     convG = 0.0001,
                     maxiter = 100,
                     posfix,
                     partialH = FALSE,
                     verbose = FALSE,
                     nproc = 1){
  
  ptm <- proc.time()
  
  if(missing(model)) stop("model argument must be given")
  if(!inherits(model, c("hlme", "lcmm", "multlcmm", "Jointlcmm", "mpjlcmm"))) stop('primary model class must be either "hlme", "lcmm", "multlcmm", "Jointlcmm" or "mpjlcmm"')
  if(model$conv == 2) warning("primary model did not fully converge")
  if(sum(c(!missing(fixed), !missing(classmb), !missing(survival))) != 1) stop("One and only one in survival, fixed or classmb must be given")
  if(missing(method) | !method %in% c("twoStageJoint", "conditional")) stop('Method must be either "twoStageJoint" or "conditional"')
  if(model$ng == 1) stop("Primary model does not have latent class structure (ng=1)")
  if(!varest %in% c("none", "paramBoot", "Hessian")) stop('Variance estimation method "varest" must be either "none", "paramBoot" or "Hessian"')
  if(!is.null(link) & missing(fixed)) stop("The argument link is not to be used with external class predictor")

  
  if(missing(posfix)) posfix = c()
  
  cl = match.call()
  
  #Informations about primary model
  argumentsIn = as.list(model$call)
  funIn = as.character(argumentsIn[[1]])
  argumentsIn[[1]] = NULL
  ng = model$ng
  nIn = length(model$best)
  
  if(is.null(argumentsIn[["classmb"]])){
    oldclassmb = ~ 1
  } else {
    oldclassmb = formula(argumentsIn[["classmb"]])
  }
  if(!missing(classmb) & oldclassmb != ~1) stop("Primary model already has class predictor")
  
  #number of MB parameters in primary model
  nInMB = ncol(model.matrix(oldclassmb, data))*(ng-1)
  
  #Get subject
  if(missing(subject)){
    if (model$call$subject %in% colnames(data)){
      subject = model$call$subject
    } else {
      stop("The argument subject must be specified if different from the subject argument used in the primary model")
    }
  } 
  
  #Get longitudinal
  if(funIn == "mpjlcmm"){
    if(missing(longitudinal)) stop("The argument longitudinal is mandatory with a mpjlcmm primary model")
    
    longCall = substitute(longitudinal)
    
    K = length(longitudinal)
    for(k in 1:K){
      cl_long = as.list(longitudinal[[k]]$call)
      cl_long[["computeDiscrete"]] = FALSE
      longitudinal[[k]]$call = as.call(cl_long)
      
      assign(as.character(longCall[[k+1]]), longitudinal[[k]])
    }
  }
  
  #nVCIn
  if(funIn == "mpjlcmm"){
    nVCIn = model$N[6]
    iVCIn = sum(model$N[1:5]) + 1:nVCIn
  } else if(funIn == "Jointlcmm"){
    nVCIn = model$N[5]
    iVCIn = sum(model$N[1:4]) + 1:nVCIn
  } else if(funIn == "multlcmm"){
    nVCIn = model$N[4]
    iVCIn = sum(model$N[3]) + 1:nVCIn
  } else {
    nVCIn = model$N[3]
    iVCIn = sum(model$N[1:2]) + 1:nVCIn
  }
  
  #pprob with new data
  argumentsInEdit = argumentsIn
  argumentsInEdit[["data"]] = data
  argumentsInEdit[["maxiter"]] = 0
  argumentsInEdit[["B"]] = model$best
  argumentsInEdit[["verbose"]] = FALSE
  modelEdit = do.call(funIn, argumentsInEdit)
  pprob = modelEdit$pprob
  
  
  arguments = list()
  
  #Finding out the number of parameters is needed for all survival external outcome
  if(!missing(survival)){
    #manage inputs
    if(!is.null(argumentsIn[["survival"]])) stop('secondary survival model is not supported with "twoStageJoint" method if primary model already includes survival')
    
    #Informations about secondary outcome model
    ##number of survival parameters to estimate
    ###number of events
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
    
    Tevent <- getElement(object=data,name=nom.Tevent)
    Event <- getElement(object=data,name=nom.Event)  
    nbevt <- length(attr(do.call("Surv",list(time=Tevent,event=Event,type="mstate")),"states")) 
    if(nbevt<1) nbevt <- 1
    
    ###get number of parameters for baseline functions
    hazard <- rep(hazard,length.out=nbevt)
    hazardtype <- rep(hazardtype,length.out=nbevt)
    if(any(hazard %in% c("splines","Splines")))
    {
      hazard[which(hazard %in% c("splines","Splines"))] <- "5-quant-splines" 
    }
    if(any(hazard %in% c("piecewise","Piecewise")))
    {
      hazard[which(hazard %in% c("piecewise","Piecewise"))] <- "5-quant-piecewise" 
    }
    
    hazWhat = hazard
    #when not weibull, keep only the last word of hazard
    hazWhat[hazard != 'Weibull'] = sapply(strsplit(hazard[hazard != 'Weibull'], "-"), `[`, 3)
    hazN = rep(2, nbevt)
    hazN[hazard != 'Weibull'] = sapply(strsplit(hazard[hazard != 'Weibull'], "-"), `[`, 1)
    hazN = as.integer(hazN)
    hazN = hazN +
      (hazWhat == "Weibull")*0 +
      (hazWhat == "piecewise")*(-1) +
      (hazWhat == "splines")*(2)
    hazN = hazN + hazN *
      (hazardtype == "Specific")*(ng-1)
    
    #we extract the number of base function parameter with constraints (ie, not PH parameters)
    nSurvConstraint = hazN
    
    hazN = hazN  +
      (hazardtype == "PH")*(ng-1)
    
    #we now have all base functions parameters :
    nEstY = sum(hazN)
    
    ###get number of parameters for survival covariates
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
    
    #I extract number of variable in each formula (through model.matrix) excluding intercept
    ncols = sapply(c(form.commun, form.cause, form.mixture, form.causek), function(x, data){
      mm = model.matrix(x, data)
      mm = mm[,-1]
      if(is.null(ncol(mm))) {return(1)}
      else {return(ncol(mm))}
    }, data = data)
    
    #I also need how many parameters each kind of covariate makes
    nparam = c(1, nbevt, ng, nbevt*ng)
    
    #we now have all the survival covariates parameters :
    nEstX = sum(ncols*nparam)
    
    nEst = nEstY + nEstX
  }
  
  #A model structure is needed for all longitudinal external outcome
  if(!missing(fixed)){
    #Manage inputs
    
    if(missing(mixture)) mixture = ~1
    if(missing(random)) random = ~-1
    
    if(!inherits(fixed,"formula")) stop("The argument fixed must be a formula")
    if(!inherits(mixture,"formula")) stop("The argument mixture must be a formula")
    if(!inherits(random,"formula")) stop("The argument random must be a formula")
    
    if(length(fixed[[2]]) != 1){
      argfunctionStrMod = "multlcmm"
    } else if(is.null(link)){
      argfunctionStrMod = "hlme"
    } else {
      argfunctionStrMod = "lcmm"
    }
    
    #let's create structure for secondary model
    argumentsStrMod = list()
    argumentsStrMod[["fixed"]] = fixed
    argumentsStrMod[["random"]] = random
    argumentsStrMod[["subject"]] = subject
    argumentsStrMod[["ng"]] = 1
    argumentsStrMod[["idiag"]] = idiag
    if(argfunctionStrMod %in% c("lcmm", "multlcmm", "Jointlcmm")){
      argumentsStrMod[["link"]] = link
      argumentsStrMod[["intnodes"]] = intnodes
      argumentsStrMod[["epsY"]] = epsY
    }
    argumentsStrMod[["data"]] = data
    argumentsStrMod[["maxiter"]] = 0
    argumentsStrMod[["verbose"]] = FALSE
    
    strMod = do.call(argfunctionStrMod, c(argumentsStrMod))
    
    argumentsStrMod[["mixture"]] = mixture
    argumentsStrMod[["classmb"]] = ~1
    argumentsStrMod[["ng"]] = ng
    argumentsStrMod[["nwg"]] = nwg
    argumentsStrMod[["B"]] = as.name("strMod")
    
    strMod = do.call(argfunctionStrMod, c(argumentsStrMod))
  }
  
  #Change general model arguments
  if(method == "twoStageJoint"){
    
    #Common argument in twoStageJoint
    arguments[["data"]] = data
    arguments[["ng"]] = ng
    arguments[["subject"]] = subject
    #technical options
    arguments[["maxiter"]] = maxiter
    arguments[["verbose"]] = verbose
    arguments[["nproc"]] = nproc
    arguments[["convB"]] = convB
    arguments[["convL"]] = convL
    arguments[["convG"]] = convG
    arguments[["partialH"]] = partialH
    #primary survival
    arguments[["survival"]] =  argumentsIn[["survival"]]
    arguments[["hazard"]] =  argumentsIn[["hazard"]]
    arguments[["hazardtype"]] =  argumentsIn[["hazardtype"]]
    arguments[["hazardnodes"]] =  argumentsIn[["hazardnodes"]]
    arguments[["hazardrange"]] =  argumentsIn[["hazardrange"]]
    arguments[["TimeDepVar"]] =  argumentsIn[["TimeDepVar"]]
    arguments[["logscale"]] =  argumentsIn[["logscale"]]
    
    
    #Primary Jointlcmm needs transformation into lcmm to be put into longitudinal.
    if(funIn == "Jointlcmm"){
      
      if(is.null(argumentsIn[["link"]])){
        argfunJoint = "hlme"
      } else {
        argfunJoint = "lcmm"
      }
      
      argumentsJoint = argumentsIn
      argumentsJoint[["survival"]] = NULL
      argumentsJoint[["hazard"]] = NULL
      argumentsJoint[["hazardtype"]] = NULL
      argumentsJoint[["hazardnodes"]] = NULL
      argumentsJoint[["hazardrange"]] = NULL
      argumentsJoint[["TimeDepVar"]] = NULL
      argumentsJoint[["logscale"]] = NULL
      
      argumentsJoint[["maxiter"]] = 0
      argumentsJoint[["mixture"]] = NULL
      argumentsJoint[["classmb"]] = NULL
      argumentsJoint[["ng"]] = 1
      argumentsJoint[["nwg"]] = FALSE
      argumentsJoint[["B"]] = NULL
      argumentsJoint[["verbose"]] = FALSE
      argumentsJoint[["data"]] = data
      
      modNoSurv = do.call(argfunJoint, argumentsJoint)
      
      argumentsJoint[["mixture"]] = argumentsIn[["mixture"]]
      argumentsJoint[["classmb"]]= ~1
      argumentsJoint[["ng"]] = argumentsIn[["ng"]]
      argumentsJoint[["nwg"]] = argumentsIn[["nwg"]]
      argumentsJoint[["B"]] = modNoSurv
      
      modNoSurv = do.call(argfunJoint, argumentsJoint)
    }
    
    
    #Yextern survival
    if(!missing(survival)){
      funOut = "mpjlcmm"
      
      #nOut : nuber of total final parameters
      nOut = nIn + nEst
      
      #index
      iKeepOut = c(1:nInMB, nInMB+nEst+1:(nIn-nInMB))
      iKeepIn = 1:nIn
      iEst = nInMB+1:nEst
      
      #input VC, only for keep betaa (iKeepIn)
      iVCKeep = iVCIn
      
      #Id of varcov estimates
      iVCOut = c()
      
      #list of arguments
      arguments[["survival"]] =  survival
      arguments[["hazard"]] =  hazard
      arguments[["hazardtype"]] =  hazardtype
      arguments[["hazardnodes"]] =  hazardnodes
      arguments[["TimeDepVar"]] =  TimeDepVar
      arguments[["logscale"]] =  logscale
      arguments[["posfix"]] = unique(c(iKeepOut, posfix))
      #what is in longitudinal ?
      if(funIn == "mpjlcmm"){
        arguments[["longitudinal"]] = longitudinal
      } else {
        arguments[["longitudinal"]] = list(model)
      }
      #initial values
      if(missing(B)){
        arguments[["B"]][iEst] = rep(1, nEst)
      } else {
        arguments[["B"]][iEst] = B
      }
    }
    
    #Yextern longitudinal
    if(!missing(fixed)){
      
      funOut = "mpjlcmm"
      
      #Let's change strMod's saved call
      strMod$call$data = substitute(data)
      
      ## Now join the primary and secondary model
      
      #Informations about secondary outcome model
      #number of classmb parameters to remove
      nMB = ng-1
      #number of remaining parameters to estimate
      nStr = length(strMod$best)
      nEst = nStr - nMB
      #nOut : nuber of total final parameters
      nOut = nIn + nEst
      
      #index
      iKeepOut = 1:nIn
      iKeepIn = iKeepOut
      iEst = nIn+1:nEst
      
      #input VC, only for keep betaa (iKeepIn)
      iVCKeep = iVCIn
      
      #Id of varcov estimates
      if(inherits(strMod, "multlcmm")){
        nVCStr = strMod$N[4]
        iVCStr = sum(strMod$N[3]) + 1:nVCStr
      } else {
        nVCStr = strMod$N[3]
        iVCStr = sum(strMod$N[1:2]) + 1:nVCStr
      }
      iVCOut = nIn + iVCStr - nMB
      
      #Liste des arguments
      #on fixe nos parametres
      arguments[["posfix"]] = unique(c(iKeepOut, posfix))
      
      #On donne les modeles
      if(funIn == "mpjlcmm"){
        arguments[["longitudinal"]] = c(longitudinal, list(strMod))
      } else if(funIn == "Jointlcmm"){
        arguments[["longitudinal"]] = list(modNoSurv, strMod)
      } else {
        arguments[["longitudinal"]] = list(model, strMod)
      }
      
      #initial values
      if(missing(B)){
        arguments[["B"]][iEst] = strMod$best[(nMB+1):nStr]
      } else {
        arguments[["B"]][iEst] = B
      }
    }
    
    #X extern
    if(!missing(classmb)){
      
      if(!inherits(classmb,"formula")) stop("The argument classmb must be a formula")
      
      funOut = "mpjlcmm"
      
      #nEst : number of MB parameters in output model
      nEst1G = ncol(model.matrix(classmb, data))
      nEst = nEst1G*(ng-1)
      
      #nOut : nuber of total final parameters
      nOut = nIn - nInMB + nEst
      
      #index
      iKeepOut = (nEst+1):nOut
      iKeepIn = (nInMB+1):nIn
      iEst = 1:nEst
      
      #input VC, only for keep betas (iKeepIn)
      iVCKeep = iVCIn-nInMB
      
      #Id of varcov estimates (none)
      iVCOut = c()
      
      #On recree tous nos arguments
      if(funIn == "mpjlcmm"){
        arguments[["longitudinal"]] = longitudinal
      } else if(funIn == "Jointlcmm"){
        arguments[["longitudinal"]] = list(modNoSurv)
      } else {
        arguments[["longitudinal"]] = list(model)
      }
      arguments[["classmb"]] = classmb
      #on ajoute des valeurs de base pour nos nouveaux estimateurs
      arguments[["B"]] = rep(0, nOut)
      #initial values
      if(missing(B)){
        arguments[["B"]][1:ng-1] = model$best[1:ng-1]
      } else {
        arguments[["B"]][iEst] = B
      }
      
      #on fixe nos parametres
      arguments[["posfix"]] = unique(c(posfix, iKeepOut))
      
      # #### weird... "Input" model has been changed because we need it to build output model
      # #To explain : "hessienne" function needs mpjlcmm model. So i need mpj for paramBoot
      # #And we need to create our longitudinal (and logicall) now before changinng "input" moded.
      # if(funIn != "mpjlcmm" & varest == "Hessian"){
      #   
      #   argumentsMpj = list()
      #   if(funIn == "Jointlcmm"){
      #     argumentsMpj[["longitudinal"]] = as.call(list(as.name("list"), modNoSurv))
      #   } else {
      #     argumentsMpj[["longitudinal"]] = as.call(list(as.name("list"), substitute(model)))
      #   }
      #   argumentsMpj[["maxiter"]] = 0
      #   argumentsMpj[["B"]] = model$best
      #   argumentsMpj[["ng"]] = ng
      #   argumentsMpj[["subject"]] = subject
      #   argumentsMpj[["data"]] = argumentsIn[["data"]]
      #   argumentsMpj[["classmb"]] = oldclassmb
      #   
      #   argumentsMpj[["survival"]] = argumentsIn[["survival"]]
      #   argumentsMpj[["hazard"]] = argumentsIn[["hazard"]]
      #   argumentsMpj[["hazardtype"]] = argumentsIn[["hazardtype"]]
      #   argumentsMpj[["hazardnodes"]] = argumentsIn[["hazardnodes"]]
      #   argumentsMpj[["TimeDepVar"]] = argumentsIn[["TimeDepVar"]]
      #   argumentsMpj[["logscale"]] = argumentsIn[["logscale"]]
      #   
      #   modelMpj = do.call("mpjlcmm", argumentsMpj)
      #   modelMpj$V = model$V
      #   
      #   longitudinalOut = as.call(list(as.name("list"), substitute(model)))
      #   
      #   modelOld = model
      #   model = modelMpj
      #   model$V = modelOld$V
      # }
      
      
      
    }
  }
  
  if(method == "conditional"){
    
    #Yextern survival
    if(!missing(survival)){
      ## we need it all in a function in order to be able to use parametric bootstrap later on
      
      #remaining parameters to estimate
      iEst = 1:nEst
      
      iKeepIn = 1:nIn
      iKeepOut = 1:nIn+nEst
      
      nOut = nEst+2
      
      #Id of varcov estimates to keep in bootstrap
      iVCKeep = iVCIn
      iVCOut = c()
      
      conditional = function(model,
                               data,
                               survival,
                               hazard,
                               hazardtype,
                               hazardnodes = NULL,
                               TimeDepVar = NULL,
                               logscale,
                               subject,
                               ng,
                               B,
                               link,
                               iEst,
                               maxiter,
                               verbose,
                               nproc,
                               convB,
                               convL,
                               convG){
        #First : let's compute P(C|\tilde C) (\tilde C : A)
        pAlY = sapply(1:ng, function(g){
          return(as.numeric(predictClass(model, data)[,2] == g))
        })
        pClY = as.matrix(predictClass(model, data)[,3:(2+ng)])
        if(any(is.nan(pClY))) stop("NaN in posterior classification probability")
        
        pA = apply(pAlY, 2, mean)
        pClA = t(pAlY)%*%pClY/(model$ns*pA)
        
        if(det(pClA) == 0 | is.na(det(pClA))) stop("Computed error matrix is singular. One class might be empty")
        
        #Then : let's add it to dataset for each individual
        indivProb = pAlY%*%pClA
        indivProb = cbind(predictClass(model, data)[,1], indivProb)
        colnames(indivProb)[1] = subject
        data = merge(data, indivProb, by = subject)
        
        #We need dummy Y for the model to run
        data$dummyY = 1
        
        #Finally : we need to build the model for ppriors !
        arguments = list()
        
        arguments[["data"]] = data
        arguments[["fixed"]] = dummyY~1 
        arguments[["mixture"]] =  ~-1
        arguments[["survival"]] =  survival
        arguments[["hazard"]] =  hazard
        arguments[["hazardtype"]] =  hazardtype
        arguments[["hazardnodes"]] =  hazardnodes
        arguments[["TimeDepVar"]] =  TimeDepVar
        arguments[["logscale"]] =  logscale
        arguments[["subject"]] = subject
        arguments[["classmb"]] = ~-1
        arguments[["ng"]] = ng
        arguments[["pprior"]] = c("prob1", "prob2")
        arguments[["posfix"]] = 1:2+length(iEst)
        #technical options
        arguments[["maxiter"]] = maxiter
        arguments[["verbose"]] = verbose
        arguments[["nproc"]] = nproc
        arguments[["convB"]] = convB
        arguments[["convL"]] = convL
        arguments[["convG"]] = convG
        
        arguments[["B"]] = c(B[iEst], 1, 0.000001)
        
        res = do.call("Jointlcmm", arguments)
        res$call = match.call()
        return(res)
      }
      
      #we need to build the model
      arguments[["data"]] = data
      arguments[["survival"]] =  survival
      arguments[["hazard"]] =  hazard
      arguments[["hazardtype"]] =  hazardtype
      arguments[["hazardnodes"]] =  hazardnodes
      arguments[["TimeDepVar"]] =  TimeDepVar
      arguments[["logscale"]] =  logscale
      arguments[["model"]] = model
      arguments[["subject"]] = subject
      arguments[["ng"]] = ng
      arguments[["link"]] = link
      arguments[["iEst"]] = iEst
      #technical options
      arguments[["maxiter"]] = maxiter
      arguments[["verbose"]] = verbose
      arguments[["nproc"]] = nproc
      arguments[["convB"]] = convB
      arguments[["convL"]] = convL
      arguments[["convG"]] = convG
      
      #on ajoute des valeurs de base pour nos nouveaux estimateurs
      arguments[["B"]] = rep(1, nIn+nOut)
      #initial values
      if(missing(B)){
        B = rep(0.1, nEst)
      }
      arguments[["B"]][iEst] = B
      
      funOut = "conditional"
    }
    
    #Yextern longitudinal
    if(!missing(fixed)){
      ## we need it all in a function in order to be able to use parametric bootstrap later on
      
      #number of classmb parameters to remove
      nMB = ng-1
      #number of remaining parameters to estimate
      nStr = length(strMod$best)
      nEst = nStr - nMB
      iEst = 1:nEst
      
      iKeepIn = 1:nIn
      iKeepOut = 1:nIn+nEst
      
      nOut = nEst
      
      #Id of varcov estimates to keep in bootstrap
      iVCKeep = iVCIn
      #Id of varcov estimates
      if(inherits(strMod, "multlcmm")){
        nVCStr = strMod$N[4]
        iVCStr = sum(strMod$N[3]) + 1:nVCStr
      } else {
        nVCStr = strMod$N[3]
        iVCStr = sum(strMod$N[1:2]) + 1:nVCStr
      }
      iVCOut = nIn + iVCStr - nMB
      
      
      conditional = function(model,
                               data,
                               fixed,
                               random,
                               subject,
                               mixture,
                               ng,
                               B,
                               link,
                               iEst,
                               maxiter,
                               verbose,
                               nproc,
                               convB,
                               convL,
                               convG){
        #First : let's compute P(C|\tilde C) (\tilde C : A)
        pAlY = sapply(1:ng, function(g){
          return(as.numeric(predictClass(model, data)[,2] == g))
        })
        pClY = as.matrix(predictClass(model, data)[,3:(2+ng)])
        if(any(is.nan(pClY))) stop("NaN in posterior classification probability")
        
        pA = apply(pAlY, 2, mean)
        pClA = t(pAlY)%*%pClY/(model$ns*pA)
        
        if(det(pClA) == 0 | is.na(det(pClA))) stop("Computed error matrix is singular. One class might be empty")
        
        #Then : let's add it to dataset for each individual
        indivProb = pAlY%*%pClA
        indivProb = cbind(predictClass(model, data)[,1], indivProb)
        colnames(indivProb)[1] = subject
        data = merge(data, indivProb, by = subject)
        
        #Finally : we need to build the model for ppriors !
        arguments = list()
        
        if(length(fixed[[2]]) != 1){
          funOut = "multlcmm"
          arguments[["link"]] = link
        } else if(missing(link)){
          funOut = "hlme"
          arguments[["link"]] = NULL
        } else {
          funOut = "lcmm"
          arguments[["link"]] = link
        }
        
        arguments[["data"]] = data
        arguments[["fixed"]] = fixed
        arguments[["random"]] = random
        arguments[["subject"]] = subject
        arguments[["mixture"]] = mixture
        arguments[["classmb"]] = ~-1
        arguments[["ng"]] = ng
        arguments[["pprior"]] = c("prob1", "prob2")
        #technical options
        arguments[["maxiter"]] = maxiter
        arguments[["verbose"]] = verbose
        arguments[["nproc"]] = nproc
        arguments[["convB"]] = convB
        arguments[["convL"]] = convL
        arguments[["convG"]] = convG
        
        arguments[["B"]] = B[iEst]
        
        res = do.call(funOut, arguments)
        res$call = match.call()
        return(res)
      }
      
      #we need to build the model
      arguments[["data"]] = data
      arguments[["model"]] = model
      arguments[["fixed"]] = fixed
      arguments[["random"]] = random
      arguments[["subject"]] = subject
      arguments[["mixture"]] = mixture
      arguments[["ng"]] = ng
      arguments[["link"]] = link
      arguments[["iEst"]] = iEst
      #technical options
      arguments[["maxiter"]] = maxiter
      arguments[["verbose"]] = verbose
      arguments[["nproc"]] = nproc
      arguments[["convB"]] = convB
      arguments[["convL"]] = convL
      arguments[["convG"]] = convG
      
      
      #on ajoute des valeurs de base pour nos nouveaux estimateurs
      arguments[["B"]] = rep(0, nIn+nOut)
      #initial values
      if(missing(B)){
        arguments[["B"]][iEst] = strMod$best[(nMB+1):nStr]
      } else {
        arguments[["B"]][iEst] = B
      }
      
      funOut = "conditional"
    }
    
    #Xextern
    if(!missing(classmb)){
      #we need it all in a function in order to be able to use parametric bootstrap later
      
      #nEst : number of MB parameters in output model
      nEst1G = ncol(model.matrix(classmb, data))
      nEst = nEst1G*(ng-1)
      iEst = 1:nEst
      
      iKeepIn = 1:nIn
      iKeepOut = 1:nIn+nEst
      
      nOut = nEst
      
      #Id of varcov estimates to keep in bootstrap
      iVCKeep = iVCIn
      #Id of varcov estimates (none)
      iVCOut = c()
      
      conditional = function(data,
                               ng,
                               B,
                               iKeepIn,
                               iEst,
                               argumentsIn,
                               funIn,
                               nproc,
                               maxiter,
                               verbose){
        argumentsInEdit = argumentsIn
        argumentsInEdit[["B"]] = B[iKeepOut]
        argumentsInEdit[["maxiter"]] = 0
        argumentsInEdit[["verbose"]] = F
        argumentsInEdit[["nproc"]] = nproc
        model = do.call(funIn, argumentsInEdit)
        
        B = B[iEst]
        
        #First : let's compute P(\tilde C|C) (\tilde C : A)
        pAlY = sapply(1:ng, function(g){
          return(as.numeric(predictClass(model, data)[,2] == g))
        })
        pClY = as.matrix(predictClass(model, data)[,3:(2+ng)])
        if(any(is.nan(pClY))) stop("NaN in posterior classification probability")
        
        betas = c(model$best[1:(ng-1)], 0)
        pC = sapply(betas, function(b, betas){
          exp(b)/sum(exp(betas))
        }, betas=betas)
        
        pAlC = t(pClY)%*%pAlY/(model$ns*pC)
        
        if(det(pAlC) == 0 | is.na(det(pAlC))) stop("Computed error matrix is singular. One class might be empty")
        
        
        #Then : let's add the classification to the dataset for each individual
        indivProb = pAlY%*%t(pAlC)
        indivProb = cbind(predictClass(model, data)[,1], indivProb)
        colnames(indivProb) = c(subject, paste0("class", 1:ng))
        data = merge(data, indivProb, by = subject)
        
        for(id in unique(data[[subject]])){
          if(length(data[[subject]] > 1)) data = data[-which(data[[subject]] == id)[-1],]
        }
        
        #negative log likelihood function
        nLL = function(beta, y, X) {
          beta = matrix(ncol = ncol(y)-1, byrow = T, beta)
          denom <- apply(exp(X%*%beta),1,sum) + 1
          num_mat <- y*cbind(exp(X%*%beta), 1)
          num <- apply(num_mat,1,sum)
          vrais <- sum(log(num/denom))
          return(-vrais)
        }
        
        #frame
        argmf = list(
          formula = classmb,
          data = data
        )
        mf = do.call(model.frame, argmf)
        ns = nrow(mf)
        
        y = data[, paste0("class", 1:ng)]
        y = as.matrix(y)
        nBy = ncol(y)-1
        
        X = model.matrix(classmb, data)
        nEst = ncol(X)*nBy
        
        opt = marqLevAlg(b=B, fn=nLL, y=y, X=X, print.info = verbose, nproc = nproc, maxiter = maxiter)

        namesX = c("intercept", colnames(X)[colnames(X) != "(Intercept)"])
        names(opt$b) = c(sapply(namesX, FUN = function(i){
          return(paste0(i, " ", colnames(y)[-ng]))
        }))
        Names = list(Xnsnames = namesX,
                     ID = subject)
        
        N = c(nEst)
        
        res = list(best = opt$b,
                   V = opt$v,
                   conv = opt$istop,
                   loglik = -opt$fn.value,
                   ns = length(unique(data[[subject]])),
                   ng = ng,
                   idprob = rep(1, ncol(X)),
                   nv2 = ncol(X),
                   gconv = c(opt$ca, opt$cb, opt$rdm),
                   pprob = NULL,
                   Names = Names,
                   N = N,
                   call = match.call())
        
        return(res)
      }
      
      #Finally : we need to build the model arguments
      arguments[["data"]] = data
      arguments[["ng"]] = ng
      arguments[["iKeepIn"]] = iKeepIn
      arguments[["iEst"]] = iEst
      arguments[["argumentsIn"]] = argumentsIn
      arguments[["funIn"]] = funIn
      arguments[["nproc"]] = nproc
      arguments[["maxiter"]] = maxiter
      arguments[["verbose"]] = verbose
      
      #on ajoute des valeurs de base pour nos nouveaux estimateurs
      arguments[["B"]] = rep(0, nIn+nOut)
      #initial values
      if(!missing(B)){
        arguments[["B"]][iEst] = B
      }
      
      funOut = "conditional"
      
      
      
    }
  }
  
  if(varest != "paramBoot"){
    arguments[["B"]][iKeepOut] = model$best[iKeepIn]
    if(verbose){cat("Model estimation...\n\n")}
    #Model Estimation
    modOut = do.call(funOut, c(arguments))
  }
  
  if(varest == "Hessian"){
    if(method != "twoStageJoint") stop("Hessian variance estimation method only avaliable for 'twoStageJoint' method")
    
    if(verbose){cat("Variance estimation...\n\n")}
    nb11 = length(model$best)
    V11 = matrix(0, nb11, nb11)
    V11[upper.tri(V11, diag=T)] = model$V
    V11[lower.tri(V11, diag=F)] = t(V11)[lower.tri(V11, diag=F)]
    n1 = model$ns
    V11 = V11[iKeepIn, iKeepIn]
    
    nb22 = length(modOut$best)
    V22 = matrix(0, nb22, nb22)
    V22[upper.tri(V22, diag=T)] = modOut$V
    V22[lower.tri(V22, diag=F)] = t(V22)[lower.tri(V22, diag=F)]
    saveV22 = V22
    V22 = V22[iEst, iEst]
    n2 = modOut$ns
    
    modOut
    if(nproc == 1){
      I12 = -hessienne(modOut)
    } else {
      I12 = -hessienne(modOut, method = "deriva", nproc = nproc)
    }
    I12 = I12[iEst, iKeepOut]
    
    V = V22*n2 + (V22*n2) %*% (I12/n2) %*% ((n2/n1)*(V11*n1)) %*% t(I12/n2) %*% (V22*n2)
    V = V/n2
    saveV22[iEst, iEst] = V
    V = saveV22
    
    modOut$V = V[upper.tri(V, diag = TRUE)]
  }
  
  #Get Bootstrap Models
  if(varest == "paramBoot"){
    if(verbose){cat("Bootstrap estimation...\n\n")}
    est = estimates(model)
    
    Vin = matrix(0, length(est), length(est))
    Vin[upper.tri(Vin, diag = T)] = model$V
    Vin[lower.tri(Vin, diag = F)] = t(Vin)[lower.tri(t(Vin), diag=F)]
    
    est = est[iKeepIn]
    Vin = Vin[iKeepIn,iKeepIn]
    
    coefss = rmvnorm(M, est, Vin)
    coefss = as.data.frame(coefss)
    colnames(coefss) = names(model$best)[iKeepIn]
    #we just need to build back varcov into the coefs instead of cholesky matrix
    coefss = apply(coefss, 1, function(coefs, model, data, iVCKeep){
      if(funIn == "mpjlcmm"){
        varcovMods = longitudinal
      } else {
        varcovMods = list(model)
      }
      
      chols = coefs[iVCKeep]
      model$cholesky = chols
      
      varcov = c()
      countChol = 0
      for(varcovMod in varcovMods){
        ncolRandMod = ncol(model.matrix(formula(varcovMod$call$random), data))
        
        if(varcovMod$idiag){
          nChol = ncolRandMod
          
          vc = chols[1:nChol+countChol]
          
          varcov = c(varcov, (vc^2))
        } else {
          nChol = ncolRandMod*(ncolRandMod+1)/2
          
          #cholMatrix
          cholMatrix = matrix(0, ncolRandMod, ncolRandMod)
          cholMatrix[upper.tri(cholMatrix, diag = T)] = chols[1:nChol+countChol]
          
          vc = t(cholMatrix)%*%cholMatrix
          varcov = c(varcov, vc[upper.tri(vc, diag = T)])
        }
        
        countChol = countChol + nChol
      }
      coefs[iVCKeep] = varcov
      
      return(coefs)
    }, model = model, data = data, iVCKeep = iVCKeep)

    if(nproc > 1)
    {
      clust <- parallel::makeCluster(nproc)
      
      #load all loaded packages
      packages = loadedNamespaces()
      for(pack in packages){
        clusterExport(clust, "pack", environment())
        clusterEvalQ(clust, require(pack, character.only = T))
      }
      
      survivalMissing = missing(survival)
      fixedMissing = missing(fixed)
      modOuts <- parApply(clust, coefss, 2, function(coefs, arguments, iKeepOut, funOut, iEst, nVCIn, iVCOut, survivalMissing, fixedMissing){
        arguments[["B"]][iKeepOut] = coefs
        arguments[["nproc"]] = 1
        
        #Model Estimation
        modOut = do.call(funOut, c(arguments))
        
        if(!fixedMissing){
          #cholesky not varcov as output in best
          if(inherits(modOut, "mpjlcmm") | !all(modOut$idiag)){
            modOut$best[iVCOut] = modOut$cholesky[-(1:nVCIn)]
          } else {
            cholMatrix = matrix(NA, length(iVCOut), length(iVCOut))
            cholMatrix[upper.tri(cholMatrix, diag = T)] = modOut$cholesky[-(1:nVCIn)]
            modOut$best[iVCOut] = diag(cholMatrix)
          }
          
          #Residual Error : need to be the same sign across bootstrap iterations
          countParamBeforeErrLink = sum(modOut$N[1:3])
          sumny = 0
          for(k in 1:modOut$K){
            
            for(y in 1:modOut$ny[k]){
              sumny = sumny+1
              
              if(modOut$contrainte[k] == 2){
                countParamYBeforeErr = sum(modOut$Nprm[3+1:5*modOut$ny[k]-modOut$ny[k]+y])+countParamBeforeErrLink
                modOut$best[countParamYBeforeErr:(countParamYBeforeErr+modOut$Nprm[3+6*modOut$ny[k]-modOut$ny[k]+y])]
                
                if(modOut$linktype[sumny] == 0)
                  countParamYBeforeErrLink = sum(modOut$Nprm[3+1:8*modOut$ny[k]-modOut$ny[k]+y])+countParamBeforeErrLink
                modOut$best[countParamYBeforeErrLink]
              }
            }
            
            countParamBeforeErrLink = countParamBeforeErrLink+modOut$npmK[k]
            if(modOut$contrainte[k] == 0 | (modOut$contrainte[k] == 1 & modOut$linktype[sumny] == 0))
              modOut$best[countParamBeforeErrLink] = abs(modOut$best[countParamBeforeErrLink])
          }
        }
        
        #Survival Base Function : need to be the same sign across bootstrap iterations
        if(!survivalMissing){
          iSurvConstraint = 1:nSurvConstraint+modOut$N[1]
          modOut$best[iSurvConstraint] = abs(modOut$best[iSurvConstraint])
        }
        
        return(modOut)
      }, arguments = arguments, iKeepOut = iKeepOut, funOut = funOut, iEst = iEst, nVCIn = nVCIn, iVCOut = iVCOut, survivalMissing = survivalMissing, fixedMissing = fixedMissing)
      parallel::stopCluster(clust)
      
      #format output
      bests = as.data.frame(matrix(NA, nrow = nEst, ncol = M))
      Vs = list()
      conv = c()
      
      modOut = modOuts[[M]]
      for (i in 1:M){
        #output V and betas
        bests[,i] = modOut$best[iEst]
        
        V = matrix(NA, nOut, nOut)
        V[upper.tri(V, diag = T)] = modOut$V
        V[lower.tri(V, diag = F)] = t(V)[lower.tri(t(V), diag = F)]
        Vs = c(Vs, list(V[iEst, iEst]))
        
        conv = c(conv, modOut$conv)
      }
    } else {
      #estimate final models, extract usefull information
      bests = as.data.frame(matrix(NA, nrow = nEst, ncol = M))
      Vs = list()
      conv = c()
      for(i in 1:M){
        if(verbose) cat("==================== Bootstrap Iteration", i, "====================\n")
        
        arguments[["B"]][iKeepOut] = coefss[,i]
        
        #Model Estimation
        captured_log = capture.output({modOut = do.call(funOut, c(arguments))})
        
        if(!missing(fixed)){
          #cholesky not varcov as output in best
          modOut$best = estimates(modOut)
          
          #Residual Error : need to be the same sign across bootstrap iterations
          countParamBeforeErrLink = sum(modOut$N[1:3])
          sumny = 0
          for(k in 1:modOut$K){
            
            for(y in 1:modOut$ny[k]){
              sumny = sumny+1
              
              if(modOut$contrainte[k] == 2){
                countParamYBeforeErr = sum(modOut$Nprm[3+1:5*modOut$ny[k]-modOut$ny[k]+y])+countParamBeforeErrLink
                modOut$best[countParamYBeforeErr:(countParamYBeforeErr+modOut$Nprm[3+6*modOut$ny[k]-modOut$ny[k]+y])]
                
                if(modOut$linktype[sumny] == 0)
                  countParamYBeforeErrLink = sum(modOut$Nprm[3+1:8*modOut$ny[k]-modOut$ny[k]+y])+countParamBeforeErrLink
                modOut$best[countParamYBeforeErrLink]
              }
            }
            
            countParamBeforeErrLink = countParamBeforeErrLink+modOut$npmK[k]
            if(modOut$contrainte[k] == 0 | (modOut$contrainte[k] == 1 & modOut$linktype[sumny] == 0))
              modOut$best[countParamBeforeErrLink] = abs(modOut$best[countParamBeforeErrLink])
          }
        }
        
        #Survival Base Function : need to be the same sign across bootstrap iterations
        if(!missing(survival)){
          iSurvConstraint = 1:nSurvConstraint+modOut$N[1]
          modOut$best[iSurvConstraint] = abs(modOut$best[iSurvConstraint])
        }
        
        #output V and betas
        bests[,i] = modOut$best[iEst]
        
        V = matrix(NA, nOut, nOut)
        V[upper.tri(V, diag = T)] = modOut$V
        V[lower.tri(V, diag = F)] = t(V)[lower.tri(t(V), diag = F)]
        Vs = c(Vs, list(V[iEst, iEst]))
        
        conv = c(conv, modOut$conv)
      }
    }
    
    Mconv = sum(conv %in% c(1,3))
    if(Mconv <= 1){
      stop("No parametric boostrap iteration could converge.")
    }
    
    #compute variance
    bests = bests[,conv %in% c(1,3)]
    mb = apply(bests, 1, mean, na.rm = T)
    
    Vs = Vs[conv %in% c(1,3)]
    V2 = Reduce("+", Vs)/Mconv
    V1 = (Mconv+1)/((Mconv-1)*Mconv) * Reduce("+", lapply(bests, function(best){
      (best-mb)%*%t(best-mb)
    }))
    covar = V2 + V1
    V = matrix(0, nOut, nOut)
    V[iEst, iEst] = covar
    
    #Create output model, based on last estimated bootstrap model structure
    modOut$best[iKeepOut] = model$best[iKeepIn]
    modOut$best[iEst] = mb
    modOut$V = V[upper.tri(V, diag = T)]
    modOut$Mconv = Mconv
    
    #replace chol for varcov in best
    if(!missing(fixed)){
      modOut$cholesky[-(1:nVCIn)] = modOut$best[iVCOut]
      if(idiag){
        modOut$best[iVCOut] = modOut$cholesky[-(1:nVCIn)]^2
      } else {
        NVC = sqrt(2*(length(modOut$cholesky)-nVCIn)+1/4)-1/2
        cholMatrix = matrix(0, NVC, NVC)
        cholMatrix[upper.tri(cholMatrix, diag = T)] = modOut$best[iVCOut]
        vc = t(cholMatrix)%*%cholMatrix
        modOut$best[iVCOut] = vc[upper.tri(vc, diag = T)]
      }
    }
    
    #Computations on final model
    z = modOut$call
    z$posfix = NULL
    z$B = modOut$best
    z$maxiter = 0
    z$verbose = FALSE
    modelEdit = eval(z)
    if(method == "twoStageJoint") modOut$pprob = modelEdit$pprob
    if(!missing(survival)) modOut$predSurv = modelEdit$predSurv
    modOut$loglik = modelEdit$loglik
  }
  
  #if(method == "twoStageJoint" & !missing(fixed)) modOut$strMod = strMod
  #if(method == "twoStageJoint" & !missing(fixed) & funIn == "Jointlcmm") modOut$modNoSurv = modNoSurv
  
  modOut$call$data = substitute(data)
  
  #Il ne faut que les noms de variables dans le call sorti
  #Selon si c'est Yext ou Xext
  # if(!missing(fixed)){
  #   if(funIn == "mpjlcmm"){
  #     modOut$call$longitudinal = as.call(c(as.list(longCall), as.name("strMod")))
  #   } else if(funIn == "Jointlcmm"){
  #     modOut$call$longitudinal = as.call(list(as.name("list"), as.name("modNoSurv"), as.name("strMod")))
  #   } else {
  #     modOut$call$longitudinal = as.call(list(as.name("list"), substitute(model), as.name("strMod")))
  #   }
  # }
  # if(!missing(classmb)){
  #   if(is.null(argumentsIn[["longitudinal"]])){
  #     if(varest == "Hessian"){           #Because of hessienne mpj....
  #       modOut$call$longitudinal = longitudinalOut 
  #     } else {
  #       modOut$call$longitudinal = as.call(list(as.name("list"), substitute(model)))
  #     }
  #   } else {
  #     modOut$call$longitudinal = argumentsIn[["longitudinal"]]
  #   }
  # }
  
  #On remet longicall comme considere dans environnement present, pas dans celui de mpjlcmm
  if(inherits(modOut, "mpjlcmm")){
    longicall = vector("list", modOut$K)
    for(k in 1:modOut$K){
      longicall[[k]] = eval(modOut$call$longitudinal)[[k]]$call
    }
    modOut$longicall = longicall
  }
  
  cost = proc.time()-ptm
  
  #Object Class Creation and transformation from mpjlcmm
  if(!missing(survival)){
    if(method == "twoStageJoint"){
      best = modOut$best[iEst]
      V = matrix(NA, nOut, nOut)
      V[upper.tri(V, diag = T)] = modOut$V
      V = V[iEst, iEst]
      V = V[upper.tri(V, diag = T)]
      
      #select output
      N = modOut$N[c(2:3, 11+modOut$K+modOut$nbevt)]
      Nprm = modOut$Nprm[2:(2+modOut$nbevt)]
      Names = modOut$Names[c("Xnsnames", "ID", "Tnames", "TimeDepVar.name")]
      Names[["Xnsnames"]] = Names[["Xnsnames"]][as.logical(modOut$idspecif)]
      levels = modOut$levels[c("levelsdata", "levelssurv")]
      
      modOut = list(nbevt = modOut$nbevt, ng = modOut$ng, ns = modOut$ns, idcom = modOut$idcom,
                    idspecif = modOut$idspecif, idtdv = modOut$idtdv, loglik = modOut$loglik,
                    best = best, V = V, gconv = modOut$gconv, conv = modOut$conv, call = cl,
                    niter = modOut$niter, N = N, Nprm = Nprm, pprob = pprob, Names = Names,
                    logspecif = modOut$logspecif, predSurv = modOut$predSurv, typrisq = modOut$typrisq,
                    hazardtype = modOut$hazardtype, hazardnodes = modOut$hazardnodes, nz = modOut$nz,
                    scoretest = modOut$scoretest, na.action = modOut$na.action, levels = levels, 
                    AIC = 2*(length(best)-length(posfix)-modOut$loglik),
                    BIC = (length(best)-length(posfix))*log(modOut$ns)-2*modOut$loglik,
                    varest = varest, runtime = cost[3])
    }
    if(method == "conditional"){
      data$dummyY = 1
      subLoglik = hlme(dummyY~1,
                       data = data,
                       subject = subject,
                       maxiter = 0,
                       B = c(1, 0.000001))$loglik
      loglik = modOut$loglik - subLoglik
      
      best = modOut$best[iEst]
      V = matrix(NA, nOut, nOut)
      V[upper.tri(V, diag = T)] = modOut$V
      V = V[iEst, iEst]
      V = V[upper.tri(V, diag = T)]
      
      N = modOut$N[c(2:3, 10)]
      Nprm = modOut$N[2:3]
      Names = modOut$Names[c("Xnames2", "ID", "Tnames", "TimeDepVar.name")]
      names(Names)[names(Names) == "Xnames2"] = "Xnsnames"
      levels = modOut$levels[c("levelsdata", "levelssurv")]
      
      
      modOut = list(nbevt = 1, ng = modOut$ng, ns = modOut$ns, idcom = modOut$idcom,
                    idspecif = modOut$idspecif, idtdv = modOut$idtdv, loglik = loglik,
                    best = best, V = V, gconv = modOut$gconv, conv = modOut$conv, call = cl,
                    niter = modOut$niter, N = N, Nprm = Nprm, pprob = pprob, Names = Names,
                    logspecif = modOut$logspecif, predSurv = modOut$predSurv, typrisq = modOut$hazard[[1]],
                    hazardtype = modOut$hazard[[2]], hazardnodes = modOut$hazard[[3]], nz = modOut$hazard[[4]],
                    scoretest = modOut$scoretest, na.action = modOut$na.action, levels = levels,
                    AIC = 2*(length(best)-length(posfix)-loglik),
                    BIC = (length(best)-length(posfix))*log(modOut$ns)-2*loglik,
                    varest = varest, runtime = cost[3])
    }
    
    class(modOut) = c("externSurv", "externVar")
    
  }
  if(!missing(fixed)){
    if(method == "twoStageJoint"){
      
      #Get info
      gconv = modOut$gconv
      conv = modOut$conv
      niter = modOut$niter
      
      #With exclusively longitudinal external outcome
      modUpdate = update(modOut)
      modOut = modUpdate[[length(modUpdate)]]
      
      modOut$gconv = gconv
      modOut$conv = conv
      modOut$niter = niter
      modOut$pprob = pprob
      if(inherits(modOut, "multlcmm")) modOut$N[3] = modOut$N[3]-modOut$N[1]
      modOut$N[1] = 0
      modOut$idprob0 = rep(0, length(modOut$idprob0))
      modOut$call = cl
      modOut$varest = varest
      modOut$runtime = cost[3]
      
      V = matrix(NA, length(modOut$best), length(modOut$best))
      V[upper.tri(V, diag = T)] = modOut$V
      V[lower.tri(V, diag = F)] = t(V)[lower.tri(t(V), diag = F)]
      V = V[-c(1:(modOut$ng-1)), -c(1:(modOut$ng-1))]
      V = V[upper.tri(V, diag = T)]
      
      modOut$V = V
      modOut$best = modOut$best[-c(1:nMB)]
      
      class(modOut) = c(class(modOut), "externVar")
    }
    if(method == "conditional"){
      #Get info
      modOut$pprob = pprob
      if(inherits(modOut, "multlcmm")) modOut$N[3] = modOut$N[3]-modOut$N[1]
      modOut$call = cl
      modOut$varest = varest
      modOut$runtime = cost[3]
      
      class(modOut) = c(class(modOut), "externVar")
    }
  }
  if(!missing(classmb)){
    V = matrix(NA, nOut, nOut)
    V[upper.tri(V, diag = T)] = modOut$V
    V[lower.tri(V, diag = F)] = t(V)[lower.tri(t(V), diag = F)]
    V = V[iEst, iEst]
    V = V[upper.tri(V, diag = T)]
    
    best = modOut$best[iEst]
    
    #Select output
    N = modOut$N[1]
    Names = modOut$Names[c("Xnsnames", "ID")]
    levels = modOut$levels[c("levelsdata", "levelsclassmb")]
    
    modOut = list(ng = modOut$ng, ns = modOut$ns, idprob = modOut$idprob, nv2 = modOut$nv2,
                  loglik = modOut$loglik, best = best, V = V, gconv = modOut$gconv,
                  conv = modOut$conv, call = cl, niter = modOut$niter, N = N,
                  pprob = modOut$pprob, Names = Names, na.action = modOut$na.action,
                  AIC = 2*(length(best)-length(posfix)-modOut$loglik),
                  BIC = (length(best)-length(posfix))*log(modOut$ns)-2*modOut$loglik,
                  varest = varest, method = method, runtime = cost[3])
    
    class(modOut) = c("externX", "externVar")
  }
  
  if(verbose){cat(paste("The externVar program took", round(cost[3],2), "seconds\n"))}
  
  return(modOut)
  
}
