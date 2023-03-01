#' Estimation of external variable latent class models
#' 
#' This function allows to model external variables in relationship with previously modeled
#' latent class structures, either external outcomes or external class predictors.
#' The joint likelihood, for both external class predictors and external outcome
#' of the class. It returns a model object frome one of the lcmm package model classes.
#' 
#' A. DATA STRUCTURE
#' 
#' The \code{data} argument must follow specific structure for individual variables,
#' aka variables with a unique constant value for each subject. For an individual variable
#' given as external outcome, data value must be present only once per subject.
#' For an individual variable given as external class predictor,
#' data values must be present for each
#' 
#' B. VARIANCE ESTIMATION
#' 
#' Not taking into account first stage variance with specifing \code{"none"} may lead to
#' underestimation of the final variance, even though it is much quicker.
#' With big models it is recommended used variance estimation method \code{"none"} first before using
#' method \code{"paramBoot"} with parameter \code{"B"} being given the resulting estimates values
#' 
#' @param model original latent class structure model from which external variables
#' must be modeled.
#' @param fixed two sided linear formula object for specifying the
#' fixed-effects on external outcome variable.
#' The response outcome is on the left of \code{~} and the covariates are separated
#' by \code{+} on the right of the \code{~}. By default, an intercept is included.
#' @param mixture one-sided formula object for the class-specific fixed effects
#' in the model for the external outcome. Among the list of covariates included in fixed,
#' the covariates with class-specific regression parameters are entered in
#' mixture separated by \code{+}. By default, an intercept is included.
#' If no intercept, \code{-1} should be the first term included.
#' @param random an optional one sided linear formula object for specifying the
#' random-effects on external outcome. By default, no random effect is included.
#' @param subject name of the covariate representing the grouping structure.
#' Even without random effect. By default, the funciton will try to retrieve it
#' from input model
#' @param classmb optional one-sided formula describing the covariates in the
#' class-membership multinomial logistic model. These are external covariates or
#' external class predictors in this function.
#' @param varest character string indicating the method used to account for step one
#' variability when computing the variance estimation.
#' either "none", "paramBoot" or "calc" ("calc" is implemented for "twoStageJoint"
#' method only) Defaults to \code{"paramBoot"}
#' @param M integer number of parametrical boostrap iterations when varest is "paramBoot".
#' Default to 200.
#' @param idiag optional logical for the structure of the variance-covariance
#' matrix of the random-effects. If \code{FALSE}, a non structured matrix of
#' variance-covariance is considered (by default). If \code{TRUE} a diagonal
#' matrix of variance-covariance is considered.
#' @param nwg optional logical indicating if the variance-covariance of the
#' random-effects is class-specific. If \code{FALSE} the variance-covariance
#' matrix is common over latent classes (by default). If \code{TRUE} a
#' class-specific proportional parameter multiplies the variance-covariance
#' matrix in each class (the proportional parameter in the last latent class
#' equals 1 to ensure identifiability).
#' @param link link optional family of link functions to estimate for the external outcome.
#' Defaults to NULL, corresponding to hlme function
#' @param epsY optional definite positive real used to rescale the marker in (0,1)
#' when the beta link function is used. By default, epsY=0.5.
#' @param data Data frame containing the variables named in
#' \code{fixed}, \code{mixture}, \code{random}, \code{classmb} and \code{subject},
#' for both the current function arguments and the input models arguments
#' Check \code{details} to get information on the data structure, especially with
#' individual external outcomes
#' @param method character string representing the method to be used :
#' only joint likelihood two stage estimation "twoStageJoint" implemented at the moment
#' @param B optional but recommanded specification for the initial values for the
#' parameters. Two options are allowed: (1) a vector of initial values is entered.
#' (2) nothing is specified. A preliminary analysis involving the
#' estimation of a standard model is performed to choose initial values.
#' @param convB optional threshold for the convergence criterion based on the
#' parameter stability. By default, convB=0.0001.
#' @param convL optional threshold for the convergence criterion based on the
#' log-likelihood stability. By default, convL=0.0001.
#' @param convG optional threshold for the convergence criterion based on the
#' derivatives. By default, convG=0.0001.
#' @param longitudinal optional list of longitudinal models of type hlme,
#' lcmm or multlcmm used to build mpjlcmm input model for "twoStageJoint" method.
#' By default, the function will try to retrieve it from the input model
#' @param maxiter optional maximum number of iterations for the Marquardt
#' iterative algorithm. Defaults to 100
#' @param posfix Optional vector specifying the indices in vector B of the
#' parameters that should not be estimated. Default to NULL, all external parameters are
#' estimated.
#' @param verbose logical indicating whether information about computation should be
#' reported. Default to FALSE.
#' @param nproc the number cores for parallel computation. Default to 1 (sequential mode).
#' 
#' 
#' @examples
#' 
#' \dontrun{
#' # Example with hlme input model
#' mMMSE1_1 <- hlme(MMSE~age65+I(age65^2)+CEP,
#'                  random=~age65+I(age65^2),
#'                  subject="ID",
#'                  data=paquid)
#' mMMSE2_1 <- hlme(MMSE~age65+I(age65^2)+CEP,
#'                  random=~age65+I(age65^2),
#'                  subject="ID",
#'                  data=paquid,
#'                  ng=2,
#'                  mixture=~age65+I(age65^2),
#'                  B=random(mMMSE1_1))
#' 
#' modPlusCESDnone = externVar(mMMSE2_1,
#'                             fixed = CESD~age65+I(age65^2)+male,
#'                             random = ~age65+I(age65^2),
#'                             mixture = ~age65+I(age65^2),
#'                             subject="ID",
#'                             data=paquid,
#'                             verbose = T,
#'                             varest = "none",
#'                             method = "twoStageJoint",
#'                             B = c(7.8, 6.5, -0.3, 6, 0.8, -0.5, -2.5, 74, -64, 106, 16, -29, 9, 5.4))
#' 
#' modPlusCESDboot = externVar(mMMSE2_1,
#'                             fixed = CESD~age65+I(age65^2)+male,
#'                             random = ~age65+I(age65^2),
#'                             mixture = ~age65+I(age65^2),
#'                             subject="ID",
#'                             data=paquid,
#'                             verbose = T,
#'                             varest = "paramBoot",
#'                             M = 10,
#'                             method = "twoStageJoint",
#'                             B = c(7.8, 6.5, -0.3, 6, 0.8, -0.5, -2.5, 74, -64, 106, 16, -29, 9, 5.4))
#' 
#' modPlusCESDcalc = externVar(mMMSE2_1,
#'                             fixed = CESD~age65+I(age65^2)+male,
#'                             random = ~age65+I(age65^2),
#'                             mixture = ~age65+I(age65^2),
#'                             subject="ID",
#'                             data=paquid,
#'                             verbose = T,
#'                             varest = "calc",
#'                             method = "twoStageJoint",
#'                             B = c(7.8, 6.5, -0.3, 6, 0.8, -0.5, -2.5, 74, -64, 106, 16, -29, 9, 5.4))
#'                             
#' summary(modPlusCESDnone)
#' summary(modPlusCESDboot)
#' summary(modPlusCESDcalc)
#' 
#' # Example with Jointlcmm input model
#' m1 <- Jointlcmm(fixed= Ydep1~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#' random=~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#' subject='ID',
#' survival = Surv(Tevent,Event)~ X1+X2,
#' hazard="Weibull",
#' hazardtype="Specific",
#' ng=1,
#' data=data_lcmm)
#' m2 <- Jointlcmm(fixed= Ydep1~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                 mixture=~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                 random=~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                 subject='ID',
#'                 survival = Surv(Tevent,Event)~X1+mixture(X2),
#'                 hazard=c("Weibull", "Weibull"),
#'                 hazardtype=c("Specific", "Specific"),
#'                 ng=2,
#'                 data=data_lcmm,
#'                 B=m1)
#' 
#' m2plusY2none = externVar(m2,
#'                          method = "twoStageJoint",
#'                          fixed = Ydep2~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                          random = ~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                          mixture=~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                          subject='ID',
#'                          idiag = TRUE,
#'                          varest = "none",
#'                          data=data_lcmm)
#' m2plusY2boot = externVar(m2,
#'                          method = "twoStageJoint",
#'                          fixed = Ydep2~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                          random = ~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                          mixture=~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                          subject='ID',
#'                          idiag = TRUE,
#'                          varest = "paramBoot",
#'                          M = 10,
#'                          data=data_lcmm,
#'                          B = m2plusY2none$best[-m2plusY2none$call$posfix])
#' m2plusY2calc = externVar(m2,
#'                          method = "twoStageJoint",
#'                          fixed = Ydep2~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                          random = ~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                          mixture=~ns(Time, knots = c(1), Boundary.knots = c(0.2, 2.1)),
#'                          subject='ID',
#'                          idiag = TRUE,
#'                          varest = "calc",
#'                          data=data_lcmm,
#'                          B = m2plusY2none$best[-m2plusY2none$call$posfix])
#' summary(m2plusY2none)
#' summary(m2plusY2boot)
#' summary(m2plusY2calc)
#' }
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
                     idiag = FALSE,
                     nwg = FALSE,
                     link = NULL,
                     intnodes=NULL,
                     epsY = 0.5,
                     data,
                     longitudinal,
                     method,
                     varest = "paramBoot",
                     M = 200,
                     B,
                     convB = 0.0001,
                     convL = 0.0001,
                     convG = 0.0001,
                     maxiter = 100,
                     posfix,
                     verbose = FALSE,
                     nproc = 1){
  
  ptm <- proc.time()
  
  if(missing(model)) stop("model argument must be given")
  if(!class(model) %in% c("hlme", "lcmm", "multlcmm", "Jointlcmm", "mpjlcmm")) stop('input models class must be either "hlme", "lcmm", "multlcmm", "Jointlcmm" or "mpjlcmm"')
  if(model$conv == 2) warning("input model did not fully converge")
  if(missing(fixed) & missing(classmb)) stop("Either external outcome in fixed or external class predictor in classmb must be given")
  if(!missing(fixed) & !missing(classmb)) stop("Either external outcome in fixed or external class predictor in classmb must be given")
  if(missing(method) | !method %in% c("twoStageJoint")) stop('Method must be either "twoStageJoint"')
  if(model$ng == 1) stop("Input model does not have latent class structure (ng=1)")
  if(!varest %in% c("none", "paramBoot", "calc")) stop('Variance estimation method "varest" must be either "none", "paramBoot" or "calc"')
  if(!is.null(link) & missing(fixed)) stop("The argument link is not to be used with external class predictor")
  

  if(missing(posfix)) posfix = c()
  
  cl = match.call()
  
  #Informations about input model
  argumentsIn = as.list(model$call)
  funIn = as.character(argumentsIn[[1]])
  argumentsIn[[1]] = NULL
  ng = model$ng
  nIn = length(model$best)
  
  #Get subject
  if(missing(subject)){
    if (model$call$subject %in% colnames(data)){
      subject = model$call$subject
    } else {
      stop("The argument subject must be specified if different from the subject argument used in the input model")
    }
  } 
  
  #Get longitudinal
  if(funIn == "mpjlcmm"){
    if(missing(longitudinal)){
      if(class(try(eval(as.call(model$call$longitudinal)))) == "try-error"){
        stop("longitudinal not found in environment, please specify longitudinal argument with a list of longitudinal models used to build mpjlcmm input model")
      }
      longCall = model$call$longitudinal
      longitudinal = eval(longCall)
    } else if (!missing(longitudinal)) {
      longCall = substitute(longitudinal)
    }
    
    K = length(longitudinal)
    for(k in 1:K){
      cl = as.list(longitudinal[[k]]$call)
      cl[["computeDiscrete"]] = FALSE
      longitudinal[[k]]$call = as.call(cl)
      
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
    iVCIn = sum(model$N[1:3]) + 1:nVCIn
  } else {
    nVCIn = model$N[3]
    iVCIn = sum(model$N[1:2]) + 1:nVCIn
  }
  
  
  #Change general model arguments
  if(method == "twoStageJoint"){
    
    #Arguments communs
    arguments = list()
    
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
    #survival
    arguments[["survival"]] =  argumentsIn[["survival"]]
    arguments[["hazard"]] =  argumentsIn[["hazard"]]
    arguments[["hazardtype"]] =  argumentsIn[["hazardtype"]]
    arguments[["hazardnodes"]] =  argumentsIn[["hazardnodes"]]
    arguments[["hazardrange"]] =  argumentsIn[["hazardrange"]]
    arguments[["TimeDepVar"]] =  argumentsIn[["TimeDepVar"]]
    arguments[["logscale"]] =  argumentsIn[["logscale"]]
    
    #Yextern
    if(!missing(fixed)){
      
      funOut = "mpjlcmm"
      
      #Erreurs
      
      if(missing(mixture)) mixture = ~1
      if(missing(random)) random = ~-1
      
      if(!inherits(fixed,"formula")) stop("The argument fixed must be a formula")
      if(!inherits(mixture,"formula")) stop("The argument mixture must be a formula")
      if(!inherits(random,"formula")) stop("The argument random must be a formula")
      
      if(!class(model) %in% c("hlme", "lcmm", "multlcmm", "Jointlcmm", "mpjlcmm")){
        stop('The input model is not a supported class for the "twoStageJoint" method')
      }
      
      if(length(fixed[[2]]) != 1){
        argfunctionStrMod = "multlcmm"
      } else if(is.null(link)){
        argfunctionStrMod = "hlme"
      } else {
        argfunctionStrMod = "lcmm"
      }
      
      #Jointlcmm needs transformation into lcmm to be put into longitudinal.
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
      argumentsStrMod[["maxiter"]] = 0
      argumentsStrMod[["B"]] = as.name("strMod")
      
      strMod = do.call(argfunctionStrMod, c(argumentsStrMod))
      
      #Let's change strMod's saved call
      strMod$call$data = substitute(data)
      
      #Now join the initial and the external outcome models
      
      #Informations about external outcome model
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
      
      #input VC, only for keep betaa
      iVCKeep = iVCIn
      
      #Id of varcov estimates
      if(inherits(strMod, "multlcmm")){
        nVCStr = strMod$N[4]
        iVCStr = sum(strMod$N[1:3]) + 1:nVCStr
      } else {
        nVCStr = strMod$N[3]
        iVCStr = sum(strMod$N[1:2]) + 1:nVCStr
      }
      iVCOut = length(iKeepIn) + iVCStr - nMB
      
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
      
      #changement de ma fonction en mpjlcmm (pour avec la variance)
      if(funIn != "mpjlcmm"){
        argumentsMpj = list()
        argumentsMpj[["longitudinal"]] = list(model)
        argumentsMpj[["maxiter"]] = 0
        argumentsMpj[["ng"]] = ng
        argumentsMpj[["subject"]] = subject
        argumentsMpj[["data"]] = argumentsIn[["data"]]
        if (!is.null(model$call$classmb)){
          argumentsMpj[["classmb"]] = model$call$classmb
        }
        
        modelMpj = do.call("mpjlcmm", argumentsMpj)
        modelMpj$V = model$V
      }
      
      funOut = "mpjlcmm"
      
      if(is.null(arguments[["classmb"]])){
        oldclassmb = ~ 1
      } else {
        oldclassmb = formula(arguments[["classmb"]])
      }
      
      #nInMB : number of MB parameters in input model
      nInMB = ncol(model.matrix(oldclassmb, data))*(ng-1)
      
      #nEst : number of MB parameters in output model
      nEst1G = ncol(model.matrix(classmb, data))
      nEst = nEst1G*(ng-1)
      
      #nOut : nuber of total final parameters
      nOut = nIn - nInMB + nEst
      
      #index
      iKeepOut = (nEst+1):nOut
      iKeepIn = (nInMB+1):nIn
      iEst = 1:nEst
      
      #input VC, only for keep betas
      iVCKeep = iVCIn-nInMB
      
      #Id of varcov estimates (none)
      iVCOut = c()
      
      #On recree tous nos arguments
      if(is.null(argumentsIn[["longitudinal"]])){
        arguments[["longitudinal"]] = list(model)
      } else {
        arguments[["longitudinal"]] = longitudinal
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
      
      #### weird... "Input" model has been changed because we need it to build output model
      #To explain : "hessienne" function needs mpjlcmm model. So i need mpj for paramBoot
      #And we need to create our longitudinal (and logicall) now before changinng "input" moded.
      #I really hate that hessian is only for mpj......
      if(funIn != "mpjlcmm" & varest == "calc"){
        longitudinalOut = as.call(list(as.name("list"), substitute(model)))
        
        modelOld = model
        model = modelMpj
        model$V = modelOld$V
      }
      
      
      
    }
  }
  
  if(varest != "paramBoot"){
    arguments[["B"]][iKeepOut] = model$best[iKeepIn]
    if(verbose){cat("Model estimation...\n\n")}
    #Model Estimation
    modOut = do.call(funOut, c(arguments))
  }
  
  if(varest == "calc"){
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
    est = model$best
    
    #cholesky not varcov
    if(funIn %in% c("mpjlcmm", "Jointlcmm") | !all(as.logical(model$idiag))){
      est[iVCIn] = model$cholesky
    } else {
      cholMatrix = matrix(NA, nVCIn, nVCIn)
      cholMatrix[upper.tri(cholMatrix, diag = T)] = model$cholesky
      est[iVCIn] = diag(cholMatrix)
    }
    
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
      cl <- parallel::makeCluster(nproc)
      
      #load all loaded packages
      packages = loadedNamespaces()
      for(pack in packages){
        clusterExport(cl, "pack", environment())
        clusterEvalQ(cl, require(pack, character.only = T))
      }
      
      modOuts <- parApply(cl, coefss, 2, function(coefs, arguments, iKeepOut, funOut, iEst, nVCIn, iVCOut){
        arguments[["B"]][iKeepOut] = coefs
        arguments[["nproc"]] = 1
        
        #Model Estimation
        modOut = do.call(funOut, c(arguments))
        
        #cholesky not varcov as output in best
        if(class(modOut) == "mpjlcmm" | !all(modOut$idiag)){
          modOut$best[iVCOut] = modOut$cholesky[-(1:nVCIn)]
        } else {
          cholMatrix = matrix(NA, length(iVCOut), length(iVCOut))
          cholMatrix[upper.tri(cholMatrix, diag = T)] = modOut$cholesky[-(1:nVCIn)]
          modOut$best[iVCOut] = diag(cholMatrix)
        }
        
        return(modOut)
      }, arguments = arguments, iKeepOut = iKeepOut, funOut = funOut, iEst = iEst, nVCIn = nVCIn, iVCOut = iVCOut)
      parallel::stopCluster(cl)
      
      #format output
      bests = as.data.frame(matrix(NA, nrow = nEst, ncol = M))
      Vs = list()
      conv = c()
      for (i in 1:M){
        modOut = modOuts[[i]]
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
        modOut = do.call(funOut, c(arguments))
        
        #cholesky not varcov as output in best
        if(class(modOut) == "mpjlcmm" | !all(modOut$idiag)){
          modOut$best[iVCOut] = modOut$cholesky[-(1:nVCIn)]
        } else {
          cholMatrix = matrix(NA, length(iVCOut), length(iVCOut))
          cholMatrix[upper.tri(cholMatrix, diag = T)] = modOut$cholesky[-(1:nVCIn)]
          modOut$best[iVCOut] = diag(cholMatrix)
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
    if(Mconv == 0){
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
  
  if(method == "twoStageJoint" & !missing(fixed)) modOut$strMod = strMod
  if(method == "twoStageJoint" & !missing(fixed) & funIn == "Jointlcmm") modOut$modNoSurv = modNoSurv
  
  modOut$call$data = substitute(data)
  
  #Il ne faut que les noms de variables dans le call sorti
  #Selon si c'est Yext ou Xext
  if(!missing(fixed)){
    if(funIn == "mpjlcmm"){
      modOut$call$longitudinal = as.call(c(as.list(longCall), as.name("strMod")))
    } else if(funIn == "Jointlcmm"){
      modOut$call$longitudinal = as.call(list(as.name("list"), as.name("modNoSurv"), as.name("strMod")))
    } else {
      modOut$call$longitudinal = as.call(list(as.name("list"), substitute(model), as.name("strMod")))
    }
  }
  if(!missing(classmb)){
    if(is.null(argumentsIn[["longitudinal"]])){
      if(varest == "calc"){           #Because of hessienne mpj....
        modOut$call$longitudinal = longitudinalOut 
      } else {
        modOut$call$longitudinal = as.call(list(as.name("list"), substitute(model)))
      }
    } else {
      modOut$call$longitudinal = argumentsIn[["longitudinal"]]
    }
  }
  
  #On remet longicall comme considéré dans l'environnement présent, pas dans celui de mpjlcmm
  if(inherits(modOut, "mpjlcmm")){
    longicall = vector("list", modOut$K)
    for(k in 1:modOut$K){
      longicall[[k]] = eval(modOut$call$longitudinal)[[k]]$call
    }
    modOut$longicall = longicall
  }
  
  cost = proc.time()-ptm
  if(verbose){cat(paste("The externVar program took", round(cost[3],2), "seconds\n"))}
  modOut$runtime = cost[3]
  
  return(modOut)
  
}
