#' Estimation of external variable latent class models
#' 
#' This function allows to model external variables in relationship with previously modeled
#' latent class structures. The joint likelihood, for both external class
#' predictors and external outcome of the class.
#' 
#' @param model original model with latent class structure with which external variables
#' must be modeled.
#' @param method character string representing the method to be used :
#' either "twoStageJoint"
#' @param fixed two sided linear formula object for specifying the fixed-effects on external variable
#' @param mixture one-sided formula object for the class-specific fixed effects in the linear mixed model
#' By default, an intercept is included. If no intercept, \code{-1} should be the first term included.
#' @param random optional one sided linear formula object for specifying the random-effects on external variable
#' By default, an intercept is included. If no intercept, \code{-1} should be the first term included.
#' @param subject name of the covariate representing the grouping structure.
#' Even without random effect, either the subject variable used in model or
#' the subject variable specified as parameter must be present in data.
#' @param classmb optional one-sided formula describing the covariates in the
#' class-membership multinomial logistic model
#' @param varest character string indicating the method used to account for step one variability
#' when computing the variance estimation.
#' either "None", "ParamBoot" or "calc" ("calc" is for "twoStageJoint" method only)
#' @param idiag optional logical for the structure of the variance-covariance
#' matrix of the random-effects. If \code{FALSE}, a non structured matrix of
#' variance-covariance is considered (by default).  If \code{TRUE} a diagonal
#' matrix of variance-covariance is considered.
#' @param nwg optional logical indicating if the variance-covariance of the
#' random-effects is class-specific. If \code{FALSE} the variance-covariance
#' matrix is common over latent classes (by default). If \code{TRUE} a
#' class-specific proportional parameter multiplies the variance-covariance
#' matrix in each class (the proportional parameter in the last latent class
#' equals 1 to ensure identifiability).
#' @param link link optional family of link functions to estimate for the external outcome. Defaults to NULL, corresponding to hlme function
#' @param epsY 
#' @param data : optional data frame containing the variables named used elsewhere.
#' IMPORTANT : unique individual variables (such as gender or educational level)
#' MUST only appear once per subject AND at the same visit time.
#' @param B : optional specification for the initial values for the parameters.
#' Two options are allowed: (1) a vector of initial values is entered.
#' (2) nothing is specified. A preliminary analysis involving the
#' estimation of a standard model is performed to choose initial
#' values.
#' @param longitudinal optional list of longitudinal models of type hlme,
#' lcmm or multlcmm used to build mpjlcmm input model for "twoStageJoint" method.
#' By default, the function will try to retrieve it from input model call.
#' @param maxiter optional maximum number of iterations for the Marquardt
#' iterative algorithm
#' @param posfix Optional vector specifying the indices in vector B of the
#' parameters that should not be estimated. Default to NULL, all external parameters are
#' estimated.
#' @param verbose logical indicating whether information about computation should be
#' reported. Default to TRUE.
#' @param M integer number of parametrical boostrap iterations when varest is "paramBoot". Default to 200.
#' @param nproc the number cores for parallel computation. Default to 1 (sequential mode).
#' 
#' @export
#' 
#' 
#' 


externVar = function(model,
                     method,
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
                     B,
                     convB = 0.0001,
                     convL = 0.0001,
                     convG = 0.0001,
                     longitudinal,
                     maxiter = 100,
                     posfix,
                     verbose = TRUE,
                     varest = "none",
                     M = 200,
                     nproc = 1){
  
  if(missing(fixed) & missing(classmb)) stop("Either external outcome in fixed or external class predictor in classmb must be given")
  if(!missing(fixed) & !missing(classmb)) stop("Both external outcome in fixed and external class predictor in classmb is not supported by this function at the moment")
  if(missing(method) | !method %in% c("twoStageJoint")) stop('Method must be either "twoStageJoint"')
  if(missing(subject)) stop("The argument subject must be specified in any model even without random-effects")
  if(model$ng == 1) stop("Input model does not have latent class structure (ng=1)")
  if(!varest %in% c("none", "paramBoot", "calc")) stop('Variance estimation method "varest" must be either "none", "paramBoot" or "calc"')
  if(!is.null(link) & missing(fixed)) stop("The argument link is not to be used with external class predictor")
  
  if(missing(random)) random = ~-1
  if(missing(posfix)) posfix = c()
    
  cl = match.call()
  
  #Informations about input model
  argumentsIn = as.list(model$call)
  funIn = as.character(argumentsIn[[1]])
  argumentsIn[[1]] = NULL
  ng = model$ng
  nIn = length(model$best)
  
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
  
  
  #Change general model arguments
  if(method == "twoStageJoint"){
    
    #Arguments communs
    arguments = list()
    
    arguments[["data"]] = substitute(data)
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
      argumentsStrMod[["data"]] = substitute(data)
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
      
      #Liste des arguments
      #on fixe nos paramètres
      arguments[["posfix"]] = unique(c(iKeepOut, posfix))
      #On donne les modèles
      if(funIn == "mpjlcmm"){
        arguments[["longitudinal"]] = as.call(c(as.list(longCall), as.name("strMod")))
      } else if(funIn == "Jointlcmm"){
        arguments[["longitudinal"]] = as.call(list(as.name("list"), as.name("modNoSurv"), as.name("strMod")))
      } else {
        arguments[["longitudinal"]] = as.call(list(as.name("list"), substitute(model), as.name("strMod")))
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
      
      #changement de ma fonction en mpjlcmm (pour avec la variance)
      if(funIn != "mpjlcmm"){
        argumentsMpj[["longitudinal"]] = as.call(list(as.name("list"), substitute(model)))
        argumentsMpj[["maxiter"]] = 0
        argumentsMpj[["ng"]] = ng
        argumentsMpj[["subject"]] = subject
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
      
      #On recrée tous nos arguments
      if(is.null(argumentsIn[["longitudinal"]])){
        arguments[["longitudinal"]] = as.call(list(as.name("list"), substitute(model)))
      } else {
        arguments[["longitudinal"]] = argumentsIn[["longitudinal"]]
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
      
      #on fixe nos paramètres
      arguments[["posfix"]] = unique(c(posfix, iKeepOut))
      
      #### weird... Input model has been changed because we need to output model to be
      if(funIn != "mpjlcmm"){
        modelOld = model
        model = modelMpj
        model$V = modelOld$V
      }
      
    }
  }
  
  if(varest != "paramBoot"){
    arguments[["B"]][iKeepOut] = model$best[iKeepIn]
    
    #Model Estimation
    modOut = do.call(funOut, c(arguments))
  }
  
  if(varest == "calc"){
    
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
    
    modOutModif = modOut
    modOutModif$call$longitudinal = eval(modOutModif$call$longitudinal)
    modOutModif$call$data = eval(modOutModif$call$data)
    I12 = -hessienne(modOutModif)
    I12 = I12[iEst, iKeepOut]
    
    V = V22*n2 + (V22*n2) %*% (I12/n2) %*% ((n2/n1)*(V11*n1)) %*% t(I12/n2) %*% (V22*n2)
    V = V/n2
    saveV22[iEst, iEst] = V
    V = saveV22
    
    modOut$V = V[upper.tri(V, diag = TRUE)]
  }
  
  #Get Bootstrap Models
  if(varest == "paramBoot"){
    est = model$best
    
    #cholesky not varcov
    whereRand = substr(names(est), 1, 7) == "varcov "
    if(funIn == "mpjlcmm" | !all(model$idiag)){
      est[whereRand] = model$cholesky #Warning : quickest way to do this, but maybe not the best/safest
    } else {
      chol = matrix(NA, sum(whereRand), sum(whereRand))
      chol[upper.tri(chol, diag = T)] = model$cholesky
      est[substr(names(est), 1, 7) == "varcov "] = diag(chol)
    }
    
    Vin = matrix(0, length(est), length(est))
    Vin[upper.tri(Vin, diag = T)] = model$V
    Vin[lower.tri(Vin, diag = F)] = t(Vin)[lower.tri(t(Vin), diag=F)]
    
    est = est[iKeepIn]
    Vin = Vin[iKeepIn,iKeepIn]
    
    coefss = rmvnorm(M, est, Vin)
    coefss = as.data.frame(coefss)
    colnames(coefss) = names(model$best)[iKeepIn]
    whereRand = whereRand[iKeepIn]
    #we just need to build back varcov into the coefs instead of cholesky matrix
    coefss = apply(coefss, 1, function(coefs, model, data, whereRand){
      if(funIn == "mpjlcmm"){
        varcovMods = eval(model$call$longitudinal)
      } else {
        varcovMods = list(model)
      }
      
      chols = coefs[whereRand]
      model$cholesky = chols
      
      varcov = c()
      countChol = 0
      for(varcovMod in varcovMods){
        ncolRandMod = ncol(model.matrix(formula(varcovMod$call$random), data))
        
        if(varcovMod$idiag){
          vc = chols[countChol:ncolRandMod]
          
          varcov = c(varcov, (vc^2))
        } else {
          nChol = ncolRandMod*(ncolRandMod+1)/2
          
          #cholMatrix
          cholMatrix = matrix(0, ncolRandMod, ncolRandMod)
          cholMatrix[upper.tri(cholMatrix, diag = T)] = chols[countChol:nChol]
          
          countChol = countChol + nChol
          
          vc = t(cholMatrix)%*%cholMatrix
          varcov = c(varcov, vc[upper.tri(vc, diag = T)])
        }
      }
      coefs[whereRand] = varcov
      
      return(coefs)
    }, model = model, data = data, whereRand = whereRand)
    
    #estimate final models, extract usefull information
    bests = as.data.frame(matrix(NA, nrow = nEst, ncol = M))
    Vs = list()
    conv = c()
    for(i in 1:M){
      arguments[["B"]][iKeepOut] = coefss[,i]
      
      #Model Estimation
      modOut = do.call(funOut, c(arguments))
      
      bests[,i] = modOut$best[iEst]
      
      V = matrix(NA, nOut, nOut)
      V[upper.tri(V, diag = T)] = modOut$V
      V[lower.tri(V, diag = F)] = t(V)[lower.tri(t(V), diag = F)]
      Vs = c(Vs, list(V[iEst, iEst]))
      
      conv = c(conv, modOut$conv)
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
    V1 = (Mconv+1)/((Mconv-1)*Mconv) * Reduce("+", apply(bests, 2, function(best){
      (best-mb)%*%t(best-mb)
    }, simplify = F))
    covar = V2 + V1
    V = matrix(0, nOut, nOut)
    V[iEst, iEst] = covar
    
    #Create output model, based on last estimated bootstrap model structure
    modOut$best[iKeepOut] = model$best[iKeepIn]
    modOut$best[iEst] = mb
    modOut$V = V[upper.tri(V, diag = T)]
    modOut$Mconv = Mconv
  }
  
  if(method == "twoStageJoint" & !missing(fixed)) modOut$strMod = strMod
  if(method == "twoStageJoint" & !missing(fixed) & funIn == "Jointlcmm") modOut$modNoSurv = modNoSurv
  
  return(modOut)
  
}
