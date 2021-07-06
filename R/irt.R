#' Estimation of mutlivariate mixed-effect models and multivariate latent class
#' mixed-effect models for multivariate longitudinal outcomes of possibly
#' multiple types (continuous Gaussian, continuous non-Gaussian - curvilinear)
#' that measure the same underlying latent process.
#' 
#' This function constitutes a multivariate extension of function \code{lcmm}.
#' It fits multivariate mixed models and multivariate latent class mixed models
#' for multivariate longitudinal outcomes of different types. It handles
#' continuous longitudinal outcomes (Gaussian or non-Gaussian, curvilinear) as
#' well as bounded quantitative and discrete longitudinal outcomes. Next
#' version will also handle ordinal outcomes.  The model assumes that all the
#' outcomes measure the same underlying latent process defined as their common
#' factor, and each outcome is related to this latent common factor by a
#' specific parameterized link function.  At the latent process level, the
#' model estimates a standard linear mixed model or a latent class linear mixed
#' model when heterogeneity in the population is investigated (in the same way
#' as in function \code{hlme}). Parameters of the nonlinear link functions and
#' of the latent process mixed model are estimated simultaneously using a
#' maximum likelihood method.
#' 
#' 
#' A. THE PARAMETERIZED LINK FUNCTIONS
#' 
#' \code{multlcmm} function estimates multivariate latent class mixed models
#' for different types of outcomes by assuming a parameterized link function
#' for linking each outcome Y_k(t) with the underlying latent common factor
#' L(t) they measure. To fix the latent process dimension, we chose to
#' constrain at the latent process level the (first) intercept of the latent
#' class mixed model at 0 and the standard error of the first random effect at
#' 1.
#' 
#' 1. With the "linear" link function, 2 parameters are required for the
#' following transformation (Y(t) - b1)/b2
#' 
#' 2. With the "beta" link function, 4 parameters are required for the
#' following transformation: [ h(Y(t)',b1,b2) - b3]/b4 where h is the Beta CDF
#' with canonical parameters c1 and c2 that can be derived from b1 and b2 as
#' c1=exp(b1)/[exp(b2)*(1+exp(b1))] and c2=1/[exp(b2)*(1+exp(b1))], and Y(t)'
#' is the rescaled outcome i.e. Y(t)'= [ Y(t) - min(Y(t)) + epsY ] / [
#' max(Y(t)) - min(Y(t)) +2*epsY ].
#' 
#' 3. With the "splines" link function, n+2 parameters are required for the
#' following transformation b_1 + b_2*I_1(Y(t)) + ... + b_{n+2} I_{n+1}(Y(t)),
#' where I_1,...,I_{n+1} is the basis of quadratic I-splines. To constraint the
#' parameters to be positive, except for b_1, the program estimates b_k^* (for
#' k=2,...,n+2) so that b_k=(b_k^*)^2. This parameterization may lead in some
#' cases to problems of convergence that we are currently addressing.
#' 
#' Details of these parameterized link functions can be found in the papers:
#' Proust-Lima et al. (Biometrics 2006) and Proust-Lima et al. (BJMSP 2013).
#' 
#' B. THE VECTOR OF PARAMETERS B
#' 
#' The parameters in the vector of initial values \code{B} or in the vector of
#' maximum likelihood estimates \code{best} are included in the following
#' order: (1) ng-1 parameters are required for intercepts in the latent class
#' membership model, and if covariates are included in \code{classmb}, ng-1
#' paramaters should be entered for each one; (2) for all covariates in
#' \code{fixed}, one parameter is required if the covariate is not in
#' \code{mixture}, ng paramaters are required if the covariate is also in
#' \code{mixture}; When ng=1, the intercept is not estimated and no parameter
#' should be specified in \code{B}. When ng>1, the first intercept is not
#' estimated and only ng-1 parameters should be specified in \code{B}; (3) for
#' all covariates included with \code{contrast()} in \code{fixed}, one
#' supplementary parameter per outcome is required excepted for the last
#' outcome for which the parameter is not estimated but deduced from the
#' others; (4) if \code{idiag=TRUE}, the variance of each random-effect
#' specified in \code{random} is required excepted the first one (usually the
#' intercept) which is constrained to 1. (5) if \code{idiag=FALSE}, the
#' inferior triangular variance-covariance matrix of all the random-effects is
#' required excepted the first variance (usually the intercept) which is
#' constrained to 1. (5) only if \code{nwg=TRUE} and \code{ng}>1, ng-1
#' parameters for class-specific proportional coefficients for the variance
#' covariance matrix of the random-effects; (6) if \code{cor} is specified, the
#' standard error of the Brownian motion or the standard error and the
#' correlation parameter of the autoregressive process; (7) the standard error
#' of the outcome-specific Gaussian errors (one per outcome); (8) if
#' \code{randomY=TRUE}, the standard error of the outcome-specific random
#' intercept (one per outcome); (9) the parameters of each parameterized link
#' function: 2 for "linear", 4 for "beta", n+2 for "splines" with n nodes.
#' 
#' C. CAUTIONS REGARDING THE USE OF THE PROGRAM
#' 
#' Some caution should be made when using the program. Convergence criteria are
#' very strict as they are based on the derivatives of the log-likelihood in
#' addition to the parameter and log-likelihood stability. In some cases, the
#' program may not converge and reach the maximum number of iterations fixed at
#' 100. In this case, the user should check that parameter estimates at the
#' last iteration are not on the boundaries of the parameter space.
#' 
#' If the parameters are on the boundaries of the parameter space, the
#' identifiability of the model is critical. This may happen especially with
#' splines parameters that may be too close to 0 (lower boundary) or classmb
#' parameters that are too high or low (perfect classification). When
#' identifiability of some parameters is suspected, the program can be run
#' again from the former estimates by fixing the suspected parameters to their
#' value with option posfix. This usually solves the problem. An alternative is
#' to remove the parameters of the Beta of Splines link function from the
#' inverse of the Hessian with option partialH.
#' 
#' If not, the program should be run again with other initial values, with a
#' higher maximum number of iterations or less strict convergence tolerances.
#' 
#' Specifically when investigating heterogeneity (that is with ng>1): (1) As
#' the log-likelihood of a latent class model can have multiple maxima, a
#' careful choice of the initial values is crucial for ensuring convergence
#' toward the global maximum. The program can be run without entering the
#' vector of initial values (see point 2). However, we recommend to
#' systematically enter initial values in \code{B} and try different sets of
#' initial values. (2) The automatic choice of initial values we provide
#' requires the estimation of a preliminary linear mixed model. The user should
#' be aware that first, this preliminary analysis can take time for large
#' datatsets and second, that the generated initial values can be very not
#' likely and even may converge slowly to a local maximum. This is the reason
#' why several alternatives exist. The vector of initial values can be directly
#' specified in \code{B} the initial values can be generated (automatically or
#' randomly) from a model with \code{ng=}. Finally, function \code{gridsearch}
#' performs an automatic grid search.
#' 
#' @param fixed a two-sided linear formula object for specifying the
#' fixed-effects in the linear mixed model at the latent process level. The
#' response outcomes are separated by \code{+} on the left of \code{~} and the
#' covariates are separated by \code{+} on the right of the \code{~}. For
#' identifiability purposes, the intercept specified by default should not be
#' removed by a \code{-1}. Variables on which a contrast above the different
#' outcomes should also be estimated are included with \code{contrast()}.
#' @param random an optional one-sided formula for the random-effects in the
#' latent process mixed model. At least one random effect should be included
#' for identifiability purposes. Covariates with a random-effect are separated
#' by \code{+}. By default, an intercept is included. If no intercept,
#' \code{-1} should be the first term included.
#' @param subject name of the covariate representing the grouping structure.
#' @param idiag optional logical for the variance-covariance structure of the
#' random-effects. If \code{FALSE}, a non structured matrix of
#' variance-covariance is considered (by default). If \code{TRUE} a diagonal
#' matrix of variance-covariance is considered.
#' @param randomY optional logical for including an outcome-specific random
#' intercept. If \code{FALSE} no outcome-specific random intercept is added
#' (default). If \code{TRUE} independent outcome-specific random intercepts
#' with parameterized variance are included.
#' @param link optional vector of families of parameterized link functions to
#' estimate (one by outcome). Option "linear" (by default) specifies a linear
#' link function. Other possibilities include "beta" for estimating a link
#' function from the family of Beta cumulative distribution functions and
#' "Splines" for approximating the link function by I-splines. For this latter
#' case, the number of nodes and the nodes location should be also specified.
#' The number of nodes is first entered followed by \code{-}, then the location
#' is specified with "equi", "quant" or "manual" for respectively equidistant
#' nodes, nodes at quantiles of the marker distribution or interior nodes
#' entered manually in argument \code{intnodes}. It is followed by \code{-} and
#' finally "splines" is indicated.  For example, "7-equi-splines" means
#' I-splines with 7 equidistant nodes, "6-quant-splines" means I-splines with 6
#' nodes located at the quantiles of the marker distribution and
#' "9-manual-splines" means I-splines with 9 nodes, the vector of 7 interior
#' nodes being entered in the argument \code{intnodes}.
#' @param intnodes optional vector of interior nodes. This argument is only
#' required for a I-splines link function with nodes entered manually.
#' @param epsY optional definite positive real used to rescale the marker in
#' (0,1) when the beta link function is used. By default, epsY=0.5.
#' @param cor optional indicator for inclusion of an autocorrelated Gaussian
#' process in the latent process linear (latent process) mixed model. Option
#' "BM" indicates a brownian motion with parameterized variance. Option "AR"
#' specifies an autoregressive process of order 1 with parameterized variance
#' and correlation intensity. Each option should be followed by the time
#' variable in brackets as \code{cor=BM(time)}. By default, no autocorrelated
#' Gaussian process is added.
#' @param data data frame containing the variables named in \code{fixed},
#' \code{mixture}, \code{random}, \code{classmb} and \code{subject}.
#' @param B optional specification for the initial values for the parameters.
#' Three options are allowed: (1) a vector of initial values is entered (the
#' order in which the parameters are included is detailed in \code{details}
#' section).  (2) nothing is specified. A preliminary analysis involving the
#' estimation of a standard linear mixed model is performed to choose initial
#' values.  (3) when ng>1, a multlcmm object is entered. It should correspond
#' to the exact same structure of model but with ng=1. The program will
#' automatically generate initial values from this model. This specification
#' avoids the preliminary analysis indicated in (2) Note that due to possible
#' local maxima, the \code{B} vector should be specified and several different
#' starting points should be tried.
#' @param convB optional threshold for the convergence criterion based on the
#' parameter stability. By default, convB=0.0001.
#' @param convL optional threshold for the convergence criterion based on the
#' log-likelihood stability. By default, convL=0.0001.
#' @param convG optional threshold for the convergence criterion based on the
#' derivatives. By default, convG=0.0001.
#' @param maxiter optional maximum number of iterations for the Marquardt
#' iterative algorithm. By default, maxiter=100.
#' @param nsim number of points used to plot the estimated link functions. By
#' default, nsim=100.
#' @param range optional vector indicating the range of the outcomes (that is
#' the minimum and maximum). By default, the range is defined according to the
#' minimum and maximum observed values of the outcome. The option should be
#' used only for Beta and Splines transformations.
#' @param rangeSurv optional vector indicating the range of the survival curve(s).
#' @param subset optional vector giving the subset of observations in
#' \code{data} to use. By default, all lines.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @param posfix Optional vector giving the indices in vector B of the
#' parameters that should not be estimated. Default to NULL, all parameters are
#' estimated.
#' @param partialH optional logical for Beta or Splines link functions only.
#' Indicates whether the parameters of the link functions can be dropped from
#' the Hessian matrix to define convergence criteria.
#' @param verbose logical indicating if information about computation should be
#' reported. Default to TRUE.
#' @param returndata logical indicating if data used for computation should be
#' returned. Default to FALSE, data are not returned.
#' @param methInteg character indicating the type of integration if ordinal outcomes
#' are considered. 'MCO' for ordinary Monte Carlo, 'MCA' for antithetic Monte Carlo,
#' 'QMC' for quasi MonteCarlo.
#' @param nMC integer, number of Monte Carlo simulations
#' @return The list returned is: \item{ns}{number of grouping units in the
#' dataset} \item{ng}{number of latent classes} \item{loglik}{log-likelihood of
#' the model} \item{best}{vector of parameter estimates in the same order as
#' specified in \code{B} and detailed in section \code{details}}
#' \item{V}{vector containing the upper triangle matrix of variance-covariance
#' estimates of \code{Best} with exception for variance-covariance parameters
#' of the random-effects for which \code{V} contains the variance-covariance
#' estimates of the Cholesky transformed parameters displayed in
#' \code{cholesky}} \item{gconv}{vector of convergence criteria: 1. on the
#' parameters, 2. on the likelihood, 3. on the derivatives} \item{conv}{status
#' of convergence: =1 if the convergence criteria were satisfied, =2 if the
#' maximum number of iterations was reached, =4 or 5 if a problem occured
#' during optimisation} \item{call}{the matched call} \item{niter}{number of
#' Marquardt iterations} \item{N}{internal information used in related
#' functions} \item{idiag}{internal information used in related functions}
#' \item{pred}{table of individual predictions and residuals in the underlying
#' latent process scale; it includes marginal predictions (pred_m), marginal
#' residuals (resid_m), subject-specific predictions (pred_ss) and
#' subject-specific residuals (resid_ss) averaged over classes, the transformed
#' observations in the latent process scale (obs) and finally the
#' class-specific marginal and subject-specific predictions (with the number of
#' the latent class: pred_m_1,pred_m_2,...,pred_ss_1,pred_ss_2,...).}
#' \item{pprob}{table of posterior classification and posterior individual
#' class-membership probabilities} \item{Xnames}{list of covariates included in
#' the model} \item{predRE}{table containing individual predictions of the
#' random-effects : a column per random-effect, a line per subject.}
#' \item{cholesky}{vector containing the estimates of the Cholesky transformed
#' parameters of the variance-covariance matrix of the random-effects}
#' \item{estimlink}{table containing the simulated values of each outcome and
#' the corresponding estimated link function} \item{epsY}{definite positive
#' reals used to rescale the markers in (0,1) when the beta link function is
#' used. By default, epsY=0.5.} \item{linktype}{indicators of link function
#' types: 0 for linear, 1 for beta, 2 for splines and 3 for thresholds}
#' \item{linknodes}{vector of nodes useful only for the 'splines' link
#' functions} \item{data}{the original data set (if returndata is TRUE)}
#' %% idea0,idprob0,idg0,idcontr0,idcor0,Xnames2,na.action,pred_RE_Y,Ynames,nbnodes
#' @author Viviane Philipps and Cecile Proust-Lima
#' @seealso
#' 
#' \code{\link{multlcmm}}, \code{\link{Jointlcmm}}
#' 
#' @references
#' Proust-Lima, Amieva, Jacqmin-Gadda (2013). Analysis of multivariate mixed
#' longitudinal data: A flexible latent process approach. Br J Math Stat
#' Psychol 66(3): 470-87.
#' @examples
#'
#' @export
#' 
irt <- function(fixed,random,subject,idiag=FALSE,cor=NULL,link="linear",intnodes=NULL,epsY=0.5,randomY=FALSE,
                survival=NULL,hazard="Weibull",hazardnodes=NULL,TimeDepVar=NULL,logscale=FALSE,startWeibull=0,
                methInteg="QMC",nMC=1000,data,subset=NULL,na.action=1,
                B,posfix=NULL,maxiter=100,convB=0.0001,convL=0.0001,convG=0.0001,partialH=FALSE,
                nsim=100,range=NULL,rangeSurv=NULL,verbose=TRUE,returndata=FALSE)
{
    ptm <- proc.time()
    if(verbose==TRUE) cat("Be patient, irt is running ... \n")

    cl <- match.call()

    nom.subject <- as.character(subject)

    if(missing(random)) stop("At least one random effect is required")
    if(random==~-1) stop("At least one random effect is required")
    if(missing(fixed)) stop("The argument fixed must be specified in any model")
    if(class(fixed)!="formula") stop("The argument fixed must be a formula")
    if(class(random)!="formula") stop("The argument random must be a formula")
    if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
    if(nrow(data)==0) stop("Data should not be empty")
    if(missing(subject)){ stop("The argument subject must be specified")}
    if(!is.numeric(data[,subject])) stop("The argument subject must be numeric")
    if(all(link %in% c("linear","beta","thresholds")) & !is.null(intnodes)) stop("Intnodes should only be specified with splines links")
    if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")

#    if(length(posfix) & missing(B)) stop("A set of initial parameters must be specified if some parameters are not estimated")


    ## garder data tel quel pour le renvoyer
    if(returndata==TRUE)
    {
        datareturn <- data
    }
    else
    {
        datareturn <- NULL
    }
    
    ## test de l'argument cor
    ncor <- 0
    cor.type <- cl$cor[1]
    cor.time <- cl$cor[2]
    cor <- paste(cor.type,cor.time,sep="-")
    if (!isTRUE(all.equal(cor,character(0))))
        {
            if (substr(cor,1,2)=="AR") { ncor <- 2 }
            else if (substr(cor,1,2)=="BM") { ncor <- 1  }
            else { stop("The argument cor must be of type AR or BM") }

            if(!(strsplit(cor,"-")[[1]][2] %in% colnames(data))) stop("Unable to find time variable from argument 'cor' in 'data'")
            else { cor.var.time <- strsplit(cor,"-")[[1]][2] }
        }

    ##pour acces aux attributs des formules
    afixed <- terms(fixed, specials=c("factor","contrast"))
    if(attr(afixed,"intercept")==0) stop("An intercept should appear in fixed for identifiability purposes")
                                        
    arandom <- terms(random, specials=c("factor"))
    ##fixed sans contrast
    fixed2 <- gsub("contrast","",fixed)
    fixed2 <- formula(paste(fixed2[2],fixed2[3],sep="~"))   
    afixed2 <- terms(fixed2)
                                        
    ##verifier si toutes les varialbes sont dans data
    variables <- c(attr(afixed,"variables"),attr(arandom,"variables"))
    variables <- unlist(lapply(variables,all.vars))  
    if(!all(variables %in% colnames(data))) stop(paste("Data should contain the variables",paste(unique(variables),collapse=" ")))


    ##contrast
    contr <- ~-1
    if(!is.null(attr(afixed,"specials")$contrast))
        {
            vcontr <- attr(afixed,"term.labels")[setdiff(attr(afixed,"specials")$contrast-1,untangle.specials(afixed,"contrast",2)$terms)]
            vcontr <- gsub("contrast","",vcontr)
            contr <- as.formula(paste("~-1+",paste(vcontr,collapse="+")))
        }
    acontr <- terms(contr)
    

    ##liste des outcomes
    nomsY <- as.character(attr(afixed,"variables")[2])
    nomsY <- strsplit(nomsY,split=" + ",fixed=TRUE)
    nomsY <- as.vector(nomsY[[1]])
    ny <- length(nomsY)

    ##pas de contrast ni randomY si un seul Y
    if(ny<2 & length(attr(afixed,"specials")$contrast)) stop("No contrast can be included with less than two outcomes")
    if(ny<2 & randomY==TRUE) stop("With less than 2 outcomes randomY should be FALSE")

    ##liste des variables utilisees  (sans les interactions et sans les Y)
    ttesLesVar <- colnames(get_all_vars(afixed,data=data[1,]))
    ttesLesVar <- c(ttesLesVar, colnames(get_all_vars(arandom,data=data[1,])))
    if (ncor>0) ttesLesVar <- unique(c(ttesLesVar,cor.var.time))
    else ttesLesVar <- unique(ttesLesVar)
    ttesLesVar <- setdiff(ttesLesVar, nomsY)

    ## argument subset
    form1 <- paste(c(nom.subject,nomsY,ttesLesVar),collapse="+")
    if(!isTRUE(all.equal(as.character(cl$subset),character(0))))
        {
            cc <- cl
            cc <- cc[c(1,which(names(cl)=="subset"))]
            cc[[1]] <- as.name("model.frame")
            cc$formula <- formula(paste("~",form1))
            cc$data <- data
            cc$na.action <- na.pass
            data <- eval(cc)
        }

    attributes(data)$terms <- NULL

    ## si subject est un factor
    if(is.factor(data[,nom.subject]))
        {
            data[,nom.subject] <- as.numeric(data[,nom.subject])
        }


    ## partie survie
    if(is.null(survival))
    {
        nbevt <- 0
        Tevent <- NULL
        idtrunc <- 0
        nprisq <- 0
        nrisqtot <- 0
        nvarxevt <- 0
        nvarxevt2 <- 0
        typrisq <- 0
        noms.surv <- NULL
        nom.timedepvar <- NULL
        form.commun <- ~-1
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
        
        ## nombre d'evenement concurrents
        nbevt <- length(attr(surv,"states"))
        if(nbevt<1) stop("No observed event in the data")
        
        
        ## pour la formule pour survivial, creer 3 formules : 
        ## une pour les covariables en mixture, une pour les covariables avec effet specifique a la cause, et une pour les effets communs.  
        form.surv <- cl$survival[3]
        
        noms.form.surv <- all.vars(attr(terms(formula(paste("~",form.surv))),"variables"))
        if(length(noms.form.surv)==0)
        {
            form.cause <- ~-1
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
                        
            
            ## creer la formule pour ni cause ni mixture
            asurv <- terms(formula(paste("~",form.surv)),specials=c("cause"))
            ind.commun <- setdiff(1:length(labels(asurv)),unlist(attr(asurv,"specials")))
            if(length(ind.commun))
            {
                form.commun <- paste(labels(asurv)[ind.commun],collapse="+")
                form.commun <- gsub("cause","",form.commun)   # si X1:cause(X2)
                form.commun <- formula(paste("~",form.commun))  
            }
            else
            {
                form.commun <- ~-1 
            }
        }
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

        
    ##verifier si toutes les variables sont dans data
    varSurv <- unique(all.vars(terms(survival)))
    if(!is.null(nom.timedepvar)){if(!(nom.timedepvar %in% all.vars(terms(survival)))) stop("Variable in 'TimeDepVar' should also appear as a covariate in the 'survival' argument")}  
    if(!all(varSurv %in% colnames(data))) stop(paste("Data should contain the variables",paste(varSurv,collapse=" ")))
    ttesLesVar <- unique(c(ttesLesVar,varSurv))
    
    
    ##subset de data avec les variables utilisees
    newdata <- data[,c(nom.subject,nomsY,noms.surv,ttesLesVar),drop=FALSE]


    ## remplacer les NA de TimeDepVar par Tevent
    Tint <- Tevent
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

    dataSurv <- NULL
    if((nbevt>0))
    {
        dataSurv <- data.frame(getElement(object=data,name=nom.subject),Tentry,Tevent,Event,Tint)
    }
    


    
    ##un data frame par outcome et creation Y0
    dataY <- paste("data",nomsY,sep=".")
    Y0 <- NULL
    IND <- NULL
    outcome <- NULL
    data0 <- NULL
    nayk <- vector("list",ny)
    for (k in 1:ny)
    {
        dtemp <- newdata[,c(nom.subject,nomsY[k],noms.surv,ttesLesVar)]
        ##enlever les NA
        linesNA <- apply(dtemp,2,function(v) which(is.na(v)))
        linesNA <- unique(unlist(linesNA))
        if(length(linesNA)) nayk[[k]] <- linesNA
        if(na.action==1 & length(linesNA)>0) dtemp <- dtemp[-linesNA,]
        if(na.action==2 & length(linesNA)>0) stop("Data contains missing values")
        assign(dataY[k],dtemp)
        Y0 <- c(Y0, dtemp[,nomsY[k]])
        IND <- c(IND, dtemp[,nom.subject])
        outcome <- c(outcome,rep(nomsY[k],nrow(dtemp)))
        data0 <- rbind(data0, dtemp[,setdiff(colnames(dtemp),nomsY[k]),drop=FALSE])   #dataset sans NA avec les covariables utilisees; obs ordonnees par outcome
    }

    
    ##creation de X0 (ttes les var + interactions)
    Xfixed <- model.matrix(fixed2[-2], data=data0)
    Xrandom <- model.matrix(random, data=data0)
    Xcontr <- model.matrix(contr,data=data0)
    Xsurv <- model.matrix(form.commun,data=data0)
    Xsurvcause <- model.matrix(form.cause,data=data0)
    
    
    z.fixed <- strsplit(colnames(Xfixed),split=":",fixed=TRUE)
    z.fixed <- lapply(z.fixed,sort)

    z.random <- strsplit(colnames(Xrandom),split=":",fixed=TRUE)
    z.random <- lapply(z.random,sort)
        
    if(contr != ~-1)
        {
            z.contr <- strsplit(colnames(Xcontr),split=":",fixed=TRUE)
            z.contr <- lapply(z.contr,sort)
        }
    else
        {
            z.contr <- list()
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

    if(form.cause != ~-1)
    {
        z.survcause <- strsplit(colnames(Xsurvcause),split=":",fixed=TRUE)
        z.survcause <- lapply(z.survcause,sort)
    }
    else
    {
        z.survcause <- list() 
    }
    
    
    if(!all(z.contr %in% z.fixed))  stop("The covariates in contrast should also appear in fixed")
    
    X0 <- cbind(Xfixed, Xrandom, Xsurv, Xsurvcause)
    
    nom.unique <- unique(colnames(X0))
    X0 <- X0[,nom.unique,drop=FALSE]

    form.cor <- ~-1
    if (ncor>0)
        {
            if(!(cor.var.time %in% colnames(X0)))
                {
                    X0 <- cbind(X0, data0[,cor.var.time])
                    colnames(X0) <- c(nom.unique, cor.var.time)
                    nom.unique <- c(nom.unique,cor.var.time)
                    form.cor <- formula(paste("~-1+",cor.var.time))
                }
        }

    X0 <- as.matrix(X0)
    ##X0 fini


    ##test de link
    if (length(link)!=1 & length(link)!=ny) stop("One link per outcome should be specified")
    if(any(link %in% c("splines","Splines")))
        {
            link[which(link %in% c("splines","Splines"))] <- "5-quant-splines"
        }
    if(length(link)==1 & ny>1)
        {
            link <- rep(link, ny)
        }

    idlink <- rep(2,ny)
    idlink[which(link=="linear")] <- 0
    idlink[which(link=="beta")] <- 1
    idlink[which(link=="thresholds")] <- 3

    spl <- strsplit(link[which(idlink==2)],"-")
    if(any(sapply(spl,length)!=3)) stop("Invalid argument 'link'")

    nySPL <- length(spl)
    nybeta <- sum(idlink==1)
    nyORD <- sum(idlink==3)
    ##remplir range si pas specifie
    if(!is.null(range) & length(range)!=2*(nySPL+nybeta)) stop("Length of vector range is not correct.")
    if((length(range)==2*(nySPL+nybeta)) & (nySPL+nybeta>0))
        {
            ind12 <- which(idlink==1 | idlink==2)
            for (i in 1:(nySPL+nybeta))
                {
                    rg <- range(get(dataY[ind12[i]])[,nomsY[ind12[i]]])
                    if(rg[1]<range[2*(i-1)+1] | rg[2]>range[2*(i-1)+2]) stop("The range specified do not cover the entire range of the data")
                }
        }
    if((is.null(range) & (nybeta+nySPL)>0) | length(range)!=2*(nySPL+nybeta))
        {
            range <- NULL
            for(k in which(idlink!=0))
                {
                    min1 <- min(get(dataY[k])[,nomsY[k]])
                    min2 <- round(min1,3)
                    if(min1<min2) min2 <- min2-0.001

                    max1 <- max(get(dataY[k])[,nomsY[k]])
                    max2 <- round(max1,3)
                    if(max1>max2) max2 <- max2+0.001
                    
                    range <- c(range, min2, max2)
                }
        }


    ## epsY
    if (any(idlink==1))
        {
            if (any(epsY<=0))
                {
                    stop("Argument 'epsY' should be positive.")
                }

            if(length(epsY)==1) epsY <- rep(epsY,nybeta)
            
            if(length(epsY)!=nybeta) stop(paste("Argument 'epsY' should be of length",nybeta))
            if(nybeta!=ny)
                {
                    epsY2 <- rep(0,ny)
                    epsY2[which(idlink==1)] <- epsY
                    epsY <- epsY2
                }

        } 
    

    nbzitr <- rep(2,ny) #nbzitr = nb de noeuds si splines, 2 sinon
    nbnodes <- NULL  #que pour les splines
    spltype <- NULL
    if(nySPL>0)
        {
            for (i in 1:nySPL)
                {
                    nbnodes <- c(nbnodes, spl[[i]][1])
                    spltype <- c(spltype, spl[[i]][2])
                    if(spl[[i]][3] != "splines") stop("Invalid argument link")
                }
        }
    nbnodes <- as.numeric(nbnodes)
    nbzitr[which(idlink==2)] <- nbnodes

    ##test splines
    if(!(all(spltype %in% c("equi","quant","manual")))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")

    ##tester longueur de intnodes
    if(!is.null(intnodes))
        {  
            if(length(intnodes) != sum(nbnodes[which(spltype=="manual")]-2)) stop(paste("Vector intnodes should be of length",sum(nbnodes[which(spltype=="manual")]-2)))
        }

    ##intnodes2 : contient tous les noeuds interieurs (pas seulement ceux de manual)
    intnodes2 <- rep(NA,sum(nbnodes-2))
    nb <- 0
    nbspl <- 0
    for (k in 1:ny)
        {
            if (idlink[k]!=2) next
            else
                {                                         
                    nbspl <- nbspl+1

                    if(spltype[nbspl]=="manual")
                        {
                            nodes <- intnodes[(nb+1):(nb+nbnodes[nbspl]-2)]
                            if(!length(nodes)) stop("The length of intnodes is not correct")
                            intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <-  nodes
                            nb <- nb+nbnodes[nbspl]-2

                            idrg <- length(which(idlink[1:k] != 0))
                            if(any(nodes <= range[2*(idrg-1)+1]) | any(nodes >= range[2*idrg])) stop("Interior nodes must be in the range of the outcome")
                        }

                    if(spltype[nbspl]=="equi")
                        {
                            nodes <- seq(range[2*(nbspl-1)+1], range[2*nbspl], length.out=nbnodes[nbspl])
                            nodes <- nodes[-nbnodes[nbspl]]
                            nodes <- nodes[-1]
                            intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <- nodes
                        }

                    if(spltype[nbspl]=="quant")
                        {
                            nodes <- quantile(get(dataY[k])[,nomsY[k]], probs=seq(0,1,length.out=nbnodes[nbspl]))
                            if(length(unique(nodes)) != length(nodes)) stop(paste("Some nodes are equal for link number",k,"; Please try to reduce the number of nodes or use manual location."))
                            nodes <- nodes[-nbnodes[nbspl]]
                            nodes <- nodes[-1]
                            intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <- as.vector(nodes)
                        }
                }
        }

    if(nb != length(intnodes)) stop(paste("The vector intnodes should be of length",nb))

    ##remplir zitr
    m <- 0
    if(nySPL>0) m <- max(nbnodes)
    zitr <- matrix(0,max(m,2),ny)
    nb12 <- 0
    nbspl <- 0
    for (k in 1:ny)
        {
            if((idlink[k]==0) | (idlink[k]==3)) zitr[1:2,k] <- c(min(get(dataY[k])[,nomsY[k]]),max(get(dataY[k])[,nomsY[k]]))

            if(idlink[k]==1)
                {
                    nb12 <- nb12 + 1
                    zitr[1:2,k] <- range[2*(nb12-1)+1:2]
                }

            if(idlink[k]==2)
                {
                    nb12 <- nb12+1
                    nbspl <- nbspl+1
                    zitr[2:(nbzitr[k]-1),k] <- intnodes2[ifelse(nbspl==1,0,1)*sum(nbnodes[1:(nbspl-1)]-2) + 1:(nbnodes[nbspl]-2)]
                    zitr[1,k] <- range[2*(nb12-1)+1]
                    zitr[nbnodes[nbspl],k]  <- range[2*nb12]
                    
                    ##verifier s'il y a des obs entre les noeuds
                    hcounts <- hist(get(dataY[k])[,nomsY[k]],breaks=zitr[1:nbnodes[nbspl],k],plot=FALSE,include.lowest=TRUE,right=TRUE)$counts
                    if(any(hcounts==0)) stop(paste("Link function number",k,"can not be estimated. Please try other nodes such that there are observations in each interval."))    
                }
        }

    ##uniqueY0 et indiceY0
    uniqueY0 <- NULL
    indiceY0 <- NULL
    nvalSPLORD <- rep(0,ny)
    nbmod <- rep(0,ny)
    modalites <- vector("list",ny)
    nb <- 0
    for (k in 1:ny)
        {
            if((idlink[k]!=2) & (idlink[k]!=3))
                {
                    indiceY0 <- c(indiceY0, rep(0,length(get(dataY[k])[,nomsY[k]])))
                    next
                }

            yk <- get(dataY[k])[,nomsY[k]]
            uniqueTemp <- sort(unique(yk))
            permut <- order(order(yk))  # sort(y)[order(order(y))] = y
            if(length(as.vector(table(yk)))==length(uniqueTemp))
                {
                    indice <- rep(1:length(uniqueTemp), as.vector(table(yk)))
                    if(idlink[k]==2)
                    {
                        indiceTemp <- nb + indice[permut]
                    }
                    else
                    {
                        indiceTemp <- indice[permut]
                    }
                    
                    nb <- nb + length(uniqueTemp)

                    uniqueY0 <- c(uniqueY0, uniqueTemp)
                    indiceY0 <- c(indiceY0, indiceTemp)
                    nvalSPLORD[k] <- length(uniqueTemp)
                }
            else
                {
                    uniqueY0 <- c(uniqueY0, yk)
                    indiceY0 <- c(indiceY0, ifelse(idlink[k]==2,nb,0)+c(1:length(yk)))
                    nb <- nb + length(yk)
                    nvalSPLORD[k] <- length(yk)
                }
            if(idlink[k]==3)
            {
                nbmod[k] <- length(na.omit(uniqueTemp))
                modalites[[k]] <- uniqueTemp
            }
        }


    ##ordonner les mesures par individu
    matYX <- cbind(IND,Y0,indiceY0,outcome,X0)
    matYXord <- matYX[order(IND),]
    Y0 <- as.numeric(matYXord[,2])
    X0 <- apply(matYXord[,-c(1:4),drop=FALSE],2,as.numeric)
    IND <- matYXord[,1]
    outcome <- matYXord[,4]
    indiceY0 <- as.numeric(matYXord[,3])

    dataSurv <- dataSurv[which(dataSurv[,1] %in% IND),]
    dataSurv <- dataSurv[order(dataSurv[,1]),]
    nmes <- as.vector(table(dataSurv[,1]))
    data.surv <- apply(dataSurv[cumsum(nmes),-1],2,as.numeric)
    tsurv0 <- data.surv[,1] 
    tsurv <- data.surv[,2]
    devt <- data.surv[,3]
    tsurvint <- data.surv[,4]
    ind_survint <- (tsurvint<tsurv) + 0 

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
        
    nz <- rep(2,nbevt) 
    locnodes <- NULL  
    typrisq <- rep(2,nbevt)   
    nprisq <- rep(2,nbevt) 
    
    nznodes <- 0 #longueur de hazardnodes
    ii <- 0
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
                nznodes <- nznodes + nz[i]-2
            }
            
            if(!all(locnodes %in% c("equi","quant","manual"))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")      
        }
        
        if(!is.null(hazardnodes))
        {
            if(!any(locnodes == "manual"))  stop("hazardnodes should be NULL if the nodes are not chosen manually")
            
            if(length(hazardnodes) != nznodes) stop(paste("Vector hazardnodes should be of length",nznodes)) 
        }  
    }
    else
    {
        if(!is.null(hazardnodes)) stop("hazardnodes should be NULL if Weibull baseline risk functions are chosen")
    }
    
    
    if(nbevt>1 & length(arghaz)==1 & nznodes>0)
    {
        hazardnodes <- rep(hazardnodes,length.out=nznodes*nbevt)
    }

    nrisqtot <- sum(nprisq)
    
    zi <- matrix(0,nrow=max(nz),ncol=nbevt)
    nb <- 0   
    
    minT1 <- 0
    maxT1 <- max(tsurv)
    tsurvevt <- tsurv
    
    if(idtrunc==1)
    {
        minT1 <- min(tsurv,tsurv0)
        maxT1 <- max(tsurv,tsurv0)
    }
    
    ## arrondir
    minT2 <- round(minT1,3)
    if(minT1<minT2) minT2 <- minT2-0.001
    minT <- minT2
    
    maxT2 <- round(maxT1,3)
    if(maxT1>maxT2) maxT2 <- maxT2+0.001
    maxT <- maxT2
    
    if(length(rangeSurv))
    {
        if(rangeSurv[1] > minT) stop(paste("rangeSurv[1] should be <=", minT))
        if(rangeSurv[2] < maxT) stop(paste("rangeSurv[2] should be >=", maxT))
        minT <- rangeSurv[1]
        maxT <- rangeSurv[2]
    }
    
    
    startWeib <- rep(0,nbevt)
    startWeib[which(typrisq==2)] <- rep(startWeibull, length.out=length(which(typrisq==2)))
    ii <- 0  
    for(i in 1:nbevt)
    {
        if(typrisq[i]==2)
        {
            if(minT < startWeibull[i]) stop("Some entry or event times are bellow startWeibull")
            zi[1:2,i] <- c(startWeib[i],maxT)
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
                pi <- c(1:(nz[i]-2))/(nz[i]-1)
                qi <- quantile(tsurvevt,prob=pi)
                zi[1,i] <- minT
                zi[2:(nz[i]-1),i] <- qi
                zi[nz[i],i] <- maxT
            }
        }   
    }

        

    ##parametres pour Fortran
    ns <- length(unique(IND))
    nv <- dim(X0)[2]
    nobs <- length(Y0)
    idiag0 <- ifelse(idiag==TRUE,1,0)
    nalea <- ifelse(randomY==TRUE,ny,0)
    logspecif <- as.numeric(logscale)
    loglik <- 0
    ni <- 0
    istop <- 0
    gconv <- rep(0,3)
    resid_m <- rep(0,nobs)
    resid_ss <- rep(0,nobs)
    Yobs <- rep(0,nobs)
    time <- seq(minT,maxT,length.out=nsim)
    risq_est <- matrix(0,nrow=nsim,ncol=nbevt)
    risqcum_est <- matrix(0,nrow=nsim,ncol=nbevt)                                
    predRE_Y <- rep(0,ns*nalea)
    rlindiv <- rep(0,ns)
    marker <- rep(0,nsim*ny)
    transfY <- rep(0,nsim*ny)
    
 
    
    ##nmes
    nmes <- matrix(0,ns,ny)
    for (k in 1:ny)
        {
            INDpresents <- which(unique(IND) %in% get(dataY[k])[,nom.subject])
            nmes[INDpresents,k] <- as.vector(table(get(dataY[k])[,nom.subject]))
        }
    maxmes <- max(apply(nmes,1,sum))



    ##remplir idprob, etc
    z.X0 <- strsplit(nom.unique,split=":",fixed=TRUE)
    z.X0 <- lapply(z.X0,sort)
    
    idea <- (z.X0 %in% z.random) + 0
    idg <- (z.X0 %in% z.fixed) + 0
    idcontr <- (z.X0 %in% z.contr) + 0
    
    if (ncor>0) idcor <- colnames(X0) %in% cor.var.time +0
    else idcor <- rep(0,nv)

    idsurv <- z.X0 %in% z.surv + z.X0 %in% z.survcause
    idsurv[1] <- 0 # 0 pour l'intercept

    idtdv <- z.X0 %in% nom.timedepvar + 0

    ## Si pas TimeDepVar dans formule survival
    if(length(nom.timedepvar) & all(idtdv==0))
    {
        stop("Variable in 'TimeDepVar' should also appear as a covariate in the 'survival' argument")
    }

    ## nb coef ppour survie
    nvarxevt <- sum(idsurv==1) + nbevt*sum(idsurv==2)
    
    nea <- sum(idea)
    predRE <- rep(0,ns*nea)

    ## nb parametres d'association
    nasso <- nea*nbevt

    ## prm partie long
    nef <- sum(idg==1)-1 
    ncontr <- (ny-1)*sum(idcontr)
    nvc <- ifelse(idiag0==1,nea,nea*(nea+1)/2)-1
    ntr <- rep(0,ny)
    ntr[which(idlink==0)] <- 2
    ntr[which(idlink==1)] <- 4
    ntr[which(idlink==2)] <- nbzitr[which(idlink==2)]+2
    ntr[which(idlink==3)] <- nbmod[which(idlink==3)]-1
    ntrtot <- sum(ntr)

    ##nombre total de parametres
    NPM <- nrisqtot + nvarxevt + nasso +
        nef + ncontr + nvc + ncor + ntrtot + nalea + ny
        
    
    V <- rep(0, NPM*(NPM+1)/2)  #pr variance des parametres

    ## parametres MC
    methInteg <- switch(methInteg,"MCO"=1,"MCA"=2,"QMC"=3)
    seqMC <- 0
    dimMC <- 0
    if(methInteg==3) 
    {
        dimMC <- nea+nalea
        if(ncor>0) dimMC <- dimMC+maxmes
        if(dimMC>0) seqMC <- randtoolbox::sobol(n=nMC,dim=dimMC,normal=TRUE,scrambling=1) 
    }


    
## gestion de B=random(mod)

        Brandom <- FALSE
        if(length(cl$B)==2)
            {
                if(class(eval(cl$B[[2]]))!="multlcmm") stop("The model specified in B should be of class multlcmm")
                if(as.character(cl$B[1])!="random") stop("Please use random() to specify random initial values")
                
                Brandom <- TRUE
                B <- eval(cl$B[[2]])

                if(length(posfix)) stop("Argument posfix is not compatible with random intial values")
            }
    
###valeurs initiales
    if(!(missing(B)))
        {
            if(!is.vector(B)) stop("B should be a vector")
                
            if (length(B)==NPM) b <- B
            else stop(paste("Vector B should be of length",NPM))
            
        }
    else ## B missing
        {
            b <- rep(0,NPM)
            
            if(nbevt>0)
            { #si Weibull et prms au carre pr positivite -> valeurs par defaut = 1
                if((any(hazard!="Weibull")==FALSE) & (isFALSE(logscale)==TRUE))
                {
                    for(i in 1:nbevt)
                    {
                        if(typrisq[i]==2)
                        {
                            b[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- 1
                        }
                    }
                }
            }
            

            
            if (nvc>0)
                {
                    if(idiag==1) b[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- rep(1,nvc)
                    if(idiag==0)
                        {
                            init.nvc <- diag(nea)
                            init.nvc <- init.nvc[upper.tri(init.nvc, diag=TRUE)]
                            b[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- init.nvc[-1]
                        }
                }
            
            if(ncor>0) b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor] <- 1

            if(nalea>0) b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+1:nalea] <- 1

            b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+1:ny] <-  1

            for(k in 1:ny)
            {
                if(idlink[k]==0)
                {
                    b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-1] <- mean(get(dataY[k])[,nomsY[k]])
                    b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])] <- 1
                }
                if(idlink[k]==1)
                {
                    b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-3] <- 0
                    b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-2] <- -log(2)
                    b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-1] <- 0.7
                    b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])] <- 0.1
                }
                if(idlink[k]==2)
                {
                    b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-ntr[k]+1] <- -2
                    b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-ntr[k]+2:ntr[k]] <- 0.1
                }
                if(idlink[k]==3)
                {
                    b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-ntr[k]+1] <- 0
                    if(ntr[k]>1) b[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:k])-ntr[k]+2:ntr[k]] <- 0.1
                }
            }
        }

    ## ## faire wRandom et b0Random
    ## nef2 <- sum(idg!=0)-1 + (ny-1)*sum(idcontr)
    ## NPM2 <- nef2+ nvc+ncor+nalea+ny+ntrtot
    
    ## wRandom <- rep(0,NPM)
    ## b0Random <- rep(0,ng-1)
    
    ## l <- 0
    ## t <- 0
    ## m <- 0
    ## for (i in 1:nv)
    ## {
    ##     if(idg[i]==1)
    ##     {
    ##         if(i==1) next
    ##         l <- l+1
    ##         t <- t+1
    ##         wRandom[nprob+t] <- l
    ##     }
    ##     if(idg[i]==2)
    ##     {
    ##         if(i==1)
    ##         {
    ##             t <- t+ng-1
    ##             b0Random <- c(b0Random,rep(0,ng-1))
    ##             next
    ##         }
    ##         l <- l+1
    ##         for (g in 1:ng)
    ##         {
    ##             t <- t+1
    ##             wRandom[nprob+t] <- l
    ##         }
    ##     }
    ##     if(idcontr[i]==1)
    ##     {
    ##         wRandom[nef-ncontr+m+1:(ny-1)] <- nef2-ncontr+m+1:(ny-1)
    ##         m <- m+ny-1
    ##     }
    ## }

    ## if(nvc>0)
    ## {
    ##     wRandom[nef+1:nvc] <- nef2+1:nvc
    ## }
    ## if(nw>0)
    ## {
    ##     b0Random <- c(b0Random,rep(1,ng-1))
    ## }
    
    ## if(ncor>0) wRandom[nef+nvc+nw+1:ncor] <- nef2+nvc+1:ncor
    
    ## wRandom[nef+nvc+nw+ncor+1:ny] <- nef2+nvc+ncor+1:ny

    ## if(nalea>0)
    ## {
    ##     wRandom[nef+nvc+nw+ncor+ny+1:nalea] <- nef2+nvc+ncor+ny+1:nalea
    ## }

    ## wRandom[nef+nvc+nw+ncor+ny+nalea+1:ntrtot] <- nef2+nvc+ncor+ny+nalea+1:ntrtot
    ## ## wRandom et b0Random ok.


    ##------------------------------------------
    ##------nom au vecteur best
    ##--------------------------------------------

    nom.X0 <- colnames(X0)
    nom.X0[nom.X0=="(Intercept)"] <- "intercept"

    if(nbevt>0)
    {
        ##prm fct de risque
        if(isTRUE(logscale))
        {
            for(i in 1:nbevt)
            {
                nom1 <- rep(paste("event",i,sep=""),nprisq[i])
                if(typrisq[i]==2)
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:2]," log(Weibull",1:2,")",sep="")
                }
                if(typrisq[i]==1)  
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," log(piecewise",1:(nz[i]-1),")",sep="")  
                }
                if(typrisq[i]==3)  
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," log(splines",1:(nz[i]+2),")",sep="")  
                }  
            }
        }
        else
        {
            for(i in 1:nbevt)
            {
                nom1 <- rep(paste("event",i,sep=""),nprisq[i])
                if(typrisq[i]==2)
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:2]," +/-sqrt(Weibull",1:2,")",sep="")
                }
                if(typrisq[i]==1)  
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(piecewise",1:(nz[i]-1),")",sep="")  
                }
                if(typrisq[i]==3)  
                {
                    names(b)[sum(nprisq[1:i])-nprisq[i]+1:nprisq[i]] <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(splines",1:(nz[i]+2),")",sep="")  
                }  
            }   
        }


        ##prm covariables survival
        nom1 <- NULL  
        for(j in 1:nv)
        {
            if(idsurv[j]==1) #X
            {
                if(idtdv[j]==1)
                {
                    nom1 <- c(nom1,paste("I(t>",nom.timedepvar,")",sep=""))
                }
                else
                {
                    nom1 <- c(nom1,nom.X0[j])
                }
            }
            
            if(idsurv[j]==2) #cause(X)
            {
                if(idtdv[j]==1)
                {
                    nom1 <- c(nom1,paste("I(t>",nom.timedepvar,") event",1:nbevt,sep=""))
                }
                else
                {                  
                    nom1 <- c(nom1,paste(nom.X0[j],paste("event",1:nbevt,sep="")))
                }
            }

        }
        
        if(nvarxevt>0) names(b)[nrisqtot+1:nvarxevt] <- nom1 

        for(i in 1:nbevt)
        {
            names(b)[nrisqtot+nvarxevt+(nbevt-1)*nea+1:nea] <- paste("event",i," asso",1:nea,sep="")
        }
    }
    

    names(b)[nrisqtot+nvarxevt+nasso+1:nef] <- nom.X0[-1][idg[-1]!=0]
    if(ncontr!=0) names(b)[nrisqtot+nvarxevt+nasso+nef+1:ncontr] <- paste("contrast",paste(rep(1:sum(idcontr),each=ny-1),rep(1:(ny-1),sum(idcontr)),sep=""),sep="")
    
    if(idlink[1]==0) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntr[1]]<- c("Linear 1","Linear 2")
    if(idlink[1]==1) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntr[1]]<- paste("Beta",c(1:ntr[1]),sep="")
    if(idlink[1]==2) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntr[1]]<- paste("I-splines",c(1:ntr[1]),sep="")
    if(idlink[1]==3) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+1:ntr[1]]<- paste("Thresh",c(1:ntr[1]),sep="")
    if(ny>1)
        {
            for (yk in 2:ny)
                {
                    if(idlink[yk]==0) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:(yk-1)])+1:ntr[yk]]<- c("Linear 1","Linear 2")
                    if(idlink[yk]==1) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:(yk-1)])+1:ntr[yk]]<- paste("Beta",c(1:ntr[yk]),sep="")
                    if(idlink[yk]==2) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:(yk-1)])+1:ntr[yk]]<- paste("I-splines",c(1:ntr[yk]),sep="")
                    if(idlink[yk]==3) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sum(ntr[1:(yk-1)])+1:ntr[yk]]<- paste("Thresh",c(1:ntr[yk]),sep="")
                }
        }
    if(nvc!=0)names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- paste("varcov",c(1:nvc))

    if(ncor>0) {names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+1:ncor] <- paste("cor",1:ncor,sep="")}
    if(nalea!=0) names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+1:nalea] <- paste("std.randomY",1:ny,sep="")
    
    names(b)[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+1:ny] <- paste("std.err",1:ny)



    ## prm fixes
    fix <- rep(0,NPM)
    if(length(posfix))
        {
            if(any(!(posfix %in% 1:NPM))) stop("Indexes in posfix are not correct")

            fix[posfix] <- 1
        }
    if(length(posfix)==NPM) stop("No parameter to estimate")

    
    ## pour H restreint
    Hr <- as.numeric(partialH)
    pbH <- rep(0,NPM)
    pbH[grep("I-splines",names(b))] <- 1
    pbH[posfix] <- 0
    if(sum(pbH)==0 & Hr==1) stop("No partial Hessian matrix can be defined")

    
#browser()
###estimation
        
    out <- .Fortran(C_irtsre,
                    as.double(Y0),
                    as.double(X0),
                    as.double(tsurv0),
                    as.double(tsurv),
                    as.integer(devt),
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
                    as.integer(ns),
                    as.integer(nv),
                    as.integer(nobs),
                    as.integer(nea),
                    as.integer(nmes),
                    as.integer(idiag0),
                    as.integer(ncor),
                    as.integer(nalea),
                    as.integer(NPM),
                    best=as.double(b),
                    V=as.double(V),
                    loglik=as.double(loglik),
                    niter=as.integer(ni),
                    conv=as.integer(istop),
                    gconv=as.double(gconv),
                    resid_m=as.double(resid_m),
                    resid_ss=as.double(resid_ss),
                    predRE=as.double(predRE),
                    predRE_Y=as.double(predRE_Y),
                    as.double(convB),
                    as.double(convL),
                    as.double(convG),
                    as.integer(maxiter),
                    as.double(epsY),
                    as.integer(idlink),
                    as.integer(nbzitr),
                    as.double(zitr),
                    as.double(uniqueY0),
                    as.integer(indiceY0),
                    as.integer(nvalSPLORD),
                    time=as.double(time),
                    risq_est=as.double(risq_est),
                    risqcum_est=as.double(risqcum_est),
                    marker=as.double(marker),
                    transfY=as.double(transfY),
                    as.integer(nsim),
                    Yobs=as.double(Yobs),
                    rlindiv=as.double(rlindiv),
                    as.integer(pbH),
                    as.integer(fix),
                    as.integer(methInteg),
                    as.integer(nMC),
                    as.integer(dimMC),
                    as.double(seqMC))
#}


    

    ## mettre NA pour les variances et covariances non calculees et  0 pr les prm fixes
    if(length(posfix))
        {
            if(out$conv==3)
                {
                    mr <- NPM-sum(pbH)-length(posfix)
                    Vr <- matrix(0,mr,mr)
                    Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
                    Vr <- t(Vr)
                    Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
                    V <- matrix(NA,NPM,NPM)
                    V[setdiff(1:NPM,c(which(pbH==1),posfix)),setdiff(1:NPM,c(which(pbH==1),posfix))] <- Vr
                    V[,posfix] <- 0
                    V[posfix,] <- 0
                    V <- V[upper.tri(V,diag=TRUE)]
                }
            else
                {
                    mr <- NPM-length(posfix)
                    Vr <- matrix(0,mr,mr)
                    Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
                    Vr <- t(Vr)
                    Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
                    V <- matrix(0,NPM,NPM)
                    V[setdiff(1:NPM,posfix),setdiff(1:NPM,posfix)] <- Vr
                    V <- V[upper.tri(V,diag=TRUE)]
                }
        }
    else
        {
            if(out$conv==3)
                {
                    mr <- NPM-sum(pbH)
                    Vr <- matrix(0,mr,mr)
                    Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
                    Vr <- t(Vr)
                    Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
                    V <- matrix(NA,NPM,NPM)
                    V[setdiff(1:NPM,which(pbH==1)),setdiff(1:NPM,which(pbH==1))] <- Vr
                    V <- V[upper.tri(V,diag=TRUE)]
                }
            else
                {
                    V <- out$V
                }
        }
    

    ## Creation du vecteur cholesky
    Cholesky <- rep(0,(nea*(nea+1)/2))
    if(idiag0==0 & nvc>0)
        {
            Cholesky[1:(nvc+1)] <- c(1,out$best[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc])
            ## Construction de la matrice U
            U <- matrix(0,nrow=nea,ncol=nea)
            U[upper.tri(U,diag=TRUE)] <- Cholesky[1:(nvc+1)]
            z <- t(U) %*% U
            out$best[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- z[upper.tri(z,diag=TRUE)][-1]
        }
    if(idiag0==1 & nvc>0)
        {
            id <- 1:nea
            indice <- rep(id+id*(id-1)/2)
            Cholesky[indice] <- c(1,out$best[nrisqtot+nvarxevt+nef+ncontr+1:nvc])
            out$best[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc] <- out$best[nrisqtot+nvarxevt+nasso+nef+ncontr+1:nvc]**2
        }

    ##predictions
    predRE <- matrix(out$predRE,ncol=nea,byrow=T)
    predRE <- data.frame(unique(IND),predRE)
    colnames(predRE) <- c(nom.subject,nom.X0[idea!=0])

    if (nalea!=0)
        {
            predRE_Y <- matrix(out$predRE_Y,ncol=ny,byrow=TRUE)
            predRE_Y <- data.frame(unique(IND),predRE_Y)
            colnames(predRE_Y)  <- c(nom.subject,nomsY)
        }
    else
        {
            predRE_Y <- rep(NA,nalea*ns)
        }


    ##pred
    pred_m <- out$Yobs-out$resid_m
    pred_ss <- out$Yobs - out$resid_ss
    pred <- data.frame(IND,outcome,pred_m,out$resid_m,pred_ss,out$resid_ss,out$Yobs)

    colnames(pred)<-c(nom.subject,"Yname","pred_m","resid_m","pred_ss","resid_ss","obs")
    rownames(pred) <- NULL

    
    ## risques
    if(nbevt>0)
    {
        risqcum_est <- matrix(out$risqcum_est,nrow=nsim,ncol=nbevt)
        risq_est <- matrix(out$risq_est,nrow=nsim,ncol=nbevt)
        predSurv <- cbind(time,risq_est,risqcum_est)
        
        temp <- paste("event",1:nbevt,".RiskFct",sep="")
        temp1 <- paste("event",1:nbevt,".CumRiskFct",sep="")
        colnames(predSurv) <- c("time",temp,temp1)
        rownames(predSurv) <- 1:nsim
    }
    else
    {
        predSurv <- NA
    }
    
    ##estimlink
    ysim <- matrix(out$marker,nsim,ny)
    transfo <- matrix(out$transfY,nsim,ny)
    estimlink <- as.vector(rbind(ysim,transfo))
    estimlink <- matrix(estimlink,nsim,2*ny)
    colnames(estimlink) <- paste(c("","transf"),rep(nomsY, each=2),sep="")

    if(any(idlink==3))
    {
        if(any(2*nbmod[which(idlink==3)] > nsim))
        {
            nsim2 <- maxval(2*nbmod[which(idlink==3)])

            estimlink2 <- matrix(NA, nrow=nsim2, ncol=2*ny)
            estimlink2[1:nsim,] <- estimlink
            
            colnames(estimlink2) <- colnames(estimlink)
            estimlink <- estimlink2
        }

        seuils <- function(x)
        {
            n <- length(x)
            if(n>1) res <- c(x[1],x[1]+cumsum(x[-1]^2))
            if(n==1) res <- x[1]
            return(res)
        }
        
        sumntr <- 0
        for (k in 1:ny)
        {
            if(idlink[k]==3)
            {
                estimlink[,2*(k-1)+1] <- NA
                estimlink[,2*(k-1)+2] <- NA

                nb <- nbmod[k]-1
                
                marker <- rep(modalites[[k]], each=2)
                seuilsk <- seuils(out$best[nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1:nb])
                transfY <- rep(seuilsk, each=2)
                transfY <- c(-Inf,transfY,Inf)

                estimlink[1:length(marker), 2*(k-1)+1] <- marker
                estimlink[1:length(transfY), 2*(k-1)+2] <- transfY
            }
            
            sumntr <- sumntr + ntr[k]
        }

        ## remplacer -Inf et Inf
        hy <- estimlink[,2*(1:ny)]
        maxhy <- max(hy[is.finite(hy)])
        minhy <- min(hy[is.finite(hy)])
        for (k in 1:ny)
        {
            if(idlink[k]==3)
            {
                estimlink[1, 2*(k-1)+2] <- minhy - (maxhy - minhy)/10
                estimlink[2*nbmod[k], 2*(k-1)+2] <- maxhy + (maxhy - minhy)/10
            }
        }
        
    }

    N <- rep(NA,12)
    N[1] <- 0 #nprob
    N[2] <- nrisqtot
    N[3] <- nvarxevt
    N[4] <- nasso
    N[5] <- nef
    N[6] <- ncontr
    N[7] <- nvc
    N[8] <- 0 #nw
    N[9] <- ncor
    N[10] <- ntrtot
    N[11] <- nalea
    N[12] <- ny
    N[13] <- nobs
    N[14] <- nbevt

    nevent <- rep(0,nbevt)
    for(ke in 1:nbevt)
    {
        nevent[ke] <- length(which(devt==ke))
    }

    Nprm <- c(nprisq,ntr)
    

    nom.X0[nom.X0=="(Intercept)"] <- "Intercept"

    ## noms des variables
    Names <- list(Xnames=nom.X0,Ynames=nomsY,
                  ID=nom.subject,Tnames=noms.surv,
                  TimeDepVar.name=nom.timedepvar,
                  Xvar=setdiff(ttesLesVar,noms.surv))
    
    names(modalites) <- nomsY

    form <- list(fixed=fixed2[-2], random=random, contr=contr,
                 form.commun=form.commun, form.cause=form.cause,
                 form.cor=form.cor)
    
    cost <- proc.time()-ptm
    
    res <-list(ns=ns,idg=idg,idcontr=idcontr,idea=idea,idcor=idcor,
               idsurv=idsurv,idtdv=idtdv,
               loglik=out$loglik,best=out$best,V=V,gconv=out$gconv,conv=out$conv,
               call=cl,niter=out$niter,N=N,nevent=nevent,Nprm=Nprm,
               idiag=idiag0,pred=pred,predRE=predRE,
               predRE_Y=predRE_Y,Names=Names,form=form,cholesky=Cholesky,
               logspecif=logspecif,predSurv=predSurv,typrisq=typrisq,hazardnodes=zi,nz=nz,
               estimlink=estimlink,epsY=epsY,linktype=idlink,linknodes=zitr,nbnodes=nbnodes,nbmod=nbmod,mod=modalites,
               na.action=nayk,AIC=2*(length(out$best)-length(posfix)-out$loglik),BIC=(length(out$best)-length(posfix))*log(ns)-2*out$loglik,data=datareturn,
                                        #wRandom=wRandom,b0Random=b0Random,
               CPUtime=cost[3])
    
    names(res$best) <- names(b)
    class(res) <-c("irt")

    if(verbose==TRUE) cat("The program took", round(cost[3],2), "seconds \n")

    return(res)
}
 
