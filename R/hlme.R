###################### derniere mise a jour : 2012/03/16 ############"




#' Estimation of latent class linear mixed models
#' 
#' This function fits linear mixed models and latent class linear mixed models
#' (LCLMM) also known as growth mixture models or heterogeneous linear mixed
#' models.  The LCLMM consists in assuming that the population is divided in a
#' finite number of latent classes. Each latent class is characterised by a
#' specific trajectory modelled by a class-specific linear mixed model.  Both
#' the latent class membership and the trajectory can be explained according to
#' covariates.  This function is limited to a mixture of Gaussian outcomes. For
#' other types of outcomes, please see function \code{lcmm}. For multivariate
#' longitudinal outcomes, please see \code{multlcmm}.
#' 
#' 
#' A. THE VECTOR OF PARAMETERS B
#' 
#' The parameters in the vector of initial values \code{B} or equivalently in
#' the vector of maximum likelihood estimates \code{best} are included in the
#' following order:
#' 
#' (1) ng-1 parameters are required for intercepts in the latent class
#' membership model, and when covariates are included in \code{classmb}, ng-1
#' paramaters should be entered for each covariate;
#' 
#' (2) for all covariates in \code{fixed}, one parameter is required if the
#' covariate is not in \code{mixture}, ng paramaters are required if the
#' covariate is also in \code{mixture};
#' 
#' (3) the variance of each random-effect specified in \code{random} (including
#' the intercept) when \code{idiag=TRUE}, or the inferior triangular
#' variance-covariance matrix of all the random-effects when
#' \code{idiag=FALSE};
#' 
#' (4) only when \code{nwg=TRUE}, ng-1 parameters are required for the ng-1
#' class-specific proportional coefficients in the variance covariance matrix
#' of the random-effects;
#' 
#' (5) when \code{cor} is specified, 1 parameter corresponding to the variance
#' of the Brownian motion should be entered with \code{cor=BM} and 2 parameters
#' corresponding to the correlation and the variance parameters of the
#' autoregressive process should be entered
#' 
#' (6) the standard error of the residual error.
#' 
#' B. CAUTIONS
#' 
#' Some caution should be made when using the program:
#' 
#' (1) As the log-likelihood of a latent class model can have multiple maxima,
#' a careful choice of the initial values is crucial for ensuring convergence
#' toward the global maximum.  The program can be run without entering the
#' vector of initial values (see point 2).  However, we recommend to
#' systematically enter initial values in \code{B} and try different sets of
#' initial values.
#' 
#' (2) The automatic choice of initial values we provide requires the
#' estimation of a preliminary linear mixed model. The user should be aware
#' that first, this preliminary analysis can take time for large datatsets and
#' second, that the generated initial values can be very not likely and even
#' may converge slowly to a local maximum.  This is the reason why several
#' alternatives exist. The vector of initial values can be directly specified
#' in \code{B} the initial values can be generated (automatically or randomly)
#' from a model with \code{ng=}. Finally, function \code{gridsearch} performs
#' an automatic grid search.
#' 
#' (3) Convergence criteria are very strict as they are based on the
#' derivatives of the log-likelihood in addition to the parameter stability and
#' log-likelihood stability.  In some cases, the program may not converge and
#' reach the maximum number of iterations fixed at 100.  In this case, the user
#' should check that parameter estimates at the last iteration are not on the
#' boundaries of the parameter space.  If the parameters are on the boundaries
#' of the parameter space, the identifiability of the model is critical. This
#' may happen especially with splines parameters that may be too close to 0
#' (lower boundary) or classmb parameters that are too high or low (perfect
#' classification). When identifiability of some parameters is suspected, the
#' program can be run again from the former estimates by fixing the suspected
#' parameters to their value with option posfix. This usually solves the
#' problem. An alternative is to remove the parameters of the Beta of Splines
#' link function from the inverse of the Hessian with option partialH.  If not,
#' the program should be run again with other initial values, with a higher
#' maximum number of iterations or less strict convergence tolerances.
#' 
#' @param fixed two-sided linear formula object for the fixed-effects in the
#' linear mixed model. The response outcome is on the left of \code{~} and the
#' covariates are separated by \code{+} on the right of \code{~}.  By default,
#' an intercept is included. If no intercept, \code{-1} should be the first
#' term included on the right of \code{~}.
#' @param mixture one-sided formula object for the class-specific fixed effects
#' in the linear mixed model (to specify only for a number of latent classes
#' greater than 1).  Among the list of covariates included in \code{fixed}, the
#' covariates with class-specific regression parameters are entered in
#' \code{mixture} separated by \code{+}.  By default, an intercept is included.
#' If no intercept, \code{-1} should be the first term included.
#' @param random optional one-sided formula for the random-effects in the
#' linear mixed model. Covariates with a random-effect are separated by
#' \code{+}.  By default, an intercept is included. If no intercept, \code{-1}
#' should be the first term included.
#' @param subject name of the covariate representing the grouping structure
#' specified with ''.
#' @param classmb optional one-sided formula describing the covariates in the
#' class-membership multinomial logistic model. Covariates included are
#' separated by \code{+}.  No intercept should be included in this formula.
#' @param ng optional number of latent classes considered. If \code{ng=1} (by
#' default) no \code{mixture} nor \code{classmb} should be specified. If
#' \code{ng>1}, \code{mixture} is required.
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
#' @param cor optional brownian motion or autoregressive process modeling the
#' correlation between the observations.  "BM" or "AR" should be specified,
#' followed by the time variable between brackets. By default, no correlation
#' is added.
#' @param data optional data frame containing the variables named in
#' \code{fixed}, \code{mixture}, \code{random}, \code{classmb} and
#' \code{subject}.
#' @param B optional specification for the initial values for the parameters.
#' Three options are allowed: (1) a vector of initial values is entered (the
#' order in which the parameters are included is detailed in \code{details}
#' section).  (2) nothing is specified. A preliminary analysis involving the
#' estimation of a standard linear mixed model is performed to choose initial
#' values.  (3) when ng>1, a hlme object is entered. It should correspond to
#' the exact same structure of model but with ng=1. The program will
#' automatically generate initial values from this model. This specification
#' avoids the preliminary analysis indicated in (2). Note that due to possible
#' local maxima, the \code{B} vector should be specified and several different
#' starting points should be tried.
#' @param convB optional threshold for the convergence criterion based on the
#' parameter stability. By default, convB=0.0001.
#' @param convL optional threshold for the convergence criterion based on the
#' log-likelihood stability. By default, convL=0.0001.
#' @param convG optional threshold for the convergence criterion based on the
#' derivatives. By default, convG=0.0001.
#' @param prior optional name of a covariate containing a prior information
#' about the latent class membership. The covariate should be an integer with
#' values in 0,1,...,ng. Value 0 indicates no prior for the subject while a
#' value in 1,...,ng indicates that the subject belongs to the corresponding
#' latent class.
#' @param maxiter optional maximum number of iterations for the Marquardt
#' iterative algorithm. By default, maxiter=500.
#' @param subset a specification of the rows to be used: defaults to all rows.
#' This can be any valid indexing vector for the rows of data or if that is not
#' supplied, a data frame made up of the variable used in formula.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @param posfix Optional vector specifying the indices in vector B of the
#' parameters that should not be estimated. Default to NULL, all parameters are
#' estimated.
#' @param verbose logical indicating if information about computation should be
#' reported. Default to TRUE.
#' @param returndata logical indicating if data used for computation should be
#' returned. Default to FALSE, data are not returned.
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
#' Marquardt iterations} \item{N}{internal information
#' used in related functions} \item{idiag}{internal information used in related
#' functions} \item{pred}{table of individual predictions and residuals; it
#' includes marginal predictions (pred_m), marginal residuals (resid_m),
#' subject-specific predictions (pred_ss) and subject-specific residuals
#' (resid_ss) averaged over classes, the observation (obs) and finally the
#' class-specific marginal and subject-specific predictions (with the number of
#' the latent class: pred_m_1,pred_m_2,...,pred_ss_1,pred_ss_2,...)}
#' \item{pprob}{table of posterior classification and posterior individual
#' class-membership probabilities} \item{Xnames}{list of covariates included in
#' the model} 
#' \item{predRE}{table containing individual predictions of the random-effects
#' : a column per random-effect, a line per subject} \item{cholesky}{vector
#' containing the estimates of the Cholesky transformed parameters of the
#' variance-covariance matrix of the random-effects}
#' \item{data}{the original data set (if returndata is TRUE)}
#' @author Cecile Proust-Lima, Benoit Liquet and Viviane Philipps
#' 
#' \email{cecile.proust-lima@@inserm.fr}
#' @seealso \code{\link{postprob}}, \code{\link{plot.hlme}},
#' \code{\link{summary}}, \code{\link{predictY}}
#' @references
#' 
#' Proust-Lima C, Philipps V, Liquet B (2017). Estimation of Extended Mixed 
#' Models Using Latent Classes and Latent Processes: The R Package lcmm. 
#' Journal of Statistical Software, 78(2), 1-56. doi:10.18637/jss.v078.i02
#' 
#' Verbeke G and Lesaffre E (1996). A linear mixed-effects model with
#' heterogeneity in the random-effects population. Journal of the American
#' Statistical Association 91, 217-21
#' 
#' Muthen B and Shedden K (1999). Finite mixture modeling with mixture outcomes
#' using the EM algorithm. Biometrics 55, 463-9
#' 
#' Proust C and Jacqmin-Gadda H (2005). Estimation of linear mixed models with
#' a mixture of distribution for the random-effects. Computer Methods Programs
#' Biomedicine 78, 165-73
#' @examples
#' 
#' 
#' ##### Example of a latent class model estimated for a varying number
#' # of latent classes: 
#' # The model includes a subject- (ID) and class-specific linear 
#' # trend (intercept and Time in fixed, random and mixture components)
#' # and a common effect of X1 and its interaction with time over classes 
#' # (in fixed). 
#' # The variance of the random intercept and slope are assumed to be equal 
#' # over classes (nwg=F).
#' # The covariate X3 predicts the class membership (in classmb).
#' #
#' # !CAUTION: initialization of mixed models with latent classes is 
#' # of most importance because of the problem of multimodality of the likelihood.
#' # Calls m2a-m2d illustrate the different implementations for the 
#' # initial values.
#' 
#' ### homogeneous linear mixed model (standard linear mixed model) 
#' ### with correlated random-effects
#' m1<-hlme(Y~Time*X1,random=~Time,subject='ID',ng=1,data=data_hlme)
#' summary(m1)
#' 
#' ### latent class linear mixed model with 2 classes
#' 
#' # a. automatic specification from G=1 model estimates:
#' m2a<-hlme(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,subject='ID',
#'          ng=2,data=data_hlme,B=m1)
#'          
#' # b. vector of initial values provided by the user:
#' m2b<-hlme(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,subject='ID',
#'          ng=2,data=data_hlme,B=c(0.11,-0.74,-0.07,20.71,
#'                                  29.39,-1,0.13,2.45,-0.29,4.5,0.36,0.79,0.97))
#'  
#' # c. random draws from G = 1 model estimates:
#' m2c<-hlme(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,subject='ID',
#'           ng=2,data=data_hlme,B=random(m1))
#' 
#' # d. gridsearch with 50 departures and 10 iterations of the algorithm 
#' #     (see function gridsearch for details)
#' \dontrun{
#' m2d <- gridsearch(rep = 50, maxiter = 10, minit = m1, hlme(Y ~ Time * X1, 
#' mixture =~ Time, random =~ Time, classmb =~ X2 + X3, subject = 'ID', ng = 2, 
#' data = data_hlme))
#' 
#' }  
#'           
#' 
#' 
#' # summary of the estimation process
#' summarytable(m1, m2a, m2b, m2c)
#' 
#' # summary of m2a
#' summary(m2a)
#' 
#' # posterior classification
#' postprob(m2a)
#' 
#' # plot of predicted trajectories using some newdata
#' newdata<-data.frame(Time=seq(0,5,length=100),
#' X1=rep(0,100),X2=rep(0,100),X3=rep(0,100))
#' plot(predictY(m2a,newdata,var.time="Time"),legend.loc="right",bty="l")
#' 
#' 
#' 
#' @export
#' 
#' 
#' 
hlme <-
    function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,cor=NULL,data,B,convB=0.0001,convL=0.0001,convG=0.0001,prior,maxiter=500,subset=NULL,na.action=1,posfix=NULL,verbose=TRUE,returndata=FALSE){

        ptm<-proc.time()
        if(verbose==TRUE) cat("Be patient, hlme is running ... \n")

        cl <- match.call()
        args <- as.list(match.call(hlme))[-1]

        nom.subject <- as.character(subject)
#### INCLUSION PRIOR
        nom.prior <- as.character(args$prior)
####
        if(!missing(mixture) & ng==1) stop("No mixture can be specified with ng=1")
        if(missing(mixture) & ng>1) stop("The argument mixture has to be specified for ng > 1")
        if(!missing(classmb) & ng==1) stop("No classmb can be specified with ng=1")
        if(missing(random)) random <- ~-1
        if(missing(fixed)) stop("The argument Fixed must be specified in any model")
        if(missing(classmb)) classmb <- ~-1
        if(missing(mixture)) mixture <- ~-1
        if(ng==1&nwg==TRUE) stop ("The argument nwg should be FALSE for ng=1")


        if(class(fixed)!="formula") stop("The argument fixed must be a formula")
        if(class(mixture)!="formula") stop("The argument mixture must be a formula")
        if(class(random)!="formula") stop("The argument random must be a formula")
        if(class(classmb)!="formula") stop("The argument classmb must be a formula")
        if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
        if(nrow(data)==0) stop("Data should not be empty") 
        if(missing(subject)){ stop("The argument subject must be specified in any model even without random-effects")}
        if(!is.numeric(data[,subject])) stop("The argument subject must be numeric")

        if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

        if(na.action==1){
            na.action=na.omit
        }else{
            na.action=na.fail
        }


        ## garder data tel quel pour le renvoyer
        if(returndata==TRUE)
        {
            datareturn <- data
        }
        else
        {
            datareturn <- NULL
        }
        
### test de l'argument cor
        ncor0 <- 0
        cor.type <- cl$cor[1]
        cor.time <- cl$cor[2] 
        cor <- paste(cor.type,cor.time,sep="-")
        if (all.equal(cor,character(0))!=TRUE)
            {
                if (substr(cor,1,2)=="AR") { ncor0 <- 2 }
                else if (substr(cor,1,2)=="BM") { ncor0 <- 1  }
                else { stop("The argument cor must be of type AR or BM") }
                 
                if(!(strsplit(cor,"-")[[1]][2] %in% colnames(data))) stop("Unable to find time variable from argument cor in data")
                else { cor.var.time <- strsplit(cor,"-")[[1]][2] }
            }  
### fin test argument cor 



### ad 2/04/2012
        X0.names2 <- c("intercept")
### ad 
        int.fixed <- 0
        int.mixture <- 0
        int.random <- 0
        int.classmb <- 0
                                        #7/05/2012
### Traitement des donnees manquantes
                                        # fixed
        m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]  
        m$formula <- terms(fixed)
        m$na.action=na.action 
        m[[1]] <- as.name("model.frame")	 
        m <- eval(m, sys.parent())      
        na.fixed <- attr(m,"na.action") 

                                        # mixture
        if(mixture[[2]] != "-1"){
            m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
            m$formula <- terms(mixture)
            m$na.action <- na.action
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.mixture <- attr(m,"na.action")	
        }else{
            na.mixture <- NULL
        }

                                        # random
        if(random[[2]] != "-1"){
            m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
            m$formula <- terms(random)
            m$na.action <- na.action
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.random <- attr(m,"na.action")
        }else{
            na.random <- NULL
        }

                                        # classmb
        if(classmb[[2]] != "-1"){ 
            m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]	
            m$formula <- terms(classmb)
            m$na.action <- na.action
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.classmb <- attr(m,"na.action")
        }else{
            na.classmb <- NULL
        }
        
                                        #cor     
        if(ncor0!=0)
            {
                m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
                m$formula <- as.formula(paste(cor.var.time,1,sep="~"))
                m$na.action <- na.action
                m[[1]] <- as.name("model.frame")
                m <- eval(m,sys.parent())    
                na.cor <- attr(m,"na.action") 	
            }
        else { na.cor <- NULL }

### names of covariate in intial fit     (sans les interactions)
        X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(fixed)),data=data))[-1]))
        if(mixture[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(mixture)),data=data))))
        if(random[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(mtemp <- get_all_vars(formula(terms(random)),data=data))))
                                        #7/05/2012
        if(classmb[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(classmb)),data=data))))
                                        #7/05/2012




        ## Table sans donnees manquante: newdata
	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.cor))
        ## dans na.action, on a les indices des NA dans le subset de data

                                        #prendre le subset :
        newdata <- data  
        if(!isTRUE(all.equal(as.character(cl$subset),character(0))))
            {
                cc <- cl
                cc <- cc[c(1,which(names(cl)=="subset"))]
                cc[[1]] <- as.name("model.frame")
                cc$formula <- formula(paste("~",paste(colnames(data),collapse="+")))
                cc$data <- data
                cc$na.action <- na.pass
                newdata <- eval(cc)
            }

                                        #enlever les NA
	if(!is.null(na.action)){
            newdata <- newdata[-na.action,]
	}

        attributes(newdata)$terms <- NULL

        ## Construction de nouvelles var explicatives sur la nouvelle table
        ## fixed
	X_fixed <- model.matrix(fixed,data=newdata)
	if(any(colnames(X_fixed)=="(Intercept)")){
            ii <- which(colnames(X_fixed)=="(Intercept)")
            colnames(X_fixed)[ii] <- "intercept"
            int.fixed <- 1
	}
	nom.fixed <- inddepvar.fixed <- inddepvar.fixed.nom <- colnames(X_fixed)
	if(int.fixed>0)inddepvar.fixed <- inddepvar.fixed[-ii]

        ## mixture
	if(mixture[[2]] != "-1"){
            X_mixture <- model.matrix(mixture,data=newdata)	
            if(any(colnames(X_mixture)=="(Intercept)")){
                ii <- which(colnames(X_mixture)=="(Intercept)")
                colnames(X_mixture)[ii] <- "intercept"
                int.mixture <- 1
            }
            nom.mixture <- inddepvar.mixture <- inddepvar.mixture.nom <- colnames(X_mixture)
            if(int.mixture>0)inddepvar.mixture <- inddepvar.mixture[-ii]
            id.X_mixture <- 1
	}else{
            inddepvar.mixture <- nom.mixture <- inddepvar.mixture.nom <- NULL
            id.X_mixture <- 0
	}
        ## random
	if(random[[2]] != "-1"){
            X_random <- model.matrix(random,data=newdata)	
            if(any(colnames(X_random)=="(Intercept)")){
                ii <- which(colnames(X_random)=="(Intercept)")
                colnames(X_random)[ii] <- "intercept"
                int.random <- 1
            }
            inddepvar.random <- inddepvar.random.nom <- colnames(X_random)
            if(int.random>0) inddepvar.random <- inddepvar.random[-ii]
            id.X_random <- 1
	}else{
            ## ad: add inddepvar.random.nom2 <- NULL 10/04/2012
            inddepvar.random <- inddepvar.random.nom <- NULL
            id.X_random <- 0
	}
        ## classmb
	if(classmb[[2]] != "-1"){ 
            if(attr(terms(classmb),"intercept")==0)
                {
                    classmb <- paste("~",classmb[2],"+1") 
                }
            X_classmb <- model.matrix(as.formula(classmb),data=newdata)
            if(any(colnames(X_classmb)=="(Intercept)")){
                ii <- which(colnames(X_classmb)=="(Intercept)")
                colnames(X_classmb)[ii] <- "intercept"
                int.classmb <- 1
            }
            id.X_classmb <- 1
            if(int.classmb>0) inddepvar.classmb <- colnames(X_classmb)[-ii]
            inddepvar.classmb.nom <- colnames(X_classmb)
	}
        else{
            inddepvar.classmb <- inddepvar.classmb.nom <- "intercept"
            id.X_classmb <- 0
	}	
        
        
                                        #7/05/2012
##############   COVARIATES       ##########################
                                        # intercept is always in inddepvar.classmb
        var.exp <- NULL
        var.exp <- c(var.exp,colnames(X_fixed))
        if(id.X_mixture == 1) var.exp <- c(var.exp,colnames(X_mixture))
        if(id.X_random == 1)var.exp <- c(var.exp,colnames(X_random))
        if(id.X_classmb == 1)var.exp <- c(var.exp,colnames(X_classmb))
        var.exp <- unique(var.exp)    
        if(ncor0>0) 
            { if(!(cor.var.time %in% var.exp)) 
                  {var.exp <- c(var.exp, cor.var.time)} #si la varaible de temps dans cor n'est dan sles variables expl, on l'ajoute
          }             
        
                                        #if(!(all(nom.mixture %in% nom.fixed))) stop("The covariates in mixture should be also included in the argument fixed")
                                        # controler si les variables de mixture sont toutes dans fixed : 
        z.fixed <- strsplit(nom.fixed,split=":",fixed=TRUE)
        z.fixed <- lapply(z.fixed,sort)
        
        if(id.X_mixture==1)
            {
                z.mixture <- strsplit(nom.mixture,split=":",fixed=TRUE)
                z.mixture <- lapply(z.mixture,sort)
            }
        else z.mixture <- list()

        if(!all(z.mixture %in% z.fixed))  stop("The covariates in mixture should also be included in the argument fixed")


        ## var dependante
        Y.name <- as.character(attributes(terms(fixed))$variables[2])
        Y0 <- newdata[,Y.name]

        ## var expli

        X0 <- X_fixed
        oldnames <- colnames(X0)

        z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
        z.X0 <- lapply(z.X0,sort)


        if(id.X_mixture == 1)
            {
                z.mixture <- strsplit(colnames(X_mixture),split=":",fixed=TRUE)
                z.mixture <- lapply(z.mixture,sort)
                for(i in 1:length(colnames(X_mixture)))
                    {
                                        #if((colnames(X_mixture)[i] %in% colnames(X0))==F){ 
                        if(!isTRUE(z.mixture[i] %in% z.X0))
                            {
                                X0 <- cbind(X0,X_mixture[,i])
                                colnames(X0) <- c(oldnames, colnames(X_mixture)[i])
                                oldnames <- colnames(X0)
                                
                                z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
                                z.X0 <- lapply(z.X0,sort)			
                            }
                    }
            }
        else
            {
                z.mixture <- list()
            }

        if(id.X_random == 1)
            {
                z.random <- strsplit(colnames(X_random),split=":",fixed=TRUE)
                z.random <- lapply(z.random,sort)
                for(i in 1:length(colnames(X_random)))
                    {
                                        #if((colnames(X_random)[i] %in% colnames(X0))==F){
                        if(!isTRUE(z.random[i] %in% z.X0))
                            {		
                                X0 <- cbind(X0,X_random[,i])
                                colnames(X0) <- c(oldnames, colnames(X_random)[i])
                                oldnames <- colnames(X0)
                                
                                z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
                                z.X0 <- lapply(z.X0,sort)				
                            }	                                                            
                    }
            }
        else
            {
                z.random <- list()
            }

        if(id.X_classmb == 1)
            {
                z.classmb <- strsplit(colnames(X_classmb),split=":",fixed=TRUE)
                z.classmb <- lapply(z.classmb,sort)
                for(i in 1:length(colnames(X_classmb)))
                    {
                                        #		if((colnames(X_classmb)[i] %in% colnames(X0))==F){
                        if(!isTRUE(z.classmb[i] %in% z.X0))
                            {			
                                X0 <- cbind(X0,X_classmb[,i])
                                colnames(X0) <- c(oldnames, colnames(X_classmb)[i])
                                oldnames <- colnames(X0)
                                
                                z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
                                z.X0 <- lapply(z.X0,sort)	         	 
                            }	
                    }
            }
        else
            {
                z.classmb <- list()
            }

        if(ncor0>0) 
            { 
                if(!(cor.var.time %in% colnames(X0)))   #cor.var.time jamais ne interaction donc ok (pas de z.cor)
                    {
                        X0 <- cbind(X0, newdata[,cor.var.time])
                        colnames(X0) <- c(oldnames, cor.var.time)
                    }
            }  

                                        #colnames(X0) <- var.exp # a remettre si on enleve les z.fixed etc

                                        #X0 <- X0[,-which(colnames(X0)=="intercept")]
### ad

        if((any(is.na(X0))==TRUE)|(any(is.na(Y0))==TRUE))stop("The data should not contain any missing value")
                                        # 
        n <- dim(data)[1]
                                        #if ((int.fixed+int.random)>0) X0<- cbind(intercept=rep(1,n),X0)
                                        #if (!((int.fixed+int.random)>0) & ("intercept" %in% colnames(X0))) X0 <- as.data.frame(X0[,-which(colnames(X0)=="intercept")])
        nom.X0 <- colnames(X0)
        nvar.exp <- length(nom.X0)

        IND <- newdata[,nom.subject]

        #IDnum <- as.numeric(IND)


#### INCLUSION PRIOR 
        #if(missing(prior)){ PRIOR <- seq(0,length=length(IDnum))} 
        if(missing(prior)){ PRIOR <- seq(0,length=length(IND))} 
        if(!missing(prior)){ 
            PRIOR <- newdata[,nom.prior]
            PRIOR[(is.na(PRIOR))] <- 0
        }
####

        ng0 <- ng
        idiag0 <- as.integer(idiag)
        nwg0 <- as.integer(nwg)

        idea0 <- rep(0,nvar.exp)
        idprob0 <- rep(0,nvar.exp)
        idg0 <- rep(0,nvar.exp)
        
        z.X0 <- strsplit(nom.X0,split=":",fixed=TRUE)
        z.X0 <- lapply(z.X0,sort)

        
        for (i in 1:nvar.exp){  
                                        #idea0[i] <- nom.X0[i]%in%inddepvar.random.nom
                                        #idprob0[i] <- nom.X0[i]%in%inddepvar.classmb.nom
                                        #if(nom.X0[i]%in%nom.fixed & !(nom.X0[i]%in%nom.mixture)) idg0[i] <- 1 
                                        #if(nom.X0[i]%in%nom.fixed & nom.X0[i]%in%nom.mixture) idg0[i] <- 2 

            idea0[i] <- z.X0[i] %in% z.random
            idprob0[i] <- z.X0[i] %in% z.classmb   
            if((z.X0[i] %in% z.fixed) & !(z.X0[i] %in% z.mixture)) idg0[i] <- 1 
            if((z.X0[i] %in% z.fixed) & (z.X0[i] %in% z.mixture)) idg0[i] <- 2 
        }

        idcor0 <- rep(0,length(z.X0))  
        if (ncor0!=0) idcor0 <- as.numeric(nom.X0 %in% cor.var.time) 

                                        #if((int.fixed+int.random)>0) idprob0[1] <- 0     
        if(("intercept" %in% nom.X0) | ("(Intercept)" %in% nom.X0))
            {
                ii <- which(nom.X0 %in% c("intercept","(Intercept)"))
                idprob0[ii] <- 0  
            }


                                        # on ordonne les donn es suivants la variable IND
        #matYX <- cbind(IDnum,IND,PRIOR,Y0,X0)
        #matYXord <- matYX[sort.list(matYX[,1]),]
        #Y0 <- matYXord[,4]  
        #X0 <- matYXord[,-c(1,2,3,4)]
        #IDnum <- matYXord[,1]
        #IND <- matYXord[,2]

        matYX <- cbind(IND,PRIOR,Y0,X0)
        matYXord <- matYX[sort.list(matYX[,1]),]
        Y0 <- as.numeric(matYXord[,3])
        X0 <- apply(matYXord[,-c(1,2,3),drop=FALSE],2,as.numeric)
        IND <- matYXord[,1]
        
#### INCLUSION PRIOR 
        PRIOR <- as.numeric(matYXord[,2])
        PRIOR <-as.integer(as.vector(PRIOR))
####

        X0<-as.numeric(as.matrix(X0))
        Y0<-as.numeric(as.matrix(Y0))
        #nmes0<-as.vector(table(IDnum))
        nmes0<-as.vector(table(IND))
        ns0<-length(nmes0)



##### INCLUSION PRIOR 
                                        # definition de prior a 0 pour l'analyse G=1
        prior2 <- as.integer(rep(0,ns0))
        prior0 <- prior2 
                                        # si prior pas missing alors mettre dedans la classe a priori. Attention tester q les valeurs sont dans 0, G
        if(!missing(prior)){ 
            prior0 <- PRIOR[cumsum(nmes0)]
        }
        INDuniq <- IND[cumsum(nmes0)]
        seqnG <- 0:ng0
        if (!(all(prior0  %in% seqnG))) stop ("The argument prior should contain integers between 0 and ng")
#####


        loglik <- as.double(0)
        ni <- 0
        istop <- 0
        gconv <-rep(0,3)
        ppi0 <- rep(0,ns0*ng0)
        nv0<-nvar.exp
        nobs0<-length(Y0)
        resid_m <- rep(0,nobs0)
        resid_ss <- rep(0,nobs0)
        pred_m_g <- rep(0,nobs0*ng0)
        pred_ss_g <- rep(0,nobs0*ng0)
        nea0 <- sum(idea0==1)
        predRE <- rep(0,nea0*ns0)

        ##---------------------------------------------------------------------------
        ##definition du vecteur de parametres + initialisation
        ##---------------------------------------------------------------------------

## gestion de B=random(mod)

        Brandom <- FALSE
        if(length(cl$B)==2)
            {
                if(class(eval(cl$B[[2]]))!="hlme") stop("The model specified in B should be of class hlme")
                if(as.character(cl$B[1])!="random") stop("Please use random() to specify random initial values")
                
                Brandom <- TRUE
                B <- eval(cl$B[[2]])

                if(length(posfix)) stop("Argument posfix is not compatible with random intial values")
            }

        

#####cas 1 : ng=1
        b <- NULL
        b1 <- NULL
        NPROB <- 0
        if(ng0==1| missing(B)){
            NEF<-sum(idg0!=0)  
            b1[1:NEF] <- 0
            if(int.fixed > 0)  b1[1] <- mean(Y0)

            if(idiag0==1){
                NVC<-sum(idea0==1)
                b1[(NEF+1):(NEF+NVC)] <- 1}

            if(idiag0==0){
                kk<-sum(idea0==1) 
                NVC<-(kk*(kk+1))/2
                indice<-cumsum(1:kk)
                bidiag<-rep(0,NVC)
                bidiag[indice] <- 1
                b1[(NEF+1):(NEF+NVC)] <- bidiag
            }
            if(ncor0==1)
                {b1[NEF+NVC+1] <- 1 }
            if(ncor0==2)
                {b1[(NEF+NVC+1):(NEF+NVC+ncor0)] <- c(0,1) }
            b1[NEF+NVC+ncor0+1] <- 1
            NPM <- length(b1)
            NW <- 0
            V <- rep(0,NPM*(NPM+1)/2) 
        }

#####cas 2 : ng>=2
        if(ng0>1){
            NPROB <- (sum(idprob0==1)+1)*(ng0-1)
            b[1:NPROB] <- 0
            NEF <- sum(idg0==1)+(sum(idg0==2))*ng0
            if(idiag0==1) NVC <- sum(idea0==1)
            if(idiag0==0){
                kk <- sum(idea0==1) 
                NVC <- (kk*(kk+1))/2}
            NW <- nwg0*(ng0-1)
            if(NW>0) b[(NPROB+NEF+NVC+1):(NPROB+NEF+NVC+NW)] <- 1   
            if(ncor0==1)
                {b[NEF+NVC+NW+1] <- 1 }
            if(ncor0==2)
                {b[(NPROB+NEF+NVC+NW+1):(NPROB+NEF+NVC+NW+ncor0)] <- c(0,1) }
            NPM <- NPROB+NEF+NVC+NW+ncor0+1
            V <- rep(0,NPM*(NPM+1)/2)
        }

        
        ## prm fixes
        fix0 <- rep(0,NPM)
        if(length(posfix))
            {
                if(any(!(posfix %in% 1:NPM))) stop("Indexes in posfix are not correct")
                
                fix0[posfix] <- 1
            }
        if(length(posfix)==NPM) stop("No parameter to estimate")
        
        
        if(missing(B)){

            if(ng0>1){
                idea2 <- idea0
                idprob2 <- rep(0,nv0)  
                idg2 <- rep(0,nv0) 
                idg2[idg0!=0] <- 1
                NEF2<-sum(idg2==1)
                NPM2<-NEF2+NVC+ncor0+1
                nwg2<-0
                ng2<-1
                ppi2<- rep(0,ns0)
                pred_m_g2 <- rep(0,nobs0)
                pred_ss_g2 <- rep(0,nobs0)
                maxiter2 <- min(75,maxiter)
                convB2 <- max(0.01,convB)
                convL2 <- max(0.01,convL)
                convG2 <- max(0.01,convG)

                V2 <- rep(0,NPM2*(NPM2+1)/2)
                best <- rep(0,NPM2)
                
                init <- .Fortran(C_hetmixlin,
                                 as.double(Y0),
                                 as.double(X0),
                                 as.integer(prior2),
                                 as.integer(idprob2),
                                 as.integer(idea2),
                                 as.integer(idg2),
                                 as.integer(idcor0),
                                 as.integer(ns0),
                                 as.integer(ng2),
                                 as.integer(nv0),
                                 as.integer(nobs0),
                                 as.integer(nea0),
                                 as.integer(nmes0),
                                 as.integer(idiag0),
                                 as.integer(nwg2),
                                 as.integer(ncor0),
                                 npm=as.integer(NPM2),
                                 best=as.double(b1),
                                 V=as.double(V2),
                                 loglik=as.double(loglik),
                                 niter=as.integer(ni),
                                 conv=as.integer(istop),
                                 gconv=as.double(gconv),
                                 ppi2=as.double(ppi2),
                                 resid_m=as.double(resid_m),
                                 resid_ss=as.double(resid_ss),
                                 pred_m_g=as.double(pred_m_g2),
                                 pred_ss_g=as.double(pred_ss_g2),
                                 predRE=as.double(predRE),
                                 as.double(convB2),
                                 as.double(convL2),
                                 as.double(convG2),
                                 as.integer(maxiter2),
                                 as.integer(fix0))

                k <- NPROB
                l <- 0
                t<- 0
                for (i in 1:nvar.exp)    {
                    if(idg0[i]==1){
                        l <- l+1
                        t <- t+1
                        b[k+t] <- init$best[l]
                    }
                    if(idg0[i]==2){
                        l <- l+1
                        for (g in 1:ng){
                            t <- t+1
                            if(init$conv==1) b[k+t] <- init$best[l]+(g-(ng+1)/2)*sqrt(init$V[l*(l+1)/2])
                            else b[k+t] <- init$best[l]+(g-(ng+1)/2)*init$best[l]
                        }
                    }
                }
                b[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- init$best[(NEF2+1):(NEF2+NVC)]
                if (ncor0>0) {b[(NPROB+NEF+NVC+NW+1):(NPROB+NEF+NVC+NW+ncor0)] <- init$best[(NPM2-ncor0):(NPM2-1)]}
                b[NPROB+NEF+NVC+NW+ncor0+1] <- init$best[NPM2]
            } 
            if(ng0==1 ){
                b <- b1
            }
        }
        else
            {
                if(is.vector(B))
                    {
                        if(length(B)!=NPM) stop(paste("Vector B should be of length",NPM))
                        else {b <-B}
                    }
                else
                    { 
                        if(class(B)!="hlme") stop("B should be either a vector or an object of class hlme")

                        ## B est le meme modele mais pr ng=1 :
                        if(ng>1 & B$ng==1)
                            {
                                NEF2 <- sum(idg0!=0)
                                NPM2 <- NEF2+NVC+ncor0+1
                                if(length(B$best)!=NPM2) stop("B is not correct")


                                if(Brandom==FALSE)
                                    {
                                        ## B deterministe
                                        l <- 0
                                        t <- 0
                                        for (i in 1:nvar.exp)
                                            {
                                                if(idg0[i]==1)
                                                    {
                                                        l <- l+1
                                                        t <- t+1
                                                        b[NPROB+t] <- B$best[l]
                                                    }
                                                if(idg0[i]==2)
                                                    {
                                                        l <- l+1
                                                        for (g in 1:ng)
                                                            {
                                                                t <- t+1
                                                                if(B$conv==1) b[NPROB+t] <- B$best[l]+(g-(ng+1)/2)*sqrt(B$V[l*(l+1)/2])
                                                                else b[NPROB+t] <- B$best[l]+(g-(ng+1)/2)*B$best[l]
                                                            }
                                                    }
                                            }
                                        if(NVC>0)
                                            {
                                                if(idiag==TRUE)
                                                    {
                                                        b[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- B$cholesky[(1:nea0)*(2:(nea0+1))/2]
                                                        
                                                    }
                                                else
                                                    {
                                                        b[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- B$cholesky
                                                    }
                                            }
                                        if (ncor0>0) {b[(NPROB+NEF+NVC+NW+1):(NPROB+NEF+NVC+NW+ncor0)] <- B$best[(NPM2-ncor0):(NPM2-1)]}
                                        b[NPROB+NEF+NVC+NW+ncor0+1] <- B$best[NPM2]
                                    }
                                else
                                    {
                                        ## B random
                                        bb <- rep(0,NPM-NPROB-NW)
                                        vbb <- matrix(0,NPM-NPROB-NW,NPM-NPROB-NW)
                                        
                                        VB <- matrix(0,NPM2,NPM2)
                                        VB[upper.tri(VB,diag=TRUE)] <- B$V
                                        VB <- t(VB)
                                        VB[upper.tri(VB,diag=TRUE)] <- B$V

                                        nbg <- idg0[which(idg0!=0)]
                                        nbg[which(nbg==2)] <- ng
                                        nbgnef <- unlist(sapply(nbg,function(k) if(k>1) rep(2,k) else k))
                                        
                                        vbb[which(nbgnef==1),setdiff(1:ncol(vbb),which(nbgnef!=1))] <- VB[which(nbg==1),setdiff(1:ncol(VB),which(nbg!=1))]
                                        vbb[(NEF+1):nrow(vbb),(NEF+1):ncol(vbb)] <- VB[(NEF2+1):nrow(VB),(NEF2+1):ncol(VB)]
                                        
                                        l <- 0
                                        t <- 0
                                        for (i in 1:nvar.exp)
                                            {
                                                if(idg0[i]==1)
                                                    {
                                                        l <- l+1
                                                        t <- t+1
                                                        bb[t] <- B$best[l]
                                                    }
                                                if(idg0[i]==2)
                                                    {
                                                        l <- l+1
                                                        for (g in 1:ng)
                                                            {
                                                                t <- t+1
                                                                bb[t] <- B$best[l]
                                                                vbb[t,t] <- VB[l,l]
                                                            }
                                                    }
                                            }
                                        
                                        if(NVC>0)
                                            {
                                                if(idiag==TRUE)
                                                    {
                                                        bb[NEF+1:NVC] <- B$cholesky[(1:nea0)*(2:(nea0+1))/2]
                                                    }
                                                else
                                                    {
                                                        bb[NEF+1:NVC] <- B$cholesky
                                                    }
                                            }
                            
                                
                                        if (ncor0>0)
                                            {
                                                bb[NEF+NVC+1:ncor0] <- B$best[(NPM2-ncor0):(NPM2-1)]
                                            }
                                
                                        bb[NEF+NVC+ncor0+1] <- B$best[NPM2]
                                                            
                                        up <- vbb[upper.tri(vbb,diag=TRUE)]
                                        vbb <- t(vbb)
                                        vbb[upper.tri(vbb,diag=TRUE)] <- up
                                        Chol <- chol(vbb)
                                        Chol <- t(Chol)
                                        
                                        b[c((NPROB+1):(NPROB+NEF+NVC),(NPROB+NEF+NVC+NW+1):NPM)] <- bb + Chol %*% rnorm(NPM-NPROB-NW)

                                        b[1:NPROB] <- 0
                                        if(NW>0) b[NPROB+NEF+NVC+1:NW] <- 1
                                   
                                    } 
                            }
                    }
            }
                
        ## faire wRandom et b0Random
        NEF2 <- sum(idg0!=0)
        NPM2 <- NEF2+NVC+ncor0+1

        wRandom <- rep(0,NPM)
        b0Random <- rep(0,ng-1) # nprob 
        
        l <- 0
        t <- 0
        for (i in 1:nvar.exp)
        {
            if(idg0[i]==1)
            {
                l <- l+1
                t <- t+1
                wRandom[NPROB+t] <- l
            }
            if(idg0[i]==2)
            {
                l <- l+1
                for (g in 1:ng)
                {
                    t <- t+1
                    wRandom[NPROB+t] <- l
                }
            }
        }

        if(NVC>0)
        {
            wRandom[NPROB+NEF+1:NVC] <- NEF2+1:NVC
        }
        if(NW>0)
        {
            b0Random <- c(b0Random,rep(1,ng-1))
        }
        
        if (ncor0>0) {wRandom[NPROB+NEF+NVC+NW+1:ncor0] <- NEF2+NVC+1:ncor0}
        wRandom[NPM] <- NPM2
        ## wRandom et b0Random ok.
        
        ##------------------------------------------
        ##------nom au vecteur best
        ##--------------------------------------------


        if(ng0>=2){
            nom <-rep(c("intercept",nom.X0[idprob0==1]),each=ng0-1)
            nom1 <- paste(nom," class",c(1:(ng0-1)),sep="")
            names(b)[1:NPROB]<-nom1
        }

        if(ng0==1) names(b)[1:NEF] <- nom.X0[idg0!=0]

        if(ng0>1){
            nom1<- NULL
            for (i in 1:nvar.exp) {
                if(idg0[i]==2){ nom <- paste(nom.X0[i]," class",c(1:ng0),sep="")
                                nom1 <- cbind(nom1,t(nom))}
                if(idg0[i]==1) nom1 <- cbind(nom1,nom.X0[i])
            }
            names(b)[(NPROB+1):(NPROB+NEF)]<- nom1
        }
        if(NVC!=0)names(b)[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- paste("varcov",c(1:(NVC)))
        if(NW!=0)names(b)[(NPROB+NEF+NVC+1):(NPROB+NEF+NVC+NW)] <- paste("varprop class",c(1:(ng0-1)))
        names(b)[NPM] <- "stderr"
        if(ncor0!=0) {names(b)[(NPROB+NEF+NVC+NW+1):(NPROB+NEF+NVC+NW+ncor0)] <- paste("cor",1:ncor0,sep="") }

        N <- NULL
        N[1] <- NPROB
        N[2] <- NEF
        N[3] <- NVC
        N[4] <- NW
        N[5] <- ncor0

        idiag <- as.integer(idiag0)
        idea <- as.integer(idea0)
        nv <- as.integer(nv0)


################ Sortie ###########################

        out <- .Fortran(C_hetmixlin,
                        as.double(Y0),
                        as.double(X0),
                        as.integer(prior0),
                        as.integer(idprob0),
                        as.integer(idea0),
                        as.integer(idg0),
                        as.integer(idcor0),
                        as.integer(ns0),
                        as.integer(ng0),
                        as.integer(nv0),
                        as.integer(nobs0),
                        as.integer(nea0),
                        as.integer(nmes0),
                        as.integer(idiag0),
                        as.integer(nwg0),
                        as.integer(ncor0),
                        as.integer(NPM),
                        best=as.double(b),
                        V=as.double(V),
                        loglik=as.double(loglik),
                        niter=as.integer(ni),
                        conv=as.integer(istop),
                        gconv=as.double(gconv),
                        ppi2=as.double(ppi0),
                        resid_m=as.double(resid_m),
                        resid_ss=as.double(resid_ss),
                        pred_m_g=as.double(pred_m_g),
                        pred_ss_g=as.double(pred_ss_g),
                        predRE=as.double(predRE),
                        as.double(convB),
                        as.double(convL),
                        as.double(convG),
                        as.integer(maxiter),
                        as.integer(fix0))
        
    ### mettre 0 pr les prm fixes
    if(length(posfix))
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
    else
        {
            V <- out$V
        }


            
### Creation du vecteur cholesky
        Cholesky <- rep(0,(nea0*(nea0+1)/2))
        if(idiag0==0 & NVC>0){
            Cholesky[1:NVC] <- out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)]
### Construction de la matrice U 
            U <- matrix(0,nrow=nea0,ncol=nea0)
            U[upper.tri(U,diag=TRUE)] <- Cholesky[1:NVC]
            z <- t(U) %*% U
            out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- z[upper.tri(z,diag=TRUE)]
        }
        if(idiag0==1 & NVC>0){
            id <- 1:nea0
            indice <- rep(id+id*(id-1)/2)
            Cholesky[indice] <- out$best[(NPROB+NEF+1):(NPROB+NEF+nea0)]
            out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)]**2 
        } 

####################################################
        if (nea0>0) {
            predRE <- matrix(out$predRE,ncol=nea0,byrow=T)
            predRE <- data.frame(INDuniq,predRE)
            colnames(predRE) <- c(nom.subject,nom.X0[idea0!=0])
        }


        if(ng0>1) {
            ppi<- matrix(out$ppi2,ncol=ng0,byrow=TRUE)
        }
        else {
            ppi <- matrix(rep(1,ns0),ncol=ng0)
        }


        classif<-apply(ppi,1,which.max)
        ppi<-data.frame(INDuniq,classif,ppi)
        temp<-paste("prob",1:ng0,sep="")
        colnames(ppi) <- c(nom.subject,"class",temp)
        rownames(ppi) <- 1:ns0

        pred_m_g <- matrix(out$pred_m_g,nrow=nobs0)
        pred_ss_g <- matrix(out$pred_ss_g,nrow=nobs0)
        pred_m <- Y0-out$resid_m
        pred_ss <- Y0-out$resid_ss
        pred <- data.frame(IND,pred_m,out$resid_m,pred_ss,out$resid_ss,Y0,pred_m_g,pred_ss_g)

        temp<-paste("pred_m",1:ng0,sep="")
        temp1<-paste("pred_ss",1:ng0,sep="")
        colnames(pred)<-c(nom.subject,"pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1) 

        names(out$best)<-names(b)
        btest <- out$best[1:length(inddepvar.fixed.nom)]
        names(btest) <-inddepvar.fixed.nom
        
### ad 2/04/2012
        if (!("intercept" %in% nom.X0)) X0.names2 <- X0.names2[-1]
### ad
        res <-list(ns=ns0,ng=ng0,idea0=idea0,idprob0=idprob0,idg0=idg0,idcor0=idcor0,loglik=out$loglik,best=out$best,V=V,gconv=out$gconv,conv=out$conv,call=cl,niter=out$niter,N=N,idiag=idiag0,pred=pred,pprob=ppi,predRE=predRE,Xnames=nom.X0,Xnames2=X0.names2,cholesky=Cholesky,na.action=na.action,AIC=2*(length(out$best)-length(posfix)-out$loglik),BIC=(length(out$best)-length(posfix))*log(ns0)-2*out$loglik,data=datareturn,wRandom=wRandom,b0Random=b0Random)
        class(res) <-c("hlme") 
        cost<-proc.time()-ptm
        if(verbose==TRUE) cat("The program took", round(cost[3],2), "seconds \n")

        
        res
    }
