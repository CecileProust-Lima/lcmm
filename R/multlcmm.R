#' Estimation of multivariate mixed-effect models and multivariate latent class
#' mixed-effect models for multivariate longitudinal outcomes of possibly
#' multiple types (continuous Gaussian, continuous non-Gaussian/curvilinear, ordinal)
#' that measure the same underlying latent process.
#' 
#' This function constitutes a multivariate extension of function \code{lcmm}.
#' It fits multivariate mixed models and multivariate latent class mixed models
#' for multivariate longitudinal outcomes of different types. It handles
#' continuous longitudinal outcomes (Gaussian or non-Gaussian, curvilinear) as
#' well as ordinal longitudinal outcomes (with cumulative probit measurement model). 
#' The model assumes that all the outcomes measure the same underlying latent process
#'  defined as their common factor, and each outcome is related to this latent common
#'  factor by a specific parameterized link function. At the latent process level, the
#' model estimates a standard linear mixed model or a latent class linear mixed
#' model when heterogeneity in the population is investigated (in the same way
#' as in functions \code{hlme} and \code{lcmm}). Parameters of the nonlinear link 
#' functions and of the latent process mixed model are estimated simultaneously 
#' using a maximum likelihood method.
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
#' 4. With the "thresholds" link function for an ordinal outcome with levels
#' 0,...,C, C-1 parameters are required for the following transformation:
#' Y(t)=c <=> b_c < L(t) <= b_{c+1} with b_0 = - infinity and b_{C+1}=+infinity.
#' To constraint the parameters to be increasing, except for the first
#' parameter b_1, the program estimates b_k^* (for k=2,...C-1) so that
#' b_{k}=b_{k-1}+(b_k^*)^2.
#' 
#' Details of these parameterized link functions can be found in the papers:
#' Proust-Lima et al. (Biometrics 2006), Proust-Lima et al. (BJMSP 2013),
#' Proust-Lima et al. (arxiv 2021 - https://arxiv.org/abs/2109.13064)
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
#' to remove the parameters of the Beta or Splines link function from the
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
#' D. NUMERICAL INTEGRATION WITH THE THRESHOLD LINK FUNCTION
#'
#' When dealing only with continuous outcomes, the computation of the likelihood does not
#' require any numerical integration over the random-effects, so that the estimation
#' procedure is relatively fast.
#' When at least one ordinal outcome is modeled, a numerical integration over the
#' random-effects is required in each computation of the individual contribution to the
#' likelihood. This achieved using a Monte-Carlo procedure. We allow three options:
#' the standard Monte-Carlo simulations, as well as antithetic Monte-Carlo and quasi
#' Monte-Carlo methods as proposed in Philipson et al (2020).
#' 
#' @param fixed a two-sided linear formula object for specifying the
#' fixed-effects in the linear mixed model at the latent process level. The
#' response outcomes are separated by \code{+} on the left of \code{~} and the
#' covariates are separated by \code{+} on the right of the \code{~}. For
#' identifiability purposes, the intercept specified by default should not be
#' removed by a \code{-1}. Variables on which a contrast above the different
#' outcomes should also be estimated are included with \code{contrast()}.
#' @param mixture a one-sided formula object for the class-specific fixed
#' effects in the latent process mixed model (to specify only for a number of
#' latent classes greater than 1). Among the list of covariates included in
#' \code{fixed}, the covariates with class-specific regression parameters are
#' entered in \code{mixture} separated by \code{+}. By default, an intercept is
#' included. If no intercept, \code{-1} should be the first term included.
#' @param random an optional one-sided formula for the random-effects in the
#' latent process mixed model. At least one random effect should be included
#' for identifiability purposes. Covariates with a random-effect are separated
#' by \code{+}. By default, an intercept is included. If no intercept,
#' \code{-1} should be the first term included.
#' @param subject name of the covariate representing the grouping structure.
#' @param classmb an optional one-sided formula describing the covariates in
#' the class-membership multinomial logistic model. Covariates included are
#' separated by \code{+}. No intercept should be included in this formula.
#' @param ng number of latent classes considered. If \code{ng=1} no
#' \code{mixture} nor \code{classmb} should be specified. If \code{ng>1},
#' \code{mixture} is required.
#' @param idiag optional logical for the variance-covariance structure of the
#' random-effects. If \code{FALSE}, a non structured matrix of
#' variance-covariance is considered (by default). If \code{TRUE} a diagonal
#' matrix of variance-covariance is considered.
#' @param nwg optional logical of class-specific variance-covariance of the
#' random-effects. If \code{FALSE} the variance-covariance matrix is common
#' over latent classes (by default). If \code{TRUE} a class-specific
#' proportional parameter multiplies the variance-covariance matrix in each
#' class (the proportional parameter in the last latent class equals 1 to
#' ensure identifiability).
#' @param randomY optional logical for including an outcome-specific random
#' intercept. If \code{FALSE} no outcome-specific random intercept is added
#' (default). If \code{TRUE} independent outcome-specific random intercepts
#' with parameterized variance are included.
#' @param link optional vector of families of parameterized link functions to
#' estimate (one by outcome). Option "linear" (by default) specifies a linear
#' link function. Other possibilities include "beta" for estimating a link
#' function from the family of Beta cumulative distribution functions,
#' "thresholds" for using a threshold model to describe the correspondence
#' between each level of an ordinal outcome and the underlying latent process and
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
#' @param prior name of the covariate containing the prior on the latent class
#' membership. The covariate should be an integer with values in 0,1,...,ng.
#' When there is no prior, the value should be 0. When there is a prior for the
#' subject, the value should be the number of the latent class (in 1,...,ng).
#' @param range optional vector indicating the range of the outcomes (that is
#' the minimum and maximum). By default, the range is defined according to the
#' minimum and maximum observed values of the outcome. The option should be
#' used only for Beta and Splines transformations.
#' @param subset optional vector giving the subset of observations in
#' \code{data} to use. By default, all lines.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @param posfix Optional vector giving the indices in vector B of the
#' parameters that should not be estimated. Default to NULL, all parameters are
#' estimated.
#' @param partialH optional logical for Splines link functions only.
#' Indicates whether the parameters of the link functions can be dropped from
#' the Hessian matrix to define convergence criteria.
#' @param verbose logical indicating if information about computation should be
#' reported. Default to TRUE.
#' @param returndata logical indicating if data used for computation should be
#' returned. Default to FALSE, data are not returned.
#' @param methInteg character indicating the type of integration if ordinal outcomes
#' are considered. 'MCO' for ordinary Monte Carlo, 'MCA' for antithetic Monte Carlo,
#' 'QMC' for quasi Monte Carlo. Default to "QMC".
#' @param nMC integer, number of Monte Carlo simulations. By default, 1000 points are used
#' if at least one threshold link is specified.
#' @param var.time optional character indicating the name of the time variable.
#' @param nproc the number cores for parallel computation.
#' Default to 1 (sequential mode).
#' @param clustertype optional character indicating the type of cluster for parallel computation.
#' @return The list returned is: \item{ns}{number of grouping units in the
#' dataset} \item{ng}{number of latent classes} \item{loglik}{log-likelihood of
#' the model} \item{best}{vector of parameter estimates in the same order as
#' specified in \code{B} and detailed in section \code{details}}
#' \item{V}{if the model converged (conv=1 or 3), vector containing the upper triangle
#' matrix of variance-covariance estimates of \code{Best} with exception for
#' variance-covariance parameters of the random-effects for which \code{V} contains the
#' variance-covariance estimates of the Cholesky transformed parameters displayed in
#' \code{cholesky}. If conv=2, \code{V} contains the second derivatives of the
#' log-likelihood.} \item{gconv}{vector of convergence criteria: 1. on the
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
#' the latent class: pred_m_1,pred_m_2,...,pred_ss_1,pred_ss_2,...). If \code{var.time}
#' is specified, the corresponding measurement time is also included.}
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
#' @author Cecile Proust-Lima and Viviane Philipps
#' 
#' \email{cecile.proust-lima@@inserm.fr}
#' @seealso
#' 
#' \code{\link{postprob}}, \code{\link{plot.multlcmm}}, \code{\link{predictL}},
#' \code{\link{predictY}} \code{\link{lcmm}}
#' @references
#' 
#' Proust-Lima C, Philipps V, Liquet B (2017). Estimation of Extended Mixed 
#' Models Using Latent Classes and Latent Processes: The R Package lcmm. 
#' Journal of Statistical Software, 78(2), 1-56. doi:10.18637/jss.v078.i02
#' 
#' Proust and Jacqmin-Gadda (2005). Estimation of linear mixed models with a
#' mixture of distribution for the random-effects. Comput Methods Programs
#' Biomed 78: 165-73.
#' 
#' Proust, Jacqmin-Gadda, Taylor, Ganiayre, and Commenges (2006). A nonlinear
#' model with latent process for cognitive evolution using multivariate
#' longitudinal data. Biometrics 62, 1014-24.
#' 
#' Proust-Lima, Dartigues and Jacqmin-Gadda (2011). Misuse of the linear mixed
#' model when evaluating risk factors of cognitive decline. Amer J Epidemiol
#' 174(9): 1077-88.
#' 
#' Proust-Lima, Amieva, Jacqmin-Gadda (2013). Analysis of multivariate mixed
#' longitudinal data: A flexible latent process approach. Br J Math Stat
#' Psychol 66(3): 470-87.
#' 
#' Commenges, Proust-Lima, Samieri, Liquet (2012). A universal approximate
#' cross-validation criterion and its asymptotic distribution, Arxiv.
#'
#' Philipson, Hickey, Crowther, Kolamunnage-Dona (2020). Faster Monte Carlo estimation
#' of semiparametric joint models of time-to-event and multivariate longitudinal data.
#' Computational Statistics & Data Analysis 151.
#'
#' Proust-Lima, Philipps, Perrot, Blanchin, Sebille (2021). Modeling repeated 
#' self-reported outcome data: a continuous-time longitudinal Item Response 
#' Theory model. https://arxiv.org/abs/2109.13064
#' 
#' 
#' @examples
#' 
#' \dontrun{
#' # Latent process mixed model for two curvilinear outcomes. Link functions are 
#' # aproximated by I-splines, the first one has 3 nodes (i.e. 1 internal node 8),
#' # the second one has 4 nodes (i.e. 2 internal nodes 12,25)
#' 
#' m1 <- multlcmm(Ydep1+Ydep2~1+Time*X2+contrast(X2),random=~1+Time,
#' subject="ID",randomY=TRUE,link=c("4-manual-splines","3-manual-splines"),
#' intnodes=c(8,12,25),data=data_lcmm)
#' 
#' # to reduce the computation time, the same model is estimated using 
#' # a vector of initial values
#' m1 <- multlcmm(Ydep1+Ydep2~1+Time*X2+contrast(X2),random=~1+Time,
#' subject="ID",randomY=TRUE,link=c("4-manual-splines","3-manual-splines"),
#' intnodes=c(8,12,25),data=data_lcmm, 
#' B=c(-1.071, -0.192,  0.106, -0.005, -0.193,  1.012,  0.870,  0.881,
#'   0.000,  0.000, -7.520,  1.401,  1.607 , 1.908,  1.431,  1.082,
#'  -7.528,  1.135 , 1.454 , 2.328, 1.052))
#' 
#' 
#' # output of the model
#' summary(m1)
#' # estimated link functions
#' plot(m1,which="linkfunction")
#' # variation percentages explained by linear mixed regression
#' VarExpl(m1,data.frame(Time=0))
#' 
#' #### Heterogeneous latent process mixed model with linear link functions 
#' #### and 2 latent classes of trajectory 
#' m2 <- multlcmm(Ydep1+Ydep2~1+Time*X2,random=~1+Time,subject="ID",
#' link="linear",ng=2,mixture=~1+Time,classmb=~1+X1,data=data_lcmm,
#' B=c( 18,-20.77,1.16,-1.41,-1.39,-0.32,0.16,-0.26,1.69,1.12,1.1,10.8,
#' 1.24,24.88,1.89))
#' # summary of the estimation
#' summary(m2)
#' # posterior classification
#' postprob(m2)
#' # longitudinal predictions in the outcomes scales for a given profile of covariates 
#' newdata <- data.frame(Time=seq(0,5,length=100),X1=0,X2=0,X3=0)
#' predGH <- predictY(m2,newdata,var.time="Time",methInteg=0,nsim=20) 
#' head(predGH)
#' }
#' 
#' @export
#' 
multlcmm <- function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,randomY=FALSE,link="linear",intnodes=NULL,epsY=0.5,cor=NULL,data,B,convB=0.0001,convL=0.0001,convG=0.0001,maxiter=100,nsim=100,prior,range=NULL,subset=NULL,na.action=1,posfix=NULL,partialH=FALSE,verbose=FALSE,returndata=FALSE,methInteg="QMC",nMC=NULL,var.time=NULL,nproc=1,clustertype=NULL)
{
    ptm<-proc.time()

    cl <- match.call()

    nom.subject <- as.character(subject)

#### INCLUSION PRIOR
    nom.prior <- NULL
    if(!missing(prior)) nom.prior <- as.character(prior)
####

    if(!missing(mixture) & ng==1) stop("No mixture can be specified with ng=1")
    if(missing(mixture) & ng>1) stop("The argument mixture has to be specified for ng > 1")
    if(!missing(classmb) & ng==1) stop("No classmb can be specified with ng=1")
    if(missing(random)) stop("At least one random effect is required")
    if(random==~-1) stop("At least one random effect is required")
    if(missing(fixed)) stop("The argument Fixed must be specified in any model")
    if(missing(classmb) & ng==1) classmb <- ~-1
    if(missing(classmb) & ng>1) classmb <- ~1
    if(missing(mixture)) mixture <- ~-1
    if(ng==1&nwg==TRUE) stop ("The argument nwg should be FALSE for ng=1")


    if(!inherits(fixed,"formula")) stop("The argument fixed must be a formula")
    if(!inherits(mixture,"formula")) stop("The argument mixture must be a formula")
    if(!inherits(random,"formula")) stop("The argument random must be a formula")
    if(!inherits(classmb,"formula")) stop("The argument classmb must be a formula")
    if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
    if(nrow(data)==0) stop("Data should not be empty")
    if(missing(subject)){ stop("The argument subject must be specified")}
    if(!is.numeric(data[[subject]])) stop("The argument subject must be numeric")
    if(all(link %in% c("linear","beta","thresholds")) & !is.null(intnodes)) stop("Intnodes should only be specified with splines links")

    if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")

    #if(length(posfix) & missing(B)) stop("A set of initial parameters must be specified if some parameters are not estimated")

    cholesky <- TRUE


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
    if (!isTRUE(all.equal(cor,character(0))))
        {
            if (substr(cor,1,2)=="AR") { ncor0 <- 2 }
            else if (substr(cor,1,2)=="BM") { ncor0 <- 1  }
            else { stop("The argument cor must be of type AR or BM") }

            if(!(strsplit(cor,"-")[[1]][2] %in% colnames(data))) stop("Unable to find time variable from argument 'cor' in 'data'")
            else { cor.var.time <- strsplit(cor,"-")[[1]][2] }
        }
### fin test argument cor

    ##pour acces aux attributs des formules
    afixed <- terms(fixed, specials=c("factor","contrast"))
    if(attr(afixed,"intercept")==0) stop("An intercept should appear in fixed for identifiability purposes")
    amixture <- terms(mixture, specials=c("factor"))
                                        
    arandom <- terms(random, specials=c("factor"))
    aclassmb <- terms(classmb, specials=c("factor"))
    ##fixed sans contrast
    fixed2 <- gsub("contrast","",fixed)
    fixed2 <- formula(paste(fixed2[2],fixed2[3],sep="~"))   
    afixed2 <- terms(fixed2)
                                        
    ##verifier si totes les varialbes sont dans data
    variables <- c(attr(afixed,"variables"),attr(arandom,"variables"),attr(amixture,"variables"),attr(aclassmb,"variables"))
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
    
    ##tjrs intercept dans classmb
    if(attr(aclassmb,"intercept")==0 & ng>1)
        {
            attr(aclassmb,"intercept") <- 1
            cat("The formula in classmb should always include an intercept. An intercept has been added.")
        }

###liste des outcomes
    nomsY <- as.character(attr(afixed,"variables")[2])
    nomsY <- strsplit(nomsY,split=" + ",fixed=TRUE)
    nomsY <- as.vector(nomsY[[1]])
    ny0 <- length(nomsY)

    ##pas de contrast ni randomY si un seul Y
    if(ny0<2 & length(attr(afixed,"specials")$contrast)) stop("No contrast can be included with less than two outcomes")
    if(ny0<2 & randomY==TRUE) stop("With less than 2 outcomes randomY should be FALSE")

###liste des variables utilisees  (sans les interactions et sans les Y)
    ttesLesVar <- colnames(get_all_vars(afixed,data=data[1,]))
    ttesLesVar <- c(ttesLesVar, colnames(get_all_vars(amixture,data=data[1,])))
    ttesLesVar <- c(ttesLesVar, colnames(get_all_vars(arandom,data=data[1,])))
    ttesLesVar <- c(ttesLesVar, colnames(get_all_vars(aclassmb,data=data[1,])))
    if (ncor0>0) ttesLesVar <- unique(c(ttesLesVar,cor.var.time))
    else ttesLesVar <- unique(ttesLesVar)
    ttesLesVar <- setdiff(ttesLesVar, nomsY)

### argument subset
    form1 <- paste(c(nom.subject,nomsY,ttesLesVar,nom.prior),collapse="+")
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

### si subject est un factor
    if(is.factor(data[,nom.subject]))
        {
            data[,nom.subject] <- as.numeric(data[,nom.subject])
        }

    
###subset de data avec les variables utilisees
    
    newdata <- data[,c(nom.subject,nomsY,ttesLesVar,nom.prior)]
    if(!is.null(nom.prior))
        {
            prior <- newdata[,nom.prior]
            newdata[which(is.na(prior)),nom.prior] <- 0
        }

###un data frame par outcome et creation Y0
    dataY <- paste("data",nomsY,sep=".")
    Y0 <- NULL
    IND <- NULL
    outcome <- NULL
    prior <- NULL
    data0 <- NULL
    nayk <- vector("list",ny0)
    for (k in 1:ny0)
        {
            dtemp <- newdata[,c(nom.subject,nomsY[k],ttesLesVar,nom.prior)]
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
            if(!is.null(nom.prior)) prior <- c(prior, dtemp[,nom.prior])
            data0 <- rbind(data0, dtemp[,setdiff(colnames(dtemp),nomsY[k]),drop=FALSE])   #dataset sans NA avec les covariables utilisees; obs ordonnees par outcome
        }

###prior=0 si pas specifie
    if(is.null(prior)) prior <- rep(0,length(Y0))

###creation de X0 (ttes les var + interactions)


    Xfixed <- model.matrix(fixed2[-2], data=data0)
    Xmixture <- model.matrix(mixture, data=data0)
    Xrandom <- model.matrix(random, data=data0)
    Xclassmb <- model.matrix(classmb, data=data0)
    Xcontr <- model.matrix(contr,data=data0)

    z.fixed <- strsplit(colnames(Xfixed),split=":",fixed=TRUE)
    z.fixed <- lapply(z.fixed,sort)

    z.random <- strsplit(colnames(Xrandom),split=":",fixed=TRUE)
    z.random <- lapply(z.random,sort)
    
    if(mixture != ~-1)
        {
            z.mixture <- strsplit(colnames(Xmixture),split=":",fixed=TRUE)
            z.mixture <- lapply(z.mixture,sort)
        }
    else
        {
            z.mixture <- list()
        }
    
    if(classmb != ~-1)
        {
            z.classmb <- strsplit(colnames(Xclassmb),split=":",fixed=TRUE)
            z.classmb <- lapply(z.classmb,sort)
        }
    else
        {
            z.classmb <- list()
        }
    
    if(contr != ~-1)
        {
            z.contr <- strsplit(colnames(Xcontr),split=":",fixed=TRUE)
            z.contr <- lapply(z.contr,sort)
        }
    else
        {
            z.contr <- list()
        }
    
    if(!all(z.mixture %in% z.fixed))  stop("The covariates in mixture should also be included in the argument fixed")
    if(!all(z.contr %in% z.fixed))  stop("The covariates in contrast should also appear in fixed")
    
    X0 <- cbind(Xfixed, Xrandom, Xclassmb)        
    nom.unique <- unique(colnames(X0))
    X0 <- X0[,nom.unique,drop=FALSE]
    
    if (ncor0>0)
        {
            if(!(cor.var.time %in% colnames(X0)))
                {
                    X0 <- cbind(X0, data0[,cor.var.time])
                    colnames(X0) <- c(nom.unique, cor.var.time)
                    nom.unique <- c(nom.unique,cor.var.time)
                }
        }

    X0 <- as.matrix(X0)
###X0 fini

    timeobs <- rep(0, nrow(data0))
    if(!is.null(var.time))
    {
        timeobs <- data0[,var.time]
        if(any(is.na(timeobs))) stop(paste("Cannot use",var.time,"as time variable because it contains missing data"))
    }

###test de link
    if (length(link)!=1 & length(link)!=ny0) stop("One link per outcome should be specified")
    if(any(link %in% c("splines","Splines")))
        {
            link[which(link %in% c("splines","Splines"))] <- "5-quant-splines"
        }
    if(length(link)==1 & ny0>1)
        {
            link <- rep(link, ny0)
        }

    idlink0 <- rep(2,ny0)
    idlink0[which(link=="linear")] <- 0
    idlink0[which(link=="beta")] <- 1
    idlink0[which(link=="thresholds")] <- 3

    spl <- strsplit(link[which(idlink0==2)],"-")
    if(any(sapply(spl,length)!=3)) stop("Invalid argument 'link'")

    nySPL <- length(spl)
    nybeta <- sum(idlink0==1)
    nyORD <- sum(idlink0==3)
    ##remplir range si pas specifie
    if(!is.null(range) & length(range)!=2*(nySPL+nybeta)) stop("Length of vector range is not correct.")
    if((length(range)==2*(nySPL+nybeta)) & (nySPL+nybeta>0))
        {
            ind12 <- which(idlink0==1 | idlink0==2)
            for (i in 1:(nySPL+nybeta))
                {
                    rg <- range(get(dataY[ind12[i]])[,nomsY[ind12[i]]])
                    if(rg[1]<range[2*(i-1)+1] | rg[2]>range[2*(i-1)+2]) stop("The range specified do not cover the entire range of the data")
                }
        }
    if((is.null(range) & (nybeta+nySPL)>0) | length(range)!=2*(nySPL+nybeta))
        {
            range <- NULL
            for(k in which(idlink0!=0))
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
    if (any(idlink0==1))
        {
            if (any(epsY<=0))
                {
                    stop("Argument 'epsY' should be positive.")
                }

            if(length(epsY)==1) epsY <- rep(epsY,nybeta)
            
            if(length(epsY)!=nybeta) stop(paste("Argument 'epsY' should be of length",nybeta))
            if(nybeta!=ny0)
                {
                    epsY2 <- rep(0,ny0)
                    epsY2[which(idlink0==1)] <- epsY
                    epsY <- epsY2
                }

        } 
    

    nbzitr0 <- rep(2,ny0) #nbzitr0 = nb de noeuds si splines, 2 sinon
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
    nbzitr0[which(idlink0==2)] <- nbnodes

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
    for (k in 1:ny0)
        {
            if (idlink0[k]!=2) next
            else
                {                                         
                    nbspl <- nbspl+1

                    if(spltype[nbspl]=="manual")
                        {
                            nodes <- intnodes[(nb+1):(nb+nbnodes[nbspl]-2)]
                            if(!length(nodes)) stop("The length of intnodes is not correct")
                            intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <-  nodes
                            nb <- nb+nbnodes[nbspl]-2

                            idrg <- length(which(idlink0[1:k] != 0))
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
    zitr <- matrix(0,max(m,2),ny0)
    nb12 <- 0
    nbspl <- 0
    for (k in 1:ny0)
        {
            if((idlink0[k]==0) | (idlink0[k]==3)) zitr[1:2,k] <- c(min(get(dataY[k])[,nomsY[k]]),max(get(dataY[k])[,nomsY[k]]))

            if(idlink0[k]==1)
                {
                    nb12 <- nb12 + 1
                    zitr[1:2,k] <- range[2*(nb12-1)+1:2]
                }

            if(idlink0[k]==2)
                {
                    nb12 <- nb12+1
                    nbspl <- nbspl+1
                    zitr[2:(nbzitr0[k]-1),k] <- intnodes2[ifelse(nbspl==1,0,1)*sum(nbnodes[1:(nbspl-1)]-2) + 1:(nbnodes[nbspl]-2)]
                    zitr[1,k] <- range[2*(nb12-1)+1]
                    zitr[nbnodes[nbspl],k]  <- range[2*nb12]
                    
                    ##verifier s'il y a des obs entre les noeuds
                    hcounts <- hist(get(dataY[k])[,nomsY[k]],breaks=zitr[1:nbnodes[nbspl],k],plot=FALSE,include.lowest=TRUE,right=TRUE)$counts
                    if(any(hcounts==0)) stop(paste("Link function number",k,"can not be estimated. Please try other nodes such that there are observations in each interval."))    
                }
        }

###uniqueY0 et indiceY0
    uniqueY0 <- NULL
    indiceY0 <- NULL
    nvalSPLORD0 <- rep(0,ny0)
    nbmod <- rep(0,ny0)
    modalites <- vector("list",ny0)
    nb <- 0
    for (k in 1:ny0)
        {
            if((idlink0[k]!=2) & (idlink0[k]!=3))
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
                    if(idlink0[k]==2)
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
                    nvalSPLORD0[k] <- length(uniqueTemp)
                }
            else
                {
                    uniqueY0 <- c(uniqueY0, yk)
                    indiceY0 <- c(indiceY0, ifelse(idlink0[k]==2,nb,0)+c(1:length(yk)))
                    nb <- nb + length(yk)
                    nvalSPLORD0[k] <- length(yk)
                }
            if(idlink0[k]==3)
            {
                nbmod[k] <- length(na.omit(uniqueTemp))
                modalites[[k]] <- uniqueTemp
            }
        }
    #if(is.null(nvalSPLORD0)) nvalSPLORD0 <- 0


###ordonner les mesures par individu
    #IDnum <- as.numeric(IND)
    matYX <- cbind(IND,timeobs,prior,Y0,indiceY0,outcome,X0)
    matYXord <- matYX[order(IND),]
    Y0 <- as.numeric(matYXord[,4])
    X0 <- apply(matYXord[,-c(1,2,3,4,5,6),drop=FALSE],2,as.numeric)
                                        #X0 <- as.matrix(X0)  a remettre si X0 <- as.data.frame(X0) remis l.211
    IND <- matYXord[,1]
    outcome <- matYXord[,6]
    indiceY0 <- as.numeric(matYXord[,5])
    prior0 <- as.numeric(unique(matYXord[,c(1,3)])[,2])
    if(length(prior0)!=length(unique(IND))) stop("Please check 'prior' argument. Subjects can not have multiple assigned classes.")
    timeobs <- matYXord[,2]
        

###parametres pour hetmixMult
    ns0 <- length(unique(IND))
    ng0 <- ng
    nv0 <- dim(X0)[2]
    nobs0 <- length(Y0)
    idiag0 <- ifelse(idiag==TRUE,1,0)
    nwg0 <- ifelse(nwg==TRUE,1,0)
    nalea0 <- ifelse(randomY==TRUE,ny0,0)
    chol <- ifelse(cholesky==TRUE,1,0)

    loglik <- 0
    ni <- 0
    istop <- 0
    gconv <- rep(0,3)
    ppi0 <- rep(0,ns0*ng0)
    resid_m <- rep(0,nobs0)
    resid_ss <- rep(0,nobs0)
    pred_m_g <- rep(0,nobs0*ng0)
    pred_ss_g <- rep(0,nobs0*ng0)
    Yobs <- rep(0,nobs0)
                                        
    predRE_Y <- rep(0,ns0*nalea0)
    rlindiv <- rep(0,ns0)
    marker <- rep(0,nsim*ny0)
    transfY <- rep(0,nsim*ny0)
    Ydiscrete <- 0
    UACV <- 0
    vraisdiscret <- 0
    
 
    
###nmes0
    nmes0 <- matrix(0,ns0,ny0)
    for (k in 1:ny0)
        {
            INDpresents <- which(unique(IND) %in% get(dataY[k])[,nom.subject])
            nmes0[INDpresents,k] <- as.vector(table(get(dataY[k])[,nom.subject]))
        }
    maxmes <- max(apply(nmes0,1,sum))



###remplir idprob, etc
    z.X0 <- strsplit(nom.unique,split=":",fixed=TRUE)
    z.X0 <- lapply(z.X0,sort)
    
    idprob0 <- z.X0 %in% z.classmb + 0
    idea0 <- z.X0 %in% z.random + 0
    idg0 <- (z.X0 %in% z.fixed) + (z.X0 %in% z.mixture)
    idcontr0 <- z.X0 %in% z.contr + 0
    

    if (ncor0>0) idcor0 <- colnames(X0) %in% cor.var.time +0
    else idcor0 <- rep(0,nv0)


    nea0 <- sum(idea0)
    predRE <- rep(0,ns0*nea0)

    ## parametres MC
    if(is.null(nMC)) nMC <- ifelse(all(idlink0 != 3), 0, 1000)
    methInteg <- switch(methInteg,"MCO"=1,"MCA"=2,"QMC"=3)
    seqMC <- 0
    dimMC <- 0
    if(methInteg==3) 
    {
        dimMC <- nea0+nalea0
        if(ncor0>0) dimMC <- dimMC+maxmes
        # dimMC <- max(nea0, nalea0, maxmes) #?? !**
        if(dimMC>0) seqMC <- randtoolbox::sobol(n=nMC,dim=dimMC,normal=TRUE,scrambling=1) 
    }
    
    ##nombre total de parametres
    NPM <- (ng0-1)*sum(idprob0) + sum(idg0==1)-1 + ng0*sum(idg0==2) + ncor0 + (ny0-1)*sum(idcontr0) +
        ifelse(idiag0==1,nea0,nea0*(nea0+1)/2)-1 + (ng0-1)*nwg0 + nalea0 + ny0 +
            2*sum(idlink0==0) + 4*sum(idlink0==1) + sum(nbnodes+2) + any(idlink0==3)*sum(nbmod[which(idlink0==3)]-1)

    V <- rep(0, NPM*(NPM+1)/2)  #pr variance des parametres

    nef <- (ng0-1)*sum(idprob0) + sum(idg0==1)-1 + ng0*sum(idg0==2) + (ny0-1)*sum(idcontr0)
    ncontr <- (ny0-1)*sum(idcontr0)
    nvc <- ifelse(idiag0==1,nea0,nea0*(nea0+1)/2)-1
    nw <- (ng0-1)*nwg0
    ntrtot0 <- nbzitr0+2
    ntrtot0[which(idlink0==0)] <- 2
    ntrtot0[which(idlink0==3)] <- nbmod[which(idlink0==3)]-1
    nprob <- sum(idprob0)*(ng0-1)

    ntr <- rep(0,ny0)
    ntr[which(idlink0==0)] <- 2
    ntr[which(idlink0==1)] <- 4
    ntr[which(idlink0==2)] <- nbzitr0[which(idlink0==2)]+2
    ntr[which(idlink0==3)] <- nbmod[which(idlink0==3)]-1

    
## gestion de B=random(mod)

        Brandom <- FALSE
        if(length(cl$B)==2)
            {
                if(!inherits(eval(cl$B[[2]]),"multlcmm")) stop("The model specified in B should be of class multlcmm")
                if(as.character(cl$B[1])!="random") stop("Please use random() to specify random initial values")
                
                Brandom <- TRUE
                B <- eval(cl$B[[2]])
                if(B$conv != 1) stop("Model in argument B did not converge properly")
                #if(length(posfix)) stop("Argument posfix is not compatible with random intial values")
            }
    
###valeurs initiales
    if(!(missing(B)))
        {
            if(is.vector(B))
                {
                    if (length(B)==NPM) b <- B
                    else stop(paste("Vector B should be of length",NPM))

                    if(nvc>0)
                    {
                        ## remplacer varcov des EA par les prm a estimer
                        
                        if(idiag==1) b[nef+1:nvc] <- sqrt(b[nef+1:nvc])

                        if(idiag==0)
                        {
                                varcov <- matrix(0,nrow=nea0,ncol=nea0)
                                varcov[upper.tri(varcov,diag=TRUE)] <- c(1,b[nef+1:nvc])
                                varcov <- t(varcov)
                                varcov[upper.tri(varcov,diag=TRUE)] <- c(1,b[nef+1:nvc])
                         
                            if(cholesky==TRUE)
                            {       
                                ch <- chol(varcov)
                                
                                b[nef+1:nvc] <- (ch[upper.tri(ch,diag=TRUE)])[-1]
                            }
                            else
                            {
                                corr <- cov2cor(varcov)
                                corr <- corr[upper.tri(corr)]

                                prmea <- matrix(0,nea0,nea0)
                                diag(prmea) <- sqrt(diag(varcov))
                                prmea[upper.tri(prmea)] <- log((1+corr)/(1-corr))                                
                                
                                b[nef+1:nvc] <- (prmea[upper.tri(prmea,diag=TRUE)])[-1]
                                
                            }
                        }
                    }
                }
            else
                {
                    if(!inherits(B,"multlcmm")) stop("B should be either a vector or an object of class multlcmm")

                    if(ng==1 & B$ng==1)
                        {
                            if(length(B$best)!=NPM) stop("B is not correct")
                            b <- B$best
                        }

                    if(ng>1 & B$ng==1)
                            {
                                nef2 <- sum(idg0!=0)-1 + (ny0-1)*sum(idcontr0)
                                NPM2 <- nef2+ nvc+ncor0+nalea0+ny0+sum(ntrtot0)
                              
                                if(length(B$best)!=NPM2) stop("B is not correct")

                                if(Brandom==FALSE)
                                    {
                                        ### B deterministe
                                        b <- rep(0,NPM)

                                        ## calcul des valeurs initiales pour les effets fixes
                                        l <- 0
                                        t <- 0
                                        for (i in 1:nv0)
                                            {
                                                if(idg0[i]==1 & i>1)
                                                    {
                                                        l <- l+1
                                                        t <- t+1
                                                        b[nprob+t] <- B$best[l]
                                                    }
                                                if(idg0[i]==2)
                                                    {
                                                        if (i==1)
                                                            {
                                                                for (g in 2:ng0)
                                                                    {
                                                                        t <- t+1
                                                                        b[nprob+t] <- -0.5*(g-1)
                                                                    }
                                                            }
                                                        if (i>1)
                                                            {
                                                                l <- l+1
                                                                for (g in 1:ng0)
                                                                    {
                                                                        t <- t+1
                                                                        if(B$conv==1) b[nprob+t] <- B$best[l]+(g-(ng0+1)/2)*sqrt(B$V[l*(l+1)/2])
                                                                        else b[nprob+t] <- B$best[l]+(g-(ng0+1)/2)*B$best[l]
                                                                    }
                                                            }
                                                    }
                                            }

                                        ## remplacer varcov par cholesky pour les effets aleatoires
                                        if(nvc>0)
                                            {
                                                if(idiag==TRUE)
                                                    {
                                                        b[nef+1:nvc] <- B$cholesky[(1:nea0)*(2:(nea0+1))/2][-1]
                                                        
                                                    }
                                                else
                                                    {
                                                        b[nef+1:nvc] <- B$cholesky[-1]
                                                    }
                                            }

                                        ## les autres parametres sont inchanges
                                        if (ncor0>0) {b[nef+nvc+nw+1:ncor0] <- B$best[nef2+nvc+1:ncor0]}
                                        b[nef+nvc+nw+ncor0+1:ny0] <- B$best[nef2+nvc+ncor0+1:ny0]
                                        b[nef+nvc+nw+ncor0+ny0+1:nalea0] <- B$best[nef2+nvc+ncor0+ny0+1:nalea0]
                                        b[(nef+nvc+nw+ncor0+ny0+nalea0+1):NPM] <-B$best[(nef2+nvc+ncor0+ny0+nalea0+1):NPM2]

                                    }
                                else
                                    {
                                        ### B random

                                        ## initialiser le vecteur bb contenant les prm (avec repetition) du modele dans B et sa variance vbb
                                        bb <- rep(0,NPM-nprob-nw)
                                        vbb <- matrix(0,NPM-nprob-nw,NPM-nprob-nw)
                                        
                                        VB <- matrix(0,NPM2,NPM2)
                                        VB[upper.tri(VB,diag=TRUE)] <- B$V
                                        VB <- t(VB)
                                        VB[upper.tri(VB,diag=TRUE)] <- B$V

                                        nbg <- idg0[which(idg0!=0)]
                                        nbg[which(nbg==2)] <- ng
                                        nbgnef <- unlist(sapply(nbg,function(k) if(k>1) rep(2,k) else k))
                                        nbgnef <- nbgnef[-1]
                                        nbg <- nbg[-1]
                                        
                                        vbb[which(nbgnef==1),setdiff(1:ncol(vbb),which(nbgnef!=1))] <- VB[which(nbg==1),setdiff(1:ncol(VB),which(nbg!=1))]
                                        vbb[(nef-nprob+1):nrow(vbb),(nef-nprob+1):ncol(vbb)] <- VB[(nef2+1):nrow(VB),(nef2+1):ncol(VB)]

                                        

                                        ## remplir les effets fixes (avec repetition si effet specifique a la classe)
                                        l <- 0
                                        t <- 0
                                        for (i in 1:nv0)
                                            {
                                                if(idg0[i]==1)
                                                    {
                                                        if(i==1) next
                                                        l <- l+1
                                                        t <- t+1
                                                        bb[t] <- B$best[l]
                                                    }
                                                if(idg0[i]==2)
                                                    {
                                                        if(i==1)
                                                            {
                                                                t <- t+ng-1
                                                                next
                                                            }
                                                        l <- l+1
                                                        for (g in 1:ng)
                                                            {
                                                                t <- t+1
                                                                bb[t] <- B$best[l]
                                                                vbb[t,t] <- VB[l,l]
                                                            }
                                                    }
                                            }

                                        ## remplacer les varcov par la cholesky
                                        if(nvc>0)
                                            {
                                                if(idiag==TRUE)
                                                    {
                                                        bb[nef-nprob+1:nvc] <- B$cholesky[(1:nea0)*(2:(nea0+1))/2][-1]
                                                    }
                                                else
                                                    {
                                                        bb[nef-nprob+1:nvc] <- B$cholesky[-1]
                                                    }
                                            }
                            
                                        ##les autres parametres sont inchanges
                                        if (ncor0>0)
                                            {
                                                bb[nef-nprob+nvc+1:ncor0] <- B$best[(NPM2-ncor0):(NPM2-1)]
                                            }

                                                                      
                                        bb[nef-nprob+nvc+ncor0+1:ny0] <- B$best[nef2+nvc+ncor0+1:ny0]

                                        if(nalea0>0)
                                            {
                                              bb[nef-nprob+nvc+ncor0+ny0+1:nalea0] <- B$best[nef2+nvc+ncor0+ny0+1:nalea0]  
                                            }

                                        bb[nef-nprob+nvc+ncor0+ny0+nalea0+1:sum(ntrtot0)] <- B$best[nef2+nvc+ncor0+ny0+nalea0+1:sum(ntrtot0)]

                                        ## on enleve les intercepts car ils seront tous initialises a 0
                                        if(idg0[1]>1)
                                            {
                                                bb <- bb[-(1:(ng-1))]
                                                vbb <- vbb[-(1:(ng-1)),-(1:(ng-1))]
                                            }
                                                            
                                        up <- vbb[upper.tri(vbb,diag=TRUE)]
                                        vbb <- t(vbb)
                                        vbb[upper.tri(vbb,diag=TRUE)] <- up
                                        #Chol <- chol(vbb)
                                        #Chol <- t(Chol)

                                        ## vecteur final b cree a partir de bb et vbb
                                        b <- rep(0,NPM)
                                        
                                        if(idg0[1]>1)
                                            {
                                                b[c((nprob+ng):(nef+nvc),(nef+nvc+nw+1):NPM)] <- rmvnorm(n=1,mean=bb,sigma = vbb) #bb + Chol %*% rnorm(length(bb))
                                                b[nprob+1:(ng-1)] <- 0
                                            } 
                                        else
                                            {                                        
                                                b[c((nprob+1):(nef+nvc),(nef+nvc+nw+1):NPM)] <- rmvnorm(n=1,mean=bb,sigma=vbb) #bb + Chol %*% rnorm(NPM-nprob-nw)
                                            }

                                        ## les prm de classmb et nwg sont toujours initalises a 0 et 1
                                        b[1:nprob] <- 0
                                        if(nw>0) b[nef+nvc+1:nw] <- 1

                                        ## if(nvc>0)
                                        ##     {
                                        ##         cholRE <- matrix(0,nea0,nea0)
                                        ##         cholRE[upper.tri(cholRE,diag=TRUE)] <- c(1,b[nef+1:nvc])
                                        ##         varcovRE <- t(cholRE) %*% cholRE
                                        ##         b[nef+1:nvc] <- varcovRE[upper.tri(varcovRE,diag=TRUE)][-1]   
                                        ##     }
                                                                                
                                    }
                            }
                                
                }
        }
    else ## B missing
        {
            b <- rep(0,NPM)
            if (nvc>0)
                {
                    if(idiag0==1) b[nef+1:nvc] <- rep(1,nvc)
                    if(idiag0==0)
                        {
                            init.nvc <- diag(nea0)
                            init.nvc <- init.nvc[upper.tri(init.nvc, diag=TRUE)]
                            b[nef+1:nvc] <- init.nvc[-1]
                        }
                }
            if(nwg0>0) b[nef+nvc+1:nw] <- 1
            if(ncor0==1) b[nef+nvc+nw+1] <- 1
            if(ncor0==2) b[nef+nvc+nw+1:2] <- c(0,1)

            b[nef+nvc+nw+ncor0+1:ny0] <-  1

            if(nalea0>0) b[nef+nvc+nw+ncor0+ny0+1:nalea0] <- 1

            for(k in 1:ny0)
                {
                    if(idlink0[k]==0)
                        {
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-1] <- mean(get(dataY[k])[,nomsY[k]])
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])] <- 1
                        }
                    if(idlink0[k]==1)
                        {
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-3] <- 0
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-2] <- -log(2)
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-1] <- 0.7
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])] <- 0.1
                        }
                    if(idlink0[k]==2)
                        {
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-ntrtot0[k]+1] <- -2
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-ntrtot0[k]+2:ntrtot0[k]] <- 0.1
                    }
                    if(idlink0[k]==3)
                    {
                        b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-ntrtot0[k]+1] <- 0
                        if(ntrtot0[k]>1) b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-ntrtot0[k]+2:ntrtot0[k]] <- 0.1
                    }
                }
        }

    ## faire wRandom et b0Random
    nef2 <- sum(idg0!=0)-1 + (ny0-1)*sum(idcontr0)
    NPM2 <- nef2+ nvc+ncor0+nalea0+ny0+sum(ntrtot0)
    
    wRandom <- rep(0,NPM)
    b0Random <- rep(0,ng-1)
    
    l <- 0
    t <- 0
    m <- 0
    for (i in 1:nv0)
    {
        if(idg0[i]==1)
        {
            if(i==1) next
            l <- l+1
            t <- t+1
            wRandom[nprob+t] <- l
        }
        if(idg0[i]==2)
        {
            if(i==1)
            {
                t <- t+ng-1
                b0Random <- c(b0Random,rep(0,ng-1))
                next
            }
            l <- l+1
            for (g in 1:ng)
            {
                t <- t+1
                wRandom[nprob+t] <- l
            }
        }
        if(idcontr0[i]==1)
        {
            wRandom[nef-ncontr+m+1:(ny0-1)] <- nef2-ncontr+m+1:(ny0-1)
            m <- m+ny0-1
        }
    }

    if(nvc>0)
    {
        wRandom[nef+1:nvc] <- nef2+1:nvc
    }
    if(nw>0)
    {
        b0Random <- c(b0Random,rep(1,ng-1))
    }
    
    if(ncor0>0) wRandom[nef+nvc+nw+1:ncor0] <- nef2+nvc+1:ncor0
    
    wRandom[nef+nvc+nw+ncor0+1:ny0] <- nef2+nvc+ncor0+1:ny0

    if(nalea0>0)
    {
        wRandom[nef+nvc+nw+ncor0+ny0+1:nalea0] <- nef2+nvc+ncor0+ny0+1:nalea0
    }

    wRandom[nef+nvc+nw+ncor0+ny0+nalea0+1:sum(ntrtot0)] <- nef2+nvc+ncor0+ny0+nalea0+1:sum(ntrtot0)
    ## wRandom et b0Random ok.


    ##------------------------------------------
    ##------nom au vecteur best
    ##--------------------------------------------

    nom.X0 <- colnames(X0)
    nom.X0[nom.X0=="(Intercept)"] <- "intercept"
    if(ng0>=2)
        {
            nom <-rep(nom.X0[idprob0==1],each=ng0-1)
            nom1 <- paste(nom," class",c(1:(ng0-1)),sep="")
            names(b)[1:nprob]<-nom1
        }


    if(ng0==1) names(b)[1:(nef-ncontr)] <- nom.X0[-1][idg0[-1]!=0]
    if(ng0>1){
 	nom1<- NULL
 	for (i in 1:nv0) {
            if(idg0[i]==2){
                if (i==1){
                    nom <- paste(nom.X0[i]," class",c(2:ng0),sep="")
                    nom1 <- cbind(nom1,t(nom))
                }
                if (i>1){
                    nom <- paste(nom.X0[i]," class",c(1:ng0),sep="")
                    nom1 <- cbind(nom1,t(nom))
                }
            }
            if(idg0[i]==1 & i>1) nom1 <- cbind(nom1,nom.X0[i])
 	}
        names(b)[(nprob+1):(nef-ncontr)]<- nom1
    }

    if(idlink0[1]==0) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+1:ntrtot0[1]]<- c("Linear 1","Linear 2")
    if(idlink0[1]==1) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+1:ntrtot0[1]]<- paste("Beta",c(1:ntrtot0[1]),sep="")
    if(idlink0[1]==2) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+1:ntrtot0[1]]<- paste("I-splines",c(1:ntrtot0[1]),sep="")
    if(idlink0[1]==3) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+1:ntrtot0[1]]<- paste("Thresh",c(1:ntrtot0[1]),sep="")
    if(ny0>1)
        {
            for (yk in 2:ny0)
                {
                    if(idlink0[yk]==0) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:(yk-1)])+1:ntrtot0[yk]]<- c("Linear 1","Linear 2")
                    if(idlink0[yk]==1) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:(yk-1)])+1:ntrtot0[yk]]<- paste("Beta",c(1:ntrtot0[yk]),sep="")
                    if(idlink0[yk]==2) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:(yk-1)])+1:ntrtot0[yk]]<- paste("I-splines",c(1:ntrtot0[yk]),sep="")
                    if(idlink0[yk]==3) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:(yk-1)])+1:ntrtot0[yk]]<- paste("Thresh",c(1:ntrtot0[yk]),sep="")
                }
        }
    if(nvc!=0)names(b)[nef+1:nvc] <- paste("varcov",c(1:nvc))
    if(nw!=0)names(b)[nef+nvc+1:nw] <- paste("varprop class",c(1:(ng0-1)))

    names(b)[nef+nvc+nw+ncor0+1:ny0] <- paste("std.err",1:ny0)

    if(ncor0>0) {names(b)[nef+nvc+nw+1:ncor0] <- paste("cor",1:ncor0,sep="")}
    if(nalea0!=0) names(b)[nef+nvc+nw+ncor0+ny0+1:nalea0] <- paste("std.randomY",1:ny0,sep="")
    if(ncontr!=0) names(b)[(nef-ncontr+1):nef] <- paste("contrast",paste(rep(1:sum(idcontr0),each=ny0-1),rep(1:(ny0-1),sum(idcontr0)),sep=""),sep="")


    ## prm fixes
    fix0 <- rep(0,NPM)
    if(length(posfix))
        {
            if(any(!(posfix %in% 1:NPM))) stop("Indexes in posfix are not correct")

            fix0[posfix] <- 1
        }
    if(length(posfix)==NPM) stop("No parameter to estimate")
    names_best <- names(b)

    
    ## pour H restreint
    Hr0 <- as.numeric(partialH)
    pbH0 <- rep(0,NPM)
    if(is.logical(partialH))
    {
        if(partialH) pbH0[grep("I-splines",names(b))] <- 1
        pbH0[posfix] <- 0
        if(sum(pbH0)==0 & Hr0==1) stop("No partial Hessian matrix can be defined")
    }
    else
    {
        if(!all(Hr0 %in% 1:NPM)) stop("Indexes in partialH are not correct")
        pbH0[Hr0] <- 1
        pbH0[posfix] <- 0
    }
    indexHr <- NULL
    if(sum(pbH0)>0) indexHr <- which(pbH0==1)
    
    ##initialisation pour ng>1
    if(missing(B) & ng0>1)
        {
            stop("Please specify initial values with argument 'B'")
        }
#browser()

    nfix <- sum(fix0)
    bfix <- 0
    if(nfix>0)
    {
        bfix <- b[which(fix0==1)]
        b <- b[which(fix0==0)]
        NPM <- NPM-nfix
    }

    if(maxiter==0)
    {
        vrais <- loglikmultlcmm(b,Y0,X0,prior0,idprob0,idea0,idg0,idcor0,idcontr0,
                                ny0,ns0,ng0,nv0,nobs0,nea0,nmes0,idiag0,nwg0,ncor0,
                                nalea0,NPM,epsY,idlink0,nbzitr0,zitr,uniqueY0,
                                indiceY0,nvalSPLORD0,fix0,nfix,bfix,methInteg,nMC,
                                dimMC,seqMC,chol)
        
        out <- list(conv=2, V=rep(NA, length(b)), best=b,
                    ppi2=rep(NA,ns0*ng0), predRE=rep(NA,ns0*nea0),
                    predRE_Y=rep(NA,ns0*nalea0), Yobs=rep(NA,nobs0),
                    resid_m=rep(NA,nobs0), resid_ss=rep(NA,nobs0),
                    marker=rep(NA,nsim*ny0), transfY=rep(NA,nsim*ny0),
                    pred_m_g=rep(NA,nobs0*ng0), pred_ss_g=rep(NA,nobs0*ng0),
                    gconv=rep(NA,3), niter=0, loglik=vrais)
    }
    else
    {   
        res <- mla(b=b, m=length(b), fn=loglikmultlcmm,
                   clustertype=clustertype,.packages=NULL,
                   epsa=convB,epsb=convL,epsd=convG,partialH=indexHr,
                   digits=8,print.info=verbose,blinding=FALSE,
                   multipleTry=25,file="",
                   nproc=nproc,maxiter=maxiter,minimize=FALSE,
                   Y0=Y0,X0=X0,prior0=prior0,idprob0=idprob0,idea0=idea0,idg0=idg0,
                   idcor0=idcor0,idcontr0=idcontr0,ny0=ny0,ns0=ns0,ng0=ng0,nv0=nv0,
                   nobs0=nobs0,nea0=nea0,nmes0=nmes0,idiag0=idiag0,nwg0=nwg0,
                   ncor0=ncor0,nalea0=nalea0,npm0=NPM,epsY0=epsY,idlink0=idlink0,
                   nbzitr0=nbzitr0,zitr0=zitr,uniqueY0=uniqueY0,indiceY0=indiceY0,
                   nvalSPLORD0=nvalSPLORD0,fix0=fix0,nfix0=nfix,bfix0=bfix,
                   methInteg0=methInteg,nMC0=nMC,dimMC0=dimMC,seqMC0=seqMC,chol0=chol)
        
        out <- list(conv=res$istop, V=res$v, best=res$b,
                    ppi2=rep(NA,ns0*ng0), predRE=rep(NA,ns0*nea0),
                    predRE_Y=rep(NA,ns0*nalea0), Yobs=rep(NA,nobs0),
                    resid_m=rep(NA,nobs0), resid_ss=rep(NA,nobs0),
                    marker=rep(NA,nsim), transfY=rep(NA,nsim),
                    pred_m_g=rep(NA,nobs0*ng0), pred_ss_g=rep(NA,nobs0*ng0),
                    gconv=c(res$ca, res$cb, res$rdm), niter=res$ni,
                    loglik=res$fn.value)
    }
    
    if(out$conv %in% c(1,2,3))
    {
        estim0 <- 0
        ll <- 0
        ppi0 <- rep(0,ns0*ng0)
        resid_m <- rep(0,nobs0)
        resid_ss <- rep(0,nobs0)
        pred_m_g <- rep(0,nobs0*ng0)
        pred_ss_g <- rep(0,nobs0*ng0)
        predRE <- rep(0,ns0*nea0)
        predRE_Y <- rep(0,ns0*nalea0)
        marker <- rep(0,nsim*ny0)
        transfY <- rep(0,nsim*ny0)
        Yobs <- rep(0,nobs0)
        Ydiscret <- 0
        vraisdiscret <- 0
        UACV <- 0
        rlindiv <- rep(0,ns0)
        post <- .Fortran(C_loglikmultlcmm,
                         as.double(Y0),
                         as.double(X0),
                         as.integer(prior0),
                         as.integer(idprob0),
                         as.integer(idea0),
                         as.integer(idg0),
                         as.integer(idcor0),
                         as.integer(idcontr0),
                         as.integer(ny0),
                         as.integer(ns0),
                         as.integer(ng0),
                         as.integer(nv0),
                         as.integer(nobs0),
                         as.integer(nea0),
                         as.integer(nmes0),
                         as.integer(idiag0),
                         as.integer(nwg0),
                         as.integer(ncor0),
                         as.integer(nalea0),
                         as.integer(NPM),
                         best=as.double(out$best),
                         ppi2=as.double(ppi0),
                         resid_m=as.double(resid_m),
                         resid_ss=as.double(resid_ss),
                         pred_m_g=as.double(pred_m_g),
                         pred_ss_g=as.double(pred_ss_g),
                         predRE=as.double(predRE),
                         predRE_Y=as.double(predRE_Y),
                         as.double(epsY),
                         as.integer(idlink0),
                         as.integer(nbzitr0),
                         as.double(zitr),
                         as.double(uniqueY0),
                         as.integer(indiceY0),
                         as.integer(nvalSPLORD0),
                         marker=as.double(marker),
                         transfY=as.double(transfY),
                         as.integer(nsim),
                         Yobs=as.double(Yobs),
                         as.integer(Ydiscret),
                         vraisdiscret=as.double(vraisdiscret),
                         UACV=as.double(UACV),
                         rlindiv=as.double(rlindiv),
                         as.integer(fix0),
                         as.integer(nfix),
                         as.double(bfix),
                         as.integer(methInteg),
                         as.integer(nMC),
                         as.integer(dimMC),
                         as.double(seqMC),
                         as.integer(chol),
                         as.integer(estim0),
                         as.double(ll))

        out$ppi2 <- post$ppi2
        out$predRE <- post$predRE
        out$predRE_Y <- post$predRE_Y
        out$Yobs <- post$Yobs
        out$resid_m <- post$resid_m
        out$resid_ss <- post$resid_ss
        out$marker <- post$marker
        out$transfY <- post$transfY
        out$pred_m_g <- post$pred_m_g
        out$pred_ss_g <- post$pred_ss_g
    }
    
    
    ## creer best a partir de b et bfix
    best <- rep(NA,length(fix0))
    best[which(fix0==0)] <- out$best
    best[which(fix0==1)] <- bfix
    names(best) <- names_best
    out$best <- best
    NPM <- NPM+nfix
    
    
    

    ## mettre NA pour les variances et covariances non calculees et  0 pr les prm fixes
    if(length(posfix))
    {
        if(out$conv==3)
        {
            mr <- NPM-sum(pbH0)-length(posfix)
            Vr <- matrix(0,mr,mr)
            Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
            Vr <- t(Vr)
            Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
            V <- matrix(NA,NPM,NPM)
            V[setdiff(1:NPM,c(which(pbH0==1),posfix)),setdiff(1:NPM,c(which(pbH0==1),posfix))] <- Vr
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
            mr <- NPM-sum(pbH0)
            Vr <- matrix(0,mr,mr)
            Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
            Vr <- t(Vr)
            Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
            V <- matrix(NA,NPM,NPM)
            V[setdiff(1:NPM,which(pbH0==1)),setdiff(1:NPM,which(pbH0==1))] <- Vr
            V <- V[upper.tri(V,diag=TRUE)]
        }
        else
        {
            V <- out$V
        }
    }
    
    
    ## Creation du vecteur cholesky
    Cholesky <- rep(0,(nea0*(nea0+1)/2))
    CorrEA <- rep(0,(nea0*(nea0+1)/2))
    if(nvc>0){
    if(cholesky==TRUE)
    {
        if(idiag0==0 & nvc>0)
        {
            Cholesky[1:(nvc+1)] <- c(1,out$best[nef+1:nvc])
            ## Construction de la matrice U
            U <- matrix(0,nrow=nea0,ncol=nea0)
            U[upper.tri(U,diag=TRUE)] <- Cholesky[1:(nvc+1)]
            z <- t(U) %*% U
            out$best[nef+1:nvc] <- z[upper.tri(z,diag=TRUE)][-1]
        }
        if(idiag0==1 & nvc>0)
        {
            id <- 1:nea0
            indice <- rep(id+id*(id-1)/2)
            Cholesky[indice] <- c(1,out$best[nef+1:nvc])
            out$best[nef+1:nvc] <- out$best[nef+1:nvc]**2
        }
    }
    else
    {
        if(idiag0==0)
        {
            CorrEA <- c(1,out$best[nef+1:nvc])
            
            prmea <- matrix(0,nea0,nea0)
            prmea[upper.tri(prmea,diag=TRUE)] <- c(1,out$best[nef+1:nvc])
            prmea <- t(prmea)
            prmea[upper.tri(prmea,diag=TRUE)] <- c(1,out$best[nef+1:nvc])
            
            sea <- abs(diag(prmea))
            
            corr <- (exp(prmea)-1)/(exp(prmea)+1)
            diag(corr) <- 1
            covea <- sweep(corr,1,sea,"*")
            covea <- sweep(covea,2,sea,"*")

            out$best[nef+1:nvc] <- (covea[upper.tri(covea,diag=TRUE)])[-1]
        }
        else
        {
            id <- 1:nea0
            indice <- rep(id+id*(id-1)/2)
            CorrEA[indice] <- c(1,out$best[nef+1:nvc])
            out$best[nef+1:nvc] <- out$best[nef+1:nvc]^2
        }
    }
    }

    
    ##predictions
    predRE <- matrix(out$predRE,ncol=nea0,byrow=T)
    predRE <- data.frame(unique(IND),predRE)
    colnames(predRE) <- c(nom.subject,nom.X0[idea0!=0])

    if (nalea0!=0)
        {
            predRE_Y <- matrix(out$predRE_Y,ncol=ny0,byrow=TRUE)
            predRE_Y <- data.frame(unique(IND),predRE_Y)
            colnames(predRE_Y)  <- c(nom.subject,nomsY)
        }
    else
        {
            predRE_Y <- rep(NA,nalea0*ns0)
        }

###ppi
    if(ng0>1) {
        ppi <- matrix(out$ppi2,ncol=ng0,byrow=TRUE)
    }
    else {
        ppi <- matrix(rep(1,ns0),ncol=ng0)
    }
    
    chooseClass <- function(ppi)
    {
        res <- which.max(ppi)
        if(!length(res)) res <- NA
        return(res)
    }
    classif <- apply(ppi,1,chooseClass)
    ppi <- data.frame(unique(IND),classif,ppi)
    temp <- paste("prob",1:ng0,sep="")
    colnames(ppi) <- c(nom.subject,"class",temp)
    rownames(ppi) <- 1:ns0

###pred
    pred_m_g <- matrix(out$pred_m_g,nrow=nobs0)
    pred_ss_g <- matrix(out$pred_ss_g,nrow=nobs0)
    pred_m <- out$Yobs-out$resid_m
    pred_ss <- out$Yobs - out$resid_ss
    pred <- data.frame(IND,outcome,pred_m,out$resid_m,pred_ss,out$resid_ss,out$Yobs,pred_m_g,pred_ss_g)

    temp<-paste("pred_m",1:ng0,sep="")
    temp1<-paste("pred_ss",1:ng0,sep="")
    colnames(pred)<-c(nom.subject,"Yname","pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1)
    rownames(pred) <- NULL

    if(!is.null(var.time))
    {
        pred <- data.frame(IND,outcome,pred_m,out$resid_m,pred_ss,out$resid_ss,out$Yobs,pred_m_g,pred_ss_g,timeobs)
        colnames(pred)<-c(nom.subject,"Yname","pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1,var.time)
    }

###estimlink
    ysim <- matrix(out$marker,nsim,ny0)
    transfo <- matrix(out$transfY,nsim,ny0)
    if(any(idlink0==3))
    {
        sumntr <- 0
        for(k in 1:ny0)
        {
            if(idlink0[k]==3)
            {               
                seuils <- out$best[nef+nvc+nw+ncor0+ny0+nalea0+sumntr+1:ntrtot0[k]]
                if(ntrtot0[k]>1) seuils <- c(seuils[1], seuils[1] + cumsum(seuils[-1]^2))

                Lmin <- min(-2*abs(seuils[1]), min(transfo[1,]))
                Lmax <- max(2*abs(seuils[length(seuils)]), max(transfo[nsim,]))
               
                ysim_k <- rep(modalites[[k]], each=2)
                transfo_k <- c(Lmin, rep(seuils, each=2), Lmax)

                n <- min(nsim, length(ysim_k))
                ysim[,k] <- max(modalites[[k]])
                ysim[1:n,k] <- rep(ysim_k, lentgh.out=n)                
                transfo[,k] <- Lmax
                transfo[1:n,k] <- rep(transfo_k, lentgh.out=n)
            }
            sumntr <- sumntr + ntrtot0[k]
        }
    }
    estimlink <- as.vector(rbind(ysim,transfo))
    estimlink <- matrix(estimlink,nsim,2*ny0)
    colnames(estimlink) <- paste(c("","transf"),rep(nomsY, each=2),sep="")


    N <- NULL
    N[1] <- (ng0-1)*sum(idprob0)
    N[2] <- (ny0-1)*sum(idcontr0)
    N[3] <- (ng0-1)*sum(idprob0) + sum(idg0==1)-1 + ng0*sum(idg0==2) + (ny0-1)*sum(idcontr0)  #nef
    N[4] <- ifelse(idiag0==1,nea0,nea0*(nea0+1)/2)-1  #nvc
    N[5] <- (ng0-1)*nwg0
    N[6] <- nalea0
    N[7] <- ncor0
    N[8] <- ny0
    N[9] <- nobs0

    nom.X0[nom.X0=="(Intercept)"] <- "Intercept"
    
    ## levels = modalites des variables dans X0 (si facteurs)
    levelsdata <- vector("list", length(ttesLesVar))
    levelsfixed <- vector("list", length(ttesLesVar))
    levelsrandom <- vector("list", length(ttesLesVar))
    levelsmixture <- vector("list", length(ttesLesVar))
    levelsclassmb <- vector("list", length(ttesLesVar))
    names(levelsdata) <- ttesLesVar
    names(levelsfixed) <- ttesLesVar
    names(levelsrandom) <- ttesLesVar
    names(levelsmixture) <- ttesLesVar
    names(levelsclassmb) <- ttesLesVar
    for(v in ttesLesVar)
    {
        if(v == "intercept") next
        
        if(is.factor(data[,v]))
        {
            levelsdata[[v]] <- levels(data[,v])
        }
                                        
        if(length(grep(paste("factor\\(",v,"\\)",sep=""), fixed)))
        {
            levelsfixed[[v]] <- levels(as.factor(data[,v]))
        }

        if(length(grep(paste("factor\\(",v,"\\)",sep=""), random)))
        {
            levelsrandom[[v]] <- levels(as.factor(data[,v]))
        }

        if(length(grep(paste("factor\\(",v,"\\)",sep=""), mixture)))
        {
            levelsmixture[[v]] <- levels(as.factor(data[,v]))
        }
                                        
        if(length(grep(paste("factor\\(",v,"\\)",sep=""), classmb)))
        {
            levelsclassmb[[v]] <- levels(as.factor(data[,v]))
        }
    }
    
    levels <- list(levelsdata=levelsdata,
                   levelsfixed=levelsfixed,
                   levelsrandom=levelsrandom,
                   levelsmixture=levelsmixture,
                   levelsclassmb=levelsclassmb)
    
    cost <- proc.time()-ptm
    
    res <-list(ns=ns0,ng=ng0,idea0=idea0,idprob0=idprob0,idg0=idg0,idcontr0=idcontr0,
               idcor0=idcor0,loglik=out$loglik,best=out$best,V=V,gconv=out$gconv,conv=out$conv,
               call=cl,niter=out$niter,N=N,idiag=idiag0,pred=pred,pprob=ppi,predRE=predRE,
               predRE_Y=predRE_Y,Ynames=nomsY,Xnames=nom.X0,Xnames2=ttesLesVar,cholesky=Cholesky,
               estimlink=estimlink,epsY=epsY,linktype=idlink0,linknodes=zitr,nbnodes=nbnodes,nbmod=nbmod,modalites=modalites,
               na.action=nayk,AIC=2*(length(out$best)-length(posfix)-out$loglik),BIC=(length(out$best)-length(posfix))*log(ns0)-2*out$loglik,data=datareturn,
               wRandom=wRandom,b0Random=b0Random,runtime=cost[3],CorrEA=CorrEA,
               levels=levels, var.time=var.time)
    
    class(res) <-c("multlcmm")

    if(verbose==TRUE) cat("The program took", round(cost[3],2), "seconds \n")

    res
}


#' @rdname multlcmm
#' @export
mlcmm <- multlcmm


