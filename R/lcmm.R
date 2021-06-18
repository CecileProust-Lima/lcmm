#' Estimation of mixed-effect models and latent class mixed-effect models for
#' different types of outcomes (continuous Gaussian, continuous non-Gaussian or
#' ordinal)
#' 
#' This function fits mixed models and latent class mixed models for different
#' types of outcomes. It handles continuous longitudinal outcomes (Gaussian or
#' non-Gaussian) as well as bounded quantitative, discrete and ordinal
#' longitudinal outcomes.  The different types of outcomes are taken into
#' account using parameterized nonlinear link functions between the observed
#' outcome and the underlying latent process of interest it measures.  At the
#' latent process level, the model estimates a standard linear mixed model or a
#' latent class linear mixed model when heterogeneity in the population is
#' investigated (in the same way as in function \code{hlme}). It should be
#' noted that the program also works when no random-effect is included.
#' Parameters of the nonlinear link function and of the latent process mixed
#' model are estimated simultaneously using a maximum likelihood method.
#' 
#' 
#' A. THE PARAMETERIZED LINK FUNCTIONS
#' 
#' \code{lcmm} function estimates mixed models and latent class mixed models
#' for different types of outcomes by assuming a parameterized link function
#' for linking the outcome Y(t) with the underlying latent process L(t) it
#' measures. To fix the latent process dimension, we chose to constrain the
#' (first) intercept of the latent class mixed model at the latent process
#' level at 0 and the standard error of the gaussian error of measurement at 1.
#' These two parameters are replaced by additional parameters in the
#' parameterized link function :
#' 
#' 1. With the "linear" link function, 2 parameters are required that
#' correspond directly to the intercept and the standard error: (Y - b1)/b2 =
#' L(t).
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
#' k=2,...,n+2) so that b_k=(b_k^*)^2.
#' 
#' 4. With the "thresholds" link function for an ordinal outcome in levels
#' 0,...,C. A maximumn of C parameters are required for the following
#' transformation: Y(t)=c <=> b_c < L(t) <= b_{c+1} with b_0 = - infinity and
#' b_{C+1}=+infinity. The number of parameters is reduced if some levels do not
#' have any information. For example, if a level c is not observed in the
#' dataset, the corresponding threshold b_{c+1} is constrained to be the same
#' as the previous one b_{c}. The number of parameters in the link function is
#' reduced by 1.
#' 
#' To constraint the parameters to be increasing, except for the first
#' parameter b_1, the program estimates b_k^* (for k=2,...C) so that
#' b_{k}=b_{k-1}+(b_k^*)^2.
#' 
#' Details of these parameterized link functions can be found in the referred
#' papers.
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
#' estimated and only ng-1 parameters should be specified in \code{B}; (3) the
#' variance of each random-effect specified in \code{random} (including the
#' intercept) if \code{idiag=TRUE} and the inferior triangular
#' variance-covariance matrix of all the random-effects if \code{idiag=FALSE};
#' (4) only if \code{nwg=TRUE}, ng-1 parameters for class-specific proportional
#' coefficients for the variance covariance matrix of the random-effects; (5)
#' In contrast with hlme, due to identifiability purposes, the standard error
#' of the Gaussian error is not estimated (fixed at 1), and should not be
#' specified in \code{B}; (6) The parameters of the link function: 2 for
#' "linear", 4 for "beta", n+2 for "splines" with n nodes and the number of
#' levels minus one for "thresholds".
#' 
#' C. CAUTIONS REGARDING THE USE OF THE PROGRAM
#' 
#' Some caution should be made when using the program.  convergence criteria
#' are very strict as they are based on derivatives of the log-likelihood in
#' addition to the parameter and log-likelihood stability.  In some cases, the
#' program may not converge and reach the maximum number of iterations fixed at
#' 100.  In this case, the user should check that parameter estimates at the
#' last iteration are not on the boundaries of the parameter space.  If the
#' parameters are on the boundaries of the parameter space, the identifiability
#' of the model is critical. This may happen especially with splines parameters
#' that may be too close to 0 (lower boundary) or classmb parameters that are
#' too high or low (perfect classification). When identifiability of some
#' parameters is suspected, the program can be run again from the former
#' estimates by fixing the suspected parameters to their value with option
#' posfix. This usually solves the problem. An alternative is to remove the
#' parameters of the Beta of Splines link function from the inverse of the
#' Hessian with option partialH.  If not, the program should be run again with
#' other initial values, with a higher maximum number of iterations or less
#' strict convergence tolerances.
#' 
#' Specifically when investigating heterogeneity (that is with ng>1): (1) As
#' the log-likelihood of a latent class model can have multiple maxima, a
#' careful choice of the initial values is crucial for ensuring convergence
#' toward the global maximum.  The program can be run without entering the
#' vector of initial values (see point 2).  However, we recommend to
#' systematically enter initial values in \code{B} and try different sets of
#' initial values.  (2) The automatic choice of initial values we provide
#' requires the estimation of a preliminary linear mixed model. The user should
#' be aware that first, this preliminary analysis can take time for large
#' datatsets and second, that the generated initial values can be very not
#' likely and even may converge slowly to a local maximum.  This is the reason
#' why several alternatives exist. The vector of initial values can be directly
#' specified in \code{B} the initial values can be generated (automatically or
#' randomly) from a model with \code{ng=}. Finally, function \code{gridsearch}
#' performs an automatic grid search.
#' 
#' D. NUMERICAL INTEGRATION WITH THE THRESHOLD LINK FUNCTION
#' 
#' With exception for the threshold link function, maximum likelihood
#' estimation implemented in lcmm does not require any numerical integration
#' over the random-effects so that the estimation procedure is relatively fast.
#' See Proust et al. (2006) for more details on the estimation procedure.
#' 
#' However, with the threshold link function and when at least one
#' random-effect is specified, a numerical integration over the random-effects
#' distribution is required in each computation of the individual contribution
#' to the likelihood which complicates greatly the estimation procedure. For
#' the moment, we do not allow any option regarding the numerical integration
#' technics used.  1. When a single random-effect is specified, we use a
#' standard non-adaptive Gaussian quadrature with 30 points. 2. When at least
#' two random-effects are specified, we use a multivariate non-adaptive
#' Gaussian quadrature implemented by Genz (1996) in HRMSYM Fortran subroutine.
#' 
#' Further developments should allow for adaptive technics and more options
#' regarding the numerical integration technic.
#' 
#' E. POSTERIOR DISCRETE LIKELIHOOD
#' 
#' Models involving nonlinear continuous link functions assume the continuous
#' data while the model with a threshold model assumes discrete data. As a
#' consequence, comparing likelihoods or criteria based on the likelihood (as
#' AIC) for these models is not possible as the former are based on a Lebesgue
#' measure and the latter on a counting measure. To make the comparison
#' possible, we compute the posterior discrete likelihood for all the models
#' with a nonlinear continuous link function. This posterior likelihood
#' considers the data as discrete; it is computed at the MLE (maximum
#' likelihood estimates) using the counting measure so that models with
#' threshold or continuous link functions become comparable. Further details
#' can be found in Proust-Lima, Amieva, Jacqmin-Gadda (2012).
#' 
#' In addition to the Akaike information criterion based on the discrete
#' posterior likelihood, we also compute a universal approximate
#' cross-validation criterion to compare models based on a different measure.
#' See Commenges, Proust-Lima, Samieri, Liquet (2015) for further details.
#' 
#' @param fixed a two-sided linear formula object for specifying the
#' fixed-effects in the linear mixed model at the latent process level. The
#' response outcome is on the left of \code{~} and the covariates are separated
#' by \code{+} on the right of the \code{~}.  Fo identifiability purposes, the
#' intercept specified by default should not be removed by a \code{-1}.
#' @param mixture a one-sided formula object for the class-specific fixed
#' effects in the latent process mixed model (to specify only for a number of
#' latent classes greater than 1).  Among the list of covariates included in
#' \code{fixed}, the covariates with class-specific regression parameters are
#' entered in \code{mixture} separated by \code{+}.  By default, an intercept
#' is included. If no intercept, \code{-1} should be the first term included.
#' @param random an optional one-sided formula for the random-effects in the
#' latent process mixed model. Covariates with a random-effect are separated by
#' \code{+}.  By default, no random effect is included.
#' @param subject name of the covariate representing the grouping structure.
#' @param classmb an optional one-sided formula describing the covariates in
#' the class-membership multinomial logistic model. Covariates included are
#' separated by \code{+}.  No intercept should be included in this formula.
#' @param ng number of latent classes considered. If \code{ng=1} no
#' \code{mixture} nor \code{classmb} should be specified. If \code{ng>1},
#' \code{mixture} is required.
#' @param idiag optional logical for the variance-covariance structure of the
#' random-effects. If \code{FALSE}, a non structured matrix of
#' variance-covariance is considered (by default).  If \code{TRUE} a diagonal
#' matrix of variance-covariance is considered.
#' @param nwg optional logical of class-specific variance-covariance of the
#' random-effects. If \code{FALSE} the variance-covariance matrix is common
#' over latent classes (by default). If \code{TRUE} a class-specific
#' proportional parameter multiplies the variance-covariance matrix in each
#' class (the proportional parameter in the last latent class equals 1 to
#' ensure identifiability).
#' @param link optional family of link functions to estimate. By default,
#' "linear" option specifies a linear link function leading to a standard
#' linear mixed model (homogeneous or heterogeneous as estimated in
#' \code{hlme}).  Other possibilities include "beta" for estimating a link
#' function from the family of Beta cumulative distribution functions,
#' "thresholds" for using a threshold model to describe the correspondence
#' between each level of an ordinal outcome and the underlying latent process,
#' and "Splines" for approximating the link function by I-splines. For this
#' latter case, the number of nodes and the nodes location should be also
#' specified. The number of nodes is first entered followed by \code{-}, then
#' the location is specified with "equi", "quant" or "manual" for respectively
#' equidistant nodes, nodes at quantiles of the marker distribution or interior
#' nodes entered manually in argument \code{intnodes}. It is followed by
#' \code{-} and finally "splines" is indicated. For example, "7-equi-splines"
#' means I-splines with 7 equidistant nodes, "6-quant-splines" means I-splines
#' with 6 nodes located at the quantiles of the marker distribution and
#' "9-manual-splines" means I-splines with 9 nodes, the vector of 7 interior
#' nodes being entered in the argument \code{intnodes}.
#' @param intnodes optional vector of interior nodes. This argument is only
#' required for a I-splines link function with nodes entered manually.
#' @param epsY optional definite positive real used to rescale the marker in
#' (0,1) when the beta link function is used. By default, epsY=0.5.
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
#' values.  (3) when ng>1, a lcmm object is entered. It should correspond to
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
#' @param maxiter optional maximum number of iterations for the Marquardt
#' iterative algorithm. By default, maxiter=100.
#' @param nsim number of points used to plot the estimated link function. By
#' default, nsim=100.
#' @param prior name of the covariate containing the prior on the latent class
#' membership. The covariate should be an integer with values in 0,1,...,ng.
#' When there is no prior, the value should be 0. When there is a prior for the
#' subject, the value should be the number of the latent class (in 1,...,ng).
#' @param range optional vector indicating the range of the outcome (that is
#' the minimum and maximum). By default, the range is defined according to the
#' minimum and maximum observed values of the outcome. The option should be
#' used only for Beta and Splines transformations.
#' @param subset optional vector giving the subset of observations in
#' \code{data} to use. By default, all lines.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @param posfix Optional vector specifying the indices in vector B of the
#' parameters that should not be estimated. Default to NULL, all parameters are
#' estimated.
#' @param partialH optional logical for Beta or Splines link functions only.
#' Indicates whether the parameters of the link functions can be dropped from
#' the Hessian matrix to define convergence criteria.
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
#' Marquardt iterations} \item{dataset}{dataset} \item{N}{internal information
#' used in related functions} \item{idiag}{internal information used in related
#' functions} \item{pred}{table of individual predictions and residuals in the
#' underlying latent process scale; it includes marginal predictions (pred_m),
#' marginal residuals (resid_m), subject-specific predictions (pred_ss) and
#' subject-specific residuals (resid_ss) averaged over classes, the transformed
#' observations in the latent process scale (obs) and finally the
#' class-specific marginal and subject-specific predictions (with the number of
#' the latent class: pred_m_1,pred_m_2,...,pred_ss_1,pred_ss_2,...). This
#' output is not available yet when specifying a thresholds transformation.}
#' \item{pprob}{table of posterior classification and posterior individual
#' class-membership probabilities} \item{Xnames}{list of covariates included in
#' the model} \item{predRE}{table containing individual predictions of the
#' random-effects : a column per random-effect, a line per subject. This output
#' is not available yet when specifying a thresholds transformation.}
#' \item{cholesky}{vector containing the estimates of the Cholesky transformed
#' parameters of the variance-covariance matrix of the random-effects}
#' \item{estimlink}{table containing the simulated values of the marker and
#' corresponding estimated link function} \item{epsY}{definite positive real
#' used to rescale the marker in (0,1) when the beta link function is used. By
#' default, epsY=0.5.} \item{linktype}{indicator of link function type: 0 for
#' linear, 1 for beta, 2 for splines and 3 for thresholds}
#' \item{linknodes}{vector of nodes useful only for the 'splines' link
#' function}\item{data}{the original data set (if returndata is TRUE)}
#' 
#' @author Cecile Proust-Lima, Amadou Diakite, Benoit Liquet and Viviane
#' Philipps
#' 
#' \email{cecile.proust-lima@@inserm.fr}
#' @seealso
#' 
#' \code{\link{postprob}}, \code{\link{plot.lcmm}}, \code{\link{plot.predict}},
#' \code{\link{hlme}}
#' @references
#' 
#' Proust-Lima C, Philipps V, Liquet B (2017). Estimation of Extended Mixed 
#' Models Using Latent Classes and Latent Processes: The R Package lcmm. 
#' Journal of Statistical Software, 78(2), 1-56. doi:10.18637/jss.v078.i02
#' 
#' Genz and Keister (1996). Fully symmetric interpolatory rules for multiple
#' integrals over infinite regions with gaussian weight. Journal of
#' Computational and Applied Mathematics 71: 299-309.
#' 
#' Proust and Jacqmin-Gadda (2005). Estimation of linear mixed models with a
#' mixture of distribution for the random-effects. Comput Methods Programs
#' Biomed 78: 165-73.
#' 
#' Proust, Jacqmin-Gadda, Taylor, Ganiayre, and Commenges (2006). A nonlinear
#' model with latent process for cognitive evolution using multivariate
#' longitudinal data. Biometrics 62: 1014-24.
#' 
#' Proust-Lima, Dartigues and Jacqmin-Gadda (2011). Misuse of the linear mixed
#' model when evaluating risk factors of cognitive decline. Amer J Epidemiol
#' 174(9): 1077-88.
#' 
#' Proust-Lima, Amieva and Jacqmin-Gadda (2013). Analysis of multivariate mixed
#' longitudinal data : a flexible latent process approach, British Journal of
#' Mathematical and Statistical Psychology 66(3): 470-87.
#' 
#' Commenges, Proust-Lima, Samieri, Liquet (2015). A universal approximate
#' cross-validation criterion for regular risk functions. Int J Biostat. 2015
#' May;11(1):51-67
#' @examples
#' 
#' \dontrun{
#' #### Estimation of homogeneous mixed models with different assumed link
#' #### functions, a quadratic mean trajectory for the latent process and 
#' #### correlated random intercept and slope (the random quadratic slope 
#' #### was removed as it did not improve the fit of the data).
#' #### -- comparison of linear, Beta and 3 different splines link functions --
#' # linear link function
#' m10<-lcmm(Ydep2~Time+I(Time^2),random=~Time,subject='ID',ng=1,
#' data=data_lcmm,link="linear")
#' summary(m10)
#' # Beta link function
#' m11<-lcmm(Ydep2~Time+I(Time^2),random=~Time,subject='ID',ng=1,
#' data=data_lcmm,link="beta")
#' summary(m11)
#' plot(m11,which="linkfunction",bty="l")
#' # I-splines with 3 equidistant nodes
#' m12<-lcmm(Ydep2~Time+I(Time^2),random=~Time,subject='ID',ng=1,
#' data=data_lcmm,link="3-equi-splines")
#' summary(m12)
#' # I-splines with 5 nodes at quantiles
#' m13<-lcmm(Ydep2~Time+I(Time^2),random=~Time,subject='ID',ng=1,
#' data=data_lcmm,link="5-quant-splines")
#' summary(m13)
#' # I-splines with 5 nodes, and interior nodes entered manually
#' m14<-lcmm(Ydep2~Time+I(Time^2),random=~Time,subject='ID',ng=1,
#' data=data_lcmm,link="5-manual-splines",intnodes=c(10,20,25))
#' summary(m14)
#' plot(m14,which="linkfunction",bty="l")
#' 
#' 
#' # Thresholds
#' # Especially for the threshold link function, we recommend to estimate 
#' # models with increasing complexity and use estimates of previous ones 
#' # to specify plausible initial values (we remind that estimation of
#' # models with threshold link function involves a computationally demanding 
#' # numerical integration -here of size 3)
#' m15<-lcmm(Ydep2~Time+I(Time^2),random=~Time,subject='ID',ng=1
#' ,data=data_lcmm,link="thresholds",maxiter=100,
#' B=c(-0.8379, -0.1103,  0.3832,  0.3788 , 0.4524, -7.3180,  0.5917,  0.7364,
#'  0.6530, 0.4038,  0.4290,  0.6099,  0.6014 , 0.5354 , 0.5029 , 0.5463,
#'  0.5310 , 0.5352, 0.6498,  0.6653,  0.5851,  0.6525,  0.6701 , 0.6670 ,
#'  0.6767 , 0.7394 , 0.7426, 0.7153,  0.7702,  0.6421))
#' summary(m15)
#' plot(m15,which="linkfunction",bty="l")
#' 
#' #### Plot of estimated different link functions:
#' #### (applicable for models that only differ in the "link function" used. 
#' ####  Otherwise, the latent process scale is different and a rescaling
#' ####  is necessary)
#' plot(m10,which="linkfunction",col=1,xlab="latent process",ylab="marker",
#' bty="l",xlim=c(-10,5),legend=NULL)
#' plot(m11,which="linkfunction",add=TRUE,col=2,legend=NULL)
#' plot(m12,which="linkfunction",add=TRUE,col=3,legend=NULL)
#' plot(m13,which="linkfunction",add=TRUE,col=4,legend=NULL)
#' plot(m14,which="linkfunction",add=TRUE,col=5,legend=NULL)
#' plot(m15,which="linkfunction",add=TRUE,col=6,legend=NULL)
#' legend(x="bottomright",legend=c("linear","beta","spl_3e","spl_5q","spl_5m","thresholds"),
#' col=1:6,lty=1,inset=.02,box.lty=0)
#' 
#' #### Estimation of 2-latent class mixed models with different assumed link 
#' #### functions with individual and class specific linear trend
#' #### for illustration, only default initial values where used but other
#' #### sets of initial values should also be tried to ensure convergence 
#' #### towards the golbal maximum
#' # Linear link function
#' m20<-lcmm(Ydep2~Time,random=~Time,subject='ID',mixture=~Time,ng=2,
#' idiag=TRUE,data=data_lcmm,link="linear",B=c(-0.98,0.79,-2.09,
#' -0.81,0.19,0.55,24.49,2.24))
#' summary(m20)
#' postprob(m20)
#' # Beta link function
#' m21<-lcmm(Ydep2~Time,random=~Time,subject='ID',mixture=~Time,ng=2,
#' idiag=TRUE,data=data_lcmm,link="beta",B=c(-0.1,-0.56,-0.4,-1.77,
#' 0.53,0.14,0.6,-0.83,0.73,0.09))
#' summary(m21)
#' postprob(m21)
#' # I-splines link function (and 5 nodes at quantiles)
#' m22<-lcmm(Ydep2~Time,random=~Time,subject='ID',mixture=~Time,ng=2,
#' idiag=TRUE,data=data_lcmm,link="5-quant-splines",B=c(0.12,0.63,
#' -1.76,-0.39,0.51,0.13,-7.37,1.05,1.28,1.96,1.3,0.93,1.05))
#' summary(m22)
#' postprob(m22)
#' 
#' data <- data_lcmm[data_lcmm$ID==193,]
#' plot(predictL(m22,var.time="Time",newdata=data,bty="l")
#' 
#' }    
#' 
#' @export
#' 
#' 
#' 
lcmm <- function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,link="linear",intnodes=NULL,epsY=0.5,cor=NULL,data,B,convB=0.0001,convL=0.0001,convG=0.0001,maxiter=100,nsim=100,prior,range=NULL,subset=NULL,na.action=1,posfix=NULL,partialH=FALSE,verbose=TRUE,returndata=FALSE)
{

mm <- match.call()
if(missing(fixed)) stop("The argument Fixed must be specified in any model")
if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
if(nrow(data)==0) stop("Data should not be empty")

### transformer l'argument cor en character
cor.type <- mm$cor[1]
cor.time <- mm$cor[2]
cor.char <- paste(cor.type,cor.time,sep="-") 
if (all.equal(cor.char,character(0))!=TRUE)
{
 if (link=="thresholds") stop("The argument cor is only available with linear, beta or splines link")
}
else
{
  cor.char <- NULL
}  
  
if (!is.null(cor.char))
{
 if(!(strsplit(cor.char,"-")[[1]][2] %in% colnames(data))) stop("Unable to find time variable from argument cor in data")
 else { cor.var.time <- strsplit(cor.char,"-")[[1]][2] }
}  
### fin cor


## gestion de B=random(mod)

        if(length(mm$B)==2)
            {
                if(class(eval(mm$B[[2]]))!="lcmm") stop("The model specified in B should be of class lcmm")
                if(as.character(mm$B[1])!="random") stop("Please use random() to specify random initial values")
                B <- eval(mm$B[[2]])   
                B$Brandom <- TRUE
                
                if(length(posfix)) stop("Argument posfix is not compatible with random intial values")
            }





if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

if(na.action==1){
	na.action=na.omit
}else{
	na.action=na.fail
}
#cat("ide :")
#cat(ide0,"\n")
#cat(zitr,"\n")
#7/05/2012
### Traitement des donnees manquantes
# fixed
if(missing(fixed)) stop("The argument Fixed must be specified in any model")
if(class(fixed)!="formula") stop("The argument fixed must be a formula")
m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
m$formula <- terms(fixed)
m$na.action <- na.action
m[[1]] <- as.name("model.frame")	
m <- eval(m, sys.parent()) 
na.fixed <- attr(m,"na.action")

# mixture
if(!missing(mixture)){
	if(class(mixture)=="formula"){	
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
	m$formula <- terms(mixture)
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
	na.mixture <- attr(m,"na.action")
	}	
}else{
	na.mixture <- NULL
}

# random
if(!missing(random)){
	if(class(random)=="formula"){	
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
	m$formula <- terms(random)
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.random <- attr(m,"na.action")
	}
}else{
	na.random <- NULL
}

# classmb
if(!missing(classmb)){ 
	if(class(classmb)=="formula"){	
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]	
	m$formula <- terms(classmb)
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.classmb <- attr(m,"na.action")
	}
}else{
	na.classmb <- NULL
}
 
#cor     
if(!is.null(cor.char))
{
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
	m$formula <- as.formula(paste(cor.var.time,1,sep="~"))
	m$na.action <- na.action
  m[[1]] <- as.name("model.frame")
  m <- eval(m,sys.parent())    
  na.cor <- attr(m,"na.action") 	
}
else { na.cor <- NULL }
 
 
	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.cor))
#7/05/2012


attr.fixed <- attributes(terms(fixed))
depvar <- as.character(attr.fixed$variables[2])
if(!isTRUE(all.equal(as.character(mm$subset),character(0))))
    {
        cc <- mm
        cc <- cc[c(1,which(names(mm)=="subset"))]
        cc[[1]] <- as.name("model.frame")
        cc$formula <- formula(paste("~",depvar))
        cc$data <- data
        cc$na.action <- na.pass
        ysubset <- eval(cc)
    }
else
    {
        ysubset <- data[,depvar,drop=FALSE]
    }

if(!is.null(na.action))
{
    Y0 <- ysubset[-na.action,]
}
else
{
    Y0 <- ysubset
}
Y0 <- unlist(Y0)

minY0 <- min(Y0)
maxY0 <- max(Y0)
if ((!missing(range)) & length(range)==2)
{
 if(minY0<range[1]|maxY0>range[2]) stop("The range specified do not cover the entire range of the data")
 if (minY0>range[1]|maxY0<range[2])
 {
  minY0 <- range[1]
  maxY0 <- range[2]
 } 
}
else
    {
        min2 <- round(minY0,3)
        if(minY0<min2) min2 <- min2-0.001
        minY0 <- min2

        max2 <- round(maxY0,3)
        if(maxY0>max2) max2 <- max2+0.001
        maxY0 <- max2
    }

if(all.equal((maxY0-minY0),0) == T){
	stop("All the values of the dependent variable are the same. No estimation can be performed in that case.")
}
#if((any(is.na(Y0))==TRUE)){
#	stop("The dependent variable should not contain any missing value")
#}

if(length(grep("-",unlist(strsplit(link,split="")))) > 2){
	stop("Please check and revise the 'link' argument according to the format given in the help.")
}



################################# pas de separateur "-", uniquement pour les splines


if(all.equal(length(grep("-",unlist(strsplit(link,split="")))),0)==T){

	if (!(link %in% c("linear","beta","thresholds","splines"))){
		stop("The only available link functions in lcmm are 'linear', 'beta', 'splines' and 'thresholds' functions.")
	}else{
		nbzitr0 <- switch(link,"linear"=2,"beta"=2,"splines"=5,"thresholds"=2)
		idlink0 <- switch(link,"linear"=0,"beta"=1,"splines"=2,"thresholds"=3)
		ntrtot0 <- switch(link,"linear"=2,"beta"=4, "splines"= (nbzitr0 + 2), "thresholds"= as.integer(maxY0-minY0)) 
		if(all.equal(link,"splines")==T | all.equal(link,"Splines")==T){
			link <- "splines"
			type <- "equi"	
		}
		
	} 	
	if ((link %in% c("thresholds"))){

		minY0 <- min(Y0)
		maxY0 <- max(Y0)
		############################## PARTIE A COMPLETER POUR ORDINAL #####################
		if(!(all.equal(minY0,as.integer(minY0))==T) | !(all.equal(maxY0,as.integer(maxY0))==T)|!all(Y0 %in% minY0:maxY0)){
			stop("With the threshold link function, the longitudinal outcome must be discrete")
		}

		IND <- sort(unique(Y0))
		IND <- IND[1:(length(IND)-1)]-minY0+1
		ide0 <- rep(0,as.integer(maxY0-minY0))
		ide0[IND] <- 1


	
	
	#######################################################################
	}
	zitr <- rep(0,nbzitr0)
	zitr[1] <- minY0
	zitr[nbzitr0] <- maxY0
	
}


################################# Avec un seul separateur "-", uniquement pour les splines
if(all.equal(length(grep("-",unlist(strsplit(link,split="")))),1)==T){
	stop("The number and location of the nodes only apply for the splines link function. For 'thresholds', 'linear' and 'beta' links, no nodes are required. For the splines link function, both the number and the type of location for the nodes should be specified (ex: 5-manual-splines for 5 manual nodes)")
}


if(all.equal(length(grep("-",unlist(strsplit(link,split="")))),2)==T){

	if(any(unlist(strsplit(link,"-")) %in% c("linear","beta","thresholds"))){
		stop("The number and location of the nodes only apply for the 'splines' link function. For 'thresholds', 'linear' and 'beta' links, no nodes are required.")
	}
	
### Verification de l'ordre de replissage du link
	if(!(unlist(strsplit(link,"-"))[3] %in% c("splines"))){
		stop("When defining a link function using splines, the third part of the 'link' argument must include only 'splines' (ex: 5-equi-splines for 5 equidistant nodes)")
	}
	
	if(!(unlist(strsplit(link,"-"))[2] %in% c("equi","manual","quant"))){
		stop("When defining a link function using splines, the second part of the 'link' argument must include only 'equi' 'manual' 'quant' for equidistant manual or quantile nodes")
	}
	
	nbzitr0 <- as.integer(unlist(strsplit(link,"-"))[1])
	ntrtot0 <- nbzitr0 + 2
	idlink0 <- 2 
	type <- unlist(strsplit(link,"-"))[2]	   
	link <- "splines" 
	if((nbzitr0-2) < 0) stop("At least 2 nodes should be specified for the splines link function.")
}


if (all.equal(idlink0,2)==T){

	zitr <- rep(0,nbzitr0)
	zitr[1] <- minY0
	zitr[nbzitr0] <- maxY0
	
	if(all.equal("manual",type)==T){
		if (is.null(intnodes)){
		stop("If 'manual' option is specified for the splines link function, intnodes argument should include the list of interior nodes")
		}else{            
		if(!(all.equal(length(intnodes),(nbzitr0-2))==T)==T) stop("Intnodes does not include the correct number of interior nodes")     
		intnodes <- sort(intnodes)
		if(intnodes[1]<=zitr[1]|intnodes[nbzitr0-2]>=zitr[nbzitr0])stop("Intnodes are not inside the boundaries of the marker")     
		zitr[2:(nbzitr0-1)] <- intnodes[1:(nbzitr0-2)]
		}
	}    
	if(all.equal("quant",type)==T){
		pas <-c(1:(nbzitr0-2))/(nbzitr0-1) 
		zitr[2:(nbzitr0-1)] <- quantile(sort(Y0),probs=pas)
		if(length(unique(zitr[1:nbzitr0]))!=nbzitr0) stop("The link function can not be estimated since some nodes are equal; Please try to reduce the number of nodes or use manual location.")
	}        	       
	if(all.equal("equi",type)==T){
		pas=as.double(maxY0-minY0)/as.double(nbzitr0-1)
		for(i in 2:(nbzitr0-1)){
		zitr[i] <- zitr[i-1]+pas
		}
	}

 #verifier s'il y a des obs entre les noeuds
 hcounts <- hist(Y0,breaks=zitr,plot=FALSE,include.lowest=TRUE,right=TRUE)$counts
 if(any(hcounts==0) & (maxiter != 0)) stop("Link function can not be estimated since some intervals defined by the nodes do not contain any observation.")  
}      


if (idlink0==1) {
if (epsY<=0) {
epsY <- 0.5
cat("Argument 'epsY' should be a definite positive real. It is changed to the default value of 0.5. \n") 
}
}


Ydiscrete <- 1
if (idlink0!=3) {
        if(!(all.equal(minY0,as.integer(minY0))==T) | !(all.equal(maxY0,as.integer(maxY0))==T)|!all(Y0 %in% minY0:maxY0)){
		Ydiscrete <- 0
	}

}






### test partialH que pour beta ou splines
if(!(idlink0 %in% c(1,2)) & partialH) stop("No partial Hessian can be define")

           






 

link <- as.character(link)
### appel des differents modeles selon la valeur de l'argument link
result <- switch(link
,"linear"=.Contlcmm(fixed=fixed,mixture=mixture,random=random,subject=subject,classmb=classmb,ng=ng,idiag=idiag,nwg=nwg,cor=cor.char,data=data,B=B,convB=convB,convL=convL,convG=convG,prior=prior,maxiter=maxiter,epsY=epsY,idlink0=idlink0,ntrtot0=ntrtot0,nbzitr0=nbzitr0,zitr=zitr,nsim=nsim,call=mm,Ydiscrete,subset=subset,na.action,posfix=posfix,partialH=partialH,verbose=verbose,returndata=returndata)

,"beta"=.Contlcmm(fixed=fixed,mixture=mixture,random=random,subject=subject,classmb=classmb,ng=ng,idiag=idiag,nwg=nwg,cor=cor.char,data=data,B=B,convB=convB,convL=convL,convG=convG,prior=prior,maxiter=maxiter,epsY=epsY,idlink0=idlink0,ntrtot0=ntrtot0,nbzitr0=nbzitr0,zitr=zitr,nsim=nsim,call=mm,Ydiscrete,subset=subset,na.action,posfix=posfix,partialH=partialH,verbose=verbose,returndata=returndata)

,"splines"=.Contlcmm(fixed=fixed,mixture=mixture,random=random,subject=subject,classmb=classmb,ng=ng,idiag=idiag,nwg=nwg,cor=cor.char,data=data,B=B,convB=convB,convL=convL,convG=convG,prior=prior,maxiter=maxiter,epsY=epsY,idlink0=idlink0,ntrtot0=ntrtot0,nbzitr0=nbzitr0,zitr=zitr,nsim=nsim,call=mm,Ydiscrete,subset=subset,na.action,posfix=posfix,partialH=partialH,verbose=verbose,returndata=returndata)
                 
,"thresholds"=.Ordlcmm(fixed=fixed,mixture=mixture,random=random,subject=subject,classmb=classmb,ng=ng,idiag=idiag,nwg=nwg,data=data,B=B,convB=convB,convL=convL,convG=convG,prior=prior,maxiter=maxiter,zitr=zitr,ide=ide0,call=mm,Ydiscrete,subset=subset,na.action=na.action,posfix=posfix,verbose=verbose,returndata=returndata))
  
return(result)
}




