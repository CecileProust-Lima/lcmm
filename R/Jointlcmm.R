#' Estimation of joint latent class models for longitudinal and time-to-event
#' data
#' 
#' This function fits joint latent class mixed models for a longitudinal
#' outcome and a right-censored (possibly left-truncated) time-to-event. The
#' function handles competing risks and Gaussian or non Gaussian (curvilinear)
#' longitudinal outcomes. For curvilinear longitudinal outcomes, normalizing
#' continuous functions (splines or Beta CDF) can be specified as in
#' \code{lcmm}.
#' 
#' 
#' A. BASELINE RISK FUNCTIONS
#' 
#' For the baseline risk functions, the following parameterizations were
#' considered. Be careful, parametrisations changed in lcmm_V1.5:
#' 
#' 1. With the "Weibull" function: 2 parameters are necessary w_1 and w_2 so
#' that the baseline risk function a_0(t) = w_1^2*w_2^2*(w_1^2*t)^(w_2^2-1) if
#' logscale=FALSE and a_0(t) = exp(w_1)*exp(w_2)(t)^(exp(w_2)-1) if
#' logscale=TRUE.
#' 
#' 2. with the "piecewise" step function and nz nodes (y_1,...y_nz), nz-1
#' parameters are necesssary p_1,...p_nz-1 so that the baseline risk function
#' a_0(t) = p_j^2 for y_j < t =< y_j+1 if logscale=FALSE and a_0(t) = exp(p_j)
#' for y_j < t =< y_j+1 if logscale=TRUE.
#' 
#' 3. with the "splines" function and nz nodes (y_1,...y_nz), nz+2 parameters
#' are necessary s_1,...s_nz+2 so that the baseline risk function a_0(t) =
#' sum_j s_j^2 M_j(t) if logscale=FALSE and a_0(t) = sum_j exp(s_j) M_j(t) if
#' logscale=TRUE where {M_j} is the basis of cubic M-splines.
#' 
#' Two parametrizations of the baseline risk function are proposed
#' (logscale=TRUE or FALSE) because in some cases, especially when the
#' instantaneous risks are very close to 0, some convergence problems may
#' appear with one parameterization or the other. As a consequence, we
#' recommend to try the alternative parameterization (changing logscale option)
#' when a joint latent class model does not converge (maximum number of
#' iterations reached) where as convergence criteria based on the parameters
#' and likelihood are small.
#' 
#' B. THE VECTOR OF PARAMETERS B
#' 
#' The parameters in the vector of initial values \code{B} or in the vector of
#' maximum likelihood estimates \code{best} are included in the following
#' order: (1) ng-1 parameters are required for intercepts in the latent class
#' membership model, and if covariates are included in \code{classmb}, ng-1
#' parameters should be entered for each one; (2) parameters for the baseline
#' risk function: 2 parameters for each Weibull, nz-1 for each piecewise
#' constant risk and nz+2 for each splines risk; this number should be
#' multiplied by ng if specific hazard is specified; otherwise, ng-1 additional
#' proportional effects are expected if PH hazard is specified; otherwise
#' nothing is added if common hazard is specified. In the presence of competing
#' events, the number of parameters should be adapted to the number of causes
#' of event; (3) for all covariates in \code{survival}, ng parameters are
#' required if the covariate is inside a \code{mixture()}, otherwise 1
#' parameter is required. Covariates parameters should be included in the same
#' order as in \code{survival}. In the presence of cause-specific effects, the
#' number of parameters should be multiplied by the number of causes; (4) for
#' all covariates in \code{fixed}, one parameter is required if the covariate
#' is not in \code{mixture}, ng parameters are required if the covariate is
#' also in \code{mixture}. Parameters should be included in the same order as
#' in \code{fixed}; (5) the variance of each random-effect specified in
#' \code{random} (including the intercept) if \code{idiag=TRUE} and the
#' inferior triangular variance-covariance matrix of all the random-effects if
#' \code{idiag=FALSE}; (6) only if \code{nwg=TRUE}, ng-1 parameters for
#' class-specific proportional coefficients for the variance covariance matrix
#' of the random-effects; (7) the variance of the residual error.
#' 
#' C. CAUTION
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
#' (2) The automatic choice of initial values that we provide requires the
#' estimation of a preliminary linear mixed model. The user should be aware
#' that first, this preliminary analysis can take time for large datatsets and
#' second, that the generated initial values can be very not likely and even
#' may converge slowly to a local maximum.  This is a reason why several
#' alternatives exist. The vector of initial values can be directly specified
#' in \code{B} the initial values can be generated (automatically or randomly)
#' from a model with \code{ng=}. Finally, function \code{gridsearch} performs
#' an automatic grid search.
#' 
#' (3) Convergence criteria are very strict as they are based on derivatives of
#' the log-likelihood in addition to the parameter and log-likelihood
#' stability.  In some cases, the program may not converge and reach the
#' maximum number of iterations fixed at 150.  In this case, the user should
#' check that parameter estimates at the last iteration are not on the
#' boundaries of the parameter space.  If the parameters are on the boundaries
#' of the parameter space, the identifiability of the model is critical. This
#' may happen especially when baseline risk functions involve splines (value
#' close to the lower boundary - 0 with logscale=F -infinity with logscale=F)
#' or classmb parameters that are too high or low (perfect classification) or
#' linkfunction parameters. When identifiability of some parameters is
#' suspected, the program can be run again from the former estimates by fixing
#' the suspected parameters to their value with option posfix. This usually
#' solves the problem. An alternative is to remove the parameters of the Beta
#' of Splines link function from the inverse of the Hessian with option
#' partialH.  If not, the program should be run again with other initial
#' values.  Some problems of convergence may happen when the instantaneous
#' risks of event are very low and "piecewise" or "splines" baseline risk
#' functions are specified. In this case, changing the parameterization of the
#' baseline risk functions with option logscale is recommended (see paragraph A
#' for details).
#' 
#' @aliases Jointlcmm jlcmm
#' @param fixed two-sided linear formula object for the fixed-effects in the
#' linear mixed model. The response outcome is on the left of \code{~} and the
#' covariates are separated by \code{+} on the right of the \code{~}.  By
#' default, an intercept is included. If no intercept, \code{-1} should be the
#' first term included on the right of \code{~}.
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
#' (called subject identifier) specified with ''.
#' @param classmb optional one-sided formula describing the covariates in the
#' class-membership multinomial logistic model. Covariates included are
#' separated by \code{+}. No intercept should be included in this formula.
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
#' @param survival two-sided formula object. The left side of the formula
#' corresponds to a \code{surv()} object of type "counting" for right-censored
#' and left-truncated data (example: \code{Surv(Time,EntryTime,Indicator)}) or
#' of type "right" for right-censored data (example:
#' \code{Surv(Time,Indicator)}). Multiple causes of event can be considered in
#' the Indicator (0 for censored, k for cause k of event).  The right side of
#' the formula specifies the names of covariates to include in the survival
#' model with \code{mixture()} when the effect is class-specific (example:
#' \code{Surv(Time,Indicator) ~} \code{ X1 + mixture(X2)} for a class-common
#' effect of X1 and a class-specific effect of X2). In the presence of
#' competing events, covariate effects are common by default. Code
#' \code{cause(X3)} specifies a cause-specific covariate effect for X3 on each
#' cause of event while \code{cause1(X3)} (or \code{cause2(X3)}, ...) specifies
#' a cause-specific effect of X3 on the first (or second, ...) cause only.
#' @param hazard optional family of hazard function assumed for the survival
#' model. By default, "Weibull" specifies a Weibull baseline risk function.
#' Other possibilities are "piecewise" for a piecewise constant risk function
#' or "splines" for a cubic M-splines baseline risk function. For these two
#' latter families, the number of nodes and the location of the nodes should be
#' specified as well, separated by \code{-}. The number of nodes is entered
#' first followed by \code{-}, then the location is specified with "equi",
#' "quant" or "manual" for respectively equidistant nodes, nodes at quantiles
#' of the times of event distribution or interior nodes entered manually in
#' argument \code{hazardnodes}. It is followed by \code{-} and finally
#' "piecewise" or "splines" indicates the family of baseline risk function
#' considered. Examples include "5-equi-splines" for M-splines with 5
#' equidistant nodes, "6-quant-piecewise" for piecewise constant risk over 5
#' intervals and nodes defined at the quantiles of the times of events
#' distribution and "9-manual-splines" for M-splines risk function with 9
#' nodes, the vector of 7 interior nodes being entered in the argument
#' \code{hazardnodes}. In the presence of competing events, a vector of hazards
#' should be provided such as \code{hazard=c("Weibull","splines"} with 2 causes
#' of event, the first one modelled by a Weibull baseline cause-specific risk
#' function and the second one by splines.
#' @param hazardtype optional indicator for the type of baseline risk function
#' when ng>1. By default "Specific" indicates a class-specific baseline risk
#' function. Other possibilities are "PH" for a baseline risk function
#' proportional in each latent class, and "Common" for a baseline risk function
#' that is common over classes. In the presence of competing events, a vector
#' of hazardtypes should be given.
#' @param hazardnodes optional vector containing interior nodes if
#' \code{splines} or \code{piecewise} is specified for the baseline hazard
#' function in \code{hazard}.
#' @param TimeDepVar optional vector containing an intermediate time
#' corresponding to a change in the risk of event. This time-dependent
#' covariate can only take the form of a time variable with the assumption that
#' there is no effect on the risk before this time and a constant effect on the
#' risk of event after this time (example: initiation of a treatment to account
#' for).
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
#' @param range optional vector indicating the range of the outcome (that is
#' the minimum and maximum). By default, the range is defined according to the
#' minimum and maximum observed values of the outcome. The option should be
#' used only for Beta and Splines transformations.
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
#' values.  (3) when ng>1, a Jointlcmm object is entered. It should correspond
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
#' iterative algorithm. By default, maxiter=150.
#' @param nsim optional number of points for the predicted survival curves and
#' predicted baseline risk curves. By default, nsim=100.
#' @param prior optional name of a covariate containing a prior information
#' about the latent class membership. The covariate should be an integer with
#' values in 0,1,...,ng. Value O indicates no prior for the subject while a
#' value in 1,...,ng indicates that the subject belongs to the corresponding
#' latent class.
#' @param logscale optional boolean indicating whether an exponential
#' (logscale=TRUE) or a square (logscale=FALSE -by default) transformation is
#' used to ensure positivity of parameters in the baseline risk functions. See
#' details section
#' @param subset a specification of the rows to be used: defaults to all rows.
#' This can be any valid indexing vector for the rows of data or if that is not
#' supplied, a data frame made up of the variable used in formula.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @param posfix Optional vector specifying the indices in vector B of the
#' parameters that should not be estimated. Default to NULL, all parameters are
#' estimated.
#' @param partialH optional logical for Piecewise and Splines baseline risk
#' functions only. Indicates whether the parameters of the baseline risk
#' functions can be dropped from the Hessian matrix to define convergence
#' criteria.
#' @param verbose logical indicating if information about computation should be
#' reported. Default to TRUE.
#' @param returndata logical indicating if data used for computation should be
#' returned. Default to FALSE, data are not returned.
#' @return The list returned is: \item{loglik}{log-likelihood of the model}
#' \item{best}{vector of parameter estimates in the same order as specified in
#' \code{B} and detailed in section \code{details}} \item{V}{vector containing
#' the upper triangle matrix of variance-covariance estimates of \code{Best}
#' with exception for variance-covariance parameters of the random-effects for
#' which \code{V} contains the variance-covariance estimates of the Cholesky
#' transformed parameters displayed in \code{cholesky}} \item{gconv}{vector of
#' convergence criteria: 1. on the parameters, 2. on the likelihood, 3. on the
#' derivatives} \item{conv}{status of convergence: =1 if the convergence
#' criteria were satisfied, =2 if the maximum number of iterations was reached,
#' =4 or 5 if a problem occured during optimisation} \item{call}{the matched
#' call} \item{niter}{number of Marquardt iterations} \item{pred}{table of
#' individual predictions and residuals; it includes marginal predictions
#' (pred_m), marginal residuals (resid_m), subject-specific predictions
#' (pred_ss) and subject-specific residuals (resid_ss) averaged over classes,
#' the observation (obs) and finally the class-specific marginal and
#' subject-specific predictions (with the number of the latent class:
#' pred_m_1,pred_m_2,...,pred_ss_1,pred_ss_2,...)} \item{pprob}{table of
#' posterior classification and posterior individual class-membership
#' probabilities based on the longitudinal data and the time-to-event data}
#' \item{pprobY}{table of posterior classification and posterior individual
#' class-membership probabilities based only on the longitudinal data}
#' \item{predRE}{table containing individual predictions of the random-effects:
#' a column per random-effect, a line per subject} \item{cholesky}{vector
#' containing the estimates of the Cholesky transformed parameters of the
#' variance-covariance matrix of the random-effects} \item{scoretest}{Statistic
#' of the Score Test for the conditional independence assumption of the
#' longitudinal and survival data given the latent class structure. Under the
#' null hypothesis, the statistics is a Chi-square with p degrees of freedom
#' where p indicates the number of random-effects in the longitudinal mixed
#' model. See Jacqmin-Gadda and Proust-Lima (2009) for more details.}
#' \item{predSurv}{table of predictions giving for the window of times to event
#' (called "time"), the predicted baseline risk function in each latent class
#' (called "RiskFct") and the predicted cumulative baseline risk function in
#' each latent class (called "CumRiskFct").} \item{hazard}{internal information
#' about the hazard specification used in related functions} \item{data}{the
#' original data set (if returndata is TRUE)}
#' %\item{specif}{internal information used in related functions}
#' %\item{Names}{internal information used in related fnctions}
#' %\item{Names2}{internal information used in related functions}
#' @author Cecile Proust Lima, Amadou Diakite and Viviane Philipps
#' 
#' \email{cecile.proust-lima@@inserm.fr}
#' @seealso \code{\link{postprob}}, \code{\link{plot.Jointlcmm}},
#' \code{\link{plot.predict}}, \code{\link{epoce}}
#' @references
#' 
#' Proust-Lima C, Philipps V, Liquet B (2017). Estimation of Extended Mixed 
#' Models Using Latent Classes and Latent Processes: The R Package lcmm. 
#' Journal of Statistical Software, 78(2), 1-56. doi:10.18637/jss.v078.i02
#'  
#' 
#' Lin, H., Turnbull, B. W., McCulloch, C. E. and Slate, E. H. (2002). Latent
#' class models for joint analysis of longitudinal biomarker and event process
#' data: application to longitudinal prostate-specific antigen readings and
#' prostate cancer. Journal of the American Statistical Association 97, 53-65.
#' 
#' Proust-Lima, C. and Taylor, J. (2009). Development and validation of a
#' dynamic prognostic tool for prostate cancer recurrence using repeated
#' measures of post-treatment PSA: a joint modelling approach. Biostatistics
#' 10, 535-49.
#' 
#' Jacqmin-Gadda, H. and Proust-Lima, C. (2010). Score test for conditional
#' independence between longitudinal outcome and time-to-event given the
#' classes in the joint latent class model. Biometrics 66(1), 11-9
#' 
#' Proust-Lima, Sene, Taylor and Jacqmin-Gadda (2014). Joint latent class
#' models of longitudinal and time-to-event data: a review. Statistical Methods
#' in Medical Research 23, 74-90.
#' @examples
#' 
#' 
#' 
#' #### Example of a joint latent class model estimated for a varying number
#' # of latent classes: 
#' # The linear mixed model includes a subject- (ID) and class-specific 
#' # linear trend (intercept and Time in fixed, random and mixture components)
#' # and a common effect of X1 and its interaction with time over classes 
#' # (in fixed).
#' # The variance of the random intercept and slopes are assumed to be equal 
#' # over classes (nwg=F).
#' # The covariate X3 predicts the class membership (in classmb). 
#' # The baseline hazard function is modelled with cubic M-splines -3 
#' # nodes at the quantiles- (in hazard) and a proportional hazard over 
#' # classes is assumed (in hazardtype). Covariates X1 and X2 predict the 
#' # risk of event (in survival) with a common effect over classes for X1
#' # and a class-specific effect of X2.
#' # !CAUTION: for illustration, only default initial values where used but 
#' # other sets of initial values should be tried to ensure convergence
#' # towards the global maximum.
#' 
#' 
#' \dontrun{
#' #### estimation with 1 latent class (ng=1): independent models for the 
#' # longitudinal outcome and the time of event
#' m1 <- Jointlcmm(fixed= Ydep1~X1*Time,random=~Time,subject='ID'
#' ,survival = Surv(Tevent,Event)~ X1+X2 ,hazard="3-quant-splines"
#' ,hazardtype="PH",ng=1,data=data_lcmm)
#' summary(m1)
#' #Goodness-of-fit statistics for m1:
#' #    maximum log-likelihood: -3944.77 ; AIC: 7919.54  ;  BIC: 7975.09  
#' }
#' 
#' #### estimation with 2 latent classes (ng=2)
#' m2 <- Jointlcmm(fixed= Ydep1~Time*X1,mixture=~Time,random=~Time,
#' classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~X1+mixture(X2),
#' hazard="3-quant-splines",hazardtype="PH",ng=2,data=data_lcmm,
#' B=c(0.64,-0.62,0,0,0.52,0.81,0.41,0.78,0.1,0.77,-0.05,10.43,11.3,-2.6,
#' -0.52,1.41,-0.05,0.91,0.05,0.21,1.5))
#' summary(m2)
#' #Goodness-of-fit statistics for m2:
#' #       maximum log-likelihood: -3921.27; AIC: 7884.54; BIC: 7962.32  
#' 
#' \dontrun{
#' #### estimation with 3 latent classes (ng=3)
#' m3 <- Jointlcmm(fixed= Ydep1~Time*X1,mixture=~Time,random=~Time,
#' classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~ X1+mixture(X2),
#' hazard="3-quant-splines",hazardtype="PH",ng=3,data=data_lcmm,
#' B=c(0.77,0.4,-0.82,-0.27,0,0,0,0.3,0.62,2.62,5.31,-0.03,1.36,0.82,
#' -13.5,10.17,10.24,11.51,-2.62,-0.43,-0.61,1.47,-0.04,0.85,0.04,0.26,1.5))
#' summary(m3)
#' #Goodness-of-fit statistics for m3:
#' #       maximum log-likelihood: -3890.26 ; AIC: 7834.53;  BIC: 7934.53  
#' 
#' #### estimation with 4 latent classes (ng=4)
#' m4 <- Jointlcmm(fixed= Ydep1~Time*X1,mixture=~Time,random=~Time,
#' classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~ X1+mixture(X2),
#' hazard="3-quant-splines",hazardtype="PH",ng=4,data=data_lcmm,
#' B=c(0.54,-0.42,0.36,-0.94,-0.64,-0.28,0,0,0,0.34,0.59,2.6,2.56,5.26,
#' -0.1,1.27,1.34,0.7,-5.72,10.54,9.02,10.2,11.58,-2.47,-2.78,-0.28,-0.57,
#' 1.48,-0.06,0.61,-0.07,0.31,1.5))
#' summary(m4)
#' #Goodness-of-fit statistics for m4:
#' #   maximum log-likelihood: -3886.93 ; AIC: 7839.86;  BIC: 7962.09  
#' 
#' 
#' ##### The model with 3 latent classes is retained according to the BIC  
#' ##### and the conditional independence assumption is not rejected at
#' ##### the 5% level. 
#' # posterior classification
#' plot(m3,which="postprob")
#' # Class-specific predicted baseline risk & survival functions in the 
#' # 3-class model retained (for the reference value of the covariates) 
#' plot(m3,which="baselinerisk",bty="l")
#' plot(m3,which="baselinerisk",ylim=c(0,5),bty="l")
#' plot(m3,which="survival",bty="l")
#' # class-specific predicted trajectories in the 3-class model retained 
#' # (with characteristics of subject ID=193)
#' data <- data_lcmm[data_lcmm$ID==193,]
#' plot(predictY(m3,var.time="Time",newdata=data,bty="l")
#' # predictive accuracy of the model evaluated with EPOCE
#' vect <- 1:15
#' cvpl <- epoce(m3,var.time="Time",pred.times=vect)
#' summary(cvpl)
#' plot(cvpl,bty="l",ylim=c(0,2))
#' ############## end of example ##############
#' }
#' 
#' @export
#' 
#' 
Jointlcmm <- function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,
                      survival,hazard="Weibull",hazardtype="Specific",
                      hazardnodes=NULL,TimeDepVar=NULL,link=NULL,intnodes=NULL,
                      epsY=0.5,range=NULL,cor=NULL,data,B,convB=0.0001,
                      convL=0.0001,convG=0.0001,maxiter=100,nsim=100,prior,
                      logscale=FALSE,subset=NULL,na.action=1,posfix=NULL,
                      partialH=FALSE,verbose=TRUE,returndata=FALSE)
    {
        
        ptm<-proc.time()
        if(verbose==TRUE) cat("Be patient, Jointlcmm is running ... \n")

        cl <- match.call()

        if(!missing(mixture) & ng==1) stop("No mixture can be specified with ng=1")
        if(missing(mixture) & ng>1) stop("The argument mixture has to be specified for ng > 1")
        if(!missing(classmb) & ng==1) stop("No classmb can be specified with ng=1")
        if(missing(fixed)) stop("The argument Fixed must be specified in any model")
        if(missing(random)) random <- ~-1  
        if(missing(classmb) & ng==1) classmb <- ~-1
        if(missing(classmb) & ng>1) classmb <- ~1
        if(missing(mixture)) mixture <- ~-1
        if(ng==1 & nwg==TRUE) stop ("The argument nwg should be FALSE for ng=1")
        if(ng==1) hazardtype <- "Specific"  #??
        if(any(!(hazardtype %in% c("Common","Specific","PH")))) stop("'hazardtype' should be either 'Common' or 'Specific' or 'PH'")
        
        if(class(fixed) != "formula") stop("The argument fixed must be a formula")
        if(class(mixture) != "formula") stop("The argument mixture must be a formula")
        if(class(random) != "formula") stop("The argument random must be a formula")
        if(class(classmb) != "formula") stop("The argument classmb must be a formula")
        if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
        if(nrow(data)==0) stop("Data should not be empty")
        if(missing(subject)){ stop("The argument subject must be specified in any model even without random-effects")}
        if(!is.numeric(data[,subject])) stop("The argument subject must be numeric")
        if(is.null(link)) link <- "NULL"
        if(any(link=="thresholds"))  stop("The link function thresholds is not available yet")
        if(link %in% c("linear","beta") & !is.null(intnodes)) stop("Intnodes should only be specified with splines links")
        
        nom.subject <- as.character(subject)
        if(!isTRUE(nom.subject %in% colnames(data))) stop(paste("Data should contain variable",nom.subject))
        
        nom.prior <- NULL
        if(!missing(prior))
            {
                nom.prior <- as.character(prior)
                if(!isTRUE(nom.prior %in% colnames(data))) stop(paste("Data should contain variable",nom.prior))
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
        
        if(!(na.action %in% c(1,2))) stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")

        if(length(posfix) & missing(B)) stop("A set of initial parameters must be specified if some parameters are not estimated")


        ## garder data tel quel pour le renvoyer
        if(returndata==TRUE)
        {
            datareturn <- data
        }
        else
        {
            datareturn <- NULL
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

        attributes(data)$terms <- NULL
        
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
        nbevt <- length(which(names(table(Event))>0))    #length(unique(Event))-1   
        if(nbevt<1) stop("No observed event in the data")
       
        
        ##pour acces aux attributs des formules
        afixed <- terms(fixed)
        if(link!="NULL" & attr(afixed,"intercept")==0) stop("An intercept should appear in fixed for identifiability purposes")
        amixture <- terms(mixture)
        arandom <- terms(random)
        aclassmb <- terms(classmb)
        ## pour la formule pour survivial, creer 3 formules : 
        ## une pour les covariables en mixture, une pour les covariables avec effet specifique a la cause, et une pour les effets communs.  
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
        
        ##verifier si toutes les variables sont dans data
        variables <- c(attr(afixed,"variables"),attr(arandom,"variables"),attr(amixture,"variables"),attr(aclassmb,"variables"),attr(asurv,"variables"))
        variables <- unlist(lapply(variables,all.vars))  
        if(!all(variables %in% colnames(data))) stop(paste("Data should contain the variables",paste(unique(variables),collapse=" ")))

 

###liste des variables utilisees  (sans les interactions et sans Y)
        if (ncor0>0) ttesLesVar <- unique(c(variables,cor.var.time))
        else ttesLesVar <- unique(c(variables))
        nomY <-  as.character(attr(afixed,"variables")[2])  
        ttesLesVar <- setdiff(ttesLesVar,nomY)
        
###subset de data avec les variables utilisees
        newdata <- data[,unique(c(nom.subject,nomY,noms.surv,ttesLesVar,nom.prior,nom.timedepvar)),drop=FALSE]

        ## remplacer les NA de prior par 0  
        if(!is.null(nom.prior))
            {
                prior <- newdata[,nom.prior]
                newdata[which(is.na(prior)),nom.prior] <- 0
                prior[which(is.na(prior))] <- 0
            }

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


        ##enlever les NA
        linesNA <- apply(newdata,2,function(v) which(is.na(v)))
        linesNA <- unique(unlist(linesNA))  
        
        if(length(linesNA))
            {
                if(na.action==1) newdata <- newdata[-linesNA,,drop=FALSE] 
                if(na.action==2) stop("Data contain missing values.")
                Tentry <- Tentry[-linesNA]  
                Tevent <- Tevent[-linesNA] 
                Event <- Event[-linesNA]
                Tint <- Tint[-linesNA]
                prior <- prior[-linesNA]
            }


### Y0
        Y0 <- newdata[,nomY]  

###creation de X0 (ttes les var + interactions)
        Xintercept <- model.matrix(~1,data=newdata) #pr etre sur d'avoir I en premier
        Xfixed <- model.matrix(fixed[-2], data=newdata)
        Xmixture <- model.matrix(mixture, data=newdata)
        Xrandom <- model.matrix(random, data=newdata)
        Xclassmb <- model.matrix(classmb, data=newdata)
        Xsurv <- model.matrix(form.commun,data=newdata)
        Xsurvmix <- model.matrix(form.mixture,data=newdata)
        Xsurvcause <- model.matrix(form.cause,data=newdata)
        for (k in 1:nbevt)
            {
                assign(paste("Xsurvcause",k,sep=""),model.matrix(form.causek[[k]],data=newdata))
            }        


        z.fixed <- strsplit(colnames(Xfixed),split=":",fixed=TRUE)
        z.fixed <- lapply(z.fixed,sort)
        
        if(random != ~-1)
            {   
                z.random <- strsplit(colnames(Xrandom),split=":",fixed=TRUE)
                z.random <- lapply(z.random,sort)
            }
        else
            {
                z.random <- list() 
            }
        
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
        
        if(form.commun != ~-1)
            {
                z.surv <- strsplit(colnames(Xsurv),split=":",fixed=TRUE)
                z.surv <- lapply(z.surv,sort)
            }
        else
            {
                z.surv <- list() 
            }
        
        if(form.mixture != ~-1)
            {
                z.survmix <- strsplit(colnames(Xsurvmix),split=":",fixed=TRUE)
                z.survmix <- lapply(z.survmix,sort)
            }
        else
            {
                z.survmix <- list() 
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
        
        for(k in 1:nbevt)
            {
                if(form.causek[[k]] != ~-1)
                    {
                        assign(paste("z.survcause",k,sep=""),strsplit(colnames(get(paste("Xsurvcause",k,sep=""))),split=":",fixed=TRUE))
                        assign(paste("z.survcause",k,sep=""),lapply(get(paste("z.survcause",k,sep="")),sort))
                    }
                else
                    {
                        assign(paste("z.survcause",k,sep=""),list())
                    }
            }

        

        if(!all(z.mixture %in% z.fixed))  stop("The covariates in mixture should also be included in the argument fixed")

        X0 <- cbind(Xintercept,Xfixed, Xrandom, Xclassmb,Xsurv,Xsurvmix,Xsurvcause)
        for(k in 1:nbevt)
            {
                X0 <- cbind(X0,get(paste("Xsurvcause",k,sep="")))
            }
        

        nom.unique <- unique(colnames(X0))
        X0 <- X0[,nom.unique,drop=FALSE]  

        if(ncor0>0)
            {
                if(!(cor.var.time %in% colnames(X0)))
                    {
                        X0 <- cbind(X0, newdata[,cor.var.time])
                        colnames(X0) <- c(nom.unique, cor.var.time)
                        nom.unique <- colnames(X0)  
                    }
            }



        X0 <- as.matrix(X0)
###X0 fini  


###test de link
        if(link %in% c("splines","Splines"))
            {
                link <- "5-quant-splines"
            }


        idlink <- 2
        if(link=="linear") idlink <- 0
        if(link=="beta") idlink <- 1
        if(link=="NULL") idlink <- -1
        
        if (idlink==1)
            {
                if (epsY<=0)
                    {
                        stop("Argument 'epsY' should be positive.")
                    }
            } 


        ##remplir range si pas specifie
        if(!is.null(range) & length(range)!=2) stop("Vector range should be of length 2.")
        if(length(range)==2)
            {
                if(max(Y0)>range[2] | min(Y0)<range[1]) stop("The range specified do not cover the entire range of the data")
            }
        
        if(is.null(range))
            {
                min1 <- min(Y0)
                min2 <- round(min1,3)
                if(min1<min2) min2 <- min2-0.001

                max1 <- max(Y0)
                max2 <- round(max1,3)
                if(max1>max2) max2 <- max2+0.001
                
                range <- c(min2,max2)
            }

        nbzitr <- 2  #nbzitr = nb de noeuds si splines, 2 sinon  
        spltype <- NULL

        if(idlink==2)
            {
                spl <- strsplit(link,"-")
                if(length(spl[[1]])!=3) stop("Invalid argument link")

                nbzitr <- as.numeric(spl[[1]][1])
                spltype <- spl[[1]][2]

                if(!(spltype %in% c("equi","quant","manual"))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")
                if(!(spl[[1]][3] %in% c("splines","Splines"))) stop("Invalid argument link")

                if(!is.null(intnodes))
                    {
                        if(spltype != "manual")  stop("intnodes should be NULL if the nodes are not chosen manually")

                        if(length(intnodes) != (nbzitr-2)) stop(paste("Vector intnodes should be of length",nbzitr-2)) 
                    }
            }

        ##remplir zitr (contient les noeuds si splines, sinon min et max)
        zitr0 <- rep(0,nbzitr)
        if(idlink %in% c(0,1)) zitr0[1:2] <- range
        
        if(idlink==2)
            {
                if(spltype=="manual")
                    {
                        zitr0[1] <- range[1]
                        zitr0[nbzitr] <- range[2]  
                        zitr0[2:(nbzitr-1)] <- intnodes  

                        ##verifier s'il y a des obs entre les noeuds
                        hcounts <- hist(Y0,breaks=zitr0,plot=FALSE,include.lowest=TRUE,right=TRUE)$counts
                        if(any(hcounts==0)) stop("Link function can not be estimated. Please try other nodes such that there are observations in each interval.")      
                    }

                if(spltype=="equi") zitr0[1:nbzitr] <- seq(range[1],range[2],length.out=nbzitr)
                if(spltype=="quant")
                    {
                        Y0bis <- Y0
                        if(range[1]<min(Y0)) Y0bis <- c(range[1],Y0)
                        if(range[2]>max(Y0)) Y0bis <- c(range[2],Y0bis)
                        zitr0[1:nbzitr] <- quantile(Y0bis,probs=seq(0,1,length.out=nbzitr))
                    }
            }

###uniqueY0 et indiceY0
        uniqueY0 <- 0
        indiceY0 <- 0
        nvalSPL0 <- 0

        if(idlink == 2)
            {
                uniqueY0 <- sort(unique(Y0))
                permut <- order(order(Y0))  # sort(y)[order(order(y))] = y
                indice <- rep(1:length(uniqueY0), as.vector(table(Y0)))
                indiceY0 <- indice[permut]
                nvalSPL0 <- length(uniqueY0)
            } 


###ordonner les mesures par individu
        IND <- newdata[,nom.subject]
        #IDnum <- as.numeric(IND)
        if(is.null(nom.prior)) prior <- rep(0,length(Y0))  
        if(!length(indiceY0)) indiceY0 <- rep(0,length(Y0))  
        matYX <- cbind(IND,prior,Y0,indiceY0,Tentry,Tevent,Event,Tint,X0)
        matYXord <- matYX[order(IND),]
        Y0 <- as.numeric(matYXord[,3])
        X0 <- apply(matYXord[,-c(1:8),drop=FALSE],2,as.numeric)
        #IDnum <- matYXord[,1]
        IND <- matYXord[,1]
        indiceY0 <- as.numeric(matYXord[,4])
        prior0 <- as.numeric(unique(matYXord[,c(1,2)])[,2])
        if(length(prior0)!=length(unique(IND))) stop("Please check 'prior' argument. Subjects can not have multiple assigned classes.")

        ## nombre de sujets
        ns0 <- length(unique(IND))  
        
### Tevent, Tentry et Event de dim ns  
        #data.surv <- unique(cbind(IND,matYX[,c(6,7,8,9)]))
                                        #if(nrow(data.surv) != ns0) stop("Subjects cannot have several times to event.")
        
        nmes <- as.vector(table(IND))
        data.surv <- apply(matYXord[cumsum(nmes),c(5,6,7,8)],2,as.numeric)
        
        tsurv0 <- data.surv[,1] 
        tsurv <- data.surv[,2]
        devt <- data.surv[,3]
        tsurvint <- data.surv[,4]
        ind_survint <- (tsurvint<tsurv) + 0 



### test de hazard
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

        ## si plusieurs splines, noeuds doivent etre identiques
        nbspl <- length(which(sapply(haz13,getElement,3) %in% c("splines","Splines")))
        if(nbspl>1)
            {
                haz3 <- haz13[which(sapply(haz13,getElement,3) %in% c("splines","Splines"))]
                if(length(unique(sapply(haz3,getElement,1)))>1) stop("The nodes location of all splines hazard functions must be identical")

                if(length(unique(sapply(haz3,getElement,2)))>1) stop("The nodes location of all splines hazard functions must be identical")
            }

        nz <- rep(2,nbevt) 
        locnodes <- NULL  
        typrisq <- rep(2,nbevt)   
        nprisq <- rep(2,nbevt) 

        nbnodes <- 0 #longueur de hazardnodes
        ii <- 0
        dejaspl <- 0
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
                                if(typrisq[i]==1 | dejaspl==0)
                                    {
                                        if(length(arghaz)>1 | i==1 )
                                            {
                                                nbnodes <- nbnodes + nz[i]-2

                                            }
                                    }
                                if(typrisq[i]==3) dejaspl <- 1
                            }
                
                if(!all(locnodes %in% c("equi","quant","manual"))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")      
                    }
                
                if(!is.null(hazardnodes))
                    {
                        if(!any(locnodes == "manual"))  stop("hazardnodes should be NULL if the nodes are not chosen manually")
                        
                        if(length(hazardnodes) != nbnodes) stop(paste("Vector hazardnodes should be of length",nbnodes)) 
                    }  
            }
        else
            {
                if(!is.null(hazardnodes)) stop("hazardnodes should be NULL if Weibull baseline risk functions are chosen")
            }


        if(nbevt>1 & length(arghaz)==1 & nbnodes>0)
            {
                hazardnodes <- rep(hazardnodes,length.out=nbnodes*nbevt)
            }

        zi <- matrix(0,nrow=max(nz),ncol=nbevt)
        nb <- 0   
 
        minT1 <- 0
        maxT1 <- max(tsurv)
        tsurvevt <- tsurv #[which(devt!=0)] 

        if(idtrunc==1)
            {
                minT1 <- min(tsurv,tsurv0)
                maxT1 <- max(tsurv,tsurv0)
            }
        ##??
        #if(!(minT %in% tsurvevt)) tsurvevt <- c(minT,tsurvevt)
        #if(!(maxT %in% tsurvevt)) tsurvevt <- c(tsurvevt,maxT)

        ## arrondir
        minT2 <- round(minT1,3)
        if(minT1<minT2) minT2 <- minT2-0.001
        minT <- minT2

        maxT2 <- round(maxT1,3)
        if(maxT1>maxT2) maxT2 <- maxT2+0.001
        maxT <- maxT2
        
        
        ii <- 0  
        for(i in 1:nbevt)
            {
                if(typrisq[i]==2)
                    {
                        zi[1:2,i] <- c(minT,maxT)
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
                                #pi <- seq(0,1,length.out=nz[i])
                                #pi <- pi[-length(pi)]
                                #pi <- pi[-1]
                                pi <- c(1:(nz[i]-2))/(nz[i]-1)
                                qi <- quantile(tsurvevt,prob=pi)
                                zi[1,i] <- minT
                                zi[2:(nz[i]-1),i] <- qi
                                zi[nz[i],i] <- maxT
                            }
                    }   
            }
        
        hazardtype <- rep(hazardtype,length.out=nbevt)  
        risqcom <- (hazardtype=="Common") + (hazardtype=="PH")*2   
        nrisq <- (risqcom==1)*nprisq + (risqcom==0)*nprisq*ng + (risqcom==2)*(nprisq+ng-1)  


        
###parametres pour fortran
        ng0 <- ng
        nv0 <- dim(X0)[2]
        nobs0 <- length(Y0)
        idiag0 <- ifelse(idiag==TRUE,1,0)
        nwg0 <- ifelse(nwg==TRUE,1,0)
        nrisqtot <- sum(nrisq)   
        
        loglik <- 0
        ni <- 0
        istop <- 0
        gconv <- rep(0,3)
        ppi0 <- rep(0,ns0*ng0)
        ppitest0 <- rep(0,ns0*ng0)
        resid_m <- rep(0,nobs0)
        resid_ss <- rep(0,nobs0)
        pred_m_g <- rep(0,nobs0*ng0)
        pred_ss_g <- rep(0,nobs0*ng0)
        Yobs <- rep(0,nobs0)
        #rlindiv <- rep(0,ns0)
        marker <- rep(0,nsim)
        transfY <- rep(0,nsim)
        #Ydiscret <- 0
        #UACV <- 0
        #vraisdiscret <- 0
        nmes0 <- as.vector(table(IND))
        logspecif <- as.numeric(logscale)
        time <- seq(minT,maxT,length.out=nsim)
        risq_est <- matrix(0,nrow=nsim*ng0,ncol=nbevt)
        risqcum_est <- matrix(0,nrow=nsim*ng0,ncol=nbevt)
        statglob <- 0
        statevt <- rep(0,nbevt)

###remplir idprob, etc
        z.X0 <- strsplit(nom.unique,split=":",fixed=TRUE)
        z.X0 <- lapply(z.X0,sort)

        idprob0 <- z.X0 %in% z.classmb + 0
        idea0 <- z.X0 %in% z.random + 0
        idg0 <- (z.X0 %in% z.fixed) + (z.X0 %in% z.mixture)


        ## indicateurs pr variables survie
        idcause <- rep(0,nv0)
        for(k in 1:nbevt)
            {
                idcause <- idcause + (z.X0 %in% get(paste("z.survcause",k,sep="")) )
            }

        idcom <- rep(0,nv0)
        idspecif <- matrix(0,ncol=nv0,nrow=nbevt)
        for(j in 1:nv0)
            {
                if((z.X0[j] %in% z.surv) &  !(z.X0[j] %in% z.survcause) & (idcause[j]==0))
                    {
                        idcom[j] <- 1
                        idspecif[,j] <- 1
                        #print("X")
                    }

                 if((z.X0[j] %in% z.survmix) & ( !(z.X0[j] %in% z.survcause) & (idcause[j]==0)))
                    {
                        idcom[j] <- 1
                        idspecif[,j] <- 2
                        #print("mixture(X)")
                    }    
                if((z.X0[j] %in% z.survmix) & (z.X0[j] %in% z.survcause))
                    {
                        idcom[j] <- 0
                        idspecif[,j] <- 2
                        #print("cause(mixture(X))")
                    }

                  if((z.X0[j] %in% z.survcause) & (!(z.X0[j] %in% z.survmix)))
                    {
                        idcom[j] <- 0
                        idspecif[,j] <- 1
                        #print("cause(X)")
                    }              

                if(idcause[j]!=0)
                    {
                        if(z.X0[j] %in% z.survmix)
                            {
                                for(k in 1:nbevt)
                                    {
                                        if(z.X0[j] %in% get(paste("z.survcause",k,sep="")))
                                            {
                                                idcom[j] <- 0
                                                idspecif[k,j] <- 2
                                                #cat("causek(mixture(X)) ,k=",k,"\n")
                                            }
                                    }
                            }
                        else
                            {
                                for(k in 1:nbevt)
                                    {
                                        if(z.X0[j] %in% get(paste("z.survcause",k,sep="")))
                                            {
                                                idcom[j] <- 0
                                                idspecif[k,j] <- 1
                                                #cat("causek(X) ,k=",k,"\n")
                                            }                                        
                                    }
                            }
                    }
                
            }
        
        
        
        ## mettre des 0 pour l'intercept
        idcom[1] <- 0
        idspecif[,1] <- 0
       
       
        ## quelle variable est TimeDepVar
        idtdv <- z.X0 %in% nom.timedepvar + 0

        
        #if(any(idcom>1) & nbevt<2) stop("No event specific effect can be estimated with less than two events")
        
        if (ncor0>0) idcor0 <- z.X0 %in% cor.var.time +0
        else idcor0 <- rep(0,nv0)


        ## Si pas TimeDepVar dans formule survival
        if(length(nom.timedepvar) & all(idtdv==0))
            {
                stop("Variable in 'TimeDepVar' should also appear as a covariate in the 'survival' argument")
                ## ou l'ajouter?
                
            }

        
        ## si on a ignore TimeDepVar #pas utile si stop
        if(length(nom.timedepvar) & nvdepsurv==0)
            {
                jj <- which(idtdv==1)
                idtdv[jj] <- 0
                idcom[jj] <- 0
                idspecif[,jj] <- 0
                idcause[jj] <- 0
                nom.timedepvar <- NULL
            }
    
        
        
        nea0 <- sum(idea0)
        predRE <- rep(0,ns0*nea0)

        ## nb coef pr survie
        nvarxevt <- 0
        nvarxevt2 <- 0 # pr valeurs initiales calculees avec Fortran
        for(j in 1:nv0)
            {
                if(idcom[j]==1)
                    {
                        if(all(idspecif[,j]==1))
                            {
                                nvarxevt <- nvarxevt + 1
                                nvarxevt2 <- nvarxevt2 + 1
                            }
                        if(all(idspecif[,j]==2))
                            {
                                nvarxevt <- nvarxevt + ng
                                nvarxevt2 <- nvarxevt2 + 1
                            }
                    }

                if(idcom[j]==0)
                    {
                        if(all(idspecif[,j]==0)) next
                        for(k in 1:nbevt)
                            {
                                if(idspecif[k,j]==1)
                                    {
                                        nvarxevt <- nvarxevt + 1
                                        nvarxevt2 <- nvarxevt2 + 1
                                    }
                                if(idspecif[k,j]==2)
                                    {
                                        nvarxevt <- nvarxevt + ng
                                        nvarxevt2 <- nvarxevt2 + 1
                                    }
                            }
                    }
            }


        ntrtot0 <- (idlink==-1) + 2*(idlink==0) + 4*(idlink==1) + (nbzitr+2)*(idlink==2)
        nef <- sum(idg0==1) + ng0*sum(idg0==2)
        if(idlink!=-1) nef <- nef-1

        nvc <- ifelse(idiag0==1,nea0,nea0*(nea0+1)/2)
        nw <- (ng0-1)*nwg0
        nprob <- sum(idprob0)*(ng0-1)

        
        ##nombre total de parametres
        NPM <- (ng0-1)*sum(idprob0) +
            nrisqtot + nvarxevt +
                nef + nvc + nw + ncor0 + ntrtot0
                     
        V <- rep(0, NPM*(NPM+1)/2)  #pr variance des parametres

        ## prm fixes
        fix0 <- rep(0,NPM)
        if(length(posfix))
            {
                if(any(!(posfix %in% 1:NPM))) stop("Indexes in posfix are not correct")
                
                fix0[posfix] <- 1
            }        
        if(length(posfix)==NPM) stop("No parameters to estimate")

        
        ## pour H restreint
        pbH0 <- rep(0,NPM)
        if(any(typrisq %in% c(1,3)))
            {
                for(k in 1:nbevt)
                    {
                        if(typrisq[k] %in% c(1,3))
                            {
                                pbH0[nprob+sum(nrisq[1:k])-nrisq[k]+1:nrisq[k]] <- 1
                                if(risqcom[k]==2)
                                    {
                                        pbH0[nprob+sum(nrisq[1:k])] <- 0
                                    }
                            }
                    }
            }
        pbH0[posfix] <- 0
        Hr0 <- as.numeric(partialH)
        if(sum(pbH0)==0 & Hr0==1) stop("No partial Hessian matrix can be defined in the absence of baseline risk function approximated by splines or piecewise linear function")
        
    
        ## gestion de B=random(mod)

        Brandom <- FALSE
        if(length(cl$B)==2)
            {
                if(class(eval(cl$B[[2]]))!="Jointlcmm") stop("The model specified in B should be of class Jointlcmm")
                if(as.character(cl$B[1])!="random") stop("Please use random() to specify random initial values")
                
                Brandom <- TRUE
                B <- eval(cl$B[[2]])

                if(length(posfix)) stop("Argument posfix is not compatible with random intial values")
            }
    
        ##valeurs initiales
        if(!(missing(B)))
            {
                if(is.vector(B))
                    {
                        if (length(B)==NPM) b <- B
                        else stop(paste("Vector B should be of length",NPM))
                    }
                else
                    {
                        if(class(B)!="Jointlcmm") stop("B should be either a vector or an object of class Jointlcmm")

                        b <- rep(0,NPM)

                        if(ng>1 & B$ng==1)
                            {
                                nef2 <- sum(idg0!=0)
                                if(idlink != -1)
                                    {
                                        nef2 <- nef2-1 #car intercept fixe a 0
                                    }
                                NPM2 <- sum(nprisq)+nvarxevt2+nef2+nvc+ncor0+ntrtot0
                                if(length(B$best)!=NPM2) stop("B is not correct")

                                
                                if(Brandom==FALSE)
                                    {
                                        ## B deterministe
                                        for(ke in 1:nbevt)
                                            {
                                                if(risqcom[ke]==0)
                                                    {
                                                        indb <- nprob+sum(nrisq[1:ke])-nrisq[ke]+1:nrisq[ke]
                                                        indinit <- sum(nprisq[1:ke])-nprisq[ke]+1:nprisq[ke]
                                                        if(B$conv==1)
                                                            {
                                                                b[indb] <- abs(rep(B$best[indinit],ng0)+rep((1:ng0)-(ng0+1)/2,each=nprisq[ke])*rep(sqrt(B$V[indinit*(indinit+1)/2]),ng0))
                                                            }
                                                        else
                                                            {
                                                                b[indb] <- abs(rep(B$best[indinit],ng0)+rep((1:ng0)-(ng0+1)/2,each=nprisq[ke])*rep(B$best[indinit],ng0))
                                                            }
                                                    }

                                                if(risqcom[ke]==1)
                                                    {
                                                        b[nprob+sum(nrisq[1:ke])-nrisq[ke]+1:nrisq[ke]] <- B$best[sum(nprisq[1:ke])-nprisq[ke]+1:nprisq[ke]]
                                                    }

                                                if(risqcom[ke]==2)
                                                    {
                                                        b[nprob+sum(nrisq[1:ke])-nrisq[ke]+1:nrisq[ke]] <- c(B$best[sum(nprisq[1:ke])-nprisq[ke]+1:nprisq[ke]],0.5+(0:(ng0-2))*0.5)
                                                    }
                                            }


                                        ## pr bevt
                                        avtj <- nprob+nrisqtot
                                        avtj2 <- sum(nprisq)
                                        for(j in 1:nv0)
                                            {
                                                if(idcom[j]==0 & all(idspecif[,j]==0)) next
                                                
                                                if(idcom[j]==1 & all(idspecif[,j]==1))
                                                    {
                                                        b[avtj+1] <- B$best[avtj2+1]
                                                        avtj <- avtj+1
                                                        avtj2 <- avtj2+1
                                                    }

                                                if(idcom[j]==1 & all(idspecif[,j]==2))
                                                    {
                                                        if(B$conv==1) b[avtj+1:ng0] <- abs(rep(B$best[avtj2+1],ng0)+c(1:ng0-(ng0+1)/2)*rep(sqrt(B$V[(avtj2+1)*(avtj2+2)/2]),ng0))
                                                        else b[avtj+1:ng0] <- abs(rep(B$best[avtj2+1],ng0)+c(1:ng0-(ng0+1)/2)*rep(B$best[avtj2+1],ng0))

                                                        avtj <- avtj+ng0
                                                        avtj2 <- avtj2+1
                                                    }

                                                if(idcom[j]==0 & any(idspecif[,j]!=0))
                                                    {
                                                        for(k in 1:nbevt)
                                                            {
                                                                if(idspecif[k,j]==0) next

                                                                if(idspecif[k,j]==1)
                                                                    {
                                                                        b[avtj+1] <- B$best[avtj2+1]
                                                                        avtj <- avtj+1
                                                                        avtj2 <- avtj2+1
                                                                    }

                                                                if(idspecif[k,j]==2)
                                                                    {
                                                                        if(B$conv==1) b[avtj+1:ng0] <- abs(rep(B$best[avtj2+1],ng0)+c(1:ng0-(ng0+1)/2)*rep(sqrt(B$V[(avtj2+1)*(avtj2+2)/2]),ng0))
                                                                        else b[avtj+1:ng0] <- abs(rep(B$best[avtj2+1],ng0)+c(1:ng0-(ng0+1)/2)*rep(B$best[avtj2+1],ng0))

                                                                        avtj <- avtj+ng0
                                                                        avtj2 <- avtj2+1 
                                                                    }
                                                            }
                                                        
                                                    }
                                            }


                                        
                                        ## pr nef
                                        avtj <- nprob+nrisqtot+nvarxevt
                                        avtj2 <- sum(nprisq)+nvarxevt2
                                        for(j in 1:nv0)
                                            {
                                                if(idg0[j]==1)
                                                    {
                                                        if(j==1 & idlink!=-1) next
                                                        else
                                                            {
                                                                if(B$conv==1) b[avtj+1] <- B$best[avtj2+1]
                                                                avtj <- avtj+1
                                                                avtj2 <- avtj2+1
                                                            }
                                                    }

                                                if(idg0[j]==2)
                                                    {
                                                        if(j==1)
                                                            {
                                                                if(idlink!=-1)
                                                                    {
                                                                        b[avtj+1:(ng0-1)] <- -0.5*(1:(ng0-1))

                                                                        avtj <- avtj+ng0-1
                                                                    }
                                                                else
                                                                    {
                                                                        b[avtj+1:ng0] <- B$best[avtj2+1]+(1:ng0-(ng0+1)/2)*B$best[avtj2+1]
                                                                        avtj <- avtj+ng0
                                                                        avtj2 <- avtj2+1
                                                                    }
                                                            }
                                                        else
                                                            {
                                                                if(B$conv==1) b[avtj+1:ng0] <- B$best[avtj2+1]+(1:ng0-(ng0+1)/2)*sqrt(B$V[(avtj2+1)*(avtj2+1+1)/2])
                                                                else b[avtj+1:ng0] <- B$best[avtj2+1]+(1:ng0-(ng0+1)/2)*B$best[avtj2+1]

                                                                avtj <- avtj+ng0
                                                                avtj2 <- avtj2+1
                                                            }
                                                        
                                                    }
                                            }


                                        ## pr nvc
                                        if(nvc>0)
                                            {
                                                b[nprob+nrisqtot+nvarxevt+nef+1:nvc] <- B$cholesky
                                            }


                                        ## cor et transfo
                                        b[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+1):NPM] <- B$best[(sum(nprisq)+nvarxevt2+nef2+nvc+1):NPM2]

                                    }
                                else
                                    {
                                        ## B random
                                        bb <- rep(0,NPM-nprob-(ng-1)*sum(length(which(risqcom==2)))-nw) # tous les prm sauf nprob, bPH et nw
                                        vbb <- matrix(0,length(bb),length(bb))
                            
                                        VB <- matrix(0,NPM2,NPM2)
                                        VB[upper.tri(VB,diag=TRUE)] <- B$V
                                        VB <- t(VB)
                                        VB[upper.tri(VB,diag=TRUE)] <- B$V

                                        
                                        copcov <- rep(0,length(bb))
                                        copcov2 <- rep(0,length(B$best))

                                        avt <- 0
                                        for(ke in 1:nbevt)
                                            {
                                                if(risqcom[ke]==0)
                                                    {
                                                        indbb <- avt+1:nrisq[ke]
                                                        indB <- sum(nprisq[1:ke])-nprisq[ke]+1:nprisq[ke]
                                                        
                                                        bb[indbb] <- rep(B$best[indB],ng)
                                                        diag(vbb[indbb,indbb]) <- rep(diag(VB[indB,indB]),ng)

                                                        avt <- avt + nrisq[ke]
                                                    }

                                                if(risqcom[ke]==1)
                                                    {
                                                        indbb <- avt+1:nrisq[ke]
                                                        indB <- sum(nprisq[1:ke])-nprisq[ke]+1:nprisq[ke]
                                                        
                                                        bb[indbb] <- B$best[indB]
                                                        copcov[indbb] <- 1
                                                        copcov2[indB] <- 1

                                                        avt <- avt + nrisq[ke]
                                                    }

                                                if(risqcom[ke]==2)
                                                    {
                                                        indbb <- avt+1:nprisq[ke]
                                                        indB <- sum(nprisq[1:ke])-nprisq[ke]+1:nprisq[ke]
                                                        
                                                        bb[indbb] <- B$best[indB]
                                                        copcov[indbb] <- 1
                                                        copcov2[indB] <- 1


                                                        avt <- avt + nprisq[ke]
                                                    }
                                            }

                                        avt2 <- sum(nprisq)
                                        for(j in 1:nv0)
                                            {
                                                if(idcom[j]==0 & all(idspecif[,j]==0)) next
                                                
                                                if(idcom[j]==1 & all(idspecif[,j]==1))
                                                    {
                                                        bb[avt+1] <- B$best[avt2+1]
                                                        copcov[avt+1] <- 1
                                                        copcov2[avt2+1] <- 1
                                                        
                                                        avt <- avt+1
                                                        avt2 <- avt2+1
                                                    }

                                                if(idcom[j]==1 & all(idspecif[,j]==2))
                                                    {
                                                        bb[avt+1:ng] <- rep(B$best[avt2+1],ng)
                                                        diag(vbb[avt+1:ng,avt+1:ng]) <- rep(VB[avt2+1,avt2+1],ng)
                                                        
                                                        avt <- avt+ng
                                                        avt2 <- avt2+1
                                                    }

                                                if(idcom[j]==0 & any(idspecif[,j]!=0))
                                                    {
                                                        for(k in 1:nbevt)
                                                            {
                                                                if(idspecif[k,j]==0) next

                                                                if(idspecif[k,j]==1)
                                                                    {
                                                                        bb[avt+1] <- B$best[avt2+1]
                                                                        copcov[avt+1] <- 1
                                                                        copcov2[avt2+1] <- 1
                                                                        
                                                                        avt <- avt+1
                                                                        avt2 <- avt2+1
                                                                    }

                                                                if(idspecif[k,j]==2)
                                                                    {
                                                                        bb[avt+1:ng] <- rep(B$best[avt2+1],ng)
                                                                        diag(vbb[avt+1:ng,avt+1:ng]) <- rep(VB[avt2+1,avt2+1],ng)
                                                                        
                                                                        avt <- avt+ng
                                                                        avt2 <- avt2+1 
                                                                    }
                                                            }
                                                        
                                                    }
                                            }

                                       
                                        for(j in 1:nv0)
                                            {
                                                if(idg0[j]==1)
                                                    {
                                                        if(j==1 & idlink!=-1) next
                                                        else
                                                            {
                                                                bb[avt+1] <- B$best[avt2+1]
                                                                copcov[avt+1] <- 1
                                                                copcov2[avt2+1] <- 1
                                                                
                                                                avt <- avt+1
                                                                avt2 <- avt2+1
                                                            }
                                                    }

                                                if(idg0[j]==2)
                                                    {
                                                        if(j==1 & idlink!=-1)
                                                            {
                                                                avt <- avt+ng-1
                                                            }
                                                        else
                                                            {
                                                                bb[avt+1:ng] <- rep(B$best[avt2+1],ng)
                                                                diag(vbb[avt+1:ng,avt+1:ng]) <- rep(VB[avt2+1,avt2+1],ng)
                                                                        
                                                                avt <- avt+ng
                                                                avt2 <- avt2+1
                                                            }
                                                        
                                                    }
                                            }
                                        

                                        if(nvc>0)
                                            {
                                                bb[avt+1:nvc] <- B$cholesky
                                                copcov[avt+1:nvc] <- 1
                                                copcov2[avt2+1:nvc] <- 1

                                                avt <- avt+nvc
                                                avt2 <- avt2+nvc
                                            }

                                        
                                        bb[(avt+1):length(bb)] <- B$best[(avt2+1):NPM2]
                                        copcov[(avt+1):length(bb)] <- 1
                                        copcov2[(avt2+1):length(B$best)] <- 1

                                        vbb[which(copcov==1),which(copcov==1)] <- VB[which(copcov2==1),which(copcov2==1)]


                                        if(idlink!=-1 & idg0[1]>1)
                                            {
                                                bb <- bb[-(sum(nrisq)-(ng-1)*sum(length(which(risqcom==2)))+nvarxevt+1:(ng-1))]
                                                vbb <- vbb[-(sum(nrisq)-(ng-1)*sum(length(which(risqcom==2)))+nvarxevt+1:(ng-1)),-(sum(nrisq)-(ng-1)*sum(length(which(risqcom==2)))+nvarxevt+1:(ng-1))]
                                            }

                                        Chol <- chol(vbb)
                                        Chol <- t(Chol)
                                        
                                        bb <- bb + Chol %*% rnorm(length(bb))

                                        
                                        b[1:nprob] <- 0

                                        avt <- 0
                                        for(ke in 1:nbevt)
                                            {
                                                if(risqcom[ke]==0)
                                                    {
                                                        b[nprob+sum(nrisq[1:ke])-nrisq[ke]+1:nrisq[ke]] <- bb[avt+1:nrisq[ke]]

                                                        avt <- avt + nrisq[ke]
                                                    }

                                                if(risqcom[ke]==1)
                                                    {
                                                        b[nprob+sum(nrisq[1:ke])-nrisq[ke]+1:nrisq[ke]] <- bb[avt+1:nrisq[ke]]
                                                       
                                                        avt <- avt + nrisq[ke]
                                                    }

                                                if(risqcom[ke]==2)
                                                    {
                                                        b[nprob+sum(nrisq[1:ke])-nrisq[ke]+1:nprisq[ke]] <- bb[avt+1:nprisq[ke]]
                                                        b[nprob+sum(nrisq[1:ke])-nrisq[ke]+nprisq[ke]+1:(ng-1)] <- 1 #bPH

                                                        avt <- avt + nprisq[ke]
                                                    }
                                            }

                                        if(nvarxevt>0) b[nprob+nrisqtot+1:nvarxevt] <- bb[avt+1:nvarxevt]

                                        if(idlink!=-1 & idg0[1]>1)
                                            {
                                                b[nprob+nrisqtot+nvarxevt+1:(ng-1)] <- 0
                                                b[(nprob+nrisqtot+nvarxevt+ng):(nprob+nrisqtot+nvarxevt+nef)] <- bb[avt+nvarxevt+1:(nef-(ng-1))]
                                                nefssI <- nef-ng+1
                                            } 
                                        else
                                            {                                        
                                                b[nprob+nrisqtot+nvarxevt+1:nef] <- bb[avt+nvarxevt+1:nef]
                                                nefssI <- nef
                                            }

                                        if(nvc>0)
                                            {
                                                b[nprob+nrisqtot+nvarxevt+nef+1:nvc] <- bb[nrisqtot-((ng-1)*length(which(risqcom==2)))+nvarxevt+nefssI+1:nvc]
                                            }
                                                
                                        if(nw>0) b[nprob+nrisqtot+nvarxevt+nef+nvc+1:nw] <- 1

                                        b[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+1):NPM] <- bb[(nrisqtot-((ng-1)*length(which(risqcom==2)))+nvarxevt+nefssI+nvc+1):length(bb)]
                                        
                                    }
                            }
                    }
   
            }
        else
            {
                b <- rep(0,NPM)
                for(i in 1:nbevt)
                    {
                        if (typrisq[i]==2)
                            {
                                if(logspecif==1)
                                    {  
                                        b[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- c(rep(c(log(sum(devt==i)/sum(tsurv[devt==i])),0),ifelse(risqcom==0,ng,1)),rep(1,(ng-1)*(risqcom[i]==2)))  
                                    }
                                else
                                    {
                                        b[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- c(rep(c(sqrt(sum(devt==i)/sum(tsurv[devt==i])),1),ifelse(risqcom==0,ng,1)),rep(1,(ng-1)*(risqcom[i]==2)))
                                    }   
                            }
                        else
                            {
                              
                              if(logspecif==1)
                              {  
                                b[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- c(rep(log(1/nprisq[i]),ifelse(risqcom[i]==0,ng*nprisq[i],nprisq[i])),rep(1,(ng-1)*(risqcom[i]==2)))
                              }
                              else
                              {
                                b[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- c(rep(sqrt(1/nprisq[i]),ifelse(risqcom[i]==0,ng*nprisq[i],nprisq[i])),rep(1,(ng-1)*(risqcom[i]==2)))
                              }   
                                                          
                              
                            }
                    } 
                if (nvc>0)
                    {
                        if(idiag==1) b[nprob+nrisqtot+nvarxevt+nef+1:nvc] <- rep(1,nvc)
                        if(idiag==0)
                            {
                                init.nvc <- diag(nea0)
                                init.nvc <- init.nvc[upper.tri(init.nvc, diag=TRUE)]
                                b[nprob+nrisqtot+nvarxevt+nef+1:nvc] <- init.nvc
                            }
                    }
                if(nwg0>0) b[nprob+nrisqtot+nvarxevt+nef+nvc+1:nw] <- 1
                if(ncor0==1) b[nprob+nrisqtot+nvarxevt+nef+nvc+nw+1] <- 1
                if(ncor0==2) b[nprob+nrisqtot+nvarxevt+nef+nvc+nw+1:2] <- c(0,1)

                if(idlink==-1) b[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+1] <- 1
                
                if(idlink==0)
                    {
                        b[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+1] <- mean(Y0)
                        b[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+2] <- 1
                    }
                if(idlink==1)
                    {
                        b[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+1] <- 0
                        b[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+2] <- -log(2)
                        b[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+3] <- 0.7
                        b[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+4] <- 0.1
                    }
                if(idlink==2)
                    {
                        b[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+1] <- -2
                        b[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+2:ntrtot0] <- 0.1
                    }
#cat("B inti pour ng=1 ", b,"\n")
                if(ng>1)
                    {
                        prior02 <- rep(0,ns0)
                        idprob02 <- rep(0,nv0)
                        idg02 <- idg0
                        idg02[idg02==2] <- 1
                        ng02 <- 1
                        nw2 <- 0
                        nef2 <- sum(idg0!=0)
                        if(idlink!=-1) nef2 <- nef2-1
                        idspecif2 <- as.vector(t(idspecif))
                        idspecif2[which(idspecif2==2)] <- 1
                        risqcom2 <- rep(0,nbevt)
                                                  
                        NPM2 <- sum(nprisq)+nvarxevt2+nef2+nvc+ncor0+ntrtot0
                        
                        b2 <- c(sapply(1:nbevt,function(k) b[nprob+sum(nrisq[1:k])-nrisq[k]+1:nprisq[k]]),rep(0,nvarxevt2+nef2),b[nprob+nrisqtot+nvarxevt+nef+1:nvc],b[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+1):length(b)])
                        
                        V2 <- rep(0,NPM2*(NPM2+1)/2)
                        loglik2 <- 0
                        ppi02 <- rep(0,ns0)
                        ppitest02 <- rep(0,nobs0)
                        pred_m_g2 <- rep(0,nobs0)
                        pred_ss_g2 <- rep(0,nobs0)
                        maxiter2 <- min(75,maxiter)
                        convBLG2 <- c(max(0.01,convB), max(0.01,convL),max(0.01,convG))
                        Hr02 <- 0
                        risq_est2 <- matrix(0,nrow=nsim,ncol=nbevt)
                        risqcum_est2 <- matrix(0,nrow=nsim,ncol=nbevt)

                        ## pour reduire le nb d'arguments:
                        int6 <- c(nw2,ncor0,idiag0,idtrunc,logspecif,maxiter2)
                         
                        init <- .Fortran(C_jointhet,
                                         as.double(Y0),
                                         as.double(X0),
                                         as.integer(prior02),
                                         as.double(tsurv0),
                                         as.double(tsurv),
                                         as.integer(devt),
                                         as.integer(ind_survint),
                                         as.integer(idprob02),
                                         as.integer(idea0),
                                         as.integer(idg02),
                                         as.integer(idcor0),
                                         as.integer(idcom),
                                         as.integer(idspecif2),
                                         as.integer(idtdv),
                                         as.integer(idlink),
                                         as.double(epsY),
                                         as.integer(nbzitr),
                                         as.double(zitr0),
                                         as.double(uniqueY0),
                                         as.integer(nvalSPL0),
                                         as.integer(indiceY0),
                                         as.integer(typrisq),
                                         as.integer(risqcom2),
                                         as.integer(nz),
                                         as.double(zi),
                                         as.integer(ns0),
                                         as.integer(ng02),
                                         as.integer(nv0),
                                         as.integer(nobs0),
                                         as.integer(nmes0),
                                         as.integer(nbevt),
                                         as.integer(nea0),
                                         as.integer(int6),
                                         as.integer(NPM2),
                                         best=as.double(b2),
                                         V=as.double(V2),
                                         as.double(loglik),
                                         as.integer(ni),
                                         conv=as.integer(istop),
                                         gconv=as.double(gconv),
                                         as.double(ppi02),
                                         as.double(ppitest02),
                                         as.double(resid_m),
                                         as.double(resid_ss),
                                         as.double(pred_m_g2),
                                         as.double(pred_ss_g2),
                                         as.double(predRE),
                                         as.double(convBLG2),
                                         as.double(time),
                                         as.double(risq_est2),
                                         as.double(risqcum_est2),
                                         as.double(marker),
                                         as.double(transfY),
                                         as.integer(nsim),
                                         as.double(Yobs),
                                         as.double(statglob),
                                         as.double(statevt),
                                         as.integer(pbH0),
                                         as.integer(fix0))

                        ## faire des valeurs initiales a partir de init
#cat("B estime pour ng=1 ",init$best, "\n")
                        ##pr risq
                        for(ke in 1:nbevt)
                            {
                                if(risqcom[ke]==0)
                                    {
                                        indb <- nprob+sum(nrisq[1:ke])-nrisq[ke]+1:nrisq[ke]
                                        indinit <- sum(nprisq[1:ke])-nprisq[ke]+1:nprisq[ke]
                                       if(init$conv==1)
                                           {
                                               b[indb] <- abs(rep(init$best[indinit],ng0)+rep((1:ng0)-(ng0+1)/2,each=nprisq[ke])*rep(sqrt(init$V[indinit*(indinit+1)/2]),ng0))
                                           }
                                       else
                                           {
                                               b[indb] <- abs(rep(init$best[indinit],ng0)+rep((1:ng0)-(ng0+1)/2,each=nprisq[ke])*rep(init$best[indinit],ng0))
                                           }
                                    }

                                if(risqcom[ke]==1)
                                    {
                                        b[nprob+sum(nrisq[1:ke])-nrisq[ke]+1:nrisq[ke]] <- init$best[sum(nrisq[1:ke])-nrisq[ke]+1:nprisq[ke]]
                                    }

                                if(risqcom[ke]==2)
                                    {
                                        b[nprob+sum(nrisq[1:ke])-nrisq[ke]+1:nrisq[ke]] <- c(init$best[sum(nprisq[1:ke])-nprisq[ke]+1:nprisq[ke]],0.5+(0:(ng0-2))*0.5)
                                    }
                            }

                        ## pr bevt
                        avtj <- nprob+nrisqtot
                        avtj2 <- sum(nprisq)
                        for(j in 1:nv0)
                            {
                                if(idcom[j]==0 & all(idspecif[,j]==0)) next
                                
                                if(idcom[j]==1 & all(idspecif[,j]==1))
                                    {
                                        b[avtj+1] <- init$best[avtj2+1]
                                        avtj <- avtj+1
                                        avtj2 <- avtj2+1
                                    }

                                if(idcom[j]==1 & all(idspecif[,j]==2))
                                    {
                                        if(init$conv==1) b[avtj+1:ng0] <- abs(rep(init$best[avtj2+1],ng0)+c(1:ng0-(ng0+1)/2)*rep(sqrt(init$V[(avtj2+1)*(avtj2+2)/2]),ng0))
                                        else b[avtj+1:ng0] <- abs(rep(init$best[avtj2+1],ng0)+c(1:ng0-(ng0+1)/2)*rep(init$best[avtj2+1],ng0))

                                        avtj <- avtj+ng0
                                        avtj2 <- avtj2+1
                                    }

                                if(idcom[j]==0 & idcause[j]!=0)
                                    {
                                        for(k in 1:nbevt)
                                            {
                                                if(idspecif[k,j]==0) next

                                                if(idspecif[k,j]==1)
                                                    {
                                                        b[avtj+1] <- init$best[avtj2+1]
                                                        avtj <- avtj+1
                                                        avtj2 <- avtj2+1
                                                    }

                                                if(idspecif[k,j]==2)
                                                    {
                                                        if(init$conv==1) b[avtj+1:ng0] <- abs(rep(init$best[avtj2+1],ng0)+c(1:ng0-(ng0+1)/2)*rep(sqrt(init$V[(avtj2+1)*(avtj2+2)/2]),ng0))
                                                        else b[avtj+1:ng0] <- abs(rep(init$best[avtj2+1],ng0)+c(1:ng0-(ng0+1)/2)*rep(init$best[avtj2+1],ng0))

                                                        avtj <- avtj+ng0
                                                        avtj2 <- avtj2+1 
                                                    }
                                            }
                                        
                                    }
                                    


                                
                                ## if(idxevt0[j]==1)
                                ##     {
                                ##         nn <- ifelse(idcausespecif[j]==1,nbevt,1)
                                ##         b[nprob+nrisqtot+avtj+1:nn] <- init$best[sum(nprisq)+avtj2+1:nn]

                                ##         avtj <- avtj+nn
                                ##         avtj2 <- avtj2+nn
                                ##     }

                                ## if(idxevt0[j]==2)
                                ##     {
                                ##         nn <- ng0*ifelse(idcausespecif[j]==1,nbevt,1)
                                ##         if(init$conv==1) b[nprob+nrisqtot+avtj+1:nn] <- abs(rep(init$best[sum(nprisq)+avtj2+1:(nn/ng0)],ng0)+rep((1:ng0)-(ng0+1)/2,nn/ng0)*rep(init$V[(sum(nprisq)+avtj2+1:(nn/ng0))*(sum(nprisq)+avtj2+1:(nn/ng0)+1)/2],ng0))
                                ##         else b[nprob+nrisqtot+avtj+1:nn] <- abs(rep(init$best[sum(nprisq)+avtj2+1:(nn/ng0)],ng0)+rep((1:ng0)-(ng0+1)/2,nn/ng0)*rep(init$best[sum(nprisq)+avtj2+1:(nn/ng0)],ng0))

                                ##         avtj <- avtj+nn
                                ##         avtj2 <- avtj2+nn/ng0
                                ##         #print(b[nprob+nrisqtot+avtj+1:nn])
                                ##         #print(rep(init$best[sum(nprisq)+avtj2+1:(nn/ng0)],ng0))
                                ##         #print(rep((1:ng0)-(ng0+1)/2,nn/ng0))
                                ##         #print(rep(init$V[(sum(nprisq)+avtj2+1:nn/ng0)*(sum(nprisq)+avtj2+1:nn/ng0+1)/2],ng0))
                                ##     }
                            }

                        
                        ## pr nef
                        avtj <- nprob+nrisqtot+nvarxevt
                        avtj2 <- sum(nprisq)+nvarxevt2
                        for(j in 1:nv0)
                            {
                                if(idg0[j]==1)
                                    {
                                        if(j==1 & idlink!=-1) next
                                        else
                                            {
                                                if(init$conv==1) b[avtj+1] <- init$best[avtj2+1]
                                                avtj <- avtj+1
                                                avtj2 <- avtj2+1
                                            }
                                    }

                                if(idg0[j]==2)
                                    {
                                        if(j==1)
                                            {
                                                if(idlink!=-1)
                                                    {
                                                        b[avtj+1:(ng0-1)] <- -0.5*(1:(ng0-1))

                                                        avtj <- avtj+ng0-1
                                                    }
                                                else
                                                    {
                                                        b[avtj+1:ng0] <- init$best[avtj2+1]+(1:ng0-(ng0+1)/2)*init$best[avtj2+1]
                                                        avtj <- avtj+ng0
                                                        avtj2 <- avtj2+1
                                                    }
                                            }
                                        else
                                            {
                                                if(init$conv==1) b[avtj+1:ng0] <- init$best[avtj2+1]+(1:ng0-(ng0+1)/2)*sqrt(init$V[(avtj2+1)*(avtj2+1+1)/2])
                                                else b[avtj+1:ng0] <- init$best[avtj2+1]+(1:ng0-(ng0+1)/2)*init$best[avtj2+1]

                                                avtj <- avtj+ng0
                                                avtj2 <- avtj2+1
                                            }
                                        
                                    }
                            }

                        ## pr nvc
                        if(nvc>0)
                            {
                                b[nprob+nrisqtot+nvarxevt+nef+1:nvc] <- init$best[sum(nprisq)+nvarxevt2+nef2+1:nvc]
                            }

                        ## cor et transfo
                        b[(nprob+nrisqtot+nvarxevt+nef+nvc+nw+1):NPM] <- init$best[(sum(nprisq)+nvarxevt2+nef2+nvc+1):NPM2]
                    }
            }
#cat("B init a partir du modele pour ng=1 ", b,"\n")

### nom au vecteur best

        nom.X0 <- colnames(X0)
        nom.X0[which(nom.X0=="(Intercept)")] <- "intercept"
        
        ##prm classmb   
        if(ng0>=2)
            {
                nom <- rep(nom.X0[idprob0==1],each=ng0-1)
                nom1 <- paste(nom," class",c(1:(ng0-1)),sep="")
                names(b)[1:nprob] <- nom1
            }

        ##prm fct de risque
        if(isTRUE(logscale))
            {
                for(i in 1:nbevt)
                    {
                        nom1 <- rep(paste("event",i,sep=""),nrisq[i])
                        if(typrisq[i]==2)
                            {
                                nom2 <- paste(nom1[1:2]," log(Weibull",1:2,")",sep="")  
                                nom1[1:(2*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
                                if(risqcom[i]==2) nom1[2+1:(ng-1)] <- paste(nom1[2+1:(ng-1)],"SurvPH") 
                                names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1 
                            }
                        if(typrisq[i]==1)  
                            {
                                nom2 <- paste(nom1[1:(nz[i]-1)]," log(piecewise",1:(nz[i]-1),")",sep="")  
                                nom1[1:((nz[i]-1)*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
                                if(risqcom[i]==2) nom1[nz[i]-1+1:(ng-1)] <- paste(nom1[nz[i]-1+1:(ng-1)],"SurvPH")  
                                names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1  
                            }
                        if(typrisq[i]==3)  
                            {
                                nom2 <- paste(nom1[1:(nz[i]-1)]," log(splines",1:(nz[i]+2),")",sep="")  
                                nom1[1:((nz[i]+2)*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
                                if(risqcom[i]==2) nom1[nz[i]+2+1:(ng-1)] <- paste(nom1[nz[i]+2+1:(ng-1)],"SurvPH")
                                names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1  
                            }  
                    }
            }
        else
            {
                for(i in 1:nbevt)
                    {
                        nom1 <- rep(paste("event",i,sep=""),nrisq[i])
                        if(typrisq[i]==2)
                            {
                                nom2 <- paste(nom1[1:2]," +/-sqrt(Weibull",1:2,")",sep="")  
                                nom1[1:(2*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
                                if(risqcom[i]==2) nom1[2+1:(ng-1)] <- paste(nom1[2+1:(ng-1)],"SurvPH")
                                names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1  
                            }
                        if(typrisq[i]==1)  
                            {
                                nom2 <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(piecewise",1:(nz[i]-1),")",sep="")  
                                nom1[1:((nz[i]-1)*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
                                if(risqcom[i]==2) nom1[nz[i]-1+1:(ng-1)] <- paste(nom1[nz[i]-1+1:(ng-1)],"SurvPH")
                                names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1  
                            }
                        if(typrisq[i]==3)  
                            {
                                nom2 <- paste(nom1[1:(nz[i]-1)]," +/-sqrt(splines",1:(nz[i]+2),")",sep="")  
                                nom1[1:((nz[i]+2)*ifelse(risqcom[i]==0,ng,1))] <- rep(nom2,ifelse(risqcom[i]==0,ng*(risqcom[i]==0),1))
                                if(risqcom[i]==2) nom1[nz[i]+2+1:(ng-1)] <- paste(nom1[nz[i]+2+1:(ng-1)],"SurvPH") 
                                names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- nom1  
                            }  
                    }   
            }
        if(ng>1)
            {
                for(i in 1:nbevt)
                    {
                        if(risqcom[i]==1) next;
                        if(risqcom[i]==0)
                            {
                                names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- paste(names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]],paste("class",rep(1:ng,each=nprisq[i]))) 
                            }
                        if(risqcom[i]==2)
                            {
                                names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]] <- paste(names(b)[nprob+sum(nrisq[1:i])-nrisq[i]+1:nrisq[i]],c(rep("",nprisq[i]),paste(" class",1:(ng-1),sep="")),sep="") 
                            }    
                    }
            }

        ##prm covariables survival
        nom1 <- NULL  
        for(j in 1:nv0)
            {
                if(idcom[j]==0 & all(idspecif[,j]==0)) next
                
                if(idcom[j]==1 & all(idspecif[,j]==1)) #X
                    {
                        if(idtdv[j]==1)
                            {
                                nom1 <- c(nom1,paste("I(t>",nom.timedepvar,")",sep=""))
                            }
                        else
                            {
                                nom1 <- c(nom1,nom.X0[j])
                            }
                        next
                    }

                if(idcom[j]==1 & all(idspecif[,j]==2)) #mixture(X)
                    {
                        if(idtdv[j]==1)
                            {
                                nom1 <- c(nom1,paste("I(t>",nom.timedepvar,") class",1:ng0,sep=""))
                            }
                        else
                            {                  
                                nom1 <- c(nom1,paste(nom.X0[j],paste("class",1:ng0,sep="")))
                            }
                        next
                    }

                if(idcom[j]==0 & all(idspecif[,j]==1)) #cause(X)
                    {
                        if(idtdv[j]==1)
                            {
                                nom1 <- c(nom1,paste("I(t>",nom.timedepvar,") event",1:nbevt,sep=""))
                            }
                        else
                            {                  
                                nom1 <- c(nom1,paste(nom.X0[j],paste("event",1:nbevt,sep="")))
                            }
                        next
                    }

                if(idcom[j]==0 & all(idspecif[,j]==2)) #cause(mixture(X))
                    {
                        if(idtdv[j]==1)
                            {
                                xevt <- paste("I(t>",nom.timedepvar,") event",1:nbevt,sep="")
                                classg <- paste("class",1:ng0,sep="")
                                nom1 <- c(nom1,paste(rep(xevt,each=ng0),rep(classg,nbevt)))
                            }
                        else
                            {
                                xevt <- paste(nom.X0[j],paste("event",1:nbevt,sep=""))
                                classg <- paste("class",1:ng0,sep="")
                                nom1 <- c(nom1,paste(rep(xevt,each=ng0),rep(classg,nbevt)))
                            }
                        next
                    }

                

                if(idcom[j]==0 & idcause[j]!=0) #causek
                    {
                        for(k in 1:nbevt)
                            {
                                if(idspecif[k,j]==0) next
                                
                                if(idspecif[k,j]==1)
                                    {
                                        if(idtdv[j]==1)
                                            {
                                                nom1 <- c(nom1,paste("I(t>",nom.timedepvar,") event",k,sep=""))
                                            }
                                        else
                                            {
                                                nom1 <- c(nom1,paste(nom.X0[j],paste("event",k,sep="")))
                                            }
                                        next
                                    }
                                
                                if(idspecif[k,j]==2)
                                    {
                                        if(idtdv[j]==1)
                                            {
                                                xevtk <- paste("I(t>",nom.timedepvar,") event",k,sep="")
                                                classg <- paste("class",1:ng0,sep="")
                                                nom1 <- c(nom1,paste(xevtk,classg))
                                            }
                                        else
                                            {
                                                xevtk <- paste(nom.X0[j],paste("event",k,sep=""))
                                                classg <- paste("class",1:ng0,sep="")
                                                nom1 <- c(nom1,paste(xevtk,classg))
                                            }
                                        next
                                    }
                            }
                        
                    }                
            }
        
        if(nvarxevt>0) names(b)[nprob+nrisqtot+1:nvarxevt] <- nom1 
        ##NB : pour chaque variable, coef ranges par evenement puis par classe  


        ##prm fixed   
        if(ng0==1)
            {
                if(idlink!=-1)
                    {
                        names(b)[nprob+nrisqtot+nvarxevt+1:nef] <- nom.X0[-1][idg0[-1]!=0]
                    }
                else
                    {
                       # print(nprob+nrisqtot+nvarxevt)
                       # print(nef)
                       # print(nom.X0)
                       # print(idg0)
                        names(b)[nprob+nrisqtot+nvarxevt+1:nef] <- nom.X0[idg0!=0]
                    }
            }
        if(ng0>1)
            {
                nom1<- NULL
                for(i in 1:nv0) 
                    {
                        if(idg0[i]==2)
                            {
                                if(i==1)
                                    {
                                        if(idlink==-1)
                                            {
                                                nom <- paste("intercept class",c(1:ng0),sep="")
                                                nom1 <- c(nom1,nom)
                                            }
                                        else
                                            {
                                                nom <- paste("intercept class",c(2:ng0),sep="")
                                                nom1 <- c(nom1,nom)
                                            }
                                    }
                                if (i>1)
                                    {
                                        nom <- paste(nom.X0[i]," class",c(1:ng0),sep="")
                                        nom1 <- c(nom1,nom)
                                    }
                            }
                        if(idg0[i]==1 & i==1 & idlink==-1) nom1 <- c(nom1,nom.X0[i])
                        if(idg0[i]==1 & i>1) nom1 <- c(nom1,nom.X0[i])
                        
                    }
                names(b)[nprob+nrisqtot+nvarxevt+1:nef]<- nom1
            }

        ##prm random et nwg   
        if(nvc!=0) names(b)[nprob+nrisqtot+nvarxevt+nef+1:nvc] <- paste("varcov",c(1:nvc))
        if(nw!=0) names(b)[nprob+nrisqtot+nvarxevt+nef+nvc+1:nw] <- paste("varprop class",c(1:(ng0-1))) 

        ##prm cor   
        if(ncor0>0) {names(b)[nprob+nrisqtot+nvarxevt+nef+nvc+nw+1:ncor0] <- paste("cor",1:ncor0,sep="")}   

        ##prm link
        if(idlink==-1) names(b)[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+1] <- "stderr"
        if(idlink==0) names(b)[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+1:2]<- c("Linear 1","Linear 2")
        if(idlink==1) names(b)[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+1:4]<- paste("Beta",1:4,sep="")
        if(idlink==2) names(b)[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor0+1:ntrtot0]<- paste("I-splines",c(1:ntrtot0),sep="")

#browser()
#print(b)
#### estimation
        idspecif <- as.vector(t(idspecif))
#return(b)
        ## pour reduire le nb d'arguments:
        int6 <- c(nw,ncor0,idiag0,idtrunc,logspecif,maxiter)
        convBLG <- c(convB,convL,convG)
        
        out <- .Fortran(C_jointhet,
                        as.double(Y0),
                        as.double(X0),
                        as.integer(prior0),
                        as.double(tsurv0),
                        as.double(tsurv),
                        as.integer(devt),
                        as.integer(ind_survint),
                        as.integer(idprob0),
                        as.integer(idea0),
                        as.integer(idg0),
                        as.integer(idcor0),
                        as.integer(idcom),
                        as.integer(idspecif),
                        as.integer(idtdv),
                        as.integer(idlink),
                        as.double(epsY),
                        as.integer(nbzitr),
                        as.double(zitr0),
                        as.double(uniqueY0),
                        as.integer(nvalSPL0),
                        as.integer(indiceY0),
                        as.integer(typrisq),
                        as.integer(risqcom),
                        as.integer(nz),
                        as.double(zi),
                        as.integer(ns0),
                        as.integer(ng0),
                        as.integer(nv0),
                        as.integer(nobs0),
                        as.integer(nmes0),
                        as.integer(nbevt),
                        as.integer(nea0),
                        as.integer(int6),
                        as.integer(NPM),
                        best=as.double(b),
                        V=as.double(V),
                        loglik=as.double(loglik),
                        niter=as.integer(ni),
                        conv=as.integer(istop),
                        gconv=as.double(gconv),
                        ppi=as.double(ppi0),
                        ppitest=as.double(ppitest0),
                        resid_m=as.double(resid_m),
                        resid_ss=as.double(resid_ss),
                        pred_m_g=as.double(pred_m_g),
                        pred_ss_g=as.double(pred_ss_g),
                        predRE=as.double(predRE),
                        as.double(convBLG),
                        time=as.double(time),
                        risq_est=as.double(risq_est),
                        risqcum_est=as.double(risqcum_est),
                        marker=as.double(marker),
                        transfY=as.double(transfY),
                        as.integer(nsim),
                        Yobs=as.double(Yobs),
                        statglob=as.double(statglob),
                        statevt=as.double(statevt),
                        as.integer(pbH0),
                        as.integer(fix0))



 ### mettre NA pour les variances et covariances non calculees et  0 pr les prm fixes
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

 
        
        
### remplacer cholesky par varcov dans best
        if(nvc>0)
            {
                ch <- out$best[nprob+nrisqtot+nvarxevt+nef+1:nvc]
                if(idiag0==0)
                    {
                        U <- matrix(0,nea0,nea0)
                        U[upper.tri(U,diag=TRUE)] <- ch
                        z <- t(U) %*% U
                        out$best[nprob+nrisqtot+nvarxevt+nef+1:nvc] <- z[upper.tri(z,diag=TRUE)]
                    }
                if(idiag0==1)
                    {
                        out$best[nprob+nrisqtot+nvarxevt+nef+1:nvc] <- out$best[nprob+nrisqtot+nvarxevt+nef+1:nvc]**2
                    }
            }
        else
            {
                ch <- NA
            }

        names(out$best) <- names(b)
        nom.unique[which(nom.unique=="(Intercept)")] <- "intercept"

### predictions des effets aleatoires
        if (nea0>0)
            {
                predRE <- matrix(out$predRE,ncol=nea0,byrow=TRUE)
                predRE <- data.frame(unique(IND),predRE)
                colnames(predRE) <- c(subject,nom.unique[idea0!=0])
            }

### proba des classes a posteroiri et classification
        if(ng0>1)
            {
                ppi <- matrix(out$ppi,ncol=ng0,nrow=ns0,byrow=TRUE)
                ppitest <- matrix(out$ppitest,ncol=ng0,byrow=TRUE)
            }
        else
            {
                ppi <- matrix(rep(1,ns0),ncol=ng0)
                ppitest <- matrix(rep(1,ns0),ncol=ng0)
            }
        
        if(!(out$conv %in% c(1,2)))
            {
                classif <- rep(NA,ns0)
            }
        else
            {
                classif <- apply(ppi,1,which.max)

                if(any(!is.finite(ppi)))
                    {
                        classif <- rep(NA,ns0)
                        iok <- which(is.finite(ppi[,1]))
                        classif[iok] <- apply(ppi[iok,,drop=FALSE],1,which.max)
                    }
            }

        ppi <- data.frame(unique(IND),classif,ppi)
        temp <- paste("probYT",1:ng0,sep="")
        colnames(ppi) <- c(subject,"class",temp)
        rownames(ppi) <- 1:ns0
  
   ## faire pareil pour ppitest

        if(!(out$conv %in% c(1,2)))
            {
                classif <- rep(NA,ns0)
            }
        else
            {
                classif <- apply(ppitest,1,which.max)

                if(any(!is.finite(ppitest)))
                    {
                        classif <- rep(NA,ns0)
                        iok <- which(is.finite(ppitest[,1]))
                        classif[iok] <- apply(ppitest[iok,,drop=FALSE],1,which.max)
                    }
            }
 
        ppitest <- data.frame(unique(IND),classif,ppitest)
        temp <- paste("probY",1:ng0,sep="")
        colnames(ppitest) <- c(subject,"class",temp)
        rownames(ppitest) <- 1:ns0  
  

### predictions marginales et subject-specifiques
        pred_m_g <- matrix(out$pred_m_g,nrow=nobs0,ncol=ng0)
        pred_ss_g <- matrix(out$pred_ss_g,nrow=nobs0,ncol=ng0)
        
        if((out$conv %in% c(1,2)))
            {
                pred_m <- out$Yobs-out$resid_m
                pred_ss <- out$Yobs-out$resid_ss
            }
        else
            {
                pred_m <- rep(NA,nobs0)
                pred_ss <- rep(NA,nobs0)
            } 

        pred <- data.frame(IND,pred_m,out$resid_m,pred_ss,out$resid_ss,out$Yobs,pred_m_g,pred_ss_g)

        temp<-paste("pred_m",1:ng0,sep="")
        temp1<-paste("pred_ss",1:ng0,sep="")
        colnames(pred)<-c(subject,"pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1)


### risques
        risqcum_est <- matrix(out$risqcum_est,nrow=nsim,ncol=ng0*nbevt)
        risq_est <- matrix(out$risq_est,nrow=nsim,ncol=ng0*nbevt)
        predSurv <- cbind(time,risq_est,risqcum_est)
        
        temp <- paste(paste("event",rep(1:nbevt,each=ng0),".RiskFct",sep=""),1:ng0,sep="")
        temp1 <- paste(paste("event",rep(1:nbevt,each=ng0),".CumRiskFct",sep=""),1:ng0,sep="")
        colnames(predSurv) <- c("time",temp,temp1)
        rownames(predSurv) <- 1:nsim

###estimlink
        estimlink <- cbind(out$marker,out$transfY)
        colnames(estimlink) <- c(nomY,paste("transf",nomY,sep=".")) 

### score test
        if(out$conv!=1) stats <- rep(NA,nbevt+1)
        else
            {
                stats <- c(out$statglob,out$statevt)
                stats[which((!is.finite(stats)) | (stats==9999))] <- NA
            }
        if(nbevt==1) stats <- stats[1]
        
### N
        N <- NULL
        N[1] <- nprob
        N[2] <- nrisqtot
        N[3] <- nvarxevt
        N[4] <- nef
        N[5] <- nvc
        N[6] <- nw
        N[7] <- ncor0
        N[8] <- ntrtot0
        N[9] <- nobs0

        nevent <- rep(0,nbevt)
        for(ke in 1:nbevt)
            {
                nevent[ke] <- length(which(devt==ke))
            }
        N <- c(N,nevent)

### noms des variables
        Names <- list(Xnames=nom.unique,Xnames2=ttesLesVar,Yname=nomY,
                      ID=nom.subject,Tnames=noms.surv,prior.name=nom.prior,
                      TimeDepVar.name=nom.timedepvar)

        res <-list(ns=ns0,ng=ng0,idprob=idprob0,idcom=idcom,
                   idspecif=idspecif,idtdv=idtdv,idg=idg0,idea=idea0,
                   idcor=idcor0,loglik=out$loglik,best=out$best,V=V,
                   gconv=out$gconv,conv=out$conv,call=cl,niter=out$niter,
                   N=N,idiag=idiag0,pred=pred,pprob=ppi,pprobY=ppitest,
                   predRE=predRE,Names=Names,cholesky=ch,logspecif=logspecif,
                   estimlink=estimlink,epsY=epsY,linktype=idlink,linknodes=zitr0,
                   predSurv=predSurv,hazard=list(typrisq,hazardtype,zi,nz),
                   scoretest=stats,na.action=linesNA,
                   AIC=2*(length(out$best)-length(posfix)-out$loglik),
                   BIC=(length(out$best)-length(posfix))*log(ns0)-2*out$loglik,
                   data=datareturn)

        class(res) <-c("Jointlcmm")

        cost<-proc.time()-ptm
        if(verbose==TRUE) cat("The program took", round(cost[3],2), "seconds \n")


        return(res)
    }


#' @rdname Jointlcmm
#' @export
jlcmm <- Jointlcmm
