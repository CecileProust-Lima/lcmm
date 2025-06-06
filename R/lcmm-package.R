



#' Estimation of extended mixed models using latent classes and latent
#' processes.
#' 
#' Functions for the estimation of latent class mixed models (LCMM), joint
#' latent class mixed models for longitudinal and survival data (JLCM)
#' and latent process mixed models (with or without latent classes of trajectory) 
#' for univariate and multivariate longitudinal outcomes of different types 
#' including curvilinear and ordinal outcomes. 
#' All the models are estimated in a maximum likelihood framework using an 
#' iterative algorithm. The package also provides various post fit functions.
#'
#' 
#' \Sexpr[stage=build,results=hide]{descr <- packageDescription("lcmm")}
#' \tabular{ll}{ Package: \tab lcmm \cr Type: \tab Package \cr 
#' Version: \tab \Sexpr[stage=build]{descr$Version} \cr 
#' Date: \tab \Sexpr[stage=build]{descr$Date}\cr 
#' License: \tab GPL (>=2.0) \cr LazyLoad: \tab yes \cr }
#' The package includes for the moment the estimation of :
#' \itemize{ \item latent class mixed models for Gaussian longitudinal outcomes
#' using \code{hlme} function, \item latent class mixed models for other
#' quantitative, bounded quantitative (curvilinear) and discrete (ordinal/binary) longitudinal
#' outcomes using \code{lcmm} function, \item mixed models (with and without latent classes) for
#' multivariate longitudinal outcomes of different nature using
#' \code{multlcmm} function (this includes a longitudinal IRT model for homogeneous and heterogeneous data),
#'  \item joint latent class mixed models for a
#' Gaussian (or curvilinear) longitudinal outcome and a right-censored
#' (potentially left-truncated and of multiple causes) time-to-event using
#' \code{Jointlcmm} function, \item joint latent class mixed models for multivariate longitudinal outcomes
#'  and a right-censored (potentially left-truncated and of multiple causes) time-to-event using
#' \code{mpjlcmm} function. }
#' 
#' Please report any bug or comment regarding the package for
#' future updates VIA GITHUB ONLY.
#' 
#' @name lcmm-package
#####@docType package
#' @author Cecile Proust-Lima, Viviane Philipps, Amadou Diakite and Benoit
#' Liquet
#' 
#' \email{cecile.proust-lima@@inserm.fr}
#' @references
#' 
#' Proust-Lima C, Philipps V, Liquet B (2017). Estimation of Extended Mixed 
#' Models Using Latent Classes and Latent Processes: The R Package lcmm. 
#' Journal of Statistical Software, 78(2), 1-56. doi:10.18637/jss.v078.i02
#' 
#' Lin, Turnbull, McCulloch and Slate (2002). Latent class models for joint
#' analysis of longitudinal biomarker and event process data: application to
#' longitudinal prostate-specific antigen readings and prostate cancer. Journal
#' of the American Statistical Association 97, 53-65.
#' 
#' Muthen and Shedden (1999). Finite mixture modeling with mixture outcomes
#' using the EM algorithm. Biometrics 55, 463-9
#' 
#' Proust and Jacqmin-Gadda (2005). Estimation of linear mixed models with a
#' mixture of distribution for the random-effects. Comput Methods Programs
#' Biomed 78:165-73
#' 
#' Proust, Jacqmin-Gadda, Taylor, Ganiayre, and Commenges (2006). A nonlinear
#' model with latent process for cognitive evolution using multivariate
#' longitudinal data. Biometrics 62, 1014-24.
#' 
#' Proust-Lima, Dartigues and Jacqmin-Gadda (2011). Misuse of the linear mixed
#' model when evaluating risk factors of cognitive decline. Amer J Epidemiol
#' 174(9), 1077-88
#' 
#' Proust-Lima and Taylor (2009). Development and validation of a dynamic
#' prognostic tool for prostate cancer recurrence using repeated measures of
#' post-treatment PSA: a joint modelling approach. Biostatistics 10, 535-49.
#' 
#' Proust-Lima, Sene, Taylor, Jacqmin-Gadda (2014). Joint latent class models
#' for longitudinal and time-to-event data: a review. Statistical Methods in
#' Medical Research 23, 74-90.
#' 
#' Proust-Lima, Amieva, Jacqmin-Gadda (2013). Analysis of multivariate mixed
#' longitudinal data: A flexible latent process approach. Br J Math Stat
#' Psychol 66(3), 470-87.
#'  
#' Proust-Lima, Philipps, Perrot, Blanchin, Sebille (2021). Modeling repeated self-reported
#' outcome data: a continuous-time longitudinal Item Response Theory model. 
#' arXiv:210913064. http://arxiv.org/abs/2109.13064
#'  
#' Proust-Lima, Dartigues, Jacqmin-Gadda (2016). Joint modeling of repeated multivariate 
#' cognitive measures and competing risks of dementia and death: a latent process and
#'  latent class approach. Stat Med;35(3):382-98
#'  
#'  Proust-Lima, Philipps, Dartigues, Bennett, Glymour, Jacqmin-Gadda, et al (2019). 
#'  Are latent variable models preferable to composite score approaches when assessing
#'  risk factors of change? Evaluation of type-I error and statistical power in longitudinal
#'  cognitive studies. Stat Methods Med Res;28(7):1942-57  
#'       
#' Verbeke and Lesaffre (1996). A linear mixed-effects model with heterogeneity
#' in the random-effects population. Journal of the American Statistical
#' Association 91, 217-21
"_PACKAGE"

#' 
#' @keywords package
#' @importFrom graphics axis hist lines matlines matplot mtext par plot points segments polygon
#' @importFrom grDevices rainbow rgb col2rgb n2mfrow
#' @importFrom stats as.formula formula get_all_vars integrate median model.frame model.matrix na.fail na.omit na.pass pchisq pnorm qnorm quantile rnorm sd terms residuals vcov fitted coef update reformulate cov2cor simulate rbinom runif uniroot rweibull
#' @importFrom survival Surv untangle.specials
#' @importFrom mvtnorm rmvnorm
#' @importFrom spacefillr generate_sobol_owen_set
#' @importFrom parallel makeCluster clusterExport stopCluster parLapply clusterEvalQ clusterSetRNGStream parSapply parApply
#' @importFrom doParallel registerDoParallel
#' @importFrom nlme ranef fixef
#' @importFrom marqLevAlg mla
#' @importFrom utils capture.output
#' @importFrom numDeriv hessian
#' @useDynLib lcmm, .registration=TRUE, .fixes="C_"
NULL









