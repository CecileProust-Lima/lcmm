# Changes in Version 2.1.0 (2023-10-06) :

* new function `externVar` to perform a secondary regression analysis after the estimation of a primary latent class model
* new argument `pprior` in hlme, lcmm, multlcmm and Jointlcmm to fix the probability to belong to each latent class
* packages survival, parallel, mvtnorm, randtoolbox, marqLevAlg, doParallel, numDeriv are now listed in Imports rather than in Depends
* no subject-specific predictions in multlcmm with ordinal outcomes
* corrections in mpjlcmm
* correction in predictL without random effects
* correction in epoce and predictY.Jointlcmm
* use of R's random number generator in Fortran codes
* use double precision rather than real(kind=8) in Fortran


# Changes in Version 2.0.2 (2023-02-17)

* all vignettes excepted the introduction vignette (now renamed lcmm.Rmd) are removed from the CRAN version because of too long check time. 
* We now provide a website at https://CecileProust-Lima.github.io/lcmm


# Changes in Version 2.0.1 (2023-02-01) :

* new vignette `Joint latent class model with Jointlcmm`
* new vignette `Multivariate latent class model with mpjlcmm`
* new argument `pprior` in the `hlme` function
* new argument `computeDiscrete` in the `lcmm` function 
* `mpjlcmm` can be used with a mix of hlme/lcmm/multlcmm objects
* `summarytable` and `summaryplot` implement two versions of ICL criterion
* new output `levels` in all estimating functions
* new output `varRE` in `hlme`
* check the convergence of the initial model when using B=random()
* random parameters are generated with rmnvorm instead of using the Cholesky transformation
* `permut`, `cuminc`, `VarCov`, `coef`, `vcov` functions are available for mpjlcmm objects
* corrections in `mpjlcmm`, especially with competing risks
* correction in residuals for Jointlcmm models
* bug fixed when using `posfix` and `partialH` simultaneously
* correction in the likelihood for mutlcmm models
* bug fixed in `predictClass` and `predictRE` when using splines
* verbose=FALSE by default


# Changes in Version 2.0.0 (2022-06-15) :

* the model's estimation is now available in parallel mode! 
* The optimization relies on the parallelized marqLevAlg R package.
* models with latent classes (ng>1) require initial values
* the `hlme` function has now a pprior argument
* the `mpjlcmm` function can be used without a time-to-event model
* the `summary` functions now shorten the parameters names
* the log-likelihood functions are now exported
* bug fixed in `mpjlcmm` when no random effect is included
* bug fixed in `Jointlcmm` with Weibull hazards and competing risks
* bug fixed in `permut` when used on `Jointlcmm` objects with competing risks
* correction of the outputs of `multlcmm` models


# Changes in Version 1.9.4 (2022-01-03) :

* the multlcmm function is now available for ordinal outcomes (link="thresholds") providing a longitudinal IRT model!
* new vignette `Dynamic IRT with multlcmm`
* new dataset simdataHADS
* new function `simulate` to simulate a dataset from a hlme, lcmm, multlcmm or Jointlcmm model
* new functions `ItemInfo` and `plot.ItemInfo` to compute and plot Fisher information for ordinal outcomes
* new argument `var.time` in the hlme, lcmm, multlcmm and Jointlcmm functions (used in plot(, which="fit"); issue #91)
* fix CRAN error with as.vector.data.frame
* correction in the `permut` function (transformation parameters were not updated)
* add envir=parent.frame() in permut and gridsearch to enable the use of these functions in a parallel setting
* fix bug in the estimation functions with infinite posterior probabilities
* the `gridsearch` function now checks that the initial model converged (ie minit$conv=1)
* the `fixef` and `ranef` function are now imported from the nlme package


# Changes in Version 1.9.3 (2021-06-17):

* new functions `predictClass`, `predictRE` and `summaryplot`
* ICL computation in `summaryplot`
* use of `rmvnorm` in `multlcmm` to generate random initial values
* `maxiter` is used in the estimation of the final model in `gridsearch`
* fix bug in `cuminc` without covariates
* fix bug in the check for numeric type for argument `subject` with tibbles
* fix bug in `predictY` with hlme object when the dataset is named "x"
* fix bug in the `update` function when the model has unestimated parameters (posfix)
* fix bug in `hlme` when posterior probabilities are NA
* fix bug in `plot` with option which="fit" (observations at the maximum time measurement where not systematically included)
* correction in the outputs (ppi and resid) of the `mpjlcmm` function

# Changes in Version 1.9.2:

* event variable in joint models can be logical
* bug fixed in `Jointlcmm` with prior when there are missing data
* bug fixed in `mpjlcmm` : initial values were badly modified (with at least 3 dimensions)
* small bugs fixed in `predictY` with median=TRUE

# Changes in Version 1.9.1:

* parallel implementation of `gridsearch` function. Thanks to Raphael Peter for his suggestion.
* add `condRE_Y` option in `predictYcond`
* add `median` options in `predictY`
* corrections in `Jointlcmm`, `multlcmm` and `mpjlcmm` when prior is specified
* bugs fixed in some prediction functions
* small bugs fixed in the summary when some parameters are not estimated
* bug fixed in `VarExpl` with models including BM or AR
* bug fixed in `update.mpjlcmm` (variance matrix was not correct)
* manage infinite ppi in `hlme`
* correction of epsY type, URL in vignettes, data statements position

# Changes in Version 1.8.1:

* new function `mpjlcmm` for estimating joint latent class models with multiple markers and/or latent processes
* various post-fit functions for `mpjlcmm` objects
* new functions `permut` and `xclass`
* creation of vignettes, thanks to Samy Youbi for his help
* variable `subject` must be numeric
* in plot(which='fit'), time intervals do not depend on subset
* add score test result in summarytable
* bug fixed in `lcmm` with prior
* bug fixed in `Jointlcmm` with infinite score test
* bug fixed in `dynpred` with TimeDepVar

# Changes in Version 1.7.9:

* bug in summary when the model did not converge
* bug in dynpred when draws=TRUE and only 1 horizon or 1 landmark, or when o covariates are included in the survival model, or when using factor
* bug in Jointlcmm when using B=m1
* bug in plot.predictY with CI
* bug in Jointlcmm when B=random(m1)


# Changes in Version 1.7.8:

* shades in plot.predictlink/L/Y
* subset in plot, which="fit"


# Changes in Version 1.7.6 (2016-12-12):

* Small bugs identified and solved in multlcmm


# Changes in Version 1.7.5 (2016-03-15):

* Small bugs identified and solved in multlcmm, predictY and predictL


# Changes in Version 1.7.4 (2015-12-26):

* The package uses lazydata to automatically load the datasets of the package.

* `jlcmm` and `mlcmm` are shortcuts for functions `Jointlcmm` and `multlcmm`, respectively.

* Function `gridsearch` provides an automatic grid of departures for reducing the odds of converging towards a local maximum. 

* Initial values can be randomly generated from a model with 1 class (called m1 in next example) with option B=random(m1) in hlme, lcmm, multlcmm and Jointlcmm.


# Changes in Version 1.7.3.0 (2015-10-23):


* Functions `hlme`, `lcmm`, `multlcmm`, `Jointlcmm` now include a posfix option to specify parameters that should not be estimated.

* Functions `lcmm`, `multlcmm`, `Jointlcmm` now include a partialH option to restrict the computation of the inverse of the Hessian matrix to a submatrix

* Functions `hlme`, `lcmm`, `multlcmm`, `Jointlcmm` now allow optional vector B to be an estimated model (with G=1) to reduce calculation time of initial values.

* Bug identified and solved in calculation of subject-specific predictions in  `hlme`, `lcmm`, `multlcmm` and `Jointlcmm` when cor is not NULL.

* Bug identified and solved in the calculation of confidence bands for individual dynamic predictions in dynpred with draws=T.

* Bug identified and solved in the calculation of the explained variance for multlcmm objects when cor is not NULL.



# Changes in Version 1.7.1 & 1.7.2 (2015-02-27):


* Function plot now includes a which="fit" option to plot observed and predicted trajectories stemming from a hlme, lcmm, Jointlcmm or multlcmm object.

* Function `predictlink` replaces deprecated function `link.confint` 

* Function `plot` gathers deprecated functions `plot.linkfunction`, `plot.baselinerisk`, `plot.survival`, `plot.fit` together 



# Changes in Version 1.7.0 (2015-02-13):

* The function `Jointlcmm` now allows competing risks data for the survival part and is also available for non-Gaussian longitudinal data. All existing methods for Jointlcmm objects (except EPOCE and Diffepoce functions) are adapted to the new framework. 

* Functions `link.confint`, `plot.linkfunction`, `predictL` are now available for Jointlcmm objects.

* The new functions `incidcum` and `plot.incidcum` respectively compute and plot the cumulative incidence associated to each competing event for Jointlcmm object.

* The new function `fitY` computes the marginal predicted values of longitudinal outcomes in their natural scale for lcmm or multlcmm objects.

* Bug identified and solved in `dynpred` function when used with a joint model assuming proportional hazards between latent classes.

* The Makevars file now allows compilation of the package with parallel make.



# Changes in Version 1.6.5 & 1.6.6 (2014-09-10):

* bug solved regarding installation problem with parallel make



# Changes in Version 1.6.4 (2014-04-11):

* The new functions `dynpred` and `plot.dynpred` respectively compute and plot individual dynamic predictions obtained from a joint latent class model estimated by Jointlcmm.

* The new function `VarCovRE` computes the standard errors of the parameters of variance-covariance of the random effects for a hlme, lcmm, Jointlcmm or multlcmm object

* The new function `WaldMult` computes multivariate Wald tests and Wald tests for combinations of parameters from hlme, lcmm, Jointlcmm or multlcmm object

* The new function `VarExpl` computes the percentages of variance explained by the linear regression for a hlme, lcmm, Jointlclmm or multlcmm object

* The new functions `estimates` and `VarCov` get respectively all parameters estimated and their variance-covariance matrix for a hlme, lcmm, Jointlcmm or multlcmm object

* Function `summary` now returns the table containing the results about the fixed effects in the longitudinal model 

* All plots consider now the ... options

* Functions plot.linkfunction and plot.predict have now an add argument

* Function multlcmm now allows "splines" or "Splines" specification for the link functions

* Functions `lcmm` and `multlcmm` now compute the transformations even if the maximum number of iterations is reached without convergence

* bug identified and solved in multlcmm when the response variables are not integers

* bug identified and solved in multlcmm when using contrast

* bug identified and solved in plot.linkfunction for the y axes positions

* bug identified and solved in hlme, lcmm, Jointlcmm and multlcmm when including interactions in `mixture`.




# Changes in Version 1.6.2 (2013-03-06):

* The new function `multlcmm` now estimates latent process mixed models for multivariate curvilinear longitudinal outcomes (with link functions: linear, beta or splines). Various post-fit computation and output functions are also available including plot.linkfunction, predictY, predictL, etc

* All the functions hlme, lcmm, Jointlcmm include a `cor` option for including a brownian motion or a first-order autoregressive error process in addition to the independent errors of measurement

* bug identified and solved in predictL, predictY and plot.predict when used with factor covariate


# Changes in Version 1.5.8 (2012-10-01):

* bug identified and solved in predictY.lcmm when used with a `splines` link function and an outcome with minimum value not at 0


# Changes in Version 1.5.7 (2012-07-24):

* The function `predictY` now computes the predicted values (possibly class-specific) of the longitudinal outcome not only from a lcmm object but also from a hlme or a Jointlcmm object for a specified profile of covariates.

* bug identified and solved in predictY.lcmm when used with a `threshold` link function and a Monte Carlo method


# Changes in Version 1.5.6 (2012-07-16):

* missing data handled in hlme, lcmm and Jointlcmm using `na.action` with attributes 1 for `na.omit` or 2 for `na.fail`

* The new function `predictY.lcmm` computes predicted values of a lcmm object in the natural outcome scale for a specified profile of covariates, and also provides confidence bands using a Monte Carlo method.

* bugs in epoce computation solved (with splines baseline risk function, and/or NaN values under solaris system)

* bug identified and solved in summary functions regarding the labels of covariate effects in peculiar cases


# Changes in Version 1.5.2 (2012-04-06):

* improved variable specification in the estimating functions Jointlcmm, lcmm and hlme with 
	- categorical variables using factor() 
	- variables entered as functions using I() 
	- interaction terms using "*" and ":"

* computation of the predictive accuracy measure EPOCE from a Jointlcmm object either on the training data or on external data (post-fit functions epoce and Diffepoce)

* for discrete outcomes, lcmm function now computates the posterior discrete log-likelihood and the universal approximate cross-validation criterion (UACV)

* Jointlcmm now includes two parameterizations of I-splines and piecewise-constant baseline risks functions to ensure positive risks: either log/exp or sqrt/square (option logscale=).  



