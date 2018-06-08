#' Simulated dataset for hlme function
#' 
#' The data were simulated from a 3-latent class linear mixed model. Repeated
#' data for 100 subjects were simulated. The three latent classes are predicted
#' by X2 and X3. In each latent class, Y follows a linear mixed model including
#' intercept and time both with correlated random-effects and class-specific
#' fixed effects. In addition, X1 and X1*time have a common impact over classes
#' on the Y trajectory.
#' 
#' 
#' @name data_hlme
#' @docType data
#' @format A data frame with 326 observations on the following 9 variables.
#' \describe{ \item{ID}{subject identification number}
#' \item{Y}{longitudinal outcome} \item{Time}{time of
#' measurement} \item{X1}{binary covariate} \item{X2}{binary
#' covariate} \item{X3}{binary covariate} }
#' @seealso \code{\link{hlme}}, \code{\link{postprob}},
#' \code{\link{summary.lcmm}}, \code{\link{plot.predict}}
#' @keywords datasets
NULL
