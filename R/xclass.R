#' Cross classifications
#'
#' This function crosses the posterior classifications of two estimated models
#'
#' @param m1 an object inheriting from classes \code{hlme}, \code{lcmm}, \code{multlcmm}, \code{Jointlcmm} or \code{mpjlcmm}
#' @param m2 an object inheriting from classes \code{hlme}, \code{lcmm}, \code{multlcmm}, \code{Jointlcmm} or \code{mpjlcmm}
#' @return the contingency table of the two classifications
#' @author Viviane Philipps and Cecile Proust-Lima
#'
#' @examples
#'
#' ## Estimation of the models
#' m2 <- hlme(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,subject='ID',ng=2,
#' data=data_hlme,B=c(0.11,-0.74,-0.07,20.71,29.39,-1,0.13,2.45,-0.29,4.5,0.36,0.79,0.97))
#' m3 <- hlme(fixed = Y ~ Time * X1, mixture = ~Time, random = ~Time,subject = "ID",
#' classmb = ~X2 + X3, ng = 3, data = data_hlme,B=c(-0.21, 0.31, -2.11, -0.81, -0.24,
#' -0.18, 25.4, 20.09, 30.18, -0.43, -1.1, 0.25, 2.37, -0.29, 2.34, 0.03, 0.74, 0.97))
#'
#' ## Compare the classifications
#' xclass(m2,m3)
#' # The 39 subjects in class 2 of m3 come from class 1 of m2.
#' # In the same way, all the subjects in class 3 come from class 2 of m2.
#' # Class 1 of m3 mixes subject from class 1 and class 2 of m2.
#' 
#' @export
#'

xclass <- function(m1,m2)
{
    if(!(class(m1) %in% c("hlme","lcmm","multlcmm","Jointlcmm","mpjlcmm"))) stop("Please use this function only with hlme, lcmm, multlcmm, Jointlcmm or mpjlcmm models")
    if(!(class(m2) %in% c("hlme","lcmm","multlcmm","Jointlcmm","mpjlcmm"))) stop("Please use this function only with hlme, lcmm, multlcmm, Jointlcmm or mpjlcmm models")

    table(m1$pprob[,2],m2$pprob[,2])
}
