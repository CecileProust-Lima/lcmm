#' Summary of models
#' 
#' This function provides a plot summarizing the results of different models
#' fitted by \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}.
#' 
#' Can be reported the usual criteria used to assess the fit and the clustering
#'  of the data:
#'  - maximum log-likelihood L (the higher the better)
#'  - number of parameters P, number of classes G, convergence criterion (1 = converged)
#'  - AIC (the lower the better) computed as -2L+2P 
#'  - BIC (the lower the better) computed as -2L+ P log(N) where N is the number of subjects
#'  - SABIC (the lower the better) computed as 
#'  - Entropy (the closer to one the better) computed as 1-sum[pi_ig*log(pi_ig)]/(N*log(G))
#'    where pi_ig is the posterior probability that subject i belongs to class g
#'  - ICL (the lower the better) computed as BIC -2*sum[c_ig*log(pi_ig)]
#'    where c_ig is the posterior class membership
#'  - %Class computed as the proportion of each class based on c_ig
#' 
#' @param m1 an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or
#' \code{Jointlcmm}
#' @param \dots  further arguments, in particular other objects of class
#' \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}
#' @param which character vector indicating which results should be returned.
#' Possible values are "loglik", "conv", "npm", "AIC", "BIC", "SABIC",
#' "entropy", "ICL", "\%class".
#' @param xaxis the abscissa of the plot, default to "G
#' @param display display of the table. By default, False.
#' @author Cecile Proust-Lima, Viviane Philipps, Sasha Cuau
#' @seealso \code{\link{summary}}, \code{\link{summarytable}} 
#' 
#' @export
#' @examples   
# paquid$age65 <- (paquid$age - 65)/10
# m1 <- hlme(MMSE ~ age65+I(age65^2)+CEP,random =~ age65+I(age65^2), subject = 'ID', data = paquid) # ng=1
# m2 <- hlme(MMSE ~ age65+I(age65^2)+CEP, random =~ age65+I(age65^2), subject = 'ID', data = paquid, ng = 2, mixture=~age65+I(age65^2), B=m1)
# m3g <- gridsearch(hlme(MMSE ~ age65+I(age65^2)+CEP,  random =~ age65+I(age65^2), subject = 'ID', data=paquid, ng = 3, mixture=~age65+I(age65^2)), rep=100, maxiter=30, minit=m1)
# summaryplot(m1,m2,m3g,which=c("G","BIC","entropy","ICL","loglik"))




summaryplot <- function(m1,...,which,width=length(which),height=1,xaxis="G",display=FALSE){
  dots <- list(...)
  names.plot <- c("lty","type","col","pch","xlab","lwd","ylab")
  dots.plot <- dots[intersect(names(dots),names.plot)]
  models <- setdiff(dots,dots.plot)
  
  if(!length(dots.plot$xlab))
  {
    dots.plot$xlab <- "#classes"
  }
  
  if(!length(dots.plot$ylab))
  {
    dots.plot$ylab <- ""
  }
  if(!length(dots.plot$col))
  {
    dots.plot$col <- "orchid4"
  }
  if(!length(dots.plot$type))
  {
    dots.plot$type<- "o"
  }
  if(!length(dots.plot$pch))
  {
    dots.plot$pch<- 10
  }
  if(!length(dots.plot$lty))
  {
    dots.plot$lty<- 1
  }
  if(!length(dots.plot$lwd))
  {
    dots.plot$lwd<- 3
  }
  if(!length(dots.plot$bty))
  {
    dots.plot$bty<- "l"
  }

  summ <- summarytable(m1,...,which=c(xaxis,which),display= FALSE)
  par(mfrow=c(height,width))
  for(x in which){
    if(x!=xaxis){
      do.call(plot,c(summ[,c(x)] ~ summ[,c(xaxis)],dots.plot,main=x,xaxt="n"))
      axis(1, at = 1:(length(models)+1))
      }
  }
  return(invisible(NULL))
}








