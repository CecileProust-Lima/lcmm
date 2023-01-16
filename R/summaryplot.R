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
#'  - SABIC (the lower the better) computed as -2L+ P log((N+2)/24)
#'  - Entropy (the closer to one the better) computed as 1-sum[pi_ig*log(pi_ig)]/(N*log(G))
#'    where pi_ig is the posterior probability that subject i belongs to class g
#'  - ICL (the lower the better) computed in two ways : ICL1 = BIC - sum[pi_ig*log(pi_ig)]
#'    or ICL2 = BIC - 2*sum(log(max(pi_ig)), where the max is taken over the classes for each subject.
#'  - %Class computed as the proportion of each class based on c_ig
#' 
#' @param m1 an object of class \code{hlme}, \code{lcmm}, \code{multlcmm},
#' \code{Jointlcmm} or \code{mpjlcmm}
#' @param \dots  further arguments, in particular other objects of class
#' \code{hlme}, \code{lcmm}, \code{multlcmm}, \code{Jointlcmm} or \code{mpjlcmm}, and
#' graphical parameters.
#' @param which character vector indicating which results should be plotted.
#' Possible values are "loglik", "conv", "npm", "AIC", "BIC", "SABIC",
#' "entropy", "ICL", "ICL1", "ICL2".
#' @param mfrow for multiple plots, number of rows and columns to split the graphical device.
#' Default to one line and length(which) columns. 
#' @param xaxis the abscissa of the plot. Default to "G", the number of latent classes.
#' @author Sasha Cuau, Viviane Philipps, Cecile Proust-Lima
#' @seealso \code{\link{summary}}, \code{\link{summarytable}} 
#' 
#' @export
#' @examples
#' \dontrun{
#' library(NormPsy)
#' paquid$normMMSE <- normMMSE(paquid$MMSE)
#' paquid$age65 <- (paquid$age - 65)/10
#' m1 <- hlme(normMMSE~age65+I(age65^2)+CEP, random=~age65+I(age65^2), subject='ID', data=paquid)
#' m2 <- hlme(normMMSE~age65+I(age65^2)+CEP, random=~age65+I(age65^2), subject='ID', data=paquid,
#' ng = 2, mixture=~age65+I(age65^2), B=m1)
#' m3g <- gridsearch(hlme(normMMSE~age65+I(age65^2)+CEP, random=~age65+I(age65^2), subject='ID',
#' data=paquid, ng=3, mixture=~age65+I(age65^2)), rep=100, maxiter=30, minit=m1)
#' summaryplot(m1, m2, m3g, which=c("BIC","entropy","ICL"),bty="l",pch=20,col=2)
#'}



summaryplot <- function(m1, ..., which=c("BIC", "entropy", "ICL"), mfrow=c(1,length(which)), xaxis="G")
{
    if(!all(which %in% c("loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "ICL", "ICL1", "ICL2"))) stop("Argument which should only contain elements among loglik, conv, npm, AIC, BIC, SABIC, entropy, ICL, ICL1, ICL2")
    
    dots <- list(...)
    names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab",
                    "cex.main","cex.sub","col","col.axis","col.lab","col.main",
                    "col.sub","crt","err","family","fig","fin","font","font.axis",
                    "font.lab","font.main","font.sub","frame.plot","lab","las","lend",
                    "lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex",
                    "mgp","mkh","oma","omd","omi","pch","pin","plt","ps","pty","smo",
                    "srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
                    "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
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
        dots.plot$col <- "black"
    }
    if(!length(dots.plot$type))
    {
        dots.plot$type <- "o"
    }
    if(!length(dots.plot$lty))
    {
        dots.plot$lty <- 1
    }

    if("ICL" %in% which) which <- unique(c(setdiff(which, "ICL"), "ICL1", "ICL2"))
    main <- which
    if(length(dots.plot$main))
    {
        main <- rep(dots.plot$main, length.out=length(which))
        dots.plot <- dots.plot[setdiff(names(dots.plot),"main")]
    }
    
    summ <- summarytable(m1, ..., which=c(xaxis,which), display=FALSE)

    par(mfrow=mfrow)
    on.exit(par(mfrow=c(1,1)))

    k <- 1
    for(x in which){
        do.call("plot",c(summ[,c(x)] ~ summ[,c(xaxis)], dots.plot, main=main[k], xaxt="n"))
        axis(1, at = 1:(length(models)+1))
        k <- k+1
    }
    return(invisible(NULL))
}








