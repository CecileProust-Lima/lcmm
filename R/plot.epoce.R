#' Plots
#' 
#' This function displays plots related to predictive accuracy functions:
#' \code{epoce} and \code{Diffepoce}.
#' 
#' These functions do not apply for the moment with multiple causes of event
#' (competing risks).
#' 
#' For \code{epoce} objects, the function displays the EPOCE estimate (either
#' MPOL or CVPOL) according to the time of prediction.  For \code{Diffepoce}
#' objects, \code{plot} displays the difference in EPOCE estimates (either MPOL
#' or CVPOL) and its 95\% tracking interval between two joint latent class
#' models
#' 
#' @name plot.pred.accuracy
#' 
#' @param x an object inheriting from classes \code{epoce} or \code{Diffepoce}
#' @param \dots other parameters to be passed through to plotting functions
#' @return Returns plots related to \code{epoce} and \code{Diffepoce}
#' @author Cecile Proust-Lima and Viviane Philipps
#' @seealso \code{\link{epoce}},\code{\link{Diffepoce}}
#' @examples
#' 
#' \dontrun{
#' # estimation of the joint latent class model
#' m3 <- Jointlcmm(fixed= Ydep1~Time*X1,mixture=~Time,random=~Time,
#' classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~X1+mixture(X2),
#' hazard="3-quant-splines",hazardtype="PH",ng=3,data=data_lcmm,
#' B=c(0.7667, 0.4020, -0.8243, -0.2726, 0.0000, 0.0000, 0.0000, 0.3020,
#' -0.6212, 2.6247, 5.3139, -0.0255, 1.3595, 0.8172, -11.6867, 10.1668,
#' 10.2355, 11.5137, -2.6209, -0.4328, -0.6062, 1.4718, -0.0378, 0.8505,
#' 0.0366, 0.2634, 1.4981))
#' # predictive accuracy of the model evaluated with EPOCE
#' VecTime <- c(1,3,5,7,9,11,13,15)
#' cvpl <- epoce(m3,var.time="Time",pred.times=VecTime)
#' summary(cvpl)
#' plot(cvpl,bty="l",ylim=c(0,2))
#' }
#' 
#' @export
#' 
plot.epoce <- function(x,...)
{
	if (!inherits(x, "epoce")) stop("use only with \"epoce\" objects")

	if (all(is.na(x$EPOCE[,4]))|all(is.na(x$EPOCE[,5]))) stop("can't produce the plot with missing EPOCE")
        
        if (all(is.infinite(x$EPOCE[,4]))|all(is.infinite(x$EPOCE[,5]))) stop("can't produce the plot with infinite EPOCE")
        
   dots <- list(...)

        if(length(list(...)$add))
            {
                add <- eval(match.call()$add)
            }
        else
            add <- FALSE
        
   if(length(list(...)$main)) 
   {
    title1 <- as.character(eval(match.call()$main))
    dots <- dots[setdiff(names(dots),"main")]
   }
   else 
   {
    if(x$new.data==FALSE) title1 <- "Cross-validated prognostic observed log-likelihood"
    else title1 <- "Mean prognostic observed log-likelihood"
   }                    

   if(length(list(...)$type))    
   {
    type1 <- eval(match.call()$type)
    dots <- dots[-which(names(dots)=="type")]
   }
   else  type1 <- "o"

   if(length(list(...)$pch))    
   {
    pch1 <- eval(match.call()$pch)
    dots <- dots[-which(names(dots)=="pch")]
   }
   else  pch1 <- 18

 
   if(length(list(...)$xlab)) 
   {
    xlab1 <- as.character(eval(match.call()$xlab))
    dots <- dots[setdiff(names(dots),"xlab")]
   }
   else xlab1 <- "prediction time"

   if(length(list(...)$ylab)) 
   {
    ylab1 <- as.character(eval(match.call()$ylab))
    dots <- dots[setdiff(names(dots),"ylab")]
   }
   else 
   {
    if(x$new.data==FALSE) ylab1 <- "CVPOL"
    else ylab1 <- "MPOL"
   } 
   
   
   if(x$new.data==FALSE) y1 <- x$EPOCE[,5]
   else y1 <- x$EPOCE[,4]

   if(length(list(...)$ylim)) 
   {
       ylim1 <- eval(match.call()$ylim)
       dots <- dots[setdiff(names(dots),"ylim")]
   }
   else
       {
           if(all(is.na(y1)) | all(is.infinite(y1)))
               {
                   ylim1 <- c(-1000,1000)
               }
           else
               {
                   ylim1 <- c(min(y1[!(is.na(y1)) & is.finite(y1)]),max(y1[!(is.na(y1)) & is.finite(y1)]))
               }
       }        

   
   names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
   "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
   "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
   "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
   "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
   dots.plot <- dots[intersect(names(dots),names.plot)]   

        if(add==FALSE)
            { 
                do.call("plot",c(dots.plot,list(x=x$EPOCE[,1],y=y1,pch=pch1,type=type1,ylab=ylab1,xlab=xlab1,main=title1,ylim=ylim1)))
            }
        else
            {
                title1 <- ""
                xlab1 <- ""
                ylab1 <- ""

                do.call("points",c(dots.plot,list(x=x$EPOCE[,1],y=y1,pch=pch1,type=type1,ylab=ylab1,xlab=xlab1,main=title1,ylim=ylim1)))
            }
}

