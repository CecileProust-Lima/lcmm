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

