

.plotbaselinerisk <- function(x,event=1,legend.loc="topleft",legend,add=FALSE,...)
{
  if(missing(x)) stop("The argument x should be specified")

  nbevt <- length(x$N)-9
  if(length(event)!=1) stop("Please specify only one event")
  if(!(event %in% c(1:nbevt))) stop("Argument event is not correct")
  if(is.na(as.logical(add))) stop("add should be TRUE or FALSE")
  
  if((x$conv==1|x$conv==2|x$conv==3)& (sum(is.na(x$predSurv)==0)))
  {
   ng <- x$ng
  
   dots <- list(...)

   if(length(list(...)$main)) 
   {
    title1 <- as.character(eval(match.call()$main))
    dots <- dots[setdiff(names(dots),"main")]
   }
   else title1 <- "Class-specific baseline risk functions"

   if(length(list(...)$col)) 
   {
    color <- as.vector(eval(match.call()$col))
    dots <- dots[-which(names(dots)=="col")]
   }
   else  color <- 1:ng                            
 
   if(length(list(...)$type))    
   {
    type1 <- eval(match.call()$type)
    dots <- dots[-which(names(dots)=="type")]
   }
   else  type1 <- "l"
      
   if(length(list(...)$xlab)) 
   {
    xlab1 <- as.character(eval(match.call()$xlab))
    dots <- dots[setdiff(names(dots),"xlab")]
   }
   else xlab1 <- "Time"
   
   if(length(list(...)$ylab)) 
   {
    ylab1 <- as.character(eval(match.call()$ylab))
    dots <- dots[setdiff(names(dots),"ylab")]
   }
   else ylab1 <- "Baseline risk function" 
   
   if(length(list(...)$lty))    
   {
    lty1 <- eval(match.call()$lty)
    dots <- dots[-which(names(dots)=="lty")]
   }
   else  lty1 <- 1:ng  
   
   if(missing(legend)) legend <- paste("class",1:ng,sep="")
   
   if(length(list(...)$box.lty)) 
   {
    box.lty1 <- as.integer(eval(match.call()$box.lty))
    dots <- dots[setdiff(names(dots),"box.lty")]
   }
   else box.lty1 <- 0
   
   if(length(list(...)$inset)) 
   {
    inset1 <- eval(match.call()$inset)
    dots <- dots[setdiff(names(dots),"inset")]
   }
   else inset1 <- c(0.05,0.05)
   

   
  names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
  "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
  "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
  "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
  "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
  dots.plot <- dots[intersect(names(dots),names.plot)]
  
  if(!isTRUE(add))
  {
   do.call("matplot",c(dots.plot,list(x=x$predSurv[,1],y=x$predSurv[,(1+(event-1)*ng+1:ng)],xlab=xlab1,ylab=ylab1,main=title1,type=type1,col=color)))
  }
  else
  {
   do.call("matlines",c(dots.plot,list(x=x$predSurv[,1],y=x$predSurv[,(1+(event-1)*ng+1:ng)],type=type1,col=color)))
  }
  
  names.legend <- c("fill","border","lty","lwd","pch","angle","density","bg","box.lwd",   
  "box.lty","box.col","pt.bg","cex","pt.cex","pt.lwd","xjust","yjust","x.intersp","y.intersp","adj","text.width",
  "text.col","text.font","merge","trace","plot","ncol","horiz","title","xpd","title.col","title.adj","seg.len")    
  dots.leg <- dots[intersect(names(dots),names.legend)]
  
                                    
   if(!is.null(legend))
    {
     do.call("legend",c(dots.leg,list(x=legend.loc,legend=legend,col=color,box.lty=box.lty1,inset=inset1,lty=lty1)))
    }
  }
  else
  {
   cat("Output can not be produced. The program stopped abnormally or there was an error in the computation of the estimated baseline risk functions and survival functions.\n")
  }
}
