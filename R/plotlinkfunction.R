.plotlinkfunction <- function(x,legend.loc="topright",legend=NULL,add=FALSE,...)
{
 if(missing(x)) stop("The argument x should be specified")
 if(x$linktype==-1) stop("The model does not define any link function")
 if(is.na(as.logical(add))) stop("add should be TRUE or FALSE")
   
 if(x$conv %in% c(1,2,3))
 {
   dots <- list(...)

   if(length(list(...)$main)) 
   {
    title1 <- as.character(eval(match.call()$main))
    dots <- dots[setdiff(names(dots),"main")]
   }
   else title1 <- "Estimated link function"                          

   if(length(list(...)$type))    
   {
    type1 <- eval(match.call()$type)
    dots <- dots[-which(names(dots)=="type")]
   }
   else  type1 <- "l" 
   
   if(length(list(...)$xlim)) 
   {
    #limx fait apres
    dots <- dots[setdiff(names(dots),"xlim")]
   }        
   if(length(list(...)$xlab)) 
   {
    xlab1 <- as.character(eval(match.call()$xlab))
    dots <- dots[setdiff(names(dots),"xlab")]
   }
   else xlab1 <- "Latent process"
   
   if(length(list(...)$ylab)) 
   {
    ylab1 <- as.character(eval(match.call()$ylab))
    dots <- dots[setdiff(names(dots),"ylab")]
   }
   else ylab1 <- "Longitudinal outcome"
  
  if (x$linktype==3 & (x$linknodes[2]-x$linknodes[1])>1)
  {
   ntrtot <- sum(x$ide==1)
   diff <- x$estimlink[(2*(x$linknodes[2]-x$linknodes[1]+1)-1),2]-x$estimlink[2,2]
   diff <- diff/ntrtot

   xlim1 <- as.vector(c(x$estimlink[2,2]-diff,x$estimlink[(2*(x$linknodes[2]-x$linknodes[1]+1)-1),2]+diff))
  } 
  else
  {
   xlim1 <- c(min(x$estimlink[,2]),max(x$estimlink[,2]))
  }

  if("xlim" %in% names(match.call())) xlim1 <- eval(match.call()$xlim) 

  
  names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
  "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
  "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
  "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
  "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
  dots.plot <- dots[intersect(names(dots),names.plot)]
 
  if(!isTRUE(add))
  {
   do.call("plot",c(dots.plot,list(x=x$estimlink[,2],y=x$estimlink[,1],xlab=xlab1,ylab=ylab1,main=title1,type=type1,xlim=xlim1)))
  }
  else
  {
   do.call("lines",c(dots.plot,list(x=x$estimlink[,2],y=x$estimlink[,1])))
  }
 }
 else
 {
  cat("Output can not be produced since the program stopped abnormally. \n")
 }
}
