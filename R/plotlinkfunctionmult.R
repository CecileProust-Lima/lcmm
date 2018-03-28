
.plotlinkfunctionmult <- function(x,legend.loc="topleft",legend=x$Ynames,add=FALSE,...)
{
  if(missing(x)) stop("The argument x should be specified")
  if (!inherits(x, "multlcmm")) stop("use only with \"multlcmm\" objects")
  if(is.na(as.logical(add))) stop("add should be TRUE or FALSE")  
  
  if(x$conv %in% c(1,2,3))
  {
   ny <- length(x$Ynames) 
   
   dots <- list(...)
   plot.axes <- list(axes=TRUE,yaxt="s",xaxt="s")
   plot.axes[names(dots[c("axes","yaxt","xaxt")])] <- dots[c("axes","yaxt","xaxt")]
   if(plot.axes$axes==FALSE) plot.axes[c("yaxt","xaxt")] <- "n"
  
   dots <- dots[setdiff(names(dots),c("ylim","ylab","yaxt","x","y","log","xaxt","axes"))]

   if(length(list(...)$main)) 
   {
    title1 <- as.character(list(...)$main)
    dots <- dots[setdiff(names(dots),"main")]
   }
   else title1 <- "Estimated link functions"

   if(length(list(...)$col)) 
   {
    color <- as.vector(list(...)$col)
    dots <- dots[-which(names(dots)=="col")]
   }
   else  color <- rainbow(ny)
   color <- rep(color,length.out=ny)                            

   if(length(list(...)$type))    
   {
    type1 <- list(...)$type
    dots <- dots[-which(names(dots)=="type")]
   }
   else  type1 <- "l" 

   if(length(list(...)$lty))    
   {
    lty1 <- dots$lty
    dots <- dots[-which(names(dots)=="lty")]
   }
   else  lty1 <- 1 
   
   if(length(list(...)$xlab)) 
   {
    xlab1 <- as.character(list(...)$xlab)
    dots <- dots[setdiff(names(dots),"xlab")]
   }
   else xlab1 <- "Latent process"
   
   if(length(list(...)$frame.plot)) 
   {
    frame.plot1 <- list(...)$frame.plot
    dots <- dots[setdiff(names(dots),"frame.plot")]
   }
   else frame.plot1 <- FALSE   
   
   if(length(list(...)$box.lty)) 
   {
    box.lty1 <- as.integer(list(...)$box.lty)
    dots <- dots[setdiff(names(dots),"box.lty")]
   }
   else box.lty1 <- 0
   
   if(length(list(...)$inset)) 
   {
    inset1 <- list(...)$inset
    dots <- dots[setdiff(names(dots),"inset")]
   }
   else inset1 <- c(0.05,0.05)

   if(length(list(...)$cex.axis)) 
   {
    cex.axis1 <- list(...)$cex.axis
    dots <- dots[setdiff(names(dots),"cex.axis")]
   }
   else cex.axis1 <- 0.8   

   if(length(list(...)$mar))
   {
    mar1 <- list(...)$mar
    dots <- dots[setdiff(names(dots),"mar")]
   }
   else
       {
           if(plot.axes$yaxt!="n" )#| isTRUE(add))
               mar1 <- c(5,ny+1,2,ny+1)+0.2 
           else
               mar1 <- c(5,4,2,4)+0.2
       }
   
   if(!isTRUE(add))
   {
    nsim <- length(x$estimlink[,1])
   }
   else
   {
    if(par("yaxs")=="r")
    {
     a <- (26*par("usr")[3])/27 + par("usr")[4]/27
     b <- (par("usr")[3]+26*par("usr")[4])/27
     
     nsim <- b-a+1
    }
    
    if(par("yaxs")=="i")
    {
     nsim <- par("usr")[4]-par("usr")[3]+1
    }
   }
 
   loc.grad <- function(y.grad,yk)
   {
    (nsim*(y.grad-min(x$linknodes[1:nbnodes[yk],yk]))-y.grad+max(x$linknodes[1:nbnodes[yk],yk]))/(max(x$linknodes[1:nbnodes[yk],yk])-min(x$linknodes[1:nbnodes[yk],yk]))
   }

   nbnodes <- rep(2,ny)
   nbnodes[which(x$linktype==2)] <- x$nbnodes
   
    
    oldmar <- par("mar")
    on.exit(par(mar=oldmar))
    par(mar=mar1)

#    list.arg <- list(x=x$estimlink[,2],y=1:nsim,type=type1,pch=pch1[1],bg=bg1[1],lty=lty1[1],lwd=lwd1[1],xlim=xlim1,ylim=c(1,nsim),xlab=xlab1,ylab="",main=title1,xaxt="n",yaxt="n",col=color[1],frame.plot=frame.plot1,mgp=mgp1)
#    do.call("plot",list.arg)
    
    names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
    "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
    "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
    "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
    "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
    dots.plot <- dots[intersect(names(dots),names.plot)]

    loc.y <- x$estimlink[,2*(1:ny)-1,drop=FALSE]
    for(yk in 1:ny)
    {
     loc.y[,yk] <- loc.grad(x$estimlink[,2*yk-1],yk) 
    }


    if(!isTRUE(add))
    {
     do.call("matplot",c(dots.plot,list(x=x$estimlink[,2*(1:ny)],y=loc.y,type=type1,col=color,axes=FALSE,ylim=c(1,nsim),xlab=xlab1,ylab="",main=title1,lty=lty1)))

    names.axis <- c("lwd","lwd.ticks","hadj","padj","cex.axis","font.axis",
    "xaxp","yaxp","tck","tcl","las","xpd","cex.axis")
    dots.axis <- dots[intersect(names(dots),names.axis)]                                                        
    
    if(plot.axes$xaxt=="s") do.call("axis",c(dots.axis,list(side=1,col=1,cex.axis=cex.axis1))) 
    #axis(1,cex.axis=0.8)
    y.grad <- pretty(min(x$linknodes[1:nbnodes[1],1]):max(x$linknodes[1:nbnodes[1],1]))
    y.grad[1] <- round(min(x$linknodes[1:nbnodes[1],1]),2)
    y.grad[length(y.grad)] <- round(max(x$linknodes[1:nbnodes[1],1]),2)
    if(plot.axes$yaxt=="s") do.call("axis",c(dots.axis,list(side=2,at=loc.grad(y.grad,1),labels=y.grad,col=color[1],col.axis=color[1],cex.axis=cex.axis1)))

    #axis(2,at=loc.grad(y.grad,1),labels=y.grad,col=color[1], col.axis=color[1], cex.axis=cex.axis1)
    if(ny>1)
    {
     for (i in 2:ny)
     {
      y.grad <- pretty(min(x$linknodes[1:nbnodes[i],i]):max(x$linknodes[1:nbnodes[i],i]))
      y.grad[1] <- round(min(x$linknodes[1:nbnodes[i],i]),2)
      y.grad[length(y.grad)] <- round(max(x$linknodes[1:nbnodes[i],i]),2)
      #list.arg <- list(x=x$estimlink[,2*i],y=1:nsim,col=color[i],type=type1,pch=pch1[i],lty=lty1[i],lwd=lwd1[i],bg=bg1[i])
      #do.call("lines",list.arg)
      if(plot.axes$yaxt=="s") do.call("axis",c(dots.axis,list(side=ifelse(i%%2==0,4,2),at=loc.grad(y.grad,i),labels=y.grad,col=color[i],col.axis=color[i],cex.axis=cex.axis1,line=(round((i+0.1)/2)-1)*2)))
      #axis(ifelse(i%%2==0,4,2),at=loc.grad(y.grad,i),labels=y.grad,col=color[i], col.axis=color[i], cex.axis=0.8,line=(round((i+0.1)/2)-1)*2)
     }
    }
   }
   else
   {
     do.call("matlines",c(dots.plot,list(x=x$estimlink[,2*(1:ny)],y=loc.y,type=type1,col=color,axes=FALSE,lty=lty1)))   
   }
   
      
   names.legend <- c("fill","border","lty","lwd","pch","angle","density","bg","box.lwd",   
   "box.lty","box.col","pt.bg","cex","pt.cex","pt.lwd","xjust","yjust","x.intersp","y.intersp","adj","text.width",
   "text.col","text.font","merge","trace","plot","ncol","horiz","title","xpd","title.col","title.adj","seg.len")    
   dots.leg <- dots[intersect(names(dots),names.legend)]
    
   if(!is.null(legend)) do.call("legend",c(dots.leg,list(x=legend.loc,legend=legend,col=color,box.lty=box.lty1,inset=inset1,lty=lty1)))
  }
  else
  {
   cat("Output can not be produced since the program stopped abnormally. \n")
  }
}
