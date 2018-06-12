#' @rdname plot.predict
#' @export
plot.predictYcond <- function(x,legend.loc="topleft",legend,add=FALSE,shades=TRUE,...)
{
  if(missing(x)) stop("The argument x should be specified")
  if(is.na(as.logical(add))) stop("add should be TRUE or FALSE")  

  if(missing(legend)) legend <- unique(x$pred[,1])
  
  if(!all(is.na(x$pred)))
  {
      ny <- max(length(x$object$Ynames),1)
   
          dots <- list(...)
          plot.axes <- list(axes=TRUE,yaxt="s",xaxt="s")
          plot.axes[names(dots[c("axes","yaxt","xaxt")])] <- dots[c("axes","yaxt","xaxt")]
          if(plot.axes$axes==FALSE) plot.axes[c("yaxt","xaxt")] <- "n"
          
          dots <- dots[setdiff(names(dots),c("ylim","ylab","yaxt","x","y","log","xaxt","axes"))]

          if(!length(dots$main)) 
              {
                  dots$main <- "Conditional predicted values" #expression(E(Y~symbol("\174")~symbol("\114")))
              }

          if(!length(dots$col)) 
              {
                  dots$col <- rainbow(ny)
              }
          color <- rep(dots$col,length.out=ny)                            

          if(!length(dots$type))    
              {
                  dots$type <- "l"
              }

          if(!length(dots$lty))    
              {
                  dots$lty <- 1
              }
          
          if(!length(dots$xlab)) 
              {
                  dots$xlab <- "Latent process"
              }
          
          if(!length(dots$frame.plot)) 
              {
                  dots$frame.plot <- FALSE
              }
          
          if(!length(dots$box.lty)) 
              {
                  dots$box.lty <- 0
              }
          
          if(!length(dots$inset)) 
              {
                  dots$inset <- c(0.05,0.05)
              }

          if(!length(dots$cex.axis)) 
              {
                  dots$ces.axis <- 0.8 
              }   

          if(length(dots$mar))
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
                  nsim <- length(x$pred[,1])/ny
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

          linknodes <- matrix(x$object$linknodes,ncol=ny)
          linktype <- x$object$linktype
          
          loc.grad <- function(y.grad,yk)
              {
                  (nsim*(y.grad-min(linknodes[1:nbnodes[yk],yk]))-y.grad+max(linknodes[1:nbnodes[yk],yk]))/(max(linknodes[1:nbnodes[yk],yk])-min(linknodes[1:nbnodes[yk],yk]))
              }

      nbnodes <- rep(2,ny)
      if(ny>1)
      {
          nbnodes[which(linktype==2)] <- x$object$nbnodes
      }
      else
      {
           nbnodes[which(linktype==2)] <- length(linknodes)
      }
          
          oldmar <- par("mar")
          on.exit(par(mar=oldmar))
          par(mar=mar1)
          
          names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
                          "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
                          "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
                          "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
                          "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
          dots.plot <- dots[intersect(names(dots),names.plot)]

          lambda <- matrix(x$pred[,2],ncol=ny)
          ypred <- matrix(x$pred[,3],ncol=ny)
          
          loc.y <- matrix(NA,nrow=nrow(x$pred)/ny,ncol=ny)
          for(yk in 1:ny)
              {
                  loc.y[,yk] <- loc.grad(ypred[,yk],yk) 
              }


          if(!isTRUE(add))
              {
                  do.call("matplot",c(dots.plot,list(x=lambda,y=loc.y,axes=FALSE,ylim=c(1,nsim),ylab="")))

                  names.axis <- c("lwd","lwd.ticks","hadj","padj","cex.axis","font.axis",
                                  "xaxp","yaxp","tck","tcl","las","xpd")
                  dots.axis <- dots[intersect(names(dots),names.axis)]                                                        
                  
                  if(plot.axes$xaxt=="s") do.call("axis",c(dots.axis,list(side=1,col=1)))
                  
                  y.grad <- pretty(min(linknodes[1:nbnodes[1],1]):max(linknodes[1:nbnodes[1],1]))
                  y.grad[1] <- round(min(linknodes[1:nbnodes[1],1]),2)
                  y.grad[length(y.grad)] <- round(max(linknodes[1:nbnodes[1],1]),2)
                  if(plot.axes$yaxt=="s") do.call("axis",c(dots.axis,list(side=2,at=loc.grad(y.grad,1),labels=y.grad,col=color[1],col.axis=color[1])))

                  if(ny>1)
                      {
                          for (i in 2:ny)
                              {
                                  y.grad <- pretty(min(linknodes[1:nbnodes[i],i]):max(linknodes[1:nbnodes[i],i]))
                                  y.grad[1] <- round(min(linknodes[1:nbnodes[i],i]),2)
                                  y.grad[length(y.grad)] <- round(max(linknodes[1:nbnodes[i],i]),2)
                                  if(plot.axes$yaxt=="s") do.call("axis",c(dots.axis,list(side=ifelse(i%%2==0,4,2),at=loc.grad(y.grad,i),labels=y.grad,col=color[i],col.axis=color[i],line=(round((i+0.1)/2)-1)*2)))
                              }
                      }
              }
          else
              {
                  do.call("matlines",c(dots.plot,list(x=lambda,y=loc.y,axes=FALSE)))
              }

          if(ncol(x$pred)==5) ## avec IC
              {
                  if(!is.na(shades))
                      {
                          if(shades==TRUE)
                              {
                                  rgbcols <- sapply(color,col2rgb)/255
                                  cols <- apply(rgbcols,2,function(x) rgb(x[1],x[2],x[3],alpha=0.15))
                                  lty <- NULL
                                  lwd <- NULL
                                  border <- NA
                                  density <- NULL
                              }
                          else
                              {
                                  cols <- color
                                  lty <- 2
                                  lwd <- dots$lwd
                                  border <- dots$col
                                  density <- 0
                              }
                          
                          binf <- matrix(x$pred[,4],ncol=ny)
                          bsup <- matrix(x$pred[,5],ncol=ny)
                          for(yk in 1:ny)
                              {
                                  binf[,yk] <- loc.grad(binf[,yk],yk) 
                                  bsup[,yk] <- loc.grad(bsup[,yk],yk) 

                                  jna <- unique(c(which(is.na(lambda[,yk])),
                                                  which(is.na(binf[,yk])),
                                                  which(is.na(bsup[,yk]))))
                                  if(length(jna))
                                      {
                                          polygon(x=c(lambda[-jna,yk],rev(lambda[-jna,yk])),y=c(binf[-jna,yk],rev(bsup[-jna,yk])),col=cols[yk],border=border,density=density,lty=lty,lwd=lwd)
                                      }
                                  else
                                      {
                                          polygon(x=c(lambda[,yk],rev(lambda[,yk])),y=c(binf[,yk],rev(bsup[,yk])),col=cols[yk],border=border,density=density,lty=lty,lwd=lwd)
                                      }
                              }
                      }
              }
          
          names.legend <- c("fill","border","lty","lwd","pch","angle","density","bg","box.lwd",   
                            "box.lty","box.col","pt.bg","cex","pt.cex","pt.lwd","xjust","yjust","x.intersp","y.intersp","adj","text.width",
                            "text.col","text.font","merge","trace","plot","ncol","horiz","title","xpd","title.col","title.adj","seg.len")    
          dots.leg <- dots[intersect(names(dots),names.legend)]
          
          if(!is.null(legend)) do.call("legend",c(dots.leg,list(x=legend.loc,legend=legend,col=color)))
      }
  else
      {
          cat("Output can not be produced since the program stopped abnormally. \n")
      }
}

