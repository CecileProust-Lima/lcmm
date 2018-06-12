#' @name plot.predict
#' @rdname plot.predict
#' @export
plot.predictL <- function(x,legend.loc="topright",legend,add=FALSE,shades=FALSE,...)
    {
        if(missing(x)) stop("The argument \'x\' is missing.")
        if(!inherits(x,"predictL")) stop("use only with \'predictL\' object")


### determiner si  draws et ng

        colx <- colnames(x$pred)
 
        if(length(grep("class2",colx))) #ng>1
            {
                if(length(grep("lower",colx)))
                    {
                        ng <- length(grep("lower.pred_class",colx))
                        Ypred <- x$pred[,1:ng]
                        lower <- x$pred[,ng+1:ng]
                        upper <- x$pred[,ng+ng+1:ng]
                    }
                else
                    {
                        ng <- length(grep("pred_class",colx))
                        Ypred <- x$pred[,1:ng]
                        lower <- NULL
                        upper <- NULL
                    }
            }
        else
            {
                if(ncol(x$pred)>1) #draws
                    {
                        ng <- 1
                        Ypred <- x$pred[,1,drop=FALSE]
                        lower <- x$pred[,2,drop=FALSE]
                        upper <- x$pred[,3,drop=FALSE]
                    }
                else
                    {
                        ng <- 1
                        Ypred <- x$pred[,1,drop=FALSE]
                        lower <- NULL
                        upper <- NULL
                    }
            }
        


        dots <- list(...)
        dots <- dots[setdiff(names(dots),c("x","y","log"))]

        if(!length(dots$main))
            {
                dots$main <- "Class-specific mean predicted trajectory"
            }

        if(!length(dots$col))
            {
                dots$col <- 1:ng
            }

        if(!length(dots$type))
            {
                dots$type <- "l"
            }

        if(!length(dots$lty))
            {
                dots$lty <- 1:ng
            }

        if(!length(dots$ylab))
            {
               dots$ylab <- "latent process"
            }

        if(!length(dots$xlab))
            {
                dots$xlab <- colnames(x$times)
            }
        
        if(missing(legend)) legend <- paste("class",1:ng,sep="")

        if(!length(dots$box.lty))
            {
                dots$box.lty <- 0
            }

        if(!length(dots$inset))
            {
                dots$inset <- c(0.02,0.02)
            }

        if(!length(dots$ylim))
            {
                dots$ylim <- range(cbind(as.matrix(Ypred),lower,upper),na.rm=TRUE)
            }

        
        
        names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab",
                        "cex.main","cex.sub","col","col.axis","col.lab","col.main",
                        "col.sub","crt","err","family","fig","fin","font","font.axis",
                        "font.lab","font.main","font.sub","frame.plot","lab","las","lend",
                        "lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex",
                        "mgp","mkh","oma","omd","omi","pch","pin","plt","ps","pty","smo",
                        "srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
                        "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
        dots.plot <- dots[intersect(names(dots),names.plot)]

        
        if(add==FALSE)
            {
                do.call("matplot",c(dots.plot,list(x=x$times,y=Ypred)))

                if(!is.null(lower))
                    {
                        if(shades==FALSE)
                            {
                                do.call("matlines",c(dots.plot[setdiff(names(dots.plot),"lty")],list(x=x$times,y=cbind(lower,upper),lty=2)))
                            }
                        else
                            {
                                rgbcols <- sapply(dots$col,col2rgb)/255
                                cols <- apply(rgbcols,2,function(x) rgb(x[1],x[2],x[3],alpha=0.15))

                                sapply(1:ng,function(k,t,yl,yu,cols) polygon(x=unlist(c(t,rev(t))),y=c(yl[,k],rev(yu[,k])),col=cols[k],border=NA),t=unlist(x$times),yl=lower,yu=upper,cols=cols)
                            }
                    }
                
                if(!is.null(legend))
                    {
                        names.legend <- c("fill","border","lty","lwd","pch","angle",
                                          "density","bg","box.lwd","box.lty",
                                          "box.col","pt.bg","cex","pt.cex","pt.lwd",
                                          "xjust","yjust","x.intersp","y.intersp",
                                          "adj","text.width","text.col","text.font",
                                          "merge","trace","plot","ncol","horiz",
                                          "title","xpd","title.col","title.adj",
                                          "seg.len","col","inset")
                        dots.leg <- dots[intersect(names(dots),names.legend)]
                        do.call("legend",c(dots.leg,list(x=legend.loc, legend=legend)))
                    }
            }
        else
            {
                do.call("matpoints",c(dots.plot,list(x=x$times,y=Ypred)))
                
                if(!is.null(lower))
                    {
                        if(shades==FALSE)
                            {
                                do.call("matlines",c(dots.plot[setdiff(names(dots.plot),"lty")],list(x=x$times,y=cbind(lower,upper),lty=2)))
                            }
                        else
                            {
                                rgbcols <- sapply(dots$col,col2rgb)/255
                                cols <- apply(rgbcols,2,function(x) rgb(x[1],x[2],x[3],alpha=0.15))
                                 
                                sapply(1:ng,function(k,t,yl,yu,cols) polygon(x=unlist(c(t,rev(t))),y=c(yl[,k],rev(yu[,k])),col=cols[k],border=NA),t=unlist(x$times),yl=lower,yu=upper,cols=cols)
                            }
                    }            
            }

        return(invisible(NULL))
    }
