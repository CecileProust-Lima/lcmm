
plot.predictlink <- function(x,legend.loc="topleft",legend,add=FALSE,shades=FALSE,...)
    {
        if(missing(x)) stop("Argument 'x' is missing")
        if(!inherits(x,"predictlink")) stop("Use only with 'predictlink' object")

        colx <- colnames(x$pred)
        if(colx[1]=="Yname") #multlcmm
            {
                ny <- length(unique(x$pred[,1]))
                nsim <- length(x$pred[,1])/ny
                
                copiex <- x$object
                ysim <- matrix(x$pred[,2],nsim,ny)
                transfo <- matrix(x$pred[,3],nsim,ny)
                estimlink <- as.vector(rbind(ysim,transfo))
                copiex$estimlink <- matrix(estimlink,nsim,2*ny)

                dots <- list(...)
                
                if(length(list(...)$xlim)==0)
                    {
                        if(all(x$object$linktype!=3))
                            {    
                                dots <- c(dots,list(xlim=range(x$pred[,3:ncol(x$pred)],na.rm=TRUE)))
                            }
                        else
                            {
                                big <- c(nsim*(which(x$object$linktype==3)-1)+1,nsim*which(x$object$linktype==3))
                                dots <- c(dots,list(xlim=range(x$pred[-big,3:ncol(x$pred)],na.rm=TRUE)))
                            }
                    }

                dots$add <- add
                if(!missing(legend)) dots$legend <- legend
                if(!missing(legend.loc)) dots$legend.loc <- legend.loc

                    
                ##tracer la transfo
                do.call("plot.multlcmm",c(dots,list(x=copiex,which="linkfunction")))
                    

                if(length(list(...)$lwd)==0)
                    {
                        dots <- c(dots,list(lwd=1))
                    }

                if(length(list(...)$col)==0)
                    {
                        dots <- c(dots,list(col=rainbow(ny)))
                    }
                    
                    
                if(ncol(x$pred)>3) #draws
                    {
                        if(shades==FALSE)
                            {
                                ##on ne garde que les options lwd et col et on impose lty=2
                                dots.bornes <- dots[intersect(names(dots),c("lwd","col"))]
                                dots.bornes <- c(dots.bornes,lty=2) 
                                
                                copiex <- x$object
                                ysim <- matrix(x$pred[,2],nsim,ny)
                                transfo <- matrix(x$pred[,4],nsim,ny)
                                estimlink <- as.vector(rbind(ysim,transfo))
                                copiex$estimlink <- matrix(estimlink,nsim,2*ny)
                                
                                ##tracer la borne inf
                                do.call("plot.multlcmm",c(dots.bornes,list(x=copiex,which="linkfunction",add=TRUE,legend=NULL)))
                                
                                
                                copiex <- x$object
                                ysim <- matrix(x$pred[,2],nsim,ny)
                                transfo <- matrix(x$pred[,5],nsim,ny)
                                estimlink <- as.vector(rbind(ysim,transfo))
                                copiex$estimlink <- matrix(estimlink,nsim,2*ny)
                                
                                ##tracer la borne sup
                                do.call("plot.multlcmm",c(dots.bornes,list(x=copiex,which="linkfunction",add=TRUE,legend=NULL)))
                            }
                        else
                            {
                                ## shades ##

                                if(!isTRUE(add))
                                    {
                                        np <- nsim
                                    }
                                else
                                    {
                                        if(par("yaxs")=="r")
                                            {
                                                a <- (26*par("usr")[3])/27 + par("usr")[4]/27
                                                b <- (par("usr")[3]+26*par("usr")[4])/27
                                                
                                                np <- b-a+1
                                            }
                                        
                                        if(par("yaxs")=="i")
                                            {
                                                np <- par("usr")[4]-par("usr")[3]+1
                                            }
                                    }
                                
                                
                                rgbcols <- sapply(dots$col,col2rgb)/255
                                cols <- apply(rgbcols,2,function(x) rgb(x[1],x[2],x[3],alpha=0.15))

                                n <- 0
                                for(k in 1:ny)
                                    {
                                        yy <- (np-1)/(nsim-1)*(c(1:nsim)-1)+1
                                        polygon(x=c(x$pred[n+1:nsim,4],x$pred[n+nsim:1,5]),y=c(yy,rev(yy)),col=cols[k],border=NA)
                                        n <- n+nsim
                                    }
                                
                            }
                    }
            }
        else #Jointlcmm ou lcmm
            {
                ny <- 1
                nsim <- length(x$pred[,1])
                
                copiex <- x$object
                ysim <- matrix(x$pred[,1],nsim,ny)
                transfo <- matrix(x$pred[,2],nsim,ny)
                estimlink <- as.vector(rbind(ysim,transfo))
                copiex$estimlink <- matrix(estimlink,nsim,2*ny)
                    
                    
                dots <- list(...)
                    
                if(length(list(...)$xlim)==0)
                    {
                        if(x$object$linktype!=3)
                            {    
                                dots <- c(dots,list(xlim=range(x$pred[,-1],na.rm=TRUE)))
                            }
                        else
                            {
                                dots <- c(dots,list(xlim=range(x$pred[-c(1,nrow(x$pred)),-1]),na.rm=TRUE))
                            }
                    }
                dots$add <- add
                if(!missing(legend)) dots$legend <- legend
                if(!missing(legend.loc)) dots$legend.loc <- legend.loc
                
                ##tracer la transfo
                if(inherits(x$object,"lcmm")) do.call("plot.lcmm",c(dots,list(x=copiex,which="link")))
                if(inherits(x$object,"Jointlcmm")) do.call("plot.Jointlcmm",c(dots,list(x=copiex,which="link")))
                
                
                if(length(list(...)$lwd)==0)
                    {
                        dots <- c(dots,list(lwd=1))
                    }
                
                if(length(list(...)$col)==0)
                    {
                        dots <- c(dots,list(col=1))
                    }
                
                    
                    if(ncol(x$pred)>2)
                        {
                            ##on ne garde que les options lwd et col et on impose lty=2
                            dots.bornes <- dots[intersect(names(dots),c("lwd","col"))]
                            dots.bornes <- c(dots.bornes,lty=2) 
                            
                            copiex <- x$object
                            ysim <- matrix(x$pred[,1],nsim,ny)
                            transfo <- matrix(x$pred[,3],nsim,ny)
                            estimlink <- as.vector(rbind(ysim,transfo))
                            copiex$estimlink <- matrix(estimlink,nsim,2*ny)

                            if(shades==FALSE)
                                {
                                    ##tracer la borne inf
                                    if(inherits(x$object,"lcmm")) do.call("plot.lcmm",c(dots.bornes,list(x=copiex,which="link",add=TRUE,legend=NULL)))
                                    if(inherits(x$object,"Jointlcmm")) do.call("plot.Jointlcmm",c(dots.bornes,list(x=copiex,which="link",add=TRUE,legend=NULL)))
                                    
                                    copiex <- x$object
                                    ysim <- matrix(x$pred[,1],nsim,ny)
                                    transfo <- matrix(x$pred[,4],nsim,ny)
                                    estimlink <- as.vector(rbind(ysim,transfo))
                                    copiex$estimlink <- matrix(estimlink,nsim,2*ny)
                                    
                                    ##tracer la borne sup
                                    if(inherits(x$object,"lcmm")) do.call("plot.lcmm",c(dots.bornes,list(x=copiex,which="link",add=TRUE,legend=NULL)))
                                    if(inherits(x$object,"Jointlcmm")) do.call("plot.Jointlcmm",c(dots.bornes,list(x=copiex,which="link",add=TRUE,legend=NULL)))
                                }
                            else
                                {
                                    ## shades ##
                                    rgbcols <- sapply(dots$col,col2rgb)/255
                                    cols <- apply(rgbcols,2,function(x) rgb(x[1],x[2],x[3],alpha=0.15))
                                    
                                    polygon(x=c(x$pred[1:nsim,3],x$pred[nsim:1,4]),y=c(x$pred[1:nsim,1],x$pred[nsim:1,1]),col=cols[1],border=NA)    
                                }
                        }
                
            }
    }

