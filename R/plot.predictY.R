#' Plot of predicted trajectories and link functions
#' 
#' This function provides the class-specific predicted trajectories stemmed
#' from a \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm} object.
#' 
#' 
#' 
#' @aliases plot.predictY.hlme plot.predictY.lcmm plot.predictY.Jointlcmm
#' plot.predictY.multlcmm plot.predictY plot.predictL.lcmm
#' plot.predictL.Jointlcmm plot.predictL.multlcmm plot.predictL
#' plot.predictlink.lcmm plot.predictlink.Jointlcmm plot.predictlink.multlcmm
#' plot.predictlink plot.predict
#' @param x an object inheriting from classes \code{predictL}, \code{predictY}
#' or \code{predictlink} representing respectively the predicted marginal mean
#' trajectory of the latent process, the predicted marginal mean trajectory of
#' the longitudinal outcome, or the predicted link function of a fitted latent
#' class model.
#' @param outcome for \code{predictY} and multivariate model fitted with
#' \code{multlcmm} only, the outcome to consider.
#' @param legend.loc keyword for the position of the legend from the list
#' \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"},
#' \code{"topleft"},\code{"top"}, \code{"topright"}, \code{"right"} and
#' \code{"center"}.
#' @param legend character or expression to appear in the legend. If no legend
#' should be added, \code{"legend"} should be NULL.
#' @param add logical indicating if the curves should be added to an existing
#' plot. Default to FALSE.
#' @param shades logical indicating if confidence intervals should be
#' represented with shades. Default to FALSE, the confidence intervals are
#' represented with dotted lines.
#' @param \dots other parameters to be passed through to plotting functions or
#' to legend
#' @author Cecile Proust-Lima, Benoit Liquet and Viviane Philipps
#' @seealso \code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{Jointlcmm}},
#' \code{\link{multlcmm}}
#' @examples
#' 
#' 
#' ################# Prediction from linear latent class model
#' ## fitted model
#' m<-lcmm(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,
#' subject='ID',ng=2,data=data_hlme,B=c(0.41,0.55,-0.18,-0.41,
#' -14.26,-0.34,1.33,13.51,24.65,2.98,1.18,26.26,0.97))
#' ## newdata for predictions plot
#' newdata<-data.frame(Time=seq(0,5,length=100),
#' X1=rep(0,100),X2=rep(0,100),X3=rep(0,100))
#' plot(predictL(m,newdata,var.time="Time"),legend.loc="right",bty="l")
#' ## data from the first subject for predictions plot
#' firstdata<-data_hlme[1:3,]
#' plot(predictL(m,firstdata,var.time="Time"),legend.loc="right",bty="l")
#' 
#'  \dontrun{
#' ################# Prediction from a joint latent class model
#' ## fitted model - see help of Jointlcmm function for details on the model
#' m3 <- Jointlcmm(fixed= Ydep1~Time*X1,mixture=~Time,random=~Time,
#' classmb=~X3,subject='ID',survival = Surv(Tevent,Event)~X1+mixture(X2),
#' hazard="3-quant-splines",hazardtype="PH",ng=3,data=data_lcmm,
#' B=c(0.7576, 0.4095, -0.8232, -0.2737, 0, 0, 0, 0.2838, -0.6338, 
#' 2.6324, 5.3963, -0.0273, 1.398, 0.8168, -15.041, 10.164, 10.2394, 
#' 11.5109, -2.6219, -0.4553, -0.6055, 1.473, -0.0383, 0.8512, 0.0389, 
#' 0.2624, 1.4982))
#' # class-specific predicted trajectories 
#' #(with characteristics of subject ID=193)
#' data <- data_lcmm[data_lcmm$ID==193,]
#' plot(predictY(m3,newdata=data,var.time="Time"),bty="l")
#' }
#' 
#' @rdname plot.predict
#' @export





plot.predictY <- function(x,outcome=1,legend.loc="topright",legend,add=FALSE,shades=FALSE,...)
    {
        if(missing(x)) stop("The argument \'x\' is missing.")
        if(!inherits(x,"predictY")) stop("use only with \'predictY\' object")
        if(all(is.na(x$times))) stop("No results can be provided for this type of object")

### determiner si mult, draws et ng

        colx <- colnames(x$pred)
        if(colx[1]=="Yname") #multlcmm
            {
                if(is.numeric(outcome)) outcome <- unique(x$pred[,1])[outcome]
                ng <- length(grep("Ypred_50_",colx))

                if(length(grep("_class2",colx))) #ng>1
                    {
                        if(length(grep("Ypred_50",colx))) #draws
                            {
                                ng <- length(grep("Ypred_50_class",colx))
                                Ypred <- x$pred[which(x$pred[,1]==outcome),1+1:ng,drop=FALSE]
                                lower <- x$pred[which(x$pred[,1]==outcome),1+ng+1:ng,drop=FALSE]
                                upper <- x$pred[which(x$pred[,1]==outcome),1+ng+ng+1:ng,drop=FALSE]
                                
                            }
                        else
                            {
                                ng <- length(grep("Ypred_class",colx))
                                Ypred <- x$pred[which(x$pred[,1]==outcome),1+1:ng,drop=FALSE]
                                lower <- NULL
                                upper <- NULL
                            }
                    }
                else
                    {
                        if(ncol(x$pred)>2)
                            {
                                ng <- 1
                                Ypred <- x$pred[which(x$pred[,1]==outcome),2,drop=FALSE]
                                lower <- x$pred[which(x$pred[,1]==outcome),3,drop=FALSE]
                                upper <- x$pred[which(x$pred[,1]==outcome),4,drop=FALSE]
                            }
                        else
                            {
                                ng <- 1
                                Ypred <- x$pred[which(x$pred[,1]==outcome),2,drop=FALSE]
                                lower <- NULL
                                upper <- NULL
                            }
                    }
            }
        else #hlme, lcmm ou Jointlcmm
            {
                if(length(grep("class2",colx))) #ng>1
                    {
                        if(length(grep("Ypred_50",colx)) | length(grep("lower",colx)))
                            {
                                ng <- length(colx)/3 #length(grep("Ypred_50_",colx))
                                Ypred <- x$pred[,1:ng]
                                lower <- x$pred[,ng+1:ng]
                                upper <- x$pred[,ng+ng+1:ng]
                            }
                        else
                            {
                                ng <- length(colx) #length(grep("Ypred_class",colx))
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
               dots$ylab <- "longitudinal outcome"
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
                                if(length(dots.plot$lwd)==3*ng | length(dots.plot$lwd)==2*ng) dots.plot$lwd <- dots.plot$lwd[(ng+1):length(dots.plot$lwd)]
                                if(length(dots.plot$lty)==3*ng | length(dots.plot$lty)==2*ng) dots.plot$lty <- dots.plot$lty[(ng+1):length(dots.plot$lty)]
                                else dots.plot$lty <- 2
                                do.call("matlines",c(dots.plot[names(dots.plot)],list(x=x$times,y=cbind(lower,upper))))
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
                                          "seg.len","inset","col")
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
                                if(length(dots.plot$lwd)==3*ng | length(dots.plot$lwd)==2*ng) dots.plot$lwd <- dots.plot$lwd[(ng+1):length(dots.plot$lwd)]
                                if(length(dots.plot$lty)==3*ng | length(dots.plot$lty)==2*ng) dots.plot$lty <- dots.plot$lty[(ng+1):length(dots.plot$lty)]
                                else dots.plot$lty <- 2
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
