#' Plot of information functions
#' 
#' This function plots the information functions stemmed
#' from a \code{lcmm} or \code{multlcmm} object with ordinal outcomes modeled via threshold links.
#' 
#' @param x an object inheriting from classes \code{ItemInfo}
#' @param which character specifying the values to plot. Should be one of 'ItemInfo' for
#' the Fisher information function of the ordinal outcomes, 'LevelInfo' for the
#' information of each item's level or 'LevelProb' for the probability of the item's
#' levels.  Default to 'ItemInfo'.
#' @param outcome character specifying the outcome to consider. Default to "all".
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
#' @author Viviane Philipps and Cecile Proust-Lima
#' @export 
plot.ItemInfo <- function(x, which="ItemInfo", outcome="all", legend.loc="topright", legend=NULL, add=FALSE, shades=TRUE, ...)
{
    if(missing(x)) stop("The argument x is missing")
    if(which=="ItemProb") stop("Please use the predictYcond function to obtain the item's expectation")
    if(!(which %in% c("ItemInfo", "LevelInfo", "LevelProb"))) stop("Argument which should be one of ItemInfo, LevelInfo or LevelProb")
 
    linktype <- x$object$linktype
    Ynames <- x$object$Ynames

    if(!(outcome %in% c("all", x$object$Ynames))) stop("Undefined outcome specified")
    if(outcome != "all")
    {
        if(length(outcome)>1) stop("Please specify only one outcome")
        if(linktype[which(Ynames==outcome)] != 3) stop("The specifed outcome is not modeled as an ordinal scale")
    }

    if(which!="ItemInfo" & outcome=="all") stop("Please specify only one outcome for level information/probability")
   
    dots <- list(...)
    dots <- dots[setdiff(names(dots),c("x","y","log"))]

    if(!length(dots$main))
    {
        if(outcome=="all")
        {
            dots$main <- switch(which, "ItemInfo"="Item information",
                            "LevelInfo"="Level information",
                            "LevelProb"="Level Probability")
        }
        else
        {
            dots$main <- switch(which, "ItemInfo"=paste("Information of item",outcome),
                                "LevelInfo"=paste("Information of",outcome,"levels"),
                                "LevelProb"=paste("Probability of", outcome, "levels"))
        }
    }
    
    if(!length(dots$col))
    {
        if(which=="ItemInfo")
        {
            if(outcome=="all") dots$col <- rainbow(length(Ynames))[which(linktype==3)]
            else dots$col <- rainbow(length(Ynames))[which(Ynames==outcome)]
        }
        else
        {
            dots$col <- 1:x$object$nbmod[which(Ynames==outcome)]
        }
    }

    if(!length(dots$type))
    {
        dots$type <- "l"
    }
    
    if(!length(dots$lty))
    {
        dots$lty <- 1
    }
    
    if(!length(dots$ylab))
    {
        dots$ylab <- ""
    }
    
    if(!length(dots$xlab))
    {
        dots$xlab <- "Latent process"
    }

    if(missing(legend))
    {
        if(which=="ItemInfo")
        {
            if(outcome=="all") legend <- Ynames[which(linktype==3)]
        }
        else
        {
            legend <- paste("Level",x$object$modalites[[which(Ynames==outcome)]])
        }
    }
    
    if(!length(dots$box.lty))
    {
        dots$box.lty <- 0
    }
    
    if(!length(dots$inset))
    {
        dots$inset <- c(0.02,0.02)
    }

    
    
    lambda <- unique(x$ItemInfo[,2])

    Ypred <- NA
    lower <- NA
    upper <- NA

    if(which=="ItemInfo")
    {
        lignes <- 1:nrow(x$ItemInfo)
        if(outcome!="all") lignes <- which(x$ItemInfo[,1]==outcome)
       
        Ypred <- matrix(x$ItemInfo[lignes,3], nrow=length(lambda))
        if(x$IC)
        {
            lower <- matrix(x$ItemInfo[lignes,4], nrow=length(lambda))
            upper <- matrix(x$ItemInfo[lignes,5], nrow=length(lambda))
        }
    }
    
    if(which=="LevelInfo")
    {
        lignes <- which(x$LevelInfo[,1]==outcome)
        
        Ypred <- matrix(x$LevelInfo[lignes,4+2*x$IC+1], nrow=length(lambda))
        if(x$IC)
        {
            lower <- matrix(x$LevelInfo[lignes,4+2*x$IC+2], nrow=length(lambda))
            upper <- matrix(x$LevelInfo[lignes,4+2*x$IC+3], nrow=length(lambda))
        }        
    }
    
    if(which=="LevelProb")
    {
        lignes <- which(x$LevelInfo[,1]==outcome)
        
        Ypred <- matrix(x$LevelInfo[lignes,4], nrow=length(lambda))
        if(x$IC)
        {
            lower <- matrix(x$LevelInfo[lignes,5], nrow=length(lambda))
            upper <- matrix(x$LevelInfo[lignes,6], nrow=length(lambda))
        }
    }

    if(!length(dots$ylim) )
    {
        dots$ylim <- range(cbind(Ypred,lower,upper), na.rm=TRUE)
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
        do.call("matplot",c(dots.plot,list(x=lambda,y=Ypred)))
        
        if(x$IC)
        {
            if(shades==FALSE)
            {
                dots.plot$lty <- 2
                do.call("matlines",c(dots.plot[names(dots.plot)],list(x=lambda,y=cbind(lower,upper))))
            }
            else
            {
                rgbcols <- sapply(dots$col,col2rgb)/255
                cols <- apply(rgbcols,2,function(x) rgb(x[1],x[2],x[3],alpha=0.15))
                
                sapply(1:ncol(lower),function(k,t,yl,yu,cols) polygon(x=unlist(c(t,rev(t))),y=c(yl[,k],rev(yu[,k])),col=cols[k],border=NA),t=lambda,yl=lower,yu=upper,cols=cols)
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
        do.call("matpoints",c(dots.plot,list(x=lambda,y=Ypred)))

        if(x$IC)
        {
            if(shades==FALSE)
            {
                dots.plot$lty <- 2
                do.call("matlines",c(dots.plot[names(dots.plot)],list(x=lambda,y=cbind(lower,upper))))
            }
            else
            {
                rgbcols <- sapply(dots$col,col2rgb)/255
                cols <- apply(rgbcols,2,function(x) rgb(x[1],x[2],x[3],alpha=0.15))
                
                sapply(1:ncol(lower),function(k,t,yl,yu,cols) polygon(x=unlist(c(t,rev(t))),y=c(yl[,k],rev(yu[,k])),col=cols[k],border=NA),t=lambda,yl=lower,yu=upper,cols=cols)
            }
        }     
    }

    return(invisible(Ypred))
}
