
plot.cuminc <- function(x,profil=1,event=1,add=FALSE,legend,legend.loc="topleft",...)
    {
        if(missing(x)) stop("The argument 'x' should be specified")
        if(!inherits(x,"cuminc")) stop("Use with 'cuminc' objects only")
        if(length(profil)>1) stop("Please specify only one profil")
        if(!(profil %in% c(1:length(x)))) stop("Wrong profil number")
        if(length(event)>1) stop("Please specify only one event")
        nbevt <- length(unique(x[[1]][,"event"]))
        if(!(event %in% c(1:nbevt))) stop(paste("event should be between 1 and",nbevt))
        
        
        xx <- x[[profil]]
        mat <- xx[which(xx[,"event"]==event),]

        if("med_class1" %in% colnames(xx) | "50_class1" %in% colnames(xx))
            {
                ic <- 1
                ng <- (ncol(xx)-2)/3
            }
        else
            {
                ic <- 0
                ng <- ncol(xx)-2
            }
        
        dots <- list(...)
        dots <- dots[setdiff(names(dots),c("x","y","log"))]

        if(!length(dots$main))
            {
                dots$main <- "Class-specific cumulative incidence"
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
                if(ic==0) dots$lty <- 1
                else dots$lty <- c(rep(1,ng),rep(2,2*ng))
            }

        if(!length(dots$ylab))
            {
                dots$ylab <- "cumulative incidence"
            }


        if(!length(dots$xlab))
            {
                dots$xlab <- "time"
            }
        
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
        else inset1 <- c(0.02,0.02)
        
        names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis",
                        "cex.lab","cex.main","cex.sub","col","col.axis",
                        "col.lab","col.main","col.sub","crt","err","family","fig",
                        "fin","font","font.axis","font.lab","font.main","font.sub",
                        "frame.plot","lab","las","lend","lheight","ljoin","lmitre",
                        "lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
                        "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub",
                        "tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
                        "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
        dots.plot <- dots[intersect(names(dots),names.plot)]

        if(!isTRUE(add))
            {
                do.call("matplot",c(dots.plot,list(x=mat[,2],y=mat[,-c(1,2)])))
            }
        else
            {
                do.call("matlines",c(dots.plot,list(x=mat[,2],y=mat[,-c(1,2)])))
            }
        
        names.legend <- c("fill","border","lty","lwd","pch","angle","density",
                          "bg","box.lwd","box.lty","box.col","pt.bg","cex","pt.cex",
                          "pt.lwd","xjust","yjust","x.intersp","y.intersp","adj",
                          "text.width","text.col","text.font","merge","trace",
                          "plot","ncol","horiz","title","xpd","title.col",
                          "title.adj","seg.len") 
        
        dots.leg <- dots[intersect(names(dots),names.legend)]
        if(!(dots$type %in% c("l","b"))) dots.leg <- dots[setdiff(names(dots),c("lty","lwd"))]
        
        if(!is.null(legend)) do.call("legend",c(dots.leg,list(x=legend.loc, legend=legend, box.lty=box.lty1, inset=inset1,col=dots$col)))

        return(invisible(NULL))
    }


#<plot.incidcum <- function(x,profil=1,event=1,add=FALSE,legend,legend.loc="topleft",...) UseMethod("plot.incidcum")
