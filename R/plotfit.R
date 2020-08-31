.plotfit <-function(x,var.time,break.times=NULL,marg=TRUE,legend.loc="bottomleft",legend,outcome=1,subset=NULL,shades=FALSE,...)
{
    if(missing(x)) stop("Argument x should be specified")
    if(!(class(x) %in% c("hlme","lcmm","Jointlcmm","multlcmm"))) stop("Use with hlme, lcmm, multlcmm or Jointlcmm only")
    if(missing(var.time)) stop("Argument var.time should be specified")
    if(class(x)=="lcmm")
        {
            if(x$linktype==3) stop("This function is not available for thresholds mixed models")
        }
    
    if(!is.null(x$data))
    {
        data <- x$data
    }
    else
    {
        data <- eval(x$call$data)
    }
    
    if(!isTRUE(all.equal(as.character(x$call$subset),character(0))))
        {
            cc <- x$call
            cc <- cc[c(1,which(names(x$call)=="subset"))]
            cc[[1]] <- as.name("model.frame")
            cc$formula <- formula(paste("~",paste(colnames(data),collapse="+")))
            cc$data <- data
            cc$na.action <- na.pass
            data <- eval(cc)
            attributes(data)$terms <- NULL
        }

    
    if(class(x) %in% c("hlme","lcmm","Jointlcmm"))
        {
            if(length(x$na.action)) data <- data[-x$na.action,]  
        }
    else
        {
            if(is.numeric(outcome))
                {
                    linesNA <- x$na.action[[outcome]]
                }
            
            if(is.character(outcome))
                {
                    numoutcome <- which(x$Ynames==outcome)
                    if(length(numoutcome)==0) stop("Argument 'outcome' is not correct")
                    linesNA <- x$na.action[[numoutcome]]
                }

            if(length(linesNA)) data <- data[-linesNA,]
        }

    timeInterv <- data[,var.time]
    
    subset2 <- try(as.numeric(subset),silent=TRUE)
    
    #if(!is.null(subset))
    if(length(subset2))
        {
            # data pprob et pred a reduire
            data <- do.call("model.frame",args=list(formula=formula(paste("~",paste(colnames(data),collapse="+"))),subset=sys.call(which=-1)$subset,data=data,na.action=na.pass))

            ids <- unique(data[,x$call$subject])
            
            xorig <- x
            x$pprob <- xorig$pprob[which(xorig$pprob[,1] %in% ids),]
            x$pred <- xorig$pred[which(xorig$pred[,1] %in% ids),]
        }

    
    

    if(is.null(break.times)) break.times <- quantile(timeInterv,prob=seq(0,1,length.out=10))
    else
        {
            if(length(break.times)==1) break.times <- quantile(timeInterv,prob=seq(0,1,length.out=break.times))
        }

    #if(min(times)<min(break.times)) stop("Values in break.times do not cover the entire range of time")
    #if(max(times)>max(break.times)) stop("Values in break.times do not cover the entire range of time")
    
    break.times <- sort(break.times)


    ng <- x$ng
    nea <- sum(x$idea)

    
    ntps <- length(break.times)-1

    numsg <- lapply(1:ng,function(g) x$pprob[which(x$pprob$class==g),1])
    times <- data[,var.time]
    maxT <- sapply(numsg,function(y) max(times[which(data[,x$call$subject] %in% y)]))




    
############  si marg=TRUE ############

    if(marg==TRUE)
        {
            if(class(x) %in% c("hlme","lcmm","Jointlcmm"))
                {
                    pred0 <- x$pred[,c(1,6,6+1:ng),drop=FALSE] #ID,obs,pred_m1...ng
                }
            else
                {
                    if(is.numeric(outcome)) indpred <- which(x$pred[,"Yname"]==x$Ynames[outcome])
                    if(is.character(outcome)) indpred <- which(x$pred[,"Yname"]==outcome)

                    if(length(indpred)==0) stop("Argument outcome is not correct")
                    
                    pred0 <- x$pred[indpred,c(1,7,7+1:ng),drop=FALSE]           
                }

            ppi <- x$pprob[,c(1,2+1:ng),drop=FALSE] #ID,prob1...ng
            name.id <- colnames(x$pred)[1]

    
            pred0 <- cbind(pred0,times)
            pred_complet0 <- merge(pred0,ppi,by=name.id)
 
            pred_times0 <- matrix(NA,length(break.times)-1,4*ng+1)
            pred_times0[,1] <- (break.times[-length(break.times)] + break.times[-1])/2
            colnames(pred_times0) <- c(var.time,paste("pred.class",1:ng,sep=""),paste("obs.class",1:ng,sep=""),
                                       paste("lower.obs.class",1:ng,sep=""),paste("upper.obs.class",1:ng,sep=""))

            
            for(i in 1:ntps)
                {
                    indg <- lapply(1:ng,function(g) which( (pred_complet0[,ng+3]>=break.times[i]) & (pred_complet0[,ng+3]<min(break.times[i+1],maxT[g])) ) )
                    if(i==(length(break.times)-1)) indg <- lapply(1:ng, function(g){ if(maxT[g]>break.times[i]) c(indg[[g]],which( (pred_complet0[,ng+3]==min(break.times[i+1],maxT[g]))  )) else indg[[g]]})
                    
                    if(any(sapply(indg,length)>0))
                        {
                            ##predictions ponderees
                            pred_times0[i,1+1:ng] <- sapply(1:ng,function(g) sum(pred_complet0[indg[[g]],2+g] * pred_complet0[indg[[g]],ng+3+g])/sum(pred_complet0[indg[[g]],ng+3+g]))
                    
                            ##observations ponderees
                            pred_times0[i,1+ng+1:ng] <- sapply(1:ng,function(g) sum(pred_complet0[indg[[g]],2] * pred_complet0[indg[[g]],ng+3+g])/sum(pred_complet0[indg[[g]],ng+3+g]))
                      
                            ##calcul des ecarts-types des moyennes ponderes
                            sd_obs <- sapply(1:ng,function(g) sd(pred_complet0[indg[[g]],2])*sqrt(sum(pred_complet0[indg[[g]],ng+3+g]**2)/sum(pred_complet0[indg[[g]],ng+3+g])**2)  )
                            ##ic des obs moyennes
                            pred_times0[i,1+2*ng+1:ng] <- pred_times0[i,1+ng+1:ng] - 1.96 * sd_obs
                            pred_times0[i,1+3*ng+1:ng] <- pred_times0[i,1+ng+1:ng] + 1.96 * sd_obs
                        }
                    
                }
        }


############ si marg=FALSE ##############
    
    pred_times <- NA
    if(marg==FALSE)
        {
            
            if(class(x) %in% c("hlme","lcmm","Jointlcmm"))
                {
                    pred <- x$pred[,c(1,6,6+ng+1:ng),drop=FALSE] #ID,obs,pred_ss1...ng
                }
            else
                {
                    if(is.numeric(outcome)) indpred <- which(x$pred[,"Yname"]==x$Ynames[outcome])
                    if(is.character(outcome)) indpred <- which(x$pred[,"Yname"]==outcome)

                    if(length(indpred)==0) stop("Argument outcome is not correct")
            
                    pred <- x$pred[indpred,c(1,7,7+ng+1:ng),drop=FALSE]           
                }   
            
            ppi <- x$pprob[,c(1,2+1:ng),drop=FALSE] #ID,prob1...ng
            name.id <- colnames(x$pred)[1]

            pred <- cbind(pred,times)
            pred_complet <- merge(pred,ppi,by=name.id)

            pred_times <- matrix(NA,length(break.times)-1,4*ng+1)
            pred_times[,1] <- (break.times[-length(break.times)] + break.times[-1])/2
            colnames(pred_times) <- c(var.time,paste("pred.class",1:ng,sep=""),paste("obs.class",1:ng,sep=""),
                                      paste("lower.obs.class",1:ng,sep=""),paste("upper.obs.class",1:ng,sep=""))

            for(i in 1:ntps)
                {
                    indg <- lapply(1:ng,function(g) which( (pred_complet[,ng+3]>=break.times[i]) & (pred_complet[,ng+3]<min(break.times[i+1],maxT[g])) ) )
                    if(i==(length(break.times)-1)) indg <- lapply(1:ng, function(g){ if(maxT[g]>break.times[i]) c(indg[[g]],which( (pred_complet[,ng+3]==min(break.times[i+1],maxT[g]))  )) else indg[[g]]})
                                        
                    if(any(sapply(indg,length)>0))
                        {
                            ##predictions ponderees
                            pred_times[i,1+1:ng] <- sapply(1:ng,function(g) sum(pred_complet[indg[[g]],2+g] * pred_complet[indg[[g]],ng+3+g])/sum(pred_complet[indg[[g]],ng+3+g]))

                            ##observations ponderees
                            pred_times[i,1+ng+1:ng] <- sapply(1:ng,function(g) sum(pred_complet[indg[[g]],2] * pred_complet[indg[[g]],ng+3+g])/sum(pred_complet[indg[[g]],ng+3+g]))


                            ##calcul des ecarts-types des moyennes ponderes
                            sd_obs <- sapply(1:ng,function(g) sd(pred_complet[indg[[g]],2])*sqrt(sum(pred_complet[indg[[g]],ng+3+g]**2)/sum(pred_complet[indg[[g]],ng+3+g])**2)  )

                            ##ic des obs moyennes
                            pred_times[i,1+2*ng+1:ng] <- pred_times[i,1+ng+1:ng] - 1.96 * sd_obs
                            pred_times[i,1+3*ng+1:ng] <- pred_times[i,1+ng+1:ng] + 1.96 * sd_obs
                        }

                }
        }
    
    ##parametres graphiques
    dots <- list(...)

    names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
                    "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
                    "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
                    "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
                    "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim")
    dots.plot <- dots[intersect(names(dots),names.plot)]
    dots.seg <- dots[intersect(names(dots),names.plot)]

    names.legend <- c("fill","border","lty","lwd","pch","angle","density","bg","box.lwd","col",
                      "box.lty","box.col","pt.bg","cex","pt.cex","pt.lwd","xjust","yjust","x.intersp","y.intersp","adj","text.width",
                      "text.col","text.font","merge","trace","plot","ncol","horiz","title","xpd","title.col","title.adj","seg.len")
    dots.leg <- dots[intersect(names(dots),names.legend)]


    if(!length(dots.plot$xlab))
        {
            dots.plot$xlab <- var.time
        }

    if(!length(dots.plot$ylab))
        {
            dots.plot$ylab <- "y"
        }

    if(!length(dots.plot$ylim))
        {
            if(marg==TRUE) dots.plot$ylim <- range(pred_times0[,-1],na.rm=TRUE)
            if(marg==FALSE) dots.plot$ylim <- range(pred_times[,-1],na.rm=TRUE)
        }


    if(!length(dots.plot$type))
        {
            dots.plot$type <- c(rep("p",ng),rep("l",2*ng))
        }
    else
        {
            dots.plot$type <- rep(dots.plot$type,length.out=3)
            dots.plot$type <- rep(dots.plot$type,each=ng)
        }

    if(!length(dots.plot$col))
        {
            dots.plot$col <- rep(1:ng,3)
            dots.leg$col <- c(rep(0,ng),rep(1:ng,3))
            dots.seg$col <- rep(1:ng,each=ntps)
        }
    else
        {
            color <- rep(dots.plot$col,length.out=ng)
            color <- rep(color,3)
            dots.plot$col <- color
            dots.leg$col <- c(rep(0,ng),color,color[ng+1:ng])
            dots.seg$col <- rep(color[ng+1:ng],each=ntps)
        }

    if(!length(dots.plot$lwd))
        {
            dots.plot$lwd <- rep(1,3*ng)
        }
    else
        {
            dots.plot$lwd <- rep(dots.plot$lwd,length.out=3)
            dots.plot$lwd <- rep(dots.plot$lwd,each=ng)
        }
    dots.plot$lwd[which(dots.plot$type=="p")] <- NA
    dots.seg$lwd <- rep(dots.plot$lwd[2*ng+1:ng],each=ntps)

    if(!length(dots.plot$lty))
        {
            dots.plot$lty <- c(rep(1,2*ng),rep(2,ng))
        }
    else
        {
            dots.plot$lty <- rep(dots.plot$lty,length.out=3)
            dots.plot$lty <- rep(dots.plot$lty,each=ng)
        }
    dots.plot$lty[which(dots.plot$type=="p")] <- 0
    dots.seg$lty <- rep(dots.plot$lty[2*ng+1:ng],each=ntps)
    

    if(!length(dots.plot$pch))
        {
            dots.plot$pch <- rep(c(20,4,4),each=ng)
        }
    else
        {
            dots.plot$pch <- rep(dots.plot$pch,length.out=3)
            dots.plot$pch <- rep(dots.plot$pch,each=ng)
        }
    dots.plot$pch[which(!(dots.plot$type %in% c("p","b")))] <- NA

    dots.leg$lty <- c(rep(NA,ng),dots.plot$lty)
    dots.leg$lwd <- c(rep(NA,ng),dots.plot$lwd)
    dots.leg$pch <- c(rep(NA,ng),dots.plot$pch)

    if(!length(dots.leg$ncol))
        {
            dots.leg$ncol <- 4
        }

    if(!length(dots.leg$seg.len))
        {
            dots.leg$seg.len <- 1.5
        }
    

    if(!length(dots.leg$box.lty))
        {
            dots.leg$box.lty <- 0
        }

    if(!length(dots.leg$inset))
        {
            dots.leg$inset <- c(0,0.02)
        }

    if(!length(dots.leg$x.intersp))
        {
            dots.leg$x.intersp <- 0.8
        }

    if(!missing(legend))
        {
            if(!is.null(legend))
                {
                    if(length(legend)==4*ng)
                        {
                            dots.leg$legend <- legend
                        }

                    if(length(legend)==ng)
                        {
                            dots.leg$legend <- c(legend,rep("pred",ng),rep("obs",ng),rep("CI (obs)",ng))
                        }

                    if(length(legend)==(ng+3))
                        {
                            dots.leg$legend <- c(legend[1:ng],rep(legend[ng+1],ng),rep(legend[ng+2],ng),rep(legend[ng+3],ng))
                        }
                }
        }

    if(missing(legend))
        {
            dots.leg$legend <- c(paste("class",1:ng,":"),rep("pred",ng),rep("obs",ng),rep("CI (obs)",ng))
            legend <- TRUE
        }





################ si marg=TRUE #############

    if(marg==TRUE)
        {

            bout <- min(diff(pred_times0[,1]))/10*ifelse(all(na.omit(dots.plot$lwd)<3),max(na.omit(dots.plot$lwd)),2)
    
            if(!length(dots.plot$main))
                {
                    dots.plot$main <- "Weighted marginal predictions"
                }
            else
                {
                    dots.plot$main <- list(...)$main[1]
                }
                                        #print(head(pred_times0))
            ##tracer les pred et les obs moyennes :
            do.call(matplot,c(dots.plot,list(x=pred_times0[,1,drop=FALSE],y=pred_times0[,1+1:(2*ng),drop=FALSE])))

            ##tracer les IC :
            if(dots.plot$type[3*ng]=="n")
                {
                    do.call(segments,c(dots.seg[setdiff(names(dots.seg),"type")],
                                       list(x0=pred_times0[,rep(1,ng),drop=FALSE],
                                            y0=pred_times0[,1+2*ng+1:ng,drop=FALSE],
                                            x=pred_times0[,rep(1,ng),drop=FALSE],
                                            y=pred_times0[,1+3*ng+1:ng,drop=FALSE])))
                    ##tracer les bouts des IC :
                    do.call(segments,c(dots.seg[setdiff(names(dots.seg),"type")],
                                      list(x0=pred_times0[,rep(1,ng),drop=FALSE]-bout,
                                            y0=pred_times0[,1+2*ng+1:ng,drop=FALSE],
                                            x=pred_times0[,rep(1,ng),drop=FALSE]+bout,
                                            y=pred_times0[,1+2*ng+1:ng,drop=FALSE])))
                    do.call(segments,c(dots.seg[setdiff(names(dots.seg),"type")],
                                      list(x0=pred_times0[,rep(1,ng),drop=FALSE]-bout,
                                            y0=pred_times0[,1+3*ng+1:ng,drop=FALSE],
                                            x=pred_times0[,rep(1,ng),drop=FALSE]+bout,
                                            y=pred_times0[,1+3*ng+1:ng,drop=FALSE])))
                }
            else
                {
                    if(shades==FALSE)
                        {
                            do.call(matlines,c(lapply(dots.plot,"[",2*ng+1:ng),
                                     list(x=pred_times0[,rep(1,2*ng),drop=FALSE],
                                        y=pred_times0[,1+2*ng+1:(2*ng),drop=FALSE])))
                        }
                    else
                        {
                            rgbcols <- sapply(dots.plot$col,col2rgb)/255
                            colors2 <- apply(rgbcols,2,function(x) rgb(x[1],x[2],x[3],alpha=0.15))
                            ##sapply(1:ng, function(k,t,binf,bsup,cols) polygon(x=unlist(c(t,rev(t))),y=c(binf[,k],rev(bsup[,k])),col=cols[k],border=NA), t=pred_times0[,1],binf=pred_times0[,1+2*ng+1:ng,drop=FALSE],bsup=pred_times0[,1+3*ng+1:ng,drop=FALSE],cols=colors2)
                            ombre <- function(k,t0,binf0,bsup0,cols)
                                {
                                    mat <- na.omit(cbind(t0,binf0[,k],bsup0[,k]))
                                    t <- mat[,1]
                                    binf <- mat[,2]
                                    bsup <- mat[,3]

                                    polygon(x=unlist(c(t,rev(t))),y=c(binf,rev(bsup)),col=cols[k],border=NA)
                                }
                            sapply(1:ng,ombre,t0=pred_times0[,1],binf0=pred_times0[,1+2*ng+1:ng,drop=FALSE],bsup0=pred_times0[,1+3*ng+1:ng,drop=FALSE],cols=colors2)
                                                        
                        }
                }
            ##ajouter la legende :
            if(!is.null(legend)) do.call("legend",c(dots.leg,list(x=legend.loc)))

            return(invisible(pred_times0))
        }


    
################ si marg=FALSE  #############
    
    if(marg==FALSE)
        {
            bout <- min(diff(pred_times[,1]))/10*ifelse(all(na.omit(dots.plot$lwd)<3),max(na.omit(dots.plot$lwd)),2)
            
            if(!length(list(...)$main))
                {
                    dots.plot$main <- "Weighted subject-specific predictions"   
                }
            else
                {
                    dots.plot$main <- list(...)$main[length(list(...)$main)]
                }
            ##tracer les pred et les obs moyennes :
            do.call(matplot,c(dots.plot,list(x=pred_times[,1,drop=FALSE],y=pred_times[,1+1:(2*ng),drop=FALSE])))

            ##tracer les IC :
            if(dots.plot$type[3*ng]=="n")
                {
                    do.call(segments,c(dots.seg[setdiff(names(dots.seg),"type")],
                                       list(x0=pred_times[,rep(1,ng),drop=FALSE],
                                            y0=pred_times[,1+2*ng+1:ng,drop=FALSE],
                                            x=pred_times[,rep(1,ng),drop=FALSE],
                                            y=pred_times[,1+3*ng+1:ng,drop=FALSE])))
                    ##tracer les bouts des IC :
                    do.call(segments,c(dots.seg[setdiff(names(dots.seg),"type")],
                                       list(x0=pred_times[,rep(1,ng),drop=FALSE]-bout,
                                            y0=pred_times[,1+2*ng+1:ng,drop=FALSE],
                                            x=pred_times[,rep(1,ng),drop=FALSE]+bout,
                                            y=pred_times[,1+2*ng+1:ng,drop=FALSE])))
                    do.call(segments,c(dots.seg[setdiff(names(dots.seg),"type")],
                                       list(x0=pred_times[,rep(1,ng),drop=FALSE]-bout,
                                            y0=pred_times[,1+3*ng+1:ng,drop=FALSE],
                                            x=pred_times[,rep(1,ng),drop=FALSE]+bout,
                                            y=pred_times[,1+3*ng+1:ng,drop=FALSE])))
                }
            else
                {
                    if(shades==FALSE)
                        {
                            do.call(matlines,c(lapply(dots.plot,"[",2*ng+1:ng),
                                       list(x=pred_times[,rep(1,2*ng),drop=FALSE],
                                            y=pred_times[,1+2*ng+1:(2*ng),drop=FALSE])))
                        }
                    else
                        {
                            rgbcols <- sapply(dots.plot$col,col2rgb)/255
                            colors2 <- apply(rgbcols,2,function(x) rgb(x[1],x[2],x[3],alpha=0.15))
                            ##sapply(1:ng, function(k,t,binf,bsup,cols) polygon(x=unlist(c(t,rev(t))),y=c(binf[,k],rev(bsup[,k])),col=cols[k],border=NA), t=pred_times[,1],binf=pred_times[,1+2*ng+1:ng,drop=FALSE],bsup=pred_times[,1+3*ng+1:ng,drop=FALSE],cols=colors2) ## pb avec les NA
                            ombre <- function(k,t0,binf0,bsup0,cols)
                                {
                                    mat <- na.omit(cbind(t0,binf0[,k],bsup0[,k]))
                                    t <- mat[,1]
                                    binf <- mat[,2]
                                    bsup <- mat[,3]

                                    polygon(x=unlist(c(t,rev(t))),y=c(binf,rev(bsup)),col=cols[k],border=NA)
                                }
                            sapply(1:ng,ombre,t0=pred_times[,1],binf0=pred_times[,1+2*ng+1:ng,drop=FALSE],bsup0=pred_times[,1+3*ng+1:ng,drop=FALSE],cols=colors2)
                            
                        }
                }
            
            ##ajouter la legende :
            if(!is.null(legend)) do.call("legend",c(dots.leg,list(x=legend.loc)))

            return(invisible(pred_times))
        }
    

}
