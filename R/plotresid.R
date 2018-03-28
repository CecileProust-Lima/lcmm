.plotresid <- function(x,...)
{

  if (any(x$linktype==3))
  {
   cat("Residual and prediction plots are not available yet for threshold mixed models \n")
  }
  else
  {
   if(x$conv==1|x$conv==2|x$conv==3)
   {
    dots <- list(...)
    
    if(length(list(...)$main)) 
    {
     main1 <- as.character(eval(match.call()$main))
     dots <- dots[setdiff(names(dots),"main")]
    }
    else main1 <- c("marginal residuals versus marginal predictions", "subject-specific residuals versus subject-specific predictions","Normal QQ Plot for marginal residuals","Normal QQ Plot for subject-specific residuals")
    main1 <- rep(main1,length.out=4)            
    
    if(length(list(...)$xlab)) 
    {
     xlab1 <- as.character(eval(match.call()$xlab))
     dots <- dots[setdiff(names(dots),"xlab")]
    }
    else xlab1 <- c("marginal predictions","subject-specific predictions","theoretical quantiles","theoretical quantiles")
    xlab1 <- rep(xlab1,length.out=4)
 
    if(length(list(...)$ylab)) 
    {
     ylab1 <- as.character(eval(match.call()$ylab))
     dots <- dots[setdiff(names(dots),"ylab")]
    }
    else ylab1 <- c("marginal residuals","subject-specific residuals","sample quantiles","sample quantiles")
    ylab1 <- rep(ylab1,length.out=4)    

     if(!length(dots$cex)) dots$cex <- 0.5
     if(!length(dots$pch)) dots$pch <- 1
     if(!length(dots$col)) dots$col <- 1
     if(!length(dots$lty)) dots$lty <- 1
    
 
    names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
    "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
    "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
    "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
    "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
    dots.plot <- dots[intersect(names(dots),names.plot)]  
  
    par(mfrow=c(2,2))
    on.exit(par(mfrow=c(1,1)))

    ## graphe des residus
    do.call("matplot",c(dots.plot,list(x=x$pred[,"pred_m"],y=x$pred[,"resid_m"],main=main1[1],xlab=xlab1[1],ylab=ylab1[1])))
    do.call("matplot",c(dots.plot,list(x=x$pred[,"pred_ss"],y=x$pred[,"resid_ss"],main=main1[2],xlab=xlab1[2],ylab=ylab1[2])))

    ## qqplot des residus

    dots.plot$lty <- c(dots$lty,2,2)
    dots.plot$type <- c("p","l","l")
 
    resid <- sort(x$pred[,"resid_m"])
    i <- 1:length(resid)
    p <- (i-0.375)/(length(resid)+0.25)
    sdp <- sqrt(p*(1-p)/length(p))
    infq <- qnorm(p)-1.96*sdp/(1/sqrt(2*pi)*exp(-(qnorm(p)**2)/2))
    supq <- qnorm(p)+1.96*sdp/(1/sqrt(2*pi)*exp(-(qnorm(p)**2)/2))

    do.call("matplot",c(dots.plot,list(x=qnorm(p),y=cbind(scale(resid),infq,supq),xlab=xlab1[3],ylab=ylab1[3],main=main1[3])))


    resid <- sort(x$pred[,"resid_ss"])
    i <- 1:length(resid)
    p <- (i-0.375)/(length(x$pred[,5])+0.25)
    sdp <- sqrt(p*(1-p)/length(p))
    infq <- qnorm(p)-1.96*sdp/(1/sqrt(2*pi)*exp(-(qnorm(p)**2)/2))
    supq <- qnorm(p)+1.96*sdp/(1/sqrt(2*pi)*exp(-(qnorm(p)**2)/2))

    do.call("matplot",c(dots.plot,list(x=qnorm(p),y=cbind(scale(resid),infq,supq),xlab=xlab1[4],ylab=ylab1[4],main=main1[4])))
    
   }
   else
   {
    cat("Output can not be produced since the program stopped abnormally.")
   }
  }
}
