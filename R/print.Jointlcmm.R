print.Jointlcmm <- function(x,...)
{

 if (!inherits(x, "Jointlcmm")) stop("use only with \"Jointlcmm\" objects")

 cat("Joint latent class model for quantitative outcome and competing risks", "\n")
 cat("     fitted by maximum likelihood method", "\n")

 cl <- x$call
 cl$B <- NULL
 if(is.data.frame(cl$data))
     {
         cl$data <- NULL
         x$call$data <- NULL    
     }
 cat(" \n")
 dput(cl)
 cat(" \n")

 posfix <- eval(cl$posfix)

 cat("Statistical Model:", "\n")
 cat(paste("     Dataset:", x$call$data),"\n")
 cat(paste("     Number of subjects:", x$ns),"\n")

 cat(paste("     Number of observations:", x$N[9]),"\n")
 cat(paste("     Number of latent classes:", x$ng), "\n")
 cat(paste("     Number of parameters:", length(x$best))," \n")
 if(length(posfix)) cat(paste("     Number of estimated parameters:", length(x$best)-length(posfix))," \n")

  nbevt <- length(x$hazard[[1]])
  nprisq <- rep(NA,nbevt)
  nrisq <- rep(NA,nbevt)

  typrisq <- x$hazard[[1]]
  hazardtype <- x$hazard[[2]]
  nz <- x$hazard[[4]] 

 for(ke in 1:nbevt)
     {
        if(typrisq[ke]==1) nprisq[ke] <- nz[ke]-1
        if(typrisq[ke]==2) nprisq[ke] <- 2
        if(typrisq[ke]==3) nprisq[ke] <- nz[ke]+2

        if(hazardtype[ke]=="Common") nrisq[ke] <- nprisq[ke]
        if(hazardtype[ke]=="PH") nrisq[ke] <- nprisq[ke]+x$ng-1
        if(hazardtype[ke]=="Specific") nrisq[ke] <- nprisq[ke]*x$ng

        cat(paste("     Event",ke,": \n"))
        cat(paste("        Number of events: ", x$N[9+ke],"\n",sep=""))
        if(x$ng>1)
            {
                if (hazardtype[ke]=="Specific") cat("        Class-specific hazards and \n")
                if (hazardtype[ke]=="PH") cat("        Proportional hazards over latent classes and \n")
                if (hazardtype[ke]=="Common") cat("        Common hazards over classes and \n")
            }
        
        if (typrisq[ke]==2)
            {
                cat("        Weibull baseline risk function \n")
            }
        if (typrisq[ke]==1)
            {
                cat("        Piecewise constant baseline risk function with nodes \n")
                cat("        ",x$hazard[[3]][1:nz[ke],ke]," \n")
            }
        if (typrisq[ke]==3)
            {
                cat("        M-splines constant baseline risk function with nodes \n")
                cat("        ",x$hazard[[3]][1:nz[ke],ke]," \n")
            }
        
        
    }

 ntrtot <- x$N[8]
 numSPL <- 0
 if(x$linktype!=-1)
     {
         cat(paste("     Link function for ",x$Names$Yname," : ",sep=""))
         if (x$linktype==0)
             {
                 cat("Linear \n")
             }
         if (x$linktype==1)
             {
                 cat("Standardised Beta CdF \n")
             }
         if (x$linktype==2) 
             {
                 cat("Quadratic I-splines with nodes ", x$linknodes ,"\n")
             }
     }

 cat(" \n")
 cat("Iteration process:", "\n")

 if(x$conv==1) cat("     Convergence criteria satisfied")
 if(x$conv==2) cat("     Maximum number of iteration reached without convergence")
 if(x$conv==3) cat("     Convergence with restrained Hessian matrix")
 if(x$conv==4|x$conv==12)
 {
  cat("     The program stopped abnormally. No results can be displayed.\n")
 }
 else
 {
  cat(" \n")
  cat("     Number of iterations: ", x$niter, "\n")
  cat("     Convergence criteria: parameters=", signif(x$gconv[1],2), "\n")
  cat("                         : likelihood=", signif(x$gconv[2],2), "\n")
  cat("                         : second derivatives=", signif(x$gconv[3],2), "\n")
  cat(" \n")
  cat("Goodness-of-fit statistics:", "\n")
  cat(paste("     maximum log-likelihood:", round(x$loglik,2))," \n")
  cat(paste("     AIC:", round(x$AIC,2))," \n")
  cat(paste("     BIC:", round(x$BIC,2))," \n")
  cat(" \n")
 }

}
