#' @export
#'
summary.externX = function(object, ...){
  x = object
  if(!inherits(x, "externX")) stop("use only with\"externX\" objects")
  
  cat("Secondary multinomial model for external class predictor", "\n")
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
  cat(paste("     Number of latent classes:", x$ng), "\n")
  cat(paste("     Number of parameters:", length(x$best))," \n")
  
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
    if(x$varest == "paramBoot") {
      cat("     Proportion of convergence on bootstrap iterations (%)=", x$Mconv, "\n")
    } else {
      cat("     Number of iterations: ", x$niter, "\n")
      cat("     Convergence criteria: parameters=", signif(x$gconv[1],2), "\n")
      cat("                         : likelihood=", signif(x$gconv[2],2), "\n")
      cat("                         : second derivatives=", signif(x$gconv[3],2), "\n")
    }
    cat(" \n")
    cat("Goodness-of-fit statistics:", "\n")
    cat(paste("     maximum log-likelihood:", round(x$loglik,2))," \n")
    cat(paste("     AIC:", round(x$AIC,2))," \n")
    cat(paste("     BIC:", round(x$BIC,2))," \n")
    cat(" \n")
    
    cat(" \n")
    cat(" \n")
    cat("Maximum Likelihood Estimates:", "\n")
    cat(" \n")
    
    nprob <- x$N[1]
    NPM <- length(x$best)
    
    ## shorten names if > 20 characters
    names_best <- names(x$best)
    if(any(sapply(names_best, nchar)>20))
    {
      islong <- which(sapply(names_best, nchar)>20)
      split_names_best <- strsplit(names_best, split=":", fixed=TRUE)
      short_names_best <- lapply(split_names_best, gsub, pattern="\\(.*\\)", replacement="(...)")
      new_names <- lapply(short_names_best, paste, collapse=":")
      names_best[islong] <- unlist(new_names)[islong]
      names(x$best) <- names_best
    }
    
    se <- rep(NA,NPM)
    if (x$conv==1 | x$conv==3)
    {
      ##recuperation des indices de V
      id <- 1:NPM
      indice <- id*(id+1)/2
      se <- sqrt(x$V[indice])
      wald <- x$best/se
      pwald <- 1-pchisq(wald**2,1)
      coef <- x$best
    }
    else
    {
      se <- NA
      wald <- NA
      pwald <- NA
      coef <- x$best
      
      sech <- rep(NA,length(coef))
      waldch <- rep(NA,length(coef))
      pwaldch <- rep(NA,length(coef))
    }
    
    ow <- options("warn")
    options(warn=-1) # to avoid warnings with conv=3
    if(x$conv!=2)
    {
      coefch <- format(as.numeric(sprintf("%.5f",coef)),nsmall=5,scientific=FALSE)
      sech <- format(as.numeric(sprintf("%.5f",se)),nsmall=5,scientific=FALSE)
      waldch <- format(as.numeric(sprintf("%.3f",wald)),nsmall=3,scientific=FALSE)
      pwaldch <- format(as.numeric(sprintf("%.5f",pwald)),nsmall=5,scientific=FALSE)
    }
    else
    {
      coefch <- format(as.numeric(sprintf("%.5f",coef)),nsmall=5,scientific=FALSE)
    }
    options(ow)
    
    if(length(posfix))
    {
      coefch[posfix] <- paste(coefch[posfix],"*",sep="")
      sech[posfix] <- ""
      waldch[posfix] <- ""
      pwaldch[posfix] <- ""
    }
    
    ## fct pr determiner la longueur max d'une chaine de caracteres
    ## (avec gestion des NA)
    maxchar <- function(x)
    {
      xx <- na.omit(x)
      if(length(xx))
      {
        res <- max(nchar(xx))
      }
      else
      {
        res <- 2
      }
      return(res)
    }
    
    if(nprob>0)
    {
      cat("Fixed effects in the class-membership model:\n" )
      cat("(the class of reference is the last class) \n")
      
      tmp <- cbind(coefch[1:nprob],sech[1:nprob],waldch[1:nprob],pwaldch[1:nprob])
      maxch <- apply(tmp,2,maxchar)
      if(any(c(1:nprob) %in% posfix)) maxch[1] <- maxch[1]-1
      dimnames(tmp) <- list(names(coef)[1:nprob],
                            c(paste(paste(rep(" ",max(maxch[1]-4,0)),collapse=""),"coef",sep=""),
                              paste(paste(rep(" ",max(maxch[2]-4,0)),collapse=""),"Se**",sep=""),
                              paste(paste(rep(" ",max(maxch[3]-4,0)),collapse=""),"Wald",sep=""),
                              paste(paste(rep(" ",max(maxch[4]-7,0)),collapse=""),"p-value",sep="")))
      
      cat("\n")
      print(tmp,quote=FALSE,na.print="")
      cat("\n")
    }
    
    tTable = tmp
    
    if(length(posfix))
    {
      cat(" *  coefficient fixed by the user \n")
    }
    if(x$varest == "none") cat(" ** total variance estimated witout correction for primary model uncertainty", "\n \n")
    if(x$varest == "Hessian") cat(" ** total variance estimated through the Hessian of the joint likelihood", "\n \n")
    if(x$varest == "paramBoot") cat(" ** total variance estimated through parametric bootstrap", "\n \n")
    
    return(invisible(tTable))
  }
}