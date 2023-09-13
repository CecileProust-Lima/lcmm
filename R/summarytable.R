#' Summary of models
#' 
#' This function provides a table summarizing the results of different models
#' fitted by \code{hlme}, \code{lcmm}, \code{multlcmm}, \code{Jointlcmm},
#' \code{mpjlcmm} or \code{externVar}.
#' 
#' Can be reported the usual criteria used to assess the fit and the clustering
#'  of the data:
#'  - maximum log-likelihood L (the higher the better)
#'  - number of parameters P, number of classes G, convergence criterion (1 = converged)
#'  - AIC (the lower the better) computed as -2L+2P 
#'  - BIC (the lower the better) computed as -2L+ P log(N) where N is the number of subjects
#'  - SABIC (the lower the better) computed as -2L+ P log((N+2)/24)
#'  - Entropy (the closer to one the better) computed in two ways : ICL1 = 1-sum[pi_ig*log(pi_ig)]/(N*log(G))
#'    where pi_ig is the posterior probability that subject i belongs to class g
#'  - ICL (the lower the better) computed in two ways : ICL1 = BIC - sum[pi_ig*log(pi_ig)]
#'    or ICL2 = BIC - 2*sum(log(max(pi_ig)), where the max is taken over the classes for each subject.
#'  - %Class computed as the proportion of each class based on c_ig
#' 
#' @param m1 an object of class \code{hlme}, \code{lcmm}, \code{multlcmm}, 
#' \code{Jointlcmm}, \code{mpjlcmm}, \code{externVar} or \code{externVar}.
#' @param \dots further arguments, in particular other objects of class
#' \code{hlme}, \code{lcmm}, \code{multlcmm}, \code{Jointlcmm} or \code{mpjlcmm}.
#' @param which character vector indicating which results should be returned.
#' Possible values are "G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC",
#' "entropy", "ICL", "ICL1", "ICL2", "\%class".
#' @param display logical indicating whether the table should be printed (the default) or not (display=FALSE)
#' @return a matrix giving for each model the values of the requested indexes.
#' By default, the number a latent classes, the
#' log-likelihood, the number of parameters, the BIC and the posterior
#' probability of the latent classes.
#' @author Cecile Proust-Lima, Viviane Philipps
#' @seealso \code{\link{summary}}, \code{\link{hlme}}, \code{\link{lcmm}},
#' \code{\link{multlcmm}}, \code{\link{Jointlcmm}}
#' 
#' @export
#'  
summarytable <- function(m1, ..., which=c("G","loglik","npm","BIC","%class"), display=TRUE)
    {
        if(missing(m1)) stop("At least one model should be specified")
        if(!inherits(m1,c("hlme","lcmm","multlcmm","Jointlcmm","mpjlcmm", "externX", "externSurv"))) stop("Use with 'hlme', 'lcmm' , 'multlcmm', 'Jointlcmm', 'mpjlcmm', 'externX' or 'externSurv' objects only")
        if(any(!(which %in% c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "scoretest","ICL", "ICL1", "ICL2","%class")))) stop(paste("which should contain elements among",paste(c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "scoretest", "ICL", "ICL1", "ICL2", "%class"),collapse=", ")))

        dots <- list(...)

        ismodel <- sapply(dots, function(m) ifelse(inherits(m,class(m1)),TRUE,FALSE))
        models <- which(ismodel==TRUE)
        nbmodels <- length(models)
        
        noms.mod <- as.character(match.call()[[2]])
        if(nbmodels>0)
            {
                noms.models <- paste("mm",2:(nbmodels+1),sep="")
                ng <- m1$ng
                for(i in 1:nbmodels)
                    {
                        assign(noms.models[i],dots[[models[i]]])
                        ng <- c(ng,dots[[models[i]]]$ng)
                        noms.mod <- c(noms.mod,as.character(match.call()[[2+models[i]]]))
                    }

                mm1 <- m1
                mm <- c("mm1",noms.models)[order(ng)]
                for(i in 1:(nbmodels+1))
                    {
                        assign(paste("m",i,sep=""),get(mm[i]))
                    }

                noms.mod <- noms.mod[order(ng)]
            }
        else
            {
                ng <- m1$ng
            }

        sabic <- function(m)
        {
            return((length(m$best)-length(m$call$posfix))*log((m$ns+2)/24)-2*m$loglik)
        }

        
        entropy <- function(x)
        {
            z <- log(as.matrix(x$pprob[,c(3:(x$ng+2))]))*as.matrix(x$pprob[,c(3:(x$ng+2))])
            if(any(!is.finite(z)))
            {
                z[which(!is.finite(z))] <- 0
            }
            res <- 1+sum(z)/(x$ns*log(x$ng))
            if(x$ng==1) res <- 1
            return(res)
        }
        
        
        ICL <- function(x)
        {
          ## ICL1 = BIC - sum(log(pprob)*pprob)
          z1 <- log(as.matrix(x$pprob[,c(3:(x$ng+2))]))*as.matrix(x$pprob[,c(3:(x$ng+2))])
          if(any(!is.finite(z1))){ z1[which(!is.finite(z1))] <- 0}
          res1 <- x$BIC - sum(z1)
          if(x$ng==1) res1 <- x$BIC
          
          ## ICL2 = BIC - 2*sum(log(pprobmax))
          z2 <- rep(0,length=length(x$pprob[,1]))
          for(g in 1:x$ng)
          {
            z2[which(x$pprob[,2]==g)] <- log(as.numeric(x$pprob[which(x$pprob[,2]==g),2+g]))
          }
          res2 <- x$BIC - 2*sum(z2)
          if(x$ng==1) res2 <- x$BIC
          
          return(c(res1, res2))
        }
        
        

        ## tous les indicateurs pour le modele m1
        pscoretest <- NA
        if(inherits(m1,"Jointlcmm")) pscoretest <- round((1-pchisq(m1$scoretest[1],sum(m1$idea))),4)
        tmp <- c(m1$ng,m1$loglik,m1$conv,length(m1$best)-length(eval(m1$call$posfix)),m1$AIC,m1$BIC,sabic(m1),entropy(m1),pscoretest,ICL(m1))
        for(g in 1:m1$ng)
            {
                tmp <- c(tmp,length(which(m1$pprob[,2]==g))/m1$ns*100)
            }
        if(m1$ng<max(ng)) tmp <- c(tmp,rep(NA,max(ng)-m1$ng))
        
        res <- matrix(tmp,nrow=1,ncol=11+max(ng))
        colnames(res) <- c("G","loglik","conv","npm","AIC","BIC","SABIC","entropy","scoretest","ICL1", "ICL2",paste("%class",1:max(ng),sep=""))

        ## tous les indicateurs pour les autres modeles
        if(nbmodels>0)
            {
                for(i in 2:(nbmodels+1))
                    {
                        m <- get(paste("m",i,sep=""))
                        pscoretest <- NA
                        if(inherits(m,"Jointlcmm")) pscoretest <- round((1-pchisq(m$scoretest[1],sum(m$idea))),4)
                        tmp <- c(m$ng,m$loglik,m$conv,length(m$best)-length(eval(m$call$posfix)),m$AIC,m$BIC,sabic(m),entropy(m),pscoretest,ICL(m))
                        
                        for(g in 1:m$ng)
                            {
                                tmp <- c(tmp,length(which(m$pprob[,2]==g))/m$ns*100)
                            }
                        if(m$ng<max(ng)) tmp <- c(tmp,rep(NA,max(ng)-m$ng))

                        res <- rbind(res,tmp)
                        
                    }
            }

        ## selection des indicateurs
        if("ICL" %in% which) which <- unique(c(setdiff(which, "ICL"), "ICL1", "ICL2"))
        which2 <- setdiff(which,"%class")
        if("%class" %in% which) which2 <- c(which2,paste("%class",1:max(ng),sep=""))
        res <- res[,which2,drop=FALSE]
        

        rownames(res) <- noms.mod
        if(display){
            prmatrix(res,na.print="")
        }
        return(invisible(res))
    }
