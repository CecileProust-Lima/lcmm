#' Summary of models
#' 
#' This function provides a table summarizing the results of different models
#' fitted by \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}.
#' 
#' 
#' @param m1 an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or
#' \code{Jointlcmm}
#' @param \dots  further arguments, in particular other objects of class
#' \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}
#' @param which character vector indicating which results should be returned.
#' Possible values are "G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC",
#' "entropy", "\%class".
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
summarytable <- function(m1,...,which=c("G","loglik","npm","BIC","%class"))
    {
        if(missing(m1)) stop("At least one model should be specified")
        if(!(class(m1) %in% c("hlme","lcmm","multlcmm","Jointlcmm","mpjlcmm"))) stop("Use with 'hlme', 'lcmm' , 'multlcmm', 'Jointlcmm' or 'mpjlcmm' objects only")
        if(any(!(which %in% c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "scoretest", "%class")))) stop(paste("which should contain elements among",paste(c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", "entropy", "scoretest", "%class"),collapse=", ")))

        dots <- list(...)

        ismodel <- sapply(dots, function(m) ifelse(class(m) %in% class(m1),TRUE,FALSE))
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
            z <- as.matrix(log(x$pprob[,c(3:(x$ng+2))])*x$pprob[,c(3:(x$ng+2))])
            if(any(!is.finite(z)))
            {
                z[which(!is.finite(z))] <- 0
            }
            res <- 1+sum(z)/(x$ns*log(x$ng))
            if(x$ng==1) res <- 1
            return(res)
        }

        ## tous les indicateurs pour le modele m1
        pscoretest <- NA
        if(class(m1)=="Jointlcmm") pscoretest <- round((1-pchisq(m1$scoretest[1],sum(m1$idea))),4)
        tmp <- c(m1$ng,m1$loglik,m1$conv,length(m1$best)-length(eval(m1$call$posfix)),m1$AIC,m1$BIC,sabic(m1),entropy(m1),pscoretest)
        for(g in 1:m1$ng)
            {
                tmp <- c(tmp,length(which(m1$pprob[,2]==g))/m1$ns*100)
            }
        if(m1$ng<max(ng)) tmp <- c(tmp,rep(NA,max(ng)-m1$ng))
        
        res <- matrix(tmp,nrow=1,ncol=9+max(ng))
        colnames(res) <- c("G","loglik","conv","npm","AIC","BIC","SABIC","entropy","scoretest",paste("%class",1:max(ng),sep=""))

        ## tous les indicateurs pour les autres modeles
        if(nbmodels>0)
            {
                for(i in 2:(nbmodels+1))
                    {
                        m <- get(paste("m",i,sep=""))
                        pscoretest <- NA
                        if(class(m)=="Jointlcmm") pscoretest <- round((1-pchisq(m$scoretest[1],sum(m$idea))),4)
                        tmp <- c(m$ng,m$loglik,m$conv,length(m$best)-length(eval(m$call$posfix)),m$AIC,m$BIC,sabic(m),entropy(m),pscoretest)
                        
                        for(g in 1:m$ng)
                            {
                                tmp <- c(tmp,length(which(m$pprob[,2]==g))/m$ns*100)
                            }
                        if(m$ng<max(ng)) tmp <- c(tmp,rep(NA,max(ng)-m$ng))

                        res <- rbind(res,tmp)
                        
                    }
            }

        ## selection des indicateurs
        which2 <- setdiff(which,"%class")
        if("%class" %in% which) which2 <- c(which2,paste("%class",1:max(ng),sep=""))
        res <- res[,which2,drop=FALSE]
        

        rownames(res) <- noms.mod
        prmatrix(res,na.print="")

        return(invisible(res))
    }
