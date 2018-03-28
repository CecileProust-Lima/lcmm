
gridsearch <- function(m,rep,maxiter,minit)
    {
        mc <- match.call()$m
        mc$maxiter <- maxiter
        
        models <- vector(mode="list",length=rep)
        assign("minit",eval(minit))

        for(k in 1:rep)
            {
                mc$B <- substitute(random(minit),environment())
                models[[k]] <- do.call(as.character(mc[[1]]),as.list(mc[-1]))
            }
        llmodels <- sapply(models,function(x){return(x$loglik)})
        kmax <- which.max(llmodels)

        mc$B <- models[[kmax]]$best
        mc$maxiter <- NULL
        
        return(do.call(as.character(mc[[1]]),as.list(mc[-1])))
    }



