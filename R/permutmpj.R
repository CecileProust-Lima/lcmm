permutmpj <- function(m,order,estim=TRUE)
{
    if(m$nbevt>0)
    {
        ## permuter les coef de la survie et classmb
        z <- m$call
        z$longitudinal <- NULL
        z$fixed <- formula(paste(m$Names$Yname[1],"~1"))
        z$random <- ~-1
        z$mixture <- ~1
        z$B <- c(m$best[1:sum(m$N[1:3])], rep(0,m$ng),1)
        z$maxiter <- 0
        z$verbose <- FALSE
        z[[1]] <- as.name("Jointlcmm")
        
        mjoint <- eval(z, envir=parent.frame())
        mjoint_perm <- permut(mjoint, order=order, estim=FALSE)

        bperm <- mjoint_perm$best[1:sum(m$N[1:3])]        
    }
    else
    {
        ## permuter les coef de classmb
        z <- call("hlme", subject=m$call$subject, data=m$call$data, ng=m$call$ng,
                  fixed=formula(paste(m$Names$Yname[1], "~1")), random=~-1, mixture=~1,
                  B=c(m$best[1:m$N[1]], rep(0,m$ng),1), maxiter=0, verbose=FALSE)
        if(length(m$call$classmb)) z$classmb <- m$call$classmb
        mclass <- eval(z, envir=parent.frame())

        mclass_perm <- permut(mclass, order=order, estim=FALSE)
        bperm <- mclass_perm$best[1:m$N[1]]
    }

    
    ## permuter les modeles longitudinaux
    #mlong <- eval(m$call$longitudinal, envir=parent.frame())
    avt <- sum(m$N[1:3])
    for(k in 1:m$K)
    {
        z <- m$longicall[[k]] #mlong[[k]]$call
        z$B <- c(rep(0,m$ng-1), m$best[avt+1:m$npmK[k]])
        z$maxiter <- 0
        z$posfix <- NULL
        z$verbose <- FALSE
        mk <- eval(z, envir=parent.frame())
        
        mk_perm <- permut(mk, order=order, estim=FALSE)        
        bperm <- c(bperm, mk_perm$best[-c(1:(m$ng-1))])

        avt <- avt + m$npmK[k]
    }

    ## modele permute
    if(estim==TRUE)
    {
        z <- m$call
        z$B <- bperm
        mnew <- eval(z, envir=parent.frame())
    }
    else
    {
        mnew <- m
        mnew$best <- bperm
    }
        
    return(mnew)
}
