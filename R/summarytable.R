summarytable <- function(m1,...)
    {
        if(missing(m1)) stop("At least one model should be specified")
        if(!(class(m1) %in% c("hlme","lcmm","multlcmm","Jointlcmm"))) stop("Use with 'hlme', 'lcmm' , 'multlcmm', or 'Jointlcmm' objects only")

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


        tmp <- c(m1$ng,m1$loglik,length(m1$best)-length(eval(m1$call$posfix)),m1$BIC)
        for(g in 1:m1$ng)
            {
                tmp <- c(tmp,length(which(m1$pprob[,2]==g))/m1$ns*100)
            }
        if(m1$ng<max(ng)) tmp <- c(tmp,rep(NA,max(ng)-m1$ng))
        
        res <- matrix(tmp,nrow=1,ncol=4+max(ng))
        colnames(res) <- c("G","loglik","npm","BIC",paste("%class",1:max(ng),sep=""))


        if(nbmodels>0)
            {
                for(i in 2:(nbmodels+1))
                    {
                        m <- get(paste("m",i,sep=""))
                        
                        tmp <- c(m$ng,m$loglik,length(m$best)-length(eval(m$call$posfix)),m$BIC)
                        
                        for(g in 1:m$ng)
                            {
                                tmp <- c(tmp,length(which(m$pprob[,2]==g))/m$ns*100)
                            }
                        if(m$ng<max(ng)) tmp <- c(tmp,rep(NA,max(ng)-m$ng))

                        res <- rbind(res,tmp)
                        
                    }
            }
        

        rownames(res) <- noms.mod
        prmatrix(res,na.print="")

        return(invisible(res))
    }
