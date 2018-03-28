
plot.hlme <- function(x,which="residuals",var.time,break.times,marg,subset,...)
    {
        if(missing(x)) stop("The model should be specified")
        if(!inherits(x,"hlme")) stop("Use with 'hlme' objects only")

        if(!(which %in% c("residuals","postprob","fit"))) stop(paste("Argument 'which' should be one of:",paste(c("residuals","postprob","fit"),collapse=" ")))

        if(which=="residuals") .plotresid(x,...)

        if(which=="postprob") .plotpostprob(x,...)

        if(which=="fit")
            {
                if(missing(var.time)) stop("Argument var.time should be specified")
                if(missing(break.times)) break.times <- NULL
                if(missing(marg)) marg <- TRUE
                if(missing(subset)) subset <- NULL
                #legend.loc?add?
                .plotfit(x,var.time=var.time,break.times,outcome=1,marg=marg,subset=subset,...)

            }

    }




plot.lcmm <- function(x,which="residuals",var.time,break.times,marg,subset,...)
    {
        if(missing(x)) stop("The model should be specified")
        if(!inherits(x,"lcmm")) stop("Use with 'lcmm' objects only")
        
        if(!(which %in% c("residuals","postprob","link","linkfunction","fit"))) stop(paste("Argument 'which' should be one of:",paste(c("residuals","postprob","link","linkfunction","fit"),collapse=" ")))

        if(which=="residuals") .plotresid(x,...)

        if(which=="postprob") .plotpostprob(x,...)
           
        if(which %in% c("link","linkfunction")) .plotlinkfunction(x,...)
        

        if(which=="fit")
            {
                if(missing(var.time)) stop("Argument var.time should be specified")
                if(missing(break.times)) break.times <- NULL
                if(missing(marg)) marg <- TRUE
                if(missing(subset)) subset <- NULL

                .plotfit(x,var.time=var.time,break.times,outcome=1,marg=marg,subset=subset,...)
            }
    }




plot.multlcmm <- function(x,which="residuals",var.time,break.times,marg,outcome,subset,...)
    {
        if(missing(x)) stop("The model should be specified")
        if(!inherits(x,"multlcmm")) stop("Use with 'multlcmm' objects only")
        
        if(!(which %in% c("residuals","postprob","link","linkfunction","fit"))) stop(paste("Argument 'which' should be one of:",paste(c("residuals","postprob","link","linkfunction","fit"),collapse=" ")))

        if(which=="residuals") .plotresid(x,...)

        if(which=="postprob") .plotpostprob(x,...)

        if(which %in% c("link","linkfunction")) .plotlinkfunctionmult(x,...)
           
        if(which=="fit")
            {
                if(missing(var.time)) stop("Argument var.time should be specified")
                if(missing(break.times)) break.times <- NULL
                if(missing(marg)) marg <- TRUE
                if(missing(subset)) subset <- NULL
                if(missing(outcome)) outcome <- 1
                
                .plotfit(x,var.time=var.time,break.times,outcome=outcome,marg=marg,subset=subset,...)
            }

    }




plot.Jointlcmm <- function(x,which="residuals",var.time,break.times,marg,event,subset,...)
    {
        if(missing(x)) stop("The model should be specified")
        if(!inherits(x,"Jointlcmm")) stop("Use with 'Jointlcmm' objects only")
        
        if(!(which %in% c("residuals","postprob","link","linkfunction","fit","hazard","baselinerisk","survival"))) stop(paste("Argument 'which' should be one of:",paste(c("residuals","postprob","link","linkfunction","fit","hazard","baselinerisk","survival"),collapse=" ")))

        if(which=="residuals") .plotresid(x,...)

        if(which=="postprob") .plotpostprob(x,...)

        if(which %in% c("link","linkfunction")) .plotlinkfunction(x,...)

        if(which %in% c("hazard","baselinerisk"))
            {
                if(missing(event)) event <- 1
                .plotbaselinerisk(x,event=event,...) #legend.loc legend add
            }
           
        if(which=="survival") .plotsurvival(x,...)

        if(which=="fit")
            {
                if(missing(var.time)) stop("Argument var.time should be specified")
                if(missing(break.times)) break.times <- NULL
                if(missing(marg)) marg <- TRUE
                if(missing(subset)) subset <- NULL
                
                .plotfit(x,var.time=var.time,break.times,outcome=1,marg=marg,subset=subset,...)
            }

    }




#plot <- function(x,...) UseMethod("plot")
