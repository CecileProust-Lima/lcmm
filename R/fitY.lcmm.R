
fitY.lcmm <- function(x)
    {
        if(missing(x)) stop("The model should be specified")
        if(!inherits(x,"lcmm")) stop("Use with 'lcmm' objects only")

        data <- eval(x$call$data)

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

        if(length(x$na.action)) data <- data[-x$na.action,] 

        id <- unique(data[,x$call$subject])

        pred <- NULL
        
        for(i in 1:length(id))
            {
                pred <- rbind(pred,predictY(x,newdata=data[which(data[,x$call$subject]==id[i]),],draws=FALSE,methInteg=1,nsim=2000)$pred)
            }

        #res <- cbind(data[,x$call$subject],pred)
        res <- data.frame(data[,x$call$subject],pred)
        colnames(res) <- c(x$call$subject,paste("Ypred_class",1:x$ng,sep=""))

        
        return(res)
    }



fitY <- function(x) UseMethod("fitY")
