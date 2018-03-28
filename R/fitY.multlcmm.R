
fitY.multlcmm <- function(x)
    {
        if(missing(x)) stop("The model should be specified")
        if(!inherits(x,"multlcmm")) stop("Use with 'multlcmm' objects only")

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

        #if(length(x$na.action)) data <- data[-x$na.action,] 

        id <- unique(data[,x$call$subject])

        ny <- length(x$Ynames)

        pred <- NULL
        
        for(i in 1:length(id))
            {
                pred <- rbind(pred,predictY(x,newdata=data[which(data[,x$call$subject]==id[i]),],draws=FALSE,methInteg=1,nsim=2000,na.action=1)$pred)
            }

        linesNA <- apply(data[,c(x$call$subject,x$Xnames2),drop=FALSE],2,function(v) which(is.na(v)))
        linesNA <- unique(unlist(linesNA))

        if(length(linesNA))
            {
                data2 <- data[-linesNA,,drop=FALSE]
                nmes <- as.vector(table(data2[,x$call$subject]))
                idres <- rep(unique(data2[,x$call$subject]),nmes*ny)
            }
        else
            {
                nmes <- as.vector(table(data[,x$call$subject]))
                idres <- rep(id,nmes*ny)           
            }
        
        
        #res <- cbind(idres,pred)
        res <- data.frame(idres,pred)
        colnames(res) <- c(x$call$subject,"Yname",paste("Ypred_class",1:x$ng,sep=""))

        
        return(res)
    }


fitY <- function(x) UseMethod("fitY")
