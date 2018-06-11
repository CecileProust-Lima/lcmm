#' @export
fitY.Jointlcmm <- function(x)
    {
        if(missing(x)) stop("The model should be specified")
        if(!inherits(x,"Jointlcmm")) stop("Use with 'Jointlcmm' objects only")

        if(!is.null(x$data))
        {
            data <- x$data
        }
        else
        {
            data <- eval(x$call$data)
        }
        
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





#' Marginal predictions of the longitudinal outcome(s) in their natural scale
#' from \code{lcmm}, \code{Jointlcmm} or \code{multlcmm} objects
#' 
#' The function computes the marginal predictions of the longitudinal
#' outcome(s) in their natural scale on the individual data used for the
#' estimation from \code{lcmm}, \code{Jointlcmm} or \code{multlcmm} objects.
#' 
#' 
#' @aliases fitY fitY.lcmm fitY.multlcmm fitY.Jointlcmm
#' @param x an object inheriting from classes \code{lcmm} or \code{multlcmm}.
#' @return For \code{lcmm} and \code{Jointlcmm} objects, returns a matrix with
#' ng+1 columns containing the subject identifier and the ng class-specific
#' marginal predicted values.
#' 
#' For \code{multlcmm} objects, returns a matrix with ng+2 columns containing
#' the subject identifier, the outcome indicator and the ng class-specific
#' predicted values.
#' @author Cecile Proust-Lima, Viviane Philipps
#' @seealso \code{\link{predictY}}, \code{\link{plot.lcmm}}
#' 
#' @export
#' 
fitY <- function(x) UseMethod("fitY")
