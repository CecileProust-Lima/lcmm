#' Marginal predictions in the natural scale of a pre-transformed outcome
#'
#' The function computes the predicted values of the longitudinal marker
#' (in each latent class if ng>1) for a specified profile of covariates, when a
#' non-parameterized pre-transformation was applied (e.g., log, square root). 
#' A Gauss-Hermite or Monte-Carlo integration is 
#' used to numerically compute the back-transformed predictions.
#' 
#' @param x an object inheriting from class \code{hlme} representing a general 
#' latent class mixed model.
#' @param newdata data frame containing the data from which predictions are to 
#' be computed. The data frame should include at least all the covariates listed
#' in x$Xnames2. Names in the data frame should be exactly x$Xnames2, i.e., 
#' the names of covariates specified in \code{hlme} calls.
#' @param var.time A character string containing the name of the variable that
#' corresponds to time in the data frame (x axis in the plot).
#' @param methInteg optional integer specifying the type of numerical
#' integration. Value 0 (by default) specifies a Gauss-Hermite integration which
#' is very rapid but neglects the correlation between the predicted values (in
#' presence of random-effects). Value 1 refers to a Monte-Carlo integration
#' which is slower but correctly accounts for the correlation between the
#' predicted values.
#' @param nsim number of points used in the numerical integration.
#' For methInteg=0, nsim should be chosen among
#' the following values: 5, 7, 9, 15, 20, 30, 40 or 50 (nsim=20 by default). If
#' methInteg=1, nsim should be relatively important (more than 200).
#' @param draws boolean specifying whether confidence bands should be computed.
#' If draws=TRUE, a Monte Carlo approximation of the posterior distribution of 
#' the predicted values is computed and the median, 2.5\% and 97.5\% percentiles
#' are given. Otherwise, the predicted values are computed at the point
#' estimate. By default, draws=FALSE.
#' @param ndraws integer. If draws=TRUE, ndraws specifies the number of draws 
#' that should be generated to approximate the posterior distribution of the 
#' predicted values. By default, ndraws=2000.
#' @param na.action Integer indicating how NAs are managed. The default is 1
#' for 'na.omit'. The alternative is 2 for 'na.fail'. Other options such as
#' 'na.pass' or 'na.exclude' are not implemented in the current version.
#' @param back function to back-transform the outcome in the original scale.
#' @param \dots further arguments to be passed to or from other methods. They
#' are ignored in this function.
#' @return An object of class \code{predictY}.
#' 
#' @examples
#' data_lcmm$transfYdep2 <- sqrt(30 - data_lcmm$Ydep2)
#' m1 <- hlme(transfYdep2 ~ Time, random=~ Time, subject="ID", data = data_lcmm)
#' pred1 <- predictYback(m1, newdata = data.frame(Time = seq(0, 3, 0.1)), 
#' var.time = "Time", back = function(x) {30 - x^2})
#' plot(pred1)
#' 
#' @export
#'
predictYback <- function(x, newdata, var.time, methInteg = 0, nsim = 20, 
                  draws = FALSE, ndraws = 2000, na.action = 1, back, ...){


    if(missing(newdata)) stop("The argument newdata should be specified")
    if(missing(x)) stop("The argument x should be specified")
    if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
    if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) {     
        cat("newdata should at least include the following covariates: ", "\n")
        cat(x$Xnames2[-1], "\n")}
    if (!all(x$Xnames2 %in% c(colnames(newdata), "intercept"))) stop("see above")
    if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
    if (!(methInteg %in% c(0,1))) stop("The integration method must be either 0 for Gauss-Hermite or 1 for Monte-Carlo")
    if ((methInteg == 0) & (!(nsim %in% c(5,7,9,15,20,30,40,50)))) stop("For Gauss-Hermite integration method, 'nsim' should be either 5,7,9,15,20,30,40 or 50")

    call_fixed <- x$call$fixed[3]
    if(is.null(x$call$random)) {call_random <- -1} else call_random <- x$call$random[2]
    if(is.null(x$call$classmb)) {call_classmb <- -1} else call_classmb <- x$call$classmb[2]
    if(is.null(x$call$mixture)) {call_mixture <- -1} else call_mixture <- x$call$mixture[2]


    if(x$conv %in% c(1, 2, 3)){

        if(x$Xnames2[1] != "intercept"){
            newdata1 <- newdata[, x$Xnames2]
            colnames(newdata1) <- x$Xnames
            newdata1 <- data.frame(newdata1)
        }else{
            newdata1 <- cbind(rep(1, length = length(newdata[, 1])), newdata[, x$Xnames2[-1]])
            colnames(newdata1) <- c("intercept", x$Xnames2[-1])
            newdata1 <- data.frame(newdata1)
        }


        if((x$conv == 2) & (draws == TRUE))
        {
            cat("No confidence interval will be provided since the program did
                not converge properly \n")
            draws <- FALSE
        }


        if((x$conv == 3) & (draws == TRUE))
        {
            cat("No confidence interval will be provided since the program 
                did not converge properly \n")
            draws <- FALSE
        }

        X1 <- NULL                                                              
        X2 <- NULL
        b1 <- NULL
        b2 <- NULL


        if(!(na.action %in% c(1, 2))) stop("only 1 for 'na.omit' or 2 for 
                            'na.fail' are required in na.action argument") 

        if(na.action == 1){
            na.action <- na.omit
        }else{
            na.action <- na.fail
        }


        ## transform to factor is the variable appears in levels$levelsdata
        for(v in colnames(newdata1))
        {
            if(v %in% names(x$levels$levelsdata))
            {
                if(!is.null(x$levels$levelsdata[[v]]))
                {
                    newdata1[, v] <- factor(newdata1[, v], 
                                            levels = x$levels$levelsdata[[v]])
                }
            }
        }        
        call_fixed <- gsub("factor","",call_fixed)
        call_random <- gsub("factor","",call_random)
        call_classmb <- gsub("factor","",call_classmb)
        call_mixture <- gsub("factor","",call_mixture)   

        call_mixture <- formula(paste("~",call_mixture,sep=""))
        call_random <- formula(paste("~",call_random,sep=""))
        call_classmb <- formula(paste("~",call_classmb,sep=""))

        ## Traitement des donnees manquantes

        mcall <- match.call()[c(1, match(c("data", "subset", "na.action"), 
                                         names(match.call()), 0))]
        mcall$na.action <- na.action
        mcall$data <- newdata1

        ## fixed
        m <- mcall
        m$formula <- formula(paste("~", call_fixed, sep=""))
        m[[1]] <- as.name("model.frame")	
        m <- eval(m, sys.parent()) 
        na.fixed <- attr(m, "na.action")

        ## mixture
        if((length(attr(terms(call_mixture),"term.labels")) + 
            attr(terms(call_mixture),"intercept")) > 0){
            id.X_mixture <- 1
            m <- mcall
            m$formula <- call_mixture
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.mixture <- attr(m,"na.action")
        }else{
            id.X_mixture <- 0
            na.mixture <- NULL
        }

        ## random
        if((length(attr(terms(call_random), "term.labels")) + 
            attr(terms(call_random),"intercept")) > 0){
            id.X_random <- 1
            m <- mcall
            m$formula <- call_random
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.random <- attr(m,"na.action")
        }else{
            id.X_random <- 0
            na.random <- NULL
        }
        
        ## classmb
        if((length(attr(terms(call_classmb),"term.labels")) + 
            attr(terms(call_classmb),"intercept")) > 0){ 
            id.X_classmb <- 1
            m <- mcall	
            m$formula <- call_classmb
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.classmb <- attr(m,"na.action")
        }else{
            id.X_classmb <- 0
            na.classmb <- NULL
        }
        
        ## cor
        na.cor <- NULL
        if(length(x$N) > 5)
        {
            if(x$N[6] > 0)
            {
                z <- which(x$idcor0 == 1)
                var.cor <- newdata1[, x$Xnames[z]]
                na.cor <- which(is.na(var.cor))
            }
        }

        ##var.time
        if(!missing(var.time))
        {
            if(!(var.time %in% colnames(newdata))) stop("'var.time' should 
                                                        be included in newdata")
            if(var.time %in% colnames(newdata1))
            {
                times <- newdata1[, var.time, drop = FALSE]
            }
            else
            {
                times <- newdata[, var.time, drop = FALSE]
            }
        }
        else
        {
            times <- newdata[, 1, drop = FALSE]
        }

        ## Table sans donnees manquante: newdata
        na.action <- unique(c(na.fixed, na.mixture, na.random, na.classmb, na.cor))
        if(length(na.action)){
            newdata1 <- newdata1[-na.action, ]
            times <- times[-na.action, , drop = FALSE]
        }

        ## create one data frame for each formula (useful with factors)
        newdata1fixed <- newdata1
        for(v in colnames(newdata1fixed))
        {
            if(v %in% names(x$levels$levelsfixed))
            {
                if(!is.null(x$levels$levelsfixed[[v]]))
                {
                    newdata1fixed[, v] <- factor(newdata1fixed[, v], 
                                            levels = x$levels$levelsfixed[[v]])
                    if(any(is.na(newdata1fixed[, v]))) stop(paste("Wrong factor 
                                                      level in variable", v))
                }
            }
        }
        newdata1mixture <- newdata1
        for(v in colnames(newdata1mixture))
        {
            if(v %in% names(x$levels$levelsmixture))
            {
                if(!is.null(x$levels$levelsmixture[[v]]))
                {
                    newdata1mixture[, v] <- factor(newdata1mixture[, v], 
                                      levels = x$levels$levelsmixture[[v]])
                    if(any(is.na(newdata1mixture[, v]))) stop(paste("Wrong 
                                                factor level in variable", v))
                }
            }
        }
        newdata1random <- newdata1
        for(v in colnames(newdata1random))
        {
            if(v %in% names(x$levels$levelsrandom))
            {
                if(!is.null(x$levels$levelsrandom[[v]]))
                {
                    newdata1random[, v] <- factor(newdata1random[, v], 
                                          levels = x$levels$levelsrandom[[v]])
                    if(any(is.na(newdata1random[, v]))) stop(paste("Wrong factor 
                                                      level in variable", v))
                }
            }
        }
        newdata1classmb <- newdata1
        for(v in colnames(newdata1classmb))
        {
            if(v %in% names(x$levels$levelsclassmb))
            {
                if(!is.null(x$levels$levelsclassmb[[v]]))
                {
                    newdata1classmb[, v] <- factor(newdata1classmb[, v], 
                                        levels = x$levels$levelsclassmb[[v]])
                    if(any(is.na(newdata1classmb[, v]))) stop(paste("Wrong 
                                              factor level in variable", v))
                }
            }
        }
        


        ## Construction de nouvelles var eplicatives sur la nouvelle table
        ## fixed
	
	X_fixed <- model.matrix(formula(paste("~", call_fixed, sep = "")), 
	                        data = newdata1fixed)
	if(colnames(X_fixed)[1] == "(Intercept)"){
            colnames(X_fixed)[1] <- "intercept"
            int.fixed <- 1
	}
        ## mixture
	if(id.X_mixture == 1){
            X_mixture <- model.matrix(call_mixture, data = newdata1mixture)	
            if(colnames(X_mixture)[1] == "(Intercept)"){
                colnames(X_mixture)[1] <- "intercept"
                int.mixture <- 1
            }
	}	
        ## random
	if(id.X_random == 1){
            X_random <- model.matrix(call_random, data = newdata1random)	
            if(colnames(X_random)[1] == "(Intercept)"){
                colnames(X_random)[1] <- "intercept"
                int.random <- 1
            }
	}	
        ## classmb
	if(id.X_classmb == 1){ 
            X_classmb <- model.matrix(call_classmb, data = newdata1classmb)
            colnames(X_classmb)[1] <- "intercept"
	}
        ##cor	
        if(x$N[5] > 0)  #on reprend la variable de temps de cor
        {
            z <- which(x$idcor0 == 1)
            var.cor <- newdata1[, x$Xnames[z]]
        }

        ## Construction des var expli
        newdata1 <- X_fixed
        colX <- colnames(X_fixed)

        if(id.X_mixture == 1){
            for(i in 1:length(colnames(X_mixture))){
		if((colnames(X_mixture)[i] %in% colnames(newdata1)) == FALSE){
                    newdata1 <- cbind(newdata1, X_mixture[, i])
                    colnames(newdata1) <- c(colX, colnames(X_mixture)[i])
                    colX <- colnames(newdata1)
		}
            }
        }
        if(id.X_random == 1){
            for(i in 1:length(colnames(X_random))){
		if((colnames(X_random)[i] %in% colnames(newdata1)) == FALSE){
                    newdata1 <- cbind(newdata1, X_random[, i])
                    colnames(newdata1) <- c(colX,colnames(X_random)[i])
                    colX <- colnames(newdata1)
		}	 
            }
        }
        if(id.X_classmb == 1){
            for(i in 1:length(colnames(X_classmb))){
		if((colnames(X_classmb)[i] %in% colnames(newdata1)) == FALSE){
                    newdata1 <- cbind(newdata1, X_classmb[,i])
                    colnames(newdata1) <- c(colX, colnames(X_classmb)[i])
                    colX <- colnames(newdata1)
		}	
            }
        }
        if(x$N[5]>0)
        { 
            if((x$idg0[z] == 0) & (x$idea0[z] == 0) & (x$idprob0[z] == 0))
            {
                newdata1 <- cbind(newdata1, var.cor)
                colnames(newdata1) <- c(colX, x$Xnames[z])
                colX <- colnames(newdata1)
            }
        }

        nv <- length(x$idg0)
        maxmes <- length(newdata1[, 1])
        npm <- length(x$best)
        best <- x$best
        if((x$idiag == 0) & (x$N[3] > 0)) best[x$N[1] + x$N[2] + 1:x$N[3]] <- x$cholesky
        if((x$idiag == 1) & (x$N[3] > 0)) best[x$N[1] + x$N[2] + 1:x$N[3]] <- sqrt(best[x$N[1] + x$N[2] + 1:x$N[3]])
        nwg <- x$N[4]
        ncor <- x$N[5]
       
    
        ## integration
        points <- rep(0, x$ng * nsim * maxmes)
        weights <- rep(0, nsim)
        #browser()
        if (!draws){ # without confidence interval
            
            out <- .Fortran(C_integ,
                            as.double(newdata1),
                            as.integer(x$idprob0),
                            as.integer(x$idea0),
                            as.integer(x$idg0),
                            as.integer(x$idcor0),
                            as.integer(x$ng),
                            as.integer(ncor),
                            as.integer(nv),
                            as.integer(maxmes),
                            as.integer(x$idiag),
                            as.integer(nwg),
                            as.integer(npm),
                            as.double(best),
                            as.integer(nsim),
                            as.integer(methInteg),
                            points=as.double(points),
                            weights = as.double(weights))
            
            out$points[which(out$points == 9999)] <- NA
            out$weights[which(out$weights == 9999)] <- NA

            backpoints <- back(out$points)
            
            Ypred <- matrix(NA, nrow = maxmes, ncol = x$ng)
            if(x$ng == 1) colnames(Ypred) <- "Ypred"
            if(x$ng > 1) colnames(Ypred) <- paste("Ypred_class", 
                                                  1:x$ng, sep = "")

            for(g in 1:x$ng)
            {
                gpoints <- matrix(backpoints[(g - 1) * maxmes * nsim 
                                             + 1:(maxmes * nsim)], maxmes, nsim)
                wgpoints <- sweep(gpoints, 2, out$weights, "*")
                gpred <- apply(wgpoints, 1, sum)
                Ypred[, g] <- gpred
            }

        } else { # with CI based on Monte Carlo draws

            ndraws <- as.integer(ndraws)
            ydraws <- NULL

            posfix <- eval(x$call$posfix)
            
            if(ndraws>0)
            {
                Mat <- matrix(0, ncol = npm, nrow = npm)
                Mat[upper.tri(Mat, diag = TRUE)] <- x$V 
                if(length(posfix))
                {
                    Mat2 <- Mat[-posfix, -posfix]
                    Chol2 <- chol(Mat2)
                    Chol <- matrix(0, npm, npm)
                    Chol[setdiff(1:npm, posfix), setdiff(1:npm, posfix)] <- Chol2
                    Chol <- t(Chol)
                }
                else
                {
                    Chol <- chol(Mat)
                    Chol <- t(Chol)
                }
            }
            
            ydraws <- matrix(NA, maxmes * x$ng, ndraws)
            for (j in 1:ndraws)
            {
                bdraw <- rnorm(npm)
                bdraw <- best + Chol %*% bdraw
                
                out <- .Fortran(C_integ,
                                as.double(newdata1),
                                as.integer(x$idprob0),
                                as.integer(x$idea0),
                                as.integer(x$idg0),
                                as.integer(x$idcor0),
                                as.integer(x$ng),
                                as.integer(ncor),
                                as.integer(nv),
                                as.integer(maxmes),
                                as.integer(x$idiag),
                                as.integer(nwg),
                                as.integer(npm),
                                as.double(bdraw),
                                as.integer(nsim),
                                as.integer(methInteg),
                                points=as.double(points),
                                weights = as.double(weights))
                
            
                out$points[which(out$points == 9999)] <- NA
                out$weights[which(out$weights == 9999)] <- NA
                
                backpoints <- back(out$points)
                
                pred <- matrix(NA, nrow = maxmes, ncol = x$ng)
                
                for(g in 1:x$ng)
                {
                    gpoints <- matrix(backpoints[(g - 1) * maxmes * nsim + 1:(maxmes * nsim)], maxmes, nsim)
                    wgpoints <- sweep(gpoints, 2, out$weights, "*")
                    gpred <- apply(wgpoints, 1, sum)
                    pred[, g] <- gpred
                }
                
                ydraws[, j] <- as.numeric(pred)
            }

            f <- function(x) {
                quantile(x[!is.na(x)], probs = c(0.025, 0.5, 0.975))
            }
            
            ydistr <- apply(ydraws, 1, FUN = f)

            Ypred_50 <- matrix(ydistr[2, ], ncol = x$ng, byrow = FALSE)
            Ypred_2.5 <- matrix(ydistr[1, ], ncol = x$ng, byrow = FALSE)
            Ypred_97.5 <- matrix(ydistr[3, ], ncol = x$ng, byrow = FALSE)

            Ypred <- cbind(Ypred_50, Ypred_2.5, Ypred_97.5)
            
            if (x$ng == 1){
                colnames(Ypred) <- c("Ypred_50","Ypred_2.5","Ypred_97.5")
            } else {
                colnames(Ypred) <- c(paste("Ypred_50_class", 1:x$ng, sep = ""), paste("Ypred_2.5_class", 1:x$ng, sep = ""), paste("Ypred_97.5_class", 1:x$ng, sep = ""))
            }
        }

        res.list <- NULL
        res.list$pred <- Ypred
        res.list$times <- times
    }
    else  #cas xconv != 1 ou 2
    { 
        cat("Predictions can not be computed since the program stopped abnormally. \n")
        res.list <- list(pred = NA, times = NA)
    }

    class(res.list) <- "predictY"
    
    return(res.list)
}
