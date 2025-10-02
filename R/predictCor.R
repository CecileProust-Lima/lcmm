#' Prediction of the Brownian motion or the autoregressive correlation
#' 
#' The function computes the predicted values of the BM or AR correlation given observed data
#' provided in input. 
#' 
#' @param model an object inheriting from class \code{hlme}, \code{lcmm}, 
#' \code{Jointlcmm} or \code{multlcmm} representing a general latent class
#' mixed model.
#' @param newdata data frame containing the data from which predictions are to be computed. 
#' The data frame should include at least all the covariates listed in model$Xnames2, and
#' the marker(s) values. Names should match exactly the names 
#' of the variables in the model.
#' @param predtimes numeric vector containing the prediction times
#' @return a matrix containing the predicted correlation in each latent class
#' @author Viviane Philipps, Cecile Proust-Lima 
#' @export 
#' @examples
#' \dontrun{
#'  m <- hlme(fixed = Y ~ Time, mixture = ~1, random = ~1, subject = "ID",
#'   ng = 2, cor = BM(Time), data = data_hlme, B = c(0,20,30,-1,5,2,0.1))
#'  predictCor(m, newdata = data_hlme[1:3, ], predtimes = seq(0, 5, 0.5))
#' }
predictCor <- function(model, newdata, predtimes){

    ## prediction = cov(W,Y) * (Var(Y))^(-1) * (Y-E(Y))

    ## type of correlation
    ncor <- NA
    if(inherits(model, "hlme"))
        ncor <- model$N[5]
    else if(inherits(model, "lcmm"))
        ncor <- model$N[6]
    else if(inherits(model, "multlcmm"))
        ncor <- model$N[7]
    else if(inherits(model, "Jointlcmm"))
        ncor <- model$N[7]
    if(is.na(ncor) | (ncor == 0)) stop("No prediction can be computed")
    
    ## name of the time variable for the correlation 
    cor.var.time <- model$Xnames[which(model$idcor == 1)]

    ## remove missing values
    fixed <- model$call$fixed
    fixed <- gsub("contrast","",fixed)
    fixed <- formula(paste(fixed[2],fixed[3],sep="~"))
    if(sum(model$idea)) random <- as.formula(model$call$random) else random <- ~ -1
    cor <- formula(paste("~-1 +", cor.var.time))
    sansNA <- removeNA(list(fixed, random, cor), newdata)
    newdata1 <- sansNA$newdata
    predtimes1 <- as.numeric(na.omit(predtimes))
    nmes <- sansNA$nmes
            
    ## compute cov(W, Y)
    covWY <- createVarianceMatrix(model, which = "cor", times1 = predtimes1, times2 = newdata1[, cor.var.time])

    ## compute Z, B, W, Valea(randomY), Verr
    Z <- model.matrix(random, data = newdata1)
    
    if(ncol(Z))
        B <- createVarianceMatrix(model, which = "random")
    else
        B <- matrix(0, 0, 0)
    
    W <- createVarianceMatrix(model, which = "cor", times1 = newdata1[, cor.var.time], times2 = newdata1[, cor.var.time])
    
    if(inherits(model, "multlcmm"))
        Valea <- createVarianceMatrix(model, which = "randomY", nmes = nmes)
    else
        Valea <- matrix(0, sum(nmes), sum(nmes))
    
    Verr <- createVarianceMatrix(model, which = "error", nmes = nmes)
    
    ## compute E(Y)
    EY <- predictY(model, newdata = newdata1)$pred

    ## Y
    Y <- newdata1[, "outcome"]


    ## compute prediction
    if(!is.list(B)){
        
        VarY <- Z %*% B %*% t(Z) + W + Valea + Verr
        
        res <- covWY %*% solve(VarY) %*% (Y - EY)
        
    } else {
      
        res <- matrix(NA, length(predtimes1), model$ng)
        for(g in 1:model$ng){

            VarYg <- Z %*% B[[g]] %*% t(Z) + W + Valea + Verr
            
            res[, g] <- covWY %*% solve(VarYg) %*% (Y - EY[, g])
        }
    }

    colnames(res) <- paste("corpred_class", 1:model$ng, sep = "")
    return(res)
}
    
