#' Compute the variance matrix
#' 
#' The function computes the variance matrix of the random effects, the correlations, or the measurement error. 
#' 
#' @param model an object inheriting from class \code{hlme}, \code{lcmm}, 
#' \code{Jointlcmm} or \code{multlcmm} representing a general latent class
#' mixed model.
#' @param which either "random", "cor", "randomY", or "error".
#' @param times1 for \code{which = "cor"}, numeric vector containing the times at which the correlation should be computed
#' @param times2 for \code{which = "cor"}, numeric vector containing the times at which the correlation should be computed
#' @param nmes for \code{which = "randomY"} and \code{which = "error"}, the number of measures (ie, the dimension of the returned variance matrix)
#' @return a matrix
#' @author Viviane Philipps 
#' @export 
#' @examples
#' \dontrun{
#'  ## The model
#'  m <- hlme(fixed = Y ~ Time, mixture = ~1, random = ~1, subject = "ID",
#'   ng = 2, cor = BM(Time), data = data_hlme, B = c(0,20,30,-1,5,2,0.1))
#'
#' ## The random effects' variance matrix (the variance of the random intercept )
#' B <- createVarianceMatrix(m, which = "random")
#'
#' ## The variance of the Brownian motion at time c(1, 2, 3, 4) 
#' W <- createVarianceMatrix(m, which = "cor", times1 = c(1, 2, 3, 4), times2 = c(1, 2, 3, 4))
#'
#' ## The variance of the measurement error at 4 visit times
#' S <- createVarianceMatrix(m, which = "error", nmes = 4)
#'
#' ## In model "m", the variance matrix of the outcome at times c(1, 2, 3, 4) is:
#' matrix(1, nrow = 4) %*% B %*% t(matrix(1, nrow = 4)) + W + S
#' }
createVarianceMatrix <- function(model, which, times1, times2, nmes){
    if(missing(model)) stop("model is missing")
    if(!inherits(model, c("hlme", "lcmm", "multlcmm", "Jointlcmm"))) stop("Use only with hlme, lcmm, multlcmm, and Jointlcmm models")
    if(!(which %in% c("random", "cor", "randomY", "error"))) stop("which should be one of: random, cor, randomY, error")

    ## random ##
    if(which == "random"){
        ## -> variance-covariance matrix of the random effects

        nea <- sum(model$idea)
        if(nea == 0) stop("No random effects in this model!")
        idiag <- model$idiag
        B <- matrix(0, nea, nea)
        
        ## hlme models ##
        if(inherits(model, "hlme")){
            if(idiag == 1)
                diag(B) <- model$best[sum(model$N[1:2]) + 1:nea]
            else {
                B[upper.tri(B, diag = TRUE)] <- model$best[sum(model$N[1:2]) + 1:model$N[3]]
                B <- t(B)
                B[upper.tri(B, diag = TRUE)] <- model$best[sum(model$N[1:2]) + 1:model$N[3]]
            }

            if(model$N[4] > 0){ # nwg
                res <- vector("list", length = model$ng)
                for(g in 1:(model$ng - 1)){
                    res[[g]] <- model$best[sum(model$N[1:3]) + g]^2 * B
                }
                res[[model$ng]] <- B
            } else {
                res <- B
            }
        }

        ## lcmm models ##
        if(inherits(model, "lcmm")){
            if(idiag == 1)
                diag(B) <- model$best[sum(model$N[1:2]) + 1:nea]
            else {
                B[upper.tri(B, diag = TRUE)] <- model$best[sum(model$N[1:2]) + 1:model$N[3]]
                B <- t(B)
                B[upper.tri(B, diag = TRUE)] <- model$best[sum(model$N[1:2]) + 1:model$N[3]]
            }

            if(model$N[4] > 0){ # nwg
                res <- vector("list", length = model$ng)
                for(g in 1:(model$ng - 1)){
                    res[[g]] <- model$best[sum(model$N[1:3]) + g]^2 * B
                }
                res[[model$ng]] <- B
            } else {
                res <- B
            }
        }

        ## multlcmm models ##
        if(inherits(model, "multlcmm")){
            if(model$N[4] > 0){
                if(idiag == 1)
                    diag(B) <- c(1, model$best[model$N[3] + 1:(nea - 1)])
                else {
                    B[upper.tri(B, diag = TRUE)] <- c(1, model$best[model$N[3] + 1:model$N[4]])
                    B <- t(B)
                    B[upper.tri(B, diag = TRUE)] <- c(1, model$best[model$N[3] + 1:model$N[4]])
                }
            } else {
                B[1, 1] <- 1
            }

            if(model$N[5] > 0){ # nwg
                res <- vector("list", length = model$ng)
                for(g in 1:(model$ng - 1)){
                    res[[g]] <- model$best[sum(model$N[3:4]) + g]^2 * B
                }
                res[[model$ng]] <- B
            } else {
                res <- B
            }
        }
        

        ## Jointlcmm models ##
        if(inherits(model, "Jointlcmm")){
            if(idiag == 1)
                diag(B) <- model$best[sum(model$N[1:4]) + 1:nea]
            else {
                B[upper.tri(B, diag = TRUE)] <- model$best[sum(model$N[1:4]) + 1:model$N[5]]
                B <- t(B)
                B[upper.tri(B, diag = TRUE)] <- model$best[sum(model$N[1:4]) + 1:model$N[5]]
            }

            if(model$N[6] > 0){ # nwg
                res <- vector("list", length = model$ng)
                for(g in 1:(model$ng - 1)){
                    res[[g]] <- model$best[sum(model$N[1:5]) + g]^2 * B
                }
                res[[model$ng]] <- B
            } else {
                res <- B
            }
        }
    }

    ## cor ##
    if(which == "cor"){
        if(missing(times1)) stop("times1 is missing")
        if(missing(times2)) stop("times2 is missing")

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
        if(is.na(ncor) | (ncor == 0)) stop("No 'cor' is this model!")

        ## correlation parameters
        bcor <- rep(NA, ncor)
        if(inherits(model, "hlme"))
            bcor <- model$best[sum(model$N[1:4]) + 1:ncor]
        if(inherits(model, "lcmm"))
            bcor <- model$best[length(model$best) - ncor + 1:ncor]
        if(inherits(model, "multlcmm"))
            bcor <- model$best[sum(model$N[3:6]) + 1:ncor]
        if(inherits(model, "Jointlcmm"))
            bcor <- model$best[sum(model$N[1:6]) + 1:ncor]
        
        ## BM ##
        if(ncor == 1){
            res <- bcor^2 * outer(times1, times2, pmin)
        } else {
            ## AR ##
            res <- bcor[2]^2 * exp(-bcor[1] * outer(times1, times2, function(x, y) {abs(x - y)}))
        }
    }

    ## randomY ##
    if(which == "randomY"){
        if(!inherits(model, "multlcmm")) stop("No randomY is this model!")

        res <- matrix(0, sum(nmes), sum(nmes))
        if(model$N[6] > 0){
            sumnmes <- 0
            for(k in 1:model$N[6]){
                res[sumnmes + 1:nmes[k], sumnmes + 1:nmes[k]] <- model$best[sum(model$N[c(3, 4, 5, 7, 8)]) + k]^2
                sumnmes <- sumnmes + nmes[k]
            }
        }
    }

    ## error ##
    if(which == "error"){
        if(!inherits(model, "multlcmm")){
            if(inherits(model, "hlme"))
                sigma <- model$best[sum(model$N[1:5]) + 1]
            if(inherits(model, "lcmm"))
                sigma <- 1
            if(inherits(model, "Jointlcmm")){
                if(model$linktype == -1)
                    sigma <- model$best[sum(model$N[1:7]) + 1]
                else
                    sigma <- 1
            }
            
            res <- sigma^2 * diag(sum(nmes))

        } else {
            res <- rep(model$best[sum(model$N[c(3,4,5,7)]) + 1:model$N[8]]^2, nmes) * diag(sum(nmes))
        }
    }

    return(res)
}

