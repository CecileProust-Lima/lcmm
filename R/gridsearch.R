#' Automatic grid search
#' 
#' This function provides an automatic grid search for latent class mixed
#' models estimated with \code{hlme}, \code{lcmm}, \code{multlcmm} and
#' \code{Jointlcmm} functions.
#' 
#' The function permits the estimation of a model from a grid of random initial
#' values to reduce the odds of a convergence towards a local maximum.
#' 
#' The function was inspired by the emEM technique described in Biernacki et
#' al. (2003). It consists in:
#' 
#' 1. randomly generating \code{rep} sets of initial values for \code{m} from
#' the estimates of \code{minit} (this is done internally using option
#' \code{B=random(minit)} \code{rep} times)
#' 
#' 2. running the optimization algorithm for the model specified in \code{m}
#' from the \code{rep} sets of initial values with a maximum number of
#' iterations of \code{maxit} each time.
#' 
#' 3. retaining the estimates of the random initialization that provides the
#' best log-likelihood after \code{maxiter} iterations.
#' 
#' 4. running the optimization algorithm from these estimates for the final
#' estimation.
#' 
#' @param m a call of \code{hlme}, \code{lcmm}, \code{multlcmm}, 
#' \code{Jointlcmm} or \code{mpjlcmm} corresponding to the model to estimate
#' @param rep the number of departures from random initial values
#' @param maxiter the number of iterations in the optimization algorithm
#' @param minit an object of class \code{hlme}, \code{lcmm}, \code{multlcmm},
#' \code{Jointlcmm} or \code{mpjlcmm} corresponding to the same model as specified
#' in m except for the number of classes (it should be one). This object is used to
#' generate random initial values
#' @param cl a cluster created by makeCluster from package parallel or an integer
#' specifying the number of cores to use for parallel computation
#' @return an object of class \code{hlme}, \code{lcmm}, \code{multlcmm},
#' \code{Jointlcmm} or \code{mpjlcmm} corresponding to the call specified in m.
#' @author Cecile Proust-Lima and Viviane Philipps
#' @references Biernacki C, Celeux G, Govaert G (2003). Choosing Starting
#' Values for the EM Algorithm for Getting the Highest Likelihood in
#' Multivariate Gaussian Mixture models. Computational Statistics and Data
#' Analysis, 41(3-4), 561-575.
#' @examples
#' 
#' \dontrun{
#' # initial model with ng=1 for the random initial values
#' m1 <- hlme(Y ~ Time * X1, random =~ Time, subject = 'ID', ng = 1, 
#'       data = data_hlme)
#' 
#' # gridsearch with 10 iterations from 50 random departures
#' m2d <- gridsearch(rep = 50, maxiter = 10, minit = m1, hlme(Y ~ Time * X1,
#'       mixture =~ Time, random =~ Time, classmb =~ X2 + X3, subject = 'ID',
#'           ng = 2, data = data_hlme))
#'         }
#' 
#' @export
#' 
gridsearch <- function(m,rep,maxiter,minit,cl=NULL)
{
    mc <- match.call()$m
    mc$maxiter <- maxiter
    
    models <- vector(mode="list",length=rep)
    assign("minit",eval(minit))

    ncl <- NULL
    
    ## parallel version
    if(!is.null(cl))
    {
        if(!inherits(cl,"cluster"))
        {
            if(!is.numeric(cl)) stop("argument cl should be either a cluster or a numeric value indicating the number of cores")

            ncl <- cl
            cl <- makeCluster(ncl)
        }
        
        ## set different seeds
        clusterSetRNGStream(cl)
        
        ## export univariate models if using mpjlcmm
        if(mc[[1]]=="mpjlcmm")
        {
            for (k in 2:length(mc[[2]]))
            {
                clusterExport(cl,list(as.character(mc[[2]][k])))
            }
        }
        
        ## export other arguments
        clusterExport(cl, list("mc", "maxiter", "minit", as.character(as.list(mc[-1])$data)), envir = environment())
        
        ## get and export loaded packages
        pck <- .packages()
        dir0 <- find.package()
        dir <- sapply(1:length(pck),function(k){gsub(pck[k],"",dir0[k])})
        clusterExport(cl,list("pck","dir"),envir=environment())
        clusterEvalQ(cl,sapply(1:length(pck),function(k){require(pck[k],lib.loc=dir[k],character.only=TRUE)}))

        ## fit models
        cat("Be patient, grid search is running ...\n")
        
        models <- parLapply(cl, 1:rep, function(X){
            mc$B <- substitute(random(minit),parent.frame(n=2))
            return(do.call(as.character(mc[[1]]),as.list(mc[-1])))
        })

        cat("Search completed, performing final estimation\n")

        if(!is.null(ncl)) stopCluster(cl)
    }
    else
    {   ## sequential version
        for(k in 1:rep)
        {
            mc$B <- substitute(random(minit),environment())
            models[[k]] <- do.call(as.character(mc[[1]]),as.list(mc[-1]))
        }
    }
    
    ## find max
    llmodels <- sapply(models,function(x){return(x$loglik)})
    kmax <- which.max(llmodels)

    mc$B <- models[[kmax]]$best
    mc$maxiter <- match.call()$m$maxiter
    
    return(do.call(as.character(mc[[1]]),as.list(mc[-1])))
}



