#' Posterior classification and class-membership probabilities
#' 
#' This function provides the posterior classification and posterior individual 
#' class-membership probabilities for external data.
#' 
#' @param model an object inheriting from class \code{hlme}, \code{lcmm}, 
#' \code{Jointlcmm} or \code{multlcmm} representing a general latent class
#' mixed model.
#' @param newdata data frame containing the data from which predictions are to be computed.
#' The data frame should include at least all the covariates listed in model$Xnames2,
#' the outcome(s) and the grouping structure. Names should match exactly.
#' @param subject character specifying the name of the grouping structure.
#' If NULL (the default), the same as in the model will be used.
#' @return a matrix with 2+ng columns: the grouping structure, the predicted class and the
#' ng posterior class-membership probabilities.
#' @author Sasha Cuau, Viviane Philipps, Cecile Proust-Lima
#' @export 
#' @examples
#' \dontrun{
#' library(NormPsy)
#' paquid$normMMSE <- normMMSE(paquid$MMSE)
#' paquid$age65 <- (paquid$age - 65)/10
#' m2b <- hlme(normMMSE ~ age65+I(age65^2)+CEP, random =~ age65+I(age65^2), subject = 'ID',
#' data = paquid, ng = 2, mixture =~ age65+I(age65^2), B = c(0, 60, 40, 0, -4, 0, -10, 10,
#' 212.869397, -216.421323,456.229910, 55.713775, -145.715516, 59.351000, 10.072221))
#' predictClass(m2b, newdata=paquid[1:6,])
#' }


predictClass <- function(model, newdata, subject=NULL){
  arguments <- as.list(model$call)
  argfunction <- as.character(arguments[[1]]) 
  arguments[[1]] <- NULL
  arguments[["data"]] <- newdata
  arguments[["B"]] <- model$best 
  arguments[["maxiter"]] <- 0
  arguments[["posfix"]] <- NULL
  arguments[["verbose"]] <- FALSE
  if(!is.null(subject)){
    arguments[['subject']] <- subject
  }
  ## specify manually the spline nodes to be sure that they are the same as initially
  if(any(model$linktype %in% c(1,2)) & !inherits(model, "mpjlcmm")){
      alink <- eval(arguments[["link"]])
      if(any(grepl("-quant-", alink)) | any(grepl("-equi-", alink))){
          arguments[["link"]] <- alink
          arguments[["link"]] <- gsub("-quant-", "-manual-", arguments[["link"]])
          arguments[["link"]] <- gsub("-equi-", "-manual-", arguments[["link"]])
      }
      if(is.matrix(model$linknodes)){ # multlcmm case
          arguments[["intnodes"]] <- as.numeric(apply(model$linknodes[,which(model$linktype == 2)],2,function(x){x[2:(which(x==max(x,na.rm=TRUE))-1)]}))
          arguments[["range"]] <- as.numeric(apply(model$linknodes[,which(model$linktype %in% c(1,2))],2,function(x){range(x,na.rm=TRUE)}))
      } else { # other models
          arguments[["intnodes"]] <- model$linknodes[-c(1, length(model$linknodes))]
          arguments[["range"]] <- model$linknodes[c(1, length(model$linknodes))]
      }
  }
  ## for Jointlcmm, specify hazard nodes
  if(any(model$hazard[[1]] %in% c(1,3))){
      ahazard <- eval(arguments[["hazard"]])
      if(any(grepl("-quant-", ahazard)) | any(grepl("-equi-", ahazard))){
          arguments[["hazard"]] <- ahazard
          arguments[["hazard"]] <- gsub("-quant-", "-manual-",  arguments[["hazard"]])
          arguments[["hazard"]] <- gsub("-equi-", "-manual-",  arguments[["hazard"]])
      }
      hnodes <- model$hazard[[3]] # matrix zi
      hcall <- eval(model$call$hazard)
      nbevt <- length(model$hazard[[1]])
      if(length(hcall)==1 | all(model$hazard[[1]]==3))
      {
          arguments[["hazardnodes"]] <- as.numeric(hnodes[-c(1, nrow(hnodes)),1])
          arguments[["hazardrange"]] <- as.numeric(hnodes[c(1, nrow(hnodes)),1])
      }
      else
      {
          kk <- order(na.omit(c(which(model$hazard[[1]]==1), which(model$hazard[[1]]==3)[1]))) # all piecewise and first spline
          arguments[["hazardnodes"]] <- as.numeric(apply(hnodes[,kk],2,function(x){x[2:(which(x==max(x,na.rm=TRUE))-1)]}))
          arguments[["hazardrange"]] <- as.numeric(apply(hnodes[,kk],2,function(x){range(x,na.rm=TRUE)}))
      }
      
  }
  w <- options()$warn
  options(warn=-1)
  on.exit(options(warn=w))
  newmodel <- do.call(argfunction , c(arguments))
  return(newmodel$pprob)
}
