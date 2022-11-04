#' Predictions of the random-effects
#' 
#' The function computes the predicted values of the random effects given observed data
#' provided in input. With multiple latent classes, these predictions are averaged over 
#' classes using the posterior class-membership probabilities.
#' 
#' @param model an object inheriting from class \code{hlme}, \code{lcmm}, 
#' \code{Jointlcmm} or \code{multlcmm} representing a general latent class
#' mixed model.
#' @param newdata data frame containing the data from which predictions are to be computed. 
#' The data frame should include at least all the covariates listed in model$Xnames2, 
#' the marker(s) values and the grouping structure. Names should match exactly the names 
#' of the variables in the model.
#' @param subject character specifying the name of the grouping structure.
#' If NULL (the default), the same as in the model will be used.
#' @return a matrix containing the grouping structure and the predicted random-effects.
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
#' predictRE(m2b,newdata=paquid[1:6,])
#' }


predictRE <- function(model, newdata, subject=NULL){
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
  w <- options()$warn
  options(warn=-1)
  on.exit(options(warn=w))
  newmodel <- do.call(argfunction, c(arguments))
  return(newmodel$predRE)
} 

