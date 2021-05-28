#' individual predictions of the random-effects
#' 
#' This function provides table containing individual predictions of the random-effects with new data : a column per random-effect, a line per subject 
#' 
#' @param model an object inheriting from class \code{lcmm}, 
#' \code{Jointlcmm} or \code{multlcmm} representing a general latent class
#' mixed model.
#' @param newdata data frame containing the data from which predictions are computed. The data frame should include at least all the covariates listed in model$Xnames2. Names in the data frame should be exactly model$Xnames2 that are the names of covariates specified in lcmm, hlme, Jointlcmm or multlcmm calls.
#' @param subject choice of subject, by default: NULL, meaning the subject identifier is the same as in the model.
#' @param verbose logical indicating if information about computation should be
#' reported. Default to FALSE.
#' @return dataframe containing the subject identifier and individual predictions of the random-effects.
#' @author Cecile Proust-Lima, Viviane Philipps, Sasha Cuau
#' @export 
#' @examples
#' paquid$normMMSE <- normMMSE(paquid$MMSE)
#' paquid$age65 <- (paquid$age - 65)/10
#' newdata <- data.frame(ID=c(0,0,0), age65=c(0.5, 1, 1.5), CEP=0, normMMSE=c(80,70,60))
#' m2b <- hlme(normMMSE ~ age65+I(age65^2)+CEP, random =~ age65+I(age65^2), subject = 'ID',data = paquid, ng = 2, mixture =~ age65+I(age65^2), B = c(0, 60, 40, 0, -4, 0, -10, 10, 212.869397, -216.421323,456.229910, 55.713775, -145.715516, 59.351000, 10.072221))
#' predictRE(m2b,newdata)


predictRE <- function(model, newdata,subject=NULL,verbose = FALSE){
  arguments<-as.list(model$call)
  argfunction <- as.character(arguments[[1]]) 
  arguments[[1]]<- NULL
  arguments[["data"]]<-newdata
  arguments[["B"]]<-model$best 
  arguments[["maxiter"]]<- 0
  if(!is.null(subject)){
    arguments[['subject']]<- subject
  }
  newmodel <- do.call(argfunction , c(arguments,verbose = verbose))
  return(newmodel$predRE)
} 

getwd()
