
applyLink <- function(model, newdata){  
  arguments <- as.list(model$call)
  argfunction <- as.character(arguments[[1]]) 
  arguments[[1]] <- NULL
  arguments[["data"]] <- newdata
  arguments[["B"]] <- model$best 
  arguments[["maxiter"]] <- 0
  arguments[["posfix"]] <- NULL
  arguments[["verbose"]] <- FALSE
  
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
  
  if(!inherits(model, "multlcmm"))
      res <- newmodel$pred[, c(1, 6)]
  else
      res <- newmodel$pred[, c(1, 2, 7)]

  return(res)
} 
