#' @export
#'
VarExpl.hlme <- function(x,values)
{
 if(missing(x)) stop("The model should be specified")
 if (!inherits(x, "hlme")) stop("use only with \"hlme\" objects")
 if(missing(values)) values <- data.frame("intercept"=1)
 if (!inherits(values, "data.frame")) stop("values should be a data.frame object")
 if(any(is.na(values))) stop("values should not contain any missing values")


 if(x$conv==1 | x$conv==2)
 {
  res <- matrix(0,nrow=1,ncol=x$ng)
 
  names.random <- NULL
  name.cor <- NULL
  if(x$N[3]>0) names.random <- x$Xnames[which(x$idea0==1)]
  if(length(x$N)>4 & x$N[5]>0) name.cor <- x$Xnames[which(x$idcor0==1)]

  if(!is.null(names.random) | !is.null(name.cor))
  {
   names.values <- unique(c(names.random,name.cor))   #contient I(T^2))

   vars <- unique(c(all.vars(x$call$random),all.vars(x$call$cor)))
   if(!all(vars %in% colnames(values))) stop(paste(c("values should give a value for each of the following covariates: ","\n",vars),collapse=" "))

   ### pour les facteurs
   for(v in colnames(values))
   {
       if(v %in% names(x$levels$levelsdata))
       {
           if(!is.null(x$levels$levelsdata[[v]]))
           {
               values[,v] <- factor(values[,v], levels=x$levels$levelsdata[[v]])
           }
       }
       if(v %in% names(x$levels$levelsrandom))
       {
           if(!is.null(x$levels$levelsrandom[[v]]))
           {
               values[,v] <- factor(values[,v], levels=x$levels$levelsrandom[[v]])
               if(any(is.na(values[,v]))) stop(paste("Wrong factor level in variable",v))
           }
       }
   }
   
   call_random <- x$call$random[2]
   call_random <- gsub("factor","",call_random)

   if(nrow(values)>1) warning("only the first line of values is used")

   values1 <- values[1,,drop=FALSE]
   var.random <- model.matrix(formula(paste("~",call_random,sep="")),data=values1)
   var.cor <- values1[,name.cor]
   if(!is.null(name.cor)) var.cor <- model.matrix(formula(paste("~-1+",name.cor,sep="")),data=values1)

   nea <- sum(x$idea0==1)
   VarU <- matrix(0,nea,nea)
   if(nea==x$N[3])
   {
    diag(VarU) <- x$best[x$N[1]+x$N[2]+1:x$N[3]]
   }
   else
   {
    VarU[lower.tri(VarU,diag=TRUE)] <- x$best[x$N[1]+x$N[2]+1:x$N[3]]
    VarU <- t(VarU)
    VarU[lower.tri(VarU,diag=TRUE)] <- x$best[x$N[1]+x$N[2]+1:x$N[3]]
   }

   numer <- var.random %*% VarU %*% t(var.random)
   if(x$ng>1)
   {
    nw <- rep(1,x$ng)
    if(x$N[4]>0) nw <- c((x$best[x$N[1]+x$N[2]+x$N[3]+1:x$N[4]])^2,1)
    numer <- rep(numer,x$ng) * nw
   }

   Corr <- 0
   if(length(x$N)>4 & x$N[5]>0)
   {
    if(x$N[5]==1)
    {
     Corr <- (x$best[x$N[1]+x$N[2]+x$N[3]+x$N[4]+x$N[5]+1])^2 * var.cor
    }
    if(x$N[5]==2)
    {
     Corr <- (x$best[x$N[1]+x$N[2]+x$N[3]+x$N[4]+x$N[5]+2])^2
    }
   }
   denom <- numer + Corr + (x$best[length(x$best)])^2

   res[1,] <- as.numeric(numer/denom *100)
  }
  
  rownames(res) <- "%Var"
  colnames(res) <- paste("class",1:x$ng,sep="")
 }
 else
 {
  cat("Output can not be produced since the program stopped abnormally. \n")
  res <- NA
 }

 
 return(res)
}



#' Percentage of variance explained by the (latent class) linear mixed model
#' regression
#' 
#' The function provides the percentage of variance explained by the (latent
#' class) linear mixed regression in a model estimated with \code{hlme},
#' \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}.
#' 
#' 
#' @aliases VarExpl VarExpl.hlme VarExpl.lcmm VarExpl.Jointlcmm
#' VarExpl.multlcmm
#' @param x an object of class \code{hlme}, \code{lcmm}, \code{multlcmm} or
#' \code{Jointlcmm}
#' @param values a data frame with a unique row that contains the values of the
#' variables in random and the time variable in the correlation process from
#' which the percentage of variance should be calculated.
#' @return For \code{hlme}, \code{lcmm}, and \code{Jointlcmm} objects, the
#' function returns a matrix with 1 row and ng (ie the number of latent
#' classes) columns containing (the class specific) percentages of variance
#' explained by the linear mixed regression.
#' 
#' For \code{multlcmm} objects, the function returns a matrix containing (the
#' class specific) percentages of variance explained by the linear mixed
#' regression for each outcome. The resulting matrix is composed of as many
#' rows as outcomes and as many columns as latent classes.
#' @author Cecile Proust-Lima, Viviane Philipps
#' @seealso \code{\link{hlme}}, \code{\link{lcmm}}, \code{\link{multlcmm}},
#' \code{\link{Jointlcmm}}
#' @examples
#' 
#' \dontrun{
#' m1 <- multlcmm(Ydep1+Ydep2~1+Time*X2+contrast(X2),random=~1+Time,
#' subject="ID",randomY=TRUE,link=c("4-manual-splines","3-manual-splines"),
#' intnodes=c(8,12,25),data=data_lcmm, 
#' B=c(-1.071, -0.192,  0.106, -0.005, -0.193,  1.012,  0.870,  0.881,
#'   0.000,  0.000, -7.520,  1.401,  1.607 , 1.908,  1.431,  1.082,
#'  -7.528,  1.135 , 1.454 , 2.328, 1.052))
#' 
#' # variation percentages explained by linear mixed regression
#' VarExpl(m1,data.frame(Time=0))
#' }
#' 
#' @export
#' 
VarExpl <- function(x,values) UseMethod("VarExpl")


