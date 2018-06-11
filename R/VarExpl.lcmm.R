#' @export
#'
VarExpl.lcmm <- function(x,values)
{
 if(missing(x)) stop("The model should be specified")
 if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")
 if(missing(values)) values <- data.frame("intercept"=1)
 if (!inherits(values, "data.frame")) stop("values should be a data.frame object")
 if(any(is.na(values))) stop("values should not contain any missing values")


 if(x$conv==1 | x$conv==2 | x$conv==3)
 {
   
  
   ntrtot <- length(x$best)- (x$N[1]+x$N[2]+x$N[3]+x$N[4]+x$N[6])
   
   
   
   
   
  res <- matrix(0,nrow=1,ncol=x$ng)

  names.random <- NULL
  name.cor <- NULL
  if(x$N[3]>0) names.random <- x$Xnames[which(x$idea0==1)]
  if(length(x$N)>5 & x$N[6]>0) name.cor <- x$Xnames[which(x$idcor0==1)]

  if(!is.null(names.random) | !is.null(name.cor))
  {
   names.values <- unique(c(names.random,name.cor))   #contient I(T^2))

   vars <- unique(c(all.vars(x$call$random),all.vars(x$call$cor)))
   if(!all(vars %in% colnames(values))) stop(paste(c("values should give a value for each of the following covariates: ","\n",vars),collapse=" "))

   ### pour les facteurs
   ##donnees de l estimation
   if(!is.null(x$data))
   {
       olddata <- x$data
   }
   else
      {
          olddata <- eval(x$call$data)
      }

   #cas ou une variable du dataset est un facteur
   for(v in setdiff(vars,"intercept"))
   {
    if(is.factor(olddata[,v]))
    {
     mod <- levels(olddata[,v])
     if (!(levels(as.factor(values[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     values[,v] <- factor(values[,v], levels=mod)
    }
   }

   #cas ou on a factor() dans l'appel
   call_random <- x$call$random
   z <- all.names(call_random)
   ind_factor <- which(z=="factor")
   if(length(ind_factor))
   {
    nom.factor <- z[ind_factor+1]
    for (v in nom.factor)
    {
     mod <- levels(as.factor(olddata[,v]))
     if (!all(levels(as.factor(values[,v])) %in% mod)) stop(paste("invalid level in factor", v))
     values[,v] <- factor(values[,v], levels=mod)
    }
   }
   call_random <- gsub("factor","",call_random)

   if (!is.null(name.cor)) values1 <- model.matrix(formula(paste("~",paste(call_random[2],name.cor,sep="+"))),data=values)
   else values1 <- model.matrix(formula(call_random),data=values)

   if(colnames(values1)[1]=="(Intercept)") colnames(values1)[1] <- "intercept"

   if(nrow(values1)>1) warning("only the first line of values is used")
   var.random <- values1[1,names.random]
   var.cor <- values1[1,name.cor]

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

   numer <- t(var.random) %*% VarU %*% var.random
   if(x$ng>0)
   {
    nw <- rep(1,x$ng)
    if(x$N[4]>0) nw <- c((x$best[x$N[1]+x$N[2]+x$N[3]+1:x$N[4]])^2,1)
    numer <- numer * nw
   }

   Corr <- 0
   if(length(x$N)>5 & x$N[6]>0)
   {
    if(x$N[6]==1)
    {
     Corr <- (x$best[x$N[1]+x$N[2]+x$N[3]+x$N[4]+ntrtot+1])^2 * var.cor
    }
    if(x$N[6]==2)
    {
     Corr <- (x$best[x$N[1]+x$N[2]+x$N[3]+x$N[4]+ntrtot+2])^2
    }
   }
  # numer <- numer + Corr
   denom <- numer + 1

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
