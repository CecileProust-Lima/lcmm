

predictY.lcmm <- function(x,newdata,var.time,methInteg=0,nsim=20,draws=FALSE,ndraws=2000,na.action=1,...){


if(missing(newdata)) stop("The argument newdata should be specified")
if(missing(x)) stop("The argument x should be specified")
if (!inherits(x, "lcmm")) stop("use only with \"lcmm\" objects")
# ad 2/04/2012 Xnames2
if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) {                               
cat("newdata should at least include the following covariates: ", "\n")
cat(x$Xnames2[-1], "\n")}
if (!all(x$Xnames2 %in% c(colnames(newdata),"intercept"))) stop("see above")
if (!inherits(newdata, "data.frame")) stop("newdata should be a data.frame object")
if (!(methInteg %in% c(0,1))) stop("The integration method must be either 0 for Gauss-Hermite or 1 for Monte-Carlo")
if ((methInteg==0)&(!(nsim %in% c(5,7,9,15,20,30,40,50)))) stop("For Gauss-Hermite integration method, 'nsim' should be either 5,7,9,15,20,30,40 or 50")
#if(missing(var.time)) stop("missing argument 'var.time'")
#if(!(var.time %in% colnames(newdata))) stop("'var.time' should be included in newdata")


call_fixed <- x$call$fixed[3]
if(is.null(x$call$random)) {call_random <- ~-1} else call_random <- x$call$random
if(is.null(x$call$classmb)) {call_classmb <- ~-1} else call_classmb <- x$call$classmb
if(is.null(x$call$mixture)) {call_mixture <- ~-1} else call_mixture <- x$call$mixture


if(x$conv==1|x$conv==2|x$conv==3) {

#------------> changement Cecile 10/04/2012
## add 12/04/2012
if(x$Xnames2[1]!="intercept"){
	newdata1 <- newdata[,x$Xnames2]
	colnames(newdata1) <- x$Xnames
	newdata1 <- data.frame(newdata1)
}else{
	newdata1 <- cbind(rep(1,length=length(newdata[,1])),newdata[,x$Xnames2[-1]])
	colnames(newdata1) <- c("intercept",x$Xnames2[-1])
	newdata1 <- data.frame(newdata1)
}


            if(x$conv==2 & draws==TRUE)
                {
                    cat("No confidence interval will be provided since the program did not converge properly \n")
                    draws <- FALSE
                }


            if(x$conv==3 & draws==TRUE)
                {
                    cat("No confidence interval will be provided since the program did not converge properly \n")
                    draws <- FALSE
                }



X1 <- NULL                                                              
X2 <- NULL
b1 <- NULL
b2 <- NULL


if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

if(na.action==1){
	na.action=na.omit
}else{
	na.action=na.fail
}

### pour les facteurs

 #cas ou une variable du dataset est un facteur
 olddata <- eval(x$call$data)
  for(v in x$Xnames2[-1])
 {
  if (is.factor(olddata[,v]))
  {
   mod <- levels(olddata[,v])
   if (!(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata1[,v] <- factor(newdata1[,v], levels=mod)
  }
 }
 
 #cas ou on a factor() dans l'appel
 z <- all.names(call_fixed)
 ind_factor <- which(z=="factor")
 if(length(ind_factor))
 {
  nom.factor <- z[ind_factor+1]  
  for (v in nom.factor)
  {
   mod <- levels(as.factor(olddata[,v]))
   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata1[,v] <- factor(newdata1[,v], levels=mod)
  }
 }
 call_fixed <- gsub("factor","",call_fixed)

 z <- all.names(call_random)
 ind_factor <- which(z=="factor")
 if(length(ind_factor))
 {
  nom.factor <- z[ind_factor+1]
  for (v in nom.factor)
  {
   mod <- levels(as.factor(olddata[,v]))
   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata1[,v] <- factor(newdata1[,v], levels=mod)
  }
 }
 call_random <- gsub("factor","",call_random)
       
 z <- all.names(call_classmb)
 ind_factor <- which(z=="factor")
 if(length(ind_factor))
 {
  nom.factor <- z[ind_factor+1]
  for (v in nom.factor)
  {
   mod <- levels(as.factor(olddata[,v]))
   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata1[,v] <- factor(newdata1[,v], levels=mod)
  }
 }
 call_classmb <- gsub("factor","",call_classmb)
        
 z <- all.names(call_mixture)
 ind_factor <- which(z=="factor")
 if(length(ind_factor))
 {
  nom.factor <- z[ind_factor+1]
  for (v in nom.factor)
  {
   mod <- levels(as.factor(olddata[,v]))
   if (!all(levels(as.factor(newdata1[,v])) %in% mod)) stop(paste("invalid level in factor", v))
   newdata1[,v] <- factor(newdata1[,v], levels=mod)
  }
 }
 call_mixture <- gsub("factor","",call_mixture)   

call_mixture <- formula(call_mixture)
call_random <- formula(call_random)
call_classmb <- formula(call_classmb)

### Traitement des donnees manquantes

# permet de conserver que data=... dans lcmm ; mcall= objet de type call
#mcall <- x$call[c(1,match(c("data"),names(x$call),0))]
mcall <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
mcall$na.action <- na.action
mcall$data <- newdata1

# fixed
m <- mcall
m$formula <- formula(paste("~",call_fixed,sep=""))
m[[1]] <- as.name("model.frame")	
m <- eval(m, sys.parent()) 
na.fixed <- attr(m,"na.action")

# mixture
if((length(attr(terms(call_mixture),"term.labels"))+attr(terms(call_mixture),"intercept"))>0){
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

# random
if((length(attr(terms(call_random),"term.labels"))+attr(terms(call_random),"intercept"))>0){
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
# classmb
if((length(attr(terms(call_classmb),"term.labels"))+attr(terms(call_classmb),"intercept"))>0){ 
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
# cor
na.cor <- NULL
if(length(x$N)>5)
{
 if(x$N[6]>0)
 {
  z <- which(x$idcor0==1)
  var.cor <- newdata1[,x$Xnames[z]]
  na.cor <- which(is.na(var.cor))
 }
}

#var.time
if(!missing( var.time))
    {
        if(!(var.time %in% colnames(newdata))) stop("'var.time' should be included in newdata")
        if(var.time %in% colnames(newdata1))
            {
                times <- newdata1[,var.time,drop=FALSE]
            }
        else
            {
                times <- newdata[,var.time,drop=FALSE]
            }
    }
else
    {
        times <- newdata[,1,drop=FALSE]
    }

## Table sans donnees manquante: newdata
na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.cor))
if(length(na.action)){
	newdata1 <- newdata1[-na.action,]
        times <- times[-na.action]
}






## Construction de nouvelles var eplicatives sur la nouvelle table
## fixed
	
	X_fixed <- model.matrix(formula(paste("~",call_fixed,sep="")),data=newdata1)
	if(colnames(X_fixed)[1]=="(Intercept)"){
		colnames(X_fixed)[1] <- "intercept"
		int.fixed <- 1
	}
## mixture
	if(id.X_mixture == 1){
            X_mixture <- model.matrix(call_mixture,data=newdata1)	
            if(colnames(X_mixture)[1]=="(Intercept)"){
                colnames(X_mixture)[1] <- "intercept"
                int.mixture <- 1
            }
	}	
## random
	if(id.X_random == 1){
            X_random <- model.matrix(call_random,data=newdata1)	
            if(colnames(X_random)[1]=="(Intercept)"){
                colnames(X_random)[1] <- "intercept"
                int.random <- 1
            }
	}	
## classmb
	if(id.X_classmb == 1){ 
            X_classmb <- model.matrix(call_classmb,data=newdata1)
            colnames(X_classmb)[1] <- "intercept"
	}
##cor	
if(length(x$N)>5)
{
 if(x$N[6]>0)  #on reprend la variable de temps de cor
 {
  z <- which(x$idcor0==1)
  var.cor <- newdata1[,x$Xnames[z]]
 }
}

## Construction des var expli
newdata1 <- X_fixed
colX <- colnames(X_fixed)

if(id.X_mixture == 1){
	for(i in 1:length(colnames(X_mixture))){
		if((colnames(X_mixture)[i] %in% colnames(newdata1))==FALSE){
			newdata1 <- cbind(newdata1,X_mixture[,i])
                        colnames(newdata1) <- c(colX,colnames(X_mixture)[i])
                        colX <- colnames(newdata1)
		}
	}
}
if(id.X_random == 1){
	for(i in 1:length(colnames(X_random))){
		if((colnames(X_random)[i] %in% colnames(newdata1))==FALSE){
			newdata1 <- cbind(newdata1,X_random[,i])
                        colnames(newdata1) <- c(colX,colnames(X_random)[i])
                        colX <- colnames(newdata1)
		}	 
	}
}
if(id.X_classmb == 1){
	for(i in 1:length(colnames(X_classmb))){
		if((colnames(X_classmb)[i] %in% colnames(newdata1))==FALSE){
			newdata1 <- cbind(newdata1,X_classmb[,i])
                        colnames(newdata1) <- c(colX,colnames(X_classmb)[i])
                        colX <- colnames(newdata1)
		}	
	}
}
if(length(x$N)>5)
{
 if(x$N[6]>0)
 { 
  if( x$idg0[z]==0 & x$idea0[z]==0 & x$idprob0[z]==0)
      {
          newdata1 <- cbind(newdata1,var.cor)
          colnames(newdata1) <- c(colX,x$Xnames[z])
          colX <- colnames(newdata1)
      }
 }
}




nv <- length(x$idg0)
maxmes <- length(newdata1[,1])
npm <- length(x$best)
best <- x$best
if(x$idiag==0 & x$N[3]>0) best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])] <- x$cholesky
if(x$idiag==1 & x$N[3]>0) best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])] <- sqrt(best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])])
nwg <- x$N[4]
ncor <- 0
if (x$linktype!=3)
{ 
 if(length(x$N)>5) {ncor <- x$N[6]}
}

 

### for linear trajectory 
                                                                                    


if (x$linktype==0){

if (!draws) {

# prediction
X1 <- NULL
X2 <- NULL
b1 <- NULL
b2 <- NULL

kk<-0
for(k in 1:length(x$idg0)){
if(x$idg0[k]==1){
X1 <- cbind(X1,newdata1[,k])
if (k==1) b1 <- c(b1,0)
if (k>1) {
place <- x$N[1]+kk
b1 <- c(b1,x$best[place+1])
kk <- kk+1
}
}

if(x$idg0[k]==2){
X2 <- cbind(X2,newdata1[,k])
if (k==1){
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng-1
b2 <- rbind(b2,c(0,x$best[place1:place2]))
kk <- kk+x$ng-1
}
if (k>1){
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng
b2 <- rbind(b2,x$best[place1:place2])
kk <- kk+x$ng
}
}
}

Ypred<-matrix(0,length(newdata1[,1]),x$ng)
colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
if (x$ng==1) colnames(Ypred) <- "Ypred"

for(g in 1:x$ng){
if(length(b1) != 0){
Ypred[,g]<- X1 %*% b1 
}
if(length(b2) != 0){
Ypred[,g]<- Ypred[,g] + X2 %*% b2[,g]
}

Ypred[,g] <- Ypred[,g]*abs(x$best[(npm-ncor)])+x$best[(npm-1-ncor)]
}
}

if (draws){


ndraws <- as.integer(ndraws)
ydraws <- NULL

            posfix <- eval(x$call$posfix)

            if(ndraws>0)
                {
                    Mat <- matrix(0,ncol=npm,nrow=npm)
                    Mat[upper.tri(Mat,diag=TRUE)]<- x$V
                    if(length(posfix))
                        {
                            Mat2 <- Mat[-posfix,-posfix]
                            Chol2 <- chol(Mat2)
                            Chol <- matrix(0,npm,npm)
                            Chol[setdiff(1:npm,posfix),setdiff(1:npm,posfix)] <- Chol2
                            Chol <- t(Chol)
                        }
                    else
                        {
                            Chol <- chol(Mat)
                            Chol <- t(Chol)
                        }
                }



for (j in 1:ndraws) {

bdraw <- rnorm(npm)
bdraw <- best + Chol %*% bdraw


# prediction
X1 <- NULL
X2 <- NULL
b1 <- NULL
b2 <- NULL


kk<-0
for(k in 1:length(x$idg0)){
if(x$idg0[k]==1){
X1 <- cbind(X1,newdata1[,k])
if (k==1) b1 <- c(b1,0)
if (k>1) {
place <- x$N[1]+kk
b1 <- c(b1,bdraw[place+1])
kk <- kk+1
}
}

if(x$idg0[k]==2){
X2 <- cbind(X2,newdata1[,k])
if (k==1){
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng-1
b2 <- rbind(b2,c(0,bdraw[place1:place2]))
kk <- kk+x$ng-1
}
if (k>1){
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng
b2 <- rbind(b2,bdraw[place1:place2])
kk <- kk+x$ng
}
}
}

Ypred<-matrix(0,length(newdata1[,1]),x$ng)
colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
if (x$ng==1) colnames(Ypred) <- "Ypred"

for(g in 1:x$ng){
if(length(b1) != 0){
Ypred[,g]<- X1 %*% b1 
}
if(length(b2) != 0){
Ypred[,g]<- Ypred[,g] + X2 %*% b2[,g]
}
Ypred[,g] <- Ypred[,g]*abs(bdraw[(npm-ncor)])+bdraw[(npm-1-ncor)]
}
pred <- as.vector(Ypred)
ydraws <- cbind(ydraws,pred)
}

f <- function(x) {
quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
}
ydistr <- apply(ydraws,1,FUN=f)
Ypred_50 <- matrix(ydistr[2,],ncol=x$ng,byrow=F)
Ypred_2.5 <- matrix(ydistr[1,],ncol=x$ng,byrow=F)
Ypred_97.5 <- matrix(ydistr[3,],ncol=x$ng,byrow=F)

Ypred <- cbind(Ypred_50,Ypred_2.5,Ypred_97.5)


if (x$ng==1){
colnames(Ypred) <- c("Ypred_50","Ypred_2.5","Ypred_97.5")
}
if (x$ng>1){
colnames(Ypred) <- c(paste("Ypred_50_class",1:x$ng,sep=""),paste("Ypred_2.5_class",1:x$ng,sep=""),paste("Ypred_97.5_class",1:x$ng,sep=""))
}




}
}



### for threshold trajectory

if (x$linktype==3){


if(!draws) {

# debut prediction
X1 <- NULL
X2 <- NULL
X3 <- NULL
b1 <- NULL
b2 <- NULL

kk<-0
for(k in 1:length(x$idg0)){
if(x$idg0[k]==1){
X1 <- cbind(X1,newdata1[,k])
if (k==1) b1 <- c(b1,0)
if (k>1) {
place <- x$N[1]+kk
b1 <- c(b1,x$best[place+1])
kk <- kk+1
}
}

if(x$idg0[k]==2){
X2 <- cbind(X2,newdata1[,k])
if (k==1){
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng-1
b2 <- rbind(b2,c(0,x$best[place1:place2]))
kk <- kk+x$ng-1
}
if (k>1){
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng
b2 <- rbind(b2,x$best[place1:place2])
kk <- kk+x$ng
}
}

if(x$idea0[k]==1){
X3 <- cbind(X3,newdata1[,k])
}
}

Ypred<-matrix(0,length(newdata1[,1]),x$ng)
colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
if (x$ng==1) colnames(Ypred) <- "Ypred"

for(g in 1:x$ng){
if(length(b1) != 0){
Ypred[,g]<- X1 %*% b1 
}
if(length(b2) != 0){
Ypred[,g]<- Ypred[,g] + X2 %*% b2[,g]
}
}


varpred <- 0
nea <- sum(x$idea0)
if(nea!=0){
if(nea==x$N[3]){
varpred <- X3 %*%  best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])]^2
}
if(nea!=x$N[3]){
U <- matrix(0,nrow=nea,ncol=nea)
U[upper.tri(U,diag=TRUE)] <- best[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])]
varpred <-  X3 %*% t(U)
varpred <- varpred %*% t(varpred)
}
if(nea>1) varpred <- diag(varpred)
}

wg <- rep(1,x$ng)
if(x$N[4]!=0&x$ng>1){
wg[1:(x$ng-1)] <- best[(x$N[1]+x$N[2]+x$N[3]+1):(x$N[1]+x$N[2]+x$N[3]+x$N[4])]
wg <- wg^2
}

ntrtot0 <- sum(x$ide==1) 
seuils <- x$ide
Nseuils <- length(x$ide)
seuils[x$ide==1] <- as.vector(best[(npm-ntrtot0+1):npm])   #ncor=0 donc ok
seuils[x$ide==0] <- 0
if (Nseuils>=2){
cumseuils <- cumsum(seuils[2:Nseuils]*seuils[2:Nseuils])
seuils[2:Nseuils] <- rep(seuils[1],(Nseuils-1))+cumseuils
}

pred <- Ypred
for(g in 1:x$ng){
Ypred[,g] <- rep(x$linknodes[2],maxmes)
for(i in 1:Nseuils)
Ypred[,g] <- Ypred[,g] - pnorm((seuils[i]-pred[,g])/sqrt(wg[g]*varpred+1))
}

if (x$ng>1) colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
if (x$ng==1) colnames(Ypred) <- "Ypred"
}

if (draws) {

ndraws <- as.integer(ndraws)
ydraws <- NULL

Mat <- matrix(0,ncol=npm,nrow=npm)
# que la partie sup utilisee donc OK si rien en bas
Mat[upper.tri(Mat,diag=TRUE)]<- x$V
Chol <- chol(Mat)
Chol <- t(Chol)

for (j in 1:ndraws) {

bdraw <- rnorm(npm)
bdraw <- best + Chol %*% bdraw

# debut prediction
X1 <- NULL
X2 <- NULL
X3 <- NULL
b1 <- NULL
b2 <- NULL

kk<-0
for(k in 1:length(x$idg0)){
if(x$idg0[k]==1){
X1 <- cbind(X1,newdata1[,k])
if (k==1) b1 <- c(b1,0)
if (k>1) {
place <- x$N[1]+kk
b1 <- c(b1,bdraw[place+1])
kk <- kk+1
}
}

if(x$idg0[k]==2){
X2 <- cbind(X2,newdata1[,k])
if (k==1){
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng-1
b2 <- rbind(b2,c(0,bdraw[place1:place2]))
kk <- kk+x$ng-1
}
if (k>1){
place1 <- x$N[1]+kk+1
place2 <- x$N[1]+kk+x$ng
b2 <- rbind(b2,bdraw[place1:place2])
kk <- kk+x$ng}
}

if(x$idea0[k]==1){
X3 <- cbind(X3,newdata1[,k])
}
}

Ypred<-matrix(0,length(newdata1[,1]),x$ng)
colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
if (x$ng==1) colnames(Ypred) <- "Ypred"

for(g in 1:x$ng){
if(length(b1) != 0){
Ypred[,g]<- X1 %*% b1 
}
if(length(b2) != 0){
Ypred[,g]<- Ypred[,g] + X2 %*% b2[,g]
}
}


varpred <- 0
nea <- sum(x$idea0)
if(nea!=0){
if(nea==x$N[3]){
varpred <- X3 %*%  bdraw[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])]^2 
}
if(nea!=x$N[3]){
U <- matrix(0,nrow=nea,ncol=nea)
U[upper.tri(U,diag=TRUE)] <- bdraw[(x$N[1]+x$N[2]+1):(x$N[1]+x$N[2]+x$N[3])]
varpred <- X3 %*% t(U) 
varpred <- varpred %*% t(varpred)
}
if(nea>1) varpred <- diag(varpred)
}

wg <- rep(1,x$ng)
if(x$N[4]!=0&x$ng>1){
wg[1:(x$ng-1)] <- bdraw[(x$N[1]+x$N[2]+x$N[3]+1):(x$N[1]+x$N[2]+x$N[3]+x$N[4])]
wg <- wg^2
}


ntrtot0 <- sum(x$ide==1) 
seuils <- x$ide
Nseuils <- length(x$ide)
seuils[x$ide==1] <- as.vector(bdraw[(npm-ntrtot0+1):npm])
seuils[x$ide==0] <- 0
if (Nseuils>=2){
cumseuils <- cumsum(seuils[2:Nseuils]*seuils[2:Nseuils])
seuils[2:Nseuils] <- rep(seuils[1],(Nseuils-1))+cumseuils
}
pred <- Ypred
for(g in 1:x$ng){
Ypred[,g] <- rep(x$linknodes[2],maxmes)
for(i in 1:Nseuils)
Ypred[,g] <- Ypred[,g] - pnorm((seuils[i]-pred[,g])/sqrt(wg[g]*varpred+1))
}
pred <- as.vector(Ypred)
ydraws <- cbind(ydraws,pred)
}

f <- function(x) {
quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
}
ydistr <- apply(ydraws,1,FUN=f)
Ypred_50 <- matrix(ydistr[2,],ncol=x$ng,byrow=F)
Ypred_2.5 <- matrix(ydistr[1,],ncol=x$ng,byrow=F)
Ypred_97.5 <- matrix(ydistr[3,],ncol=x$ng,byrow=F)

Ypred <- cbind(Ypred_50,Ypred_2.5,Ypred_97.5)

if (x$ng==1){
colnames(Ypred) <- c("Ypred_50","Ypred_2.5","Ypred_97.5")
}
if (x$ng>1){
colnames(Ypred) <- c(paste("Ypred_50_class",1:x$ng,sep=""),paste("Ypred_2.5_class",1:x$ng,sep=""),paste("Ypred_97.5_class",1:x$ng,sep=""))
}
}
}



### for splines or beta trajectory


if (x$linktype %in% c(1,2)){


nbzitr <- length(x$linknodes)
epsY <- x$epsY
Ymarg <- rep(0,maxmes*x$ng)

#cat(c(nv,x$ng,nbzitr,epsY,nwg,nsim,methInteg,x$Ydiscrete),"\n")
#cat(epsY,"\n")

if (!draws){      
    out <- .Fortran(C_predictcont,
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
                    as.double(epsY),
                    as.integer(x$linktype),
                    as.integer(nbzitr),
                    as.double(x$linknodes),
                    as.integer(nsim),
                    as.integer(methInteg),
                    as.integer(x$Ydiscrete),
                    Ymarg=as.double(Ymarg))
   
out$Ymarg[out$Ymarg==9999] <- NA

#cat(out$Ymarg)
Ypred <- matrix(out$Ymarg,ncol=x$ng,byrow=FALSE)

if (x$ng==1)colnames(Ypred) <- "Ypred"
if (x$ng>1)colnames(Ypred) <- paste("Ypred_class",1:x$ng,sep="")
}

########### ajout ndraws ###############################


if (draws) {

ndraws <- as.integer(ndraws)
ydraws <- NULL

Mat <- matrix(0,ncol=npm,nrow=npm)
# que la partie sup utilisee donc OK si rien en bas
Mat[upper.tri(Mat,diag=TRUE)]<- x$V
Chol <- chol(Mat)
Chol <- t(Chol)

for (j in 1:ndraws) {

bdraw <- rnorm(npm)
bdraw <- best + Chol %*% bdraw
 
out <- .Fortran(C_predictcont,
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
                as.double(epsY),
                as.integer(x$linktype),
                as.integer(nbzitr),
                as.double(x$linknodes),
                as.integer(nsim),
                as.integer(methInteg),
                as.integer(x$Ydiscrete),
                Ymarg=as.double(Ymarg))

out$Ymarg[out$Ymarg==9999] <- NA
ydraws <- cbind(ydraws,out$Ymarg)
}

f <- function(x) {
quantile(x[!is.na(x)],probs=c(0.025,0.5,0.975))
}
ydistr <- apply(ydraws,1,FUN=f)
Ypred_50 <- matrix(ydistr[2,],ncol=x$ng,byrow=FALSE)
Ypred_2.5 <- matrix(ydistr[1,],ncol=x$ng,byrow=FALSE)
Ypred_97.5 <- matrix(ydistr[3,],ncol=x$ng,byrow=FALSE)

Ypred <- cbind(Ypred_50,Ypred_2.5,Ypred_97.5)


if (x$ng==1){
colnames(Ypred) <- c("Ypred_50","Ypred_2.5","Ypred_97.5")
}
if (x$ng>1){
colnames(Ypred) <- c(paste("Ypred_50_class",1:x$ng,sep=""),paste("Ypred_2.5_class",1:x$ng,sep=""),paste("Ypred_97.5_class",1:x$ng,sep=""))
}
}
}


res.list <- NULL
res.list$pred <- Ypred
res.list$times <- times
}

else  #cas xconv != 1 ou 2
{ 
 cat("Predictions can not be computed since the program stopped abnormally. \n")
 res.list <- list(pred=NA,times=NA)
}

class(res.list) <- "predictY"
return(res.list)
}


predictY <- function(x,newdata,var.time,...) UseMethod("predictY")
