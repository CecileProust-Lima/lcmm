
############################### Fonction Jointlcmm ###################################

Diffepoce <- function(epoceM1,epoceM2){

cl <- match.call()
if(class(epoceM1)=="Jointlcmm") stop("Diffepoce allows only arguments of classes 'epoce'. Please run function epoce on the Jointlcmm objects before running Diffepoce.")
if(class(epoceM2)=="Jointlcmm") stop("Diffepoce allows only arguments of classes 'epoce'. Please run function epoce on the Jointlcmm objects before running Diffepoce.")
if(class(epoceM1)!="epoce") stop("Diffepoce allows only arguments of classes 'epoce'.")
if(class(epoceM2)!="epoce") stop("Diffepoce allows only arguments of classes 'epoce'.")

if(!(all.equal(epoceM1$EPOCE[,1],epoceM2$EPOCE[,1]))) stop("The two epoce objects should have the same prediction times")

if(!(all.equal(epoceM1$IndivContrib[,1],epoceM2$IndivContrib[,1]))) stop("The two epoce objects should be derived from the same dataset")


epoceM1$IndivContrib[is.na(epoceM1$IndivContrib)] <- 0
epoceM2$IndivContrib[is.na(epoceM2$IndivContrib)] <- 0

DEPOCE <- as.double(epoceM1$EPOCE[,5])-as.double(epoceM2$EPOCE[,5])
if (epoceM1$new.data==TRUE){
DEPOCE <- as.double(epoceM1$EPOCE[,4])-as.double(epoceM2$EPOCE[,4])
}

M1Contrib <- epoceM1$IndivContrib[,-1]
M2Contrib <- epoceM2$IndivContrib[,-1]	



diff <- M1Contrib - M2Contrib
diff_sq <- (M1Contrib - M2Contrib)**2


DMPOL <- apply(diff,FUN=sum,MARGIN=2)/apply(M1Contrib[,]!=0,FUN=sum,MARGIN=2)


omega <- (apply(diff_sq,FUN=sum,MARGIN=2)/apply(M1Contrib[,]!=0,FUN=sum,MARGIN=2)) - (DMPOL)**2
omega <- omega/apply(M1Contrib[,]!=0,FUN=sum,MARGIN=2) 

TIinf <- as.vector(DEPOCE) - qnorm(0.975)*sqrt(omega)
TIsup <- as.vector(DEPOCE) + qnorm(0.975)*sqrt(omega)


DiffEPOCE <- cbind(epoceM1$EPOCE[,1],DEPOCE,TIinf,TIsup)
rownames(DiffEPOCE) <- rep(" ",length(DiffEPOCE[,1]))

if (epoceM1$new.data==TRUE){
colnames(DiffEPOCE) <- c("pred. times","Diff MPOL","95%TIinf","95%TIsup")
}
else{
colnames(DiffEPOCE) <- c("pred. times","Diff CVPOL","95%TIinf","95%TIsup")
}

res <- list(call.JLCM1=epoceM1$call.Jointlcmm,call.JLCM2=epoceM2$call.Jointlcmm,call=cl,DiffEPOCE=DiffEPOCE,new.data=epoceM1$new.data)



class(res) <-c("Diffepoce")
res
}
