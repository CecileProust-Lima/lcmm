
#' Difference of expected prognostic cross-entropy (EPOCE) estimators and its
#' 95\% tracking interval between two joint latent class models estimated with
#' \code{Jointlcmm}
#' 
#' This function computes the difference of 2 EPOCE estimates (CVPOL or MPOL)
#' and its 95\% tracking interval between two joint latent class models
#' estimated using \code{Jointlcmm} and evaluated using \code{epoce} function.
#' Difference in CVPOL is computed when the EPOCE was previously estimated on
#' the same dataset as used for estimation (using an approximated
#' cross-validation), and difference in MPOL is computed when the EPOCE was
#' previously estimated on an external dataset.
#' 
#' This function does not apply for the moment with multiple causes of event
#' (competing risks).
#' 
#' From the EPOCE estimates and the individual contributions to the prognostic
#' observed log-likelihood obtained with \code{epoce} function on the same
#' dataset from two different estimated joint latent class models, the
#' difference of CVPOL (or MPOL) and its 95\% tracking interval is computed.
#' The 95\% tracking interval is:
#' 
#' Delta(MPOL) +/- qnorm(0.975)*sqrt(VARIANCE) for an external dataset
#' 
#' Delta(CVPOL) +/- qnorm(0.975)*sqrt(VARIANCE) for the dataset used in
#' \code{Jointlcmm}
#' 
#' where Delta(CVPOL) (or Delta(MPOL)) is the difference of CVPOL (or MPOL) of
#' the two joint latent class models, and VARIANCE is the empirical variance of
#' the difference of individual contributions to the prognostic observed
#' log-likelihoods of the two joint latent class models.
#' 
#' See Commenges et al. (2012) and Proust-Lima et al. (2012) for further
#' details.
#' 
#' @param epoceM1 a first object inheriting from class \code{epoce}
#' @param epoceM2 a second object inheriting from class \code{epoce}
#' @return \item{call.Jointlcmm1}{the \code{Jointlcmm} call for epoceM1 }
#' \item{call.Jointlcmm2}{the \code{Jointlcmm} call for epoceM2 }
#' \item{call}{the matched call} \item{DiffEPOCE}{Dataframe containing, for
#' each prediction time s, the difference in either MPOL or CVPOL depending on
#' the dataset used, and the 95\% tracking bands (TIinf and TIsup) }
#' \item{new.data}{a boolean for internal use only, which is FALSE if
#' computation is done on the same data as for \code{Jointlcmm} estimation, and
#' TRUE otherwise. }
#' @author Cecile Proust-Lima and Amadou Diakite
#' @seealso
#' \code{\link{Jointlcmm}}, \code{\link{epoce}}, \code{\link{summary.Diffepoce}}
#' @references Commenges, Liquet and Proust-Lima (2012). Choice of prognostic
#' estimators in joint models by estimating differences of expected conditional
#' Kullback-Leibler risks. Biometrics 68(2), 380-7.
#' 
#' Proust-Lima, Sene, Taylor, Jacqmin-Gadda (2014). Joint latent class models
#' for longitudinal and time-to-event data: a review. Statistical Methods in
#' Medical Research 23, 74-90.
#' @examples
#' 
#' \dontrun{
#' #### estimation with 2 latent classes (ng=2)
#' m2 <- Jointlcmm(fixed= Ydep1~Time*X1,random=~Time,mixture=~Time,subject='ID'
#' ,survival = Surv(Tevent,Event)~ X1+X2 ,hazard="Weibull"
#' ,hazardtype="PH",ng=2,data=data_lcmm,
#' B=c( 0.7608, -9.4974,  1.0242,  1.4331,  0.1063 , 0.6714, 10.4679, 11.3178,
#'  -2.5671, -0.5386,  1.4616, -0.0605,  0.9489,  0.1020,  0.2079,  1.5045),logscale=TRUE)
#' m1 <- Jointlcmm(fixed= Ydep1~Time*X1,random=~Time,subject='ID'
#' ,survival = Surv(Tevent,Event)~ X1+X2 ,hazard="Weibull"
#' ,hazardtype="PH",ng=1,data=data_lcmm,
#' B=c(-7.6634,  0.9136,  0.1002,  0.6641, 10.5675, -1.6589,  1.4767, -0.0806,
#'   0.9240,0.5643,  1.2277,  1.5004))
#' 
#' ## EPOCE computation for predictions times from 1 to 6 on the dataset used
#' ## for estimation of m.
#' VecTime <- c(1,3,5,7,9,11,13,15)
#' cvpol1 <- epoce(m1,var.time="Time",pred.times=VecTime)
#' cvpol1
#' cvpol2 <- epoce(m2,var.time="Time",pred.times=VecTime)
#' cvpol2
#' DeltaEPOCE <- Diffepoce(cvpol1,cvpol2)
#' summary(DeltaEPOCE)
#' plot(DeltaEPOCE,bty="l")
#' }
#' 
#' @export
#' 
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
