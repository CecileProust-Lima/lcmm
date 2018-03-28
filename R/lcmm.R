lcmm <- function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,link="linear",intnodes=NULL,epsY=0.5,cor=NULL,data,B,convB=0.0001,convL=0.0001,convG=0.0001,maxiter=100,nsim=100,prior,range=NULL,subset=NULL,na.action=1,posfix=NULL,partialH=FALSE,verbose=TRUE)
{

mm <- match.call()
if(missing(fixed)) stop("The argument Fixed must be specified in any model")
if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
if(nrow(data)==0) stop("Data should not be empty")

### transformer l'argument cor en character
cor.type <- mm$cor[1]
cor.time <- mm$cor[2]
cor.char <- paste(cor.type,cor.time,sep="-") 
if (all.equal(cor.char,character(0))!=TRUE)
{
 if (link=="thresholds") stop("The argument cor is only available with linear, beta or splines link")
}
else
{
  cor.char <- NULL
}  
  
if (!is.null(cor.char))
{
 if(!(strsplit(cor.char,"-")[[1]][2] %in% colnames(data))) stop("Unable to find time variable from argument cor in data")
 else { cor.var.time <- strsplit(cor.char,"-")[[1]][2] }
}  
### fin cor


## gestion de B=random(mod)

        if(length(mm$B)==2)
            {
                if(class(eval(mm$B[[2]]))!="lcmm") stop("The model specified in B should be of class lcmm")
                if(as.character(mm$B[1])!="random") stop("Please use random() to specify random initial values")
                
                B <- eval(mm$B[[2]])   
                B$Brandom <- TRUE
                
                if(length(posfix)) stop("Argument posfix is not compatible with random intial values")
            }





if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

if(na.action==1){
	na.action=na.omit
}else{
	na.action=na.fail
}
#cat("ide :")
#cat(ide0,"\n")
#cat(zitr,"\n")
#7/05/2012
### Traitement des donnees manquantes
# fixed
if(missing(fixed)) stop("The argument Fixed must be specified in any model")
if(class(fixed)!="formula") stop("The argument fixed must be a formula")
m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
m$formula <- terms(fixed)
m$na.action <- na.action
m[[1]] <- as.name("model.frame")	
m <- eval(m, sys.parent()) 
na.fixed <- attr(m,"na.action")

# mixture
if(!missing(mixture)){
	if(class(mixture)=="formula"){	
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
	m$formula <- terms(mixture)
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
	na.mixture <- attr(m,"na.action")
	}	
}else{
	na.mixture <- NULL
}

# random
if(!missing(random)){
	if(class(random)=="formula"){	
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
	m$formula <- terms(random)
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.random <- attr(m,"na.action")
	}
}else{
	na.random <- NULL
}

# classmb
if(!missing(classmb)){ 
	if(class(classmb)=="formula"){	
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]	
	m$formula <- terms(classmb)
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame")	
	m <- eval(m, sys.parent()) 
 	na.classmb <- attr(m,"na.action")
	}
}else{
	na.classmb <- NULL
}
 
#cor     
if(!is.null(cor.char))
{
	m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
	m$formula <- as.formula(paste(cor.var.time,1,sep="~"))
	m$na.action <- na.action
  m[[1]] <- as.name("model.frame")
  m <- eval(m,sys.parent())    
  na.cor <- attr(m,"na.action") 	
}
else { na.cor <- NULL }
 
 
	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.cor))
#7/05/2012


attr.fixed <- attributes(terms(fixed))
depvar <- as.character(attr.fixed$variables[2])
if(!isTRUE(all.equal(as.character(mm$subset),character(0))))
    {
        cc <- mm
        cc <- cc[c(1,which(names(mm)=="subset"))]
        cc[[1]] <- as.name("model.frame")
        cc$formula <- formula(paste("~",depvar))
        cc$data <- data
        cc$na.action <- na.pass
        ysubset <- eval(cc)
    }
else
    {
        ysubset <- data[,depvar,drop=FALSE]
    }

if(!is.null(na.action))
{
    Y0 <- ysubset[-na.action,]
}
else
{
    Y0 <- ysubset
}
Y0 <- unlist(Y0)

minY0 <- min(Y0)
maxY0 <- max(Y0)
if ((!missing(range)) & length(range)==2)
{
 if(minY0<range[1]|maxY0>range[2]) stop("The range specified do not cover the entire range of the data")
 if (minY0>range[1]|maxY0<range[2])
 {
  minY0 <- range[1]
  maxY0 <- range[2]
 } 
}
else
    {
        min2 <- round(minY0,3)
        if(minY0<min2) min2 <- min2-0.001
        minY0 <- min2

        max2 <- round(maxY0,3)
        if(maxY0>max2) max2 <- max2+0.001
        maxY0 <- max2
    }

if(all.equal((maxY0-minY0),0) == T){
	stop("All the values of the dependent variable are the same. No estimation can be performed in that case.")
}
#if((any(is.na(Y0))==TRUE)){
#	stop("The dependent variable should not contain any missing value")
#}

if(length(grep("-",unlist(strsplit(link,split="")))) > 2){
	stop("Please check and revise the 'link' argument according to the format given in the help.")
}



################################# pas de separateur "-", uniquement pour les splines


if(all.equal(length(grep("-",unlist(strsplit(link,split="")))),0)==T){

	if (!(link %in% c("linear","beta","thresholds","splines"))){
		stop("The only available link functions in lcmm are 'linear', 'beta', 'splines' and 'thresholds' functions.")
	}else{
		nbzitr0 <- switch(link,"linear"=2,"beta"=2,"splines"=5,"thresholds"=2)
		idlink0 <- switch(link,"linear"=0,"beta"=1,"splines"=2,"thresholds"=3)
		ntrtot0 <- switch(link,"linear"=2,"beta"=4, "splines"= (nbzitr0 + 2), "thresholds"= as.integer(maxY0-minY0)) 
		if(all.equal(link,"splines")==T | all.equal(link,"Splines")==T){
			link <- "splines"
			type <- "equi"	
		}
		
	} 	
	if ((link %in% c("thresholds"))){

		minY0 <- min(Y0)
		maxY0 <- max(Y0)
		############################## PARTIE A COMPLETER POUR ORDINAL #####################
		if(!(all.equal(minY0,as.integer(minY0))==T) | !(all.equal(maxY0,as.integer(maxY0))==T)|!all(Y0 %in% minY0:maxY0)){
			stop("With the threshold link function, the longitudinal outcome must be discrete")
		}

		IND <- sort(unique(Y0))
		IND <- IND[1:(length(IND)-1)]-minY0+1
		ide0 <- rep(0,as.integer(maxY0-minY0))
		ide0[IND] <- 1


	
	
	#######################################################################
	}
	zitr <- rep(0,nbzitr0)
	zitr[1] <- minY0
	zitr[nbzitr0] <- maxY0
	
}


################################# Avec un seul separateur "-", uniquement pour les splines
if(all.equal(length(grep("-",unlist(strsplit(link,split="")))),1)==T){
	stop("The number and location of the nodes only apply for the splines link function. For 'thresholds', 'linear' and 'beta' links, no nodes are required. For the splines link function, both the number and the type of location for the nodes should be specified (ex: 5-manual-splines for 5 manual nodes)")
}


if(all.equal(length(grep("-",unlist(strsplit(link,split="")))),2)==T){

	if(any(unlist(strsplit(link,"-")) %in% c("linear","beta","thresholds"))){
		stop("The number and location of the nodes only apply for the 'splines' link function. For 'thresholds', 'linear' and 'beta' links, no nodes are required.")
	}
	
### Verification de l'ordre de replissage du link
	if(!(unlist(strsplit(link,"-"))[3] %in% c("splines"))){
		stop("When defining a link function using splines, the third part of the 'link' argument must include only 'splines' (ex: 5-equi-splines for 5 equidistant nodes)")
	}
	
	if(!(unlist(strsplit(link,"-"))[2] %in% c("equi","manual","quant"))){
		stop("When defining a link function using splines, the second part of the 'link' argument must include only 'equi' 'manual' 'quant' for equidistant manual or quantile nodes")
	}
	
	nbzitr0 <- as.integer(unlist(strsplit(link,"-"))[1])
	ntrtot0 <- nbzitr0 + 2
	idlink0 <- 2 
	type <- unlist(strsplit(link,"-"))[2]	   
	link <- "splines" 
	if((nbzitr0-2) < 0) stop("At least 2 nodes should be specified for the splines link function.")
}


if (all.equal(idlink0,2)==T){

	zitr <- rep(0,nbzitr0)
	zitr[1] <- minY0
	zitr[nbzitr0] <- maxY0
	
	if(all.equal("manual",type)==T){
		if (is.null(intnodes)){
		stop("If 'manual' option is specified for the splines link function, intnodes argument should include the list of interior nodes")
		}else{            
		if(!(all.equal(length(intnodes),(nbzitr0-2))==T)==T) stop("Intnodes does not include the correct number of interior nodes")     
		intnodes <- sort(intnodes)
		if(intnodes[1]<=zitr[1]|intnodes[nbzitr0-2]>=zitr[nbzitr0])stop("Intnodes are not inside the boundaries of the marker")     
		zitr[2:(nbzitr0-1)] <- intnodes[1:(nbzitr0-2)]
		}
	}    
	if(all.equal("quant",type)==T){
		pas <-c(1:(nbzitr0-2))/(nbzitr0-1) 
		zitr[2:(nbzitr0-1)] <- quantile(sort(Y0),probs=pas)
		if(length(unique(zitr[1:nbzitr0]))!=nbzitr0) stop("The link function can not be estimated since some nodes are equal; Please try to reduce the number of nodes or use manual location.")
	}        	       
	if(all.equal("equi",type)==T){
		pas=as.double(maxY0-minY0)/as.double(nbzitr0-1)
		for(i in 2:(nbzitr0-1)){
		zitr[i] <- zitr[i-1]+pas
		}
	}

 #verifier s'il y a des obs entre les noeuds
 hcounts <- hist(Y0,breaks=zitr,plot=FALSE,include.lowest=TRUE,right=TRUE)$counts
 if(any(hcounts==0)) stop("Link function can not be estimated since some intervals defined by the nodes do not contain any observation.")  
}      


if (idlink0==1) {
if (epsY<=0) {
epsY <- 0.5
cat("Argument 'epsY' should be a definite positive real. It is changed to the default value of 0.5. \n") 
}
}


Ydiscrete <- 1
if (idlink0!=3) {
        if(!(all.equal(minY0,as.integer(minY0))==T) | !(all.equal(maxY0,as.integer(maxY0))==T)|!all(Y0 %in% minY0:maxY0)){
		Ydiscrete <- 0
	}

}






### test partialH que pour beta ou splines
if(!(idlink0 %in% c(1,2)) & partialH) stop("No partial Hessian can be define")

           






 

link <- as.character(link)
### appel des differents modeles selon la valeur de l'argument link
result <- switch(link
,"linear"=.Contlcmm(fixed=fixed,mixture=mixture,random=random,subject=subject,classmb=classmb,ng=ng,idiag=idiag,nwg=nwg,cor=cor.char,data=data,B=B,convB=convB,convL=convL,convG=convG,prior=prior,maxiter=maxiter,epsY=epsY,idlink0=idlink0,ntrtot0=ntrtot0,nbzitr0=nbzitr0,zitr=zitr,nsim=nsim,call=mm,Ydiscrete,subset=subset,na.action,posfix=posfix,partialH=partialH,verbose=verbose)

,"beta"=.Contlcmm(fixed=fixed,mixture=mixture,random=random,subject=subject,classmb=classmb,ng=ng,idiag=idiag,nwg=nwg,cor=cor.char,data=data,B=B,convB=convB,convL=convL,convG=convG,prior=prior,maxiter=maxiter,epsY=epsY,idlink0=idlink0,ntrtot0=ntrtot0,nbzitr0=nbzitr0,zitr=zitr,nsim=nsim,call=mm,Ydiscrete,subset=subset,na.action,posfix=posfix,partialH=partialH,verbose=verbose)

,"splines"=.Contlcmm(fixed=fixed,mixture=mixture,random=random,subject=subject,classmb=classmb,ng=ng,idiag=idiag,nwg=nwg,cor=cor.char,data=data,B=B,convB=convB,convL=convL,convG=convG,prior=prior,maxiter=maxiter,epsY=epsY,idlink0=idlink0,ntrtot0=ntrtot0,nbzitr0=nbzitr0,zitr=zitr,nsim=nsim,call=mm,Ydiscrete,subset=subset,na.action,posfix=posfix,partialH=partialH,verbose=verbose)
                 
,"thresholds"=.Ordlcmm(fixed=fixed,mixture=mixture,random=random,subject=subject,classmb=classmb,ng=ng,idiag=idiag,nwg=nwg,data=data,B=B,convB=convB,convL=convL,convG=convG,prior=prior,maxiter=maxiter,zitr=zitr,ide=ide0,call=mm,Ydiscrete,subset=subset,na.action=na.action,posfix=posfix,verbose=verbose))
  
return(result)
}




