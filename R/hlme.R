###################### derniere mise a jour : 2012/03/16 ############"


hlme <-
    function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,cor=NULL,data,B,convB=0.0001,convL=0.0001,convG=0.0001,prior,maxiter=500,subset=NULL,na.action=1,posfix=NULL,verbose=TRUE){

        ptm<-proc.time()
        if(verbose==TRUE) cat("Be patient, hlme is running ... \n")

        cl <- match.call()
        args <- as.list(match.call(hlme))[-1]

        nom.subject <- as.character(subject)
#### INCLUSION PRIOR
        nom.prior <- as.character(args$prior)
####
        if(!missing(mixture) & ng==1) stop("No mixture can be specified with ng=1")
        if(missing(mixture) & ng>1) stop("The argument mixture has to be specified for ng > 1")
        if(!missing(classmb) & ng==1) stop("No classmb can be specified with ng=1")
        if(missing(random)) random <- ~-1
        if(missing(fixed)) stop("The argument Fixed must be specified in any model")
        if(missing(classmb)) classmb <- ~-1
        if(missing(mixture)) mixture <- ~-1
        if(ng==1&nwg==TRUE) stop ("The argument nwg should be FALSE for ng=1")


        if(class(fixed)!="formula") stop("The argument fixed must be a formula")
        if(class(mixture)!="formula") stop("The argument mixture must be a formula")
        if(class(random)!="formula") stop("The argument random must be a formula")
        if(class(classmb)!="formula") stop("The argument classmb must be a formula")
        if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
        if(nrow(data)==0) stop("Data should not be empty") 
        if(missing(subject)){ stop("The argument subject must be specified in any model even without random-effects")} 

        if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument") 

        if(na.action==1){
            na.action=na.omit
        }else{
            na.action=na.fail
        }

### test de l'argument cor
        ncor0 <- 0
        cor.type <- cl$cor[1]
        cor.time <- cl$cor[2] 
        cor <- paste(cor.type,cor.time,sep="-")
        if (all.equal(cor,character(0))!=TRUE)
            {
                if (substr(cor,1,2)=="AR") { ncor0 <- 2 }
                else if (substr(cor,1,2)=="BM") { ncor0 <- 1  }
                else { stop("The argument cor must be of type AR or BM") }
                 
                if(!(strsplit(cor,"-")[[1]][2] %in% colnames(data))) stop("Unable to find time variable from argument cor in data")
                else { cor.var.time <- strsplit(cor,"-")[[1]][2] }
            }  
### fin test argument cor 



### ad 2/04/2012
        X0.names2 <- c("intercept")
### ad 
        int.fixed <- 0
        int.mixture <- 0
        int.random <- 0
        int.classmb <- 0
                                        #7/05/2012
### Traitement des donnees manquantes
                                        # fixed
        m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]  
        m$formula <- terms(fixed)
        m$na.action=na.action 
        m[[1]] <- as.name("model.frame")	 
        m <- eval(m, sys.parent())      
        na.fixed <- attr(m,"na.action") 

                                        # mixture
        if(mixture[[2]] != "-1"){
            m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
            m$formula <- terms(mixture)
            m$na.action <- na.action
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.mixture <- attr(m,"na.action")	
        }else{
            na.mixture <- NULL
        }

                                        # random
        if(random[[2]] != "-1"){
            m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
            m$formula <- terms(random)
            m$na.action <- na.action
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.random <- attr(m,"na.action")
        }else{
            na.random <- NULL
        }

                                        # classmb
        if(classmb[[2]] != "-1"){ 
            m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]	
            m$formula <- terms(classmb)
            m$na.action <- na.action
            m[[1]] <- as.name("model.frame")	
            m <- eval(m, sys.parent()) 
            na.classmb <- attr(m,"na.action")
        }else{
            na.classmb <- NULL
        }
        
                                        #cor     
        if(ncor0!=0)
            {
                m <- match.call()[c(1,match(c("data","subset","na.action"),names(match.call()),0))]
                m$formula <- as.formula(paste(cor.var.time,1,sep="~"))
                m$na.action <- na.action
                m[[1]] <- as.name("model.frame")
                m <- eval(m,sys.parent())    
                na.cor <- attr(m,"na.action") 	
            }
        else { na.cor <- NULL }

### names of covariate in intial fit     (sans les interactions)
        X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(fixed)),data=data))[-1]))
        if(mixture[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(mixture)),data=data))))
        if(random[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(mtemp <- get_all_vars(formula(terms(random)),data=data))))
                                        #7/05/2012
        if(classmb[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(classmb)),data=data))))
                                        #7/05/2012




        ## Table sans donnees manquante: newdata
	na.action <- unique(c(na.fixed,na.mixture,na.random,na.classmb,na.cor))
        ## dans na.action, on a les indices des NA dans le subset de data

                                        #prendre le subset :
        newdata <- data  
        if(!isTRUE(all.equal(as.character(cl$subset),character(0))))
            {
                cc <- cl
                cc <- cc[c(1,which(names(cl)=="subset"))]
                cc[[1]] <- as.name("model.frame")
                cc$formula <- formula(paste("~",paste(colnames(data),collapse="+")))
                cc$data <- data
                cc$na.action <- na.pass
                newdata <- eval(cc)
            }

                                        #enlever les NA
	if(!is.null(na.action)){
            newdata <- newdata[-na.action,]
	}

        attributes(newdata)$terms <- NULL

        ## Construction de nouvelles var explicatives sur la nouvelle table
        ## fixed
	X_fixed <- model.matrix(fixed,data=newdata)
	if(any(colnames(X_fixed)=="(Intercept)")){
            ii <- which(colnames(X_fixed)=="(Intercept)")
            colnames(X_fixed)[ii] <- "intercept"
            int.fixed <- 1
	}
	nom.fixed <- inddepvar.fixed <- inddepvar.fixed.nom <- colnames(X_fixed)
	if(int.fixed>0)inddepvar.fixed <- inddepvar.fixed[-ii]

        ## mixture
	if(mixture[[2]] != "-1"){
            X_mixture <- model.matrix(mixture,data=newdata)	
            if(any(colnames(X_mixture)=="(Intercept)")){
                ii <- which(colnames(X_mixture)=="(Intercept)")
                colnames(X_mixture)[ii] <- "intercept"
                int.mixture <- 1
            }
            nom.mixture <- inddepvar.mixture <- inddepvar.mixture.nom <- colnames(X_mixture)
            if(int.mixture>0)inddepvar.mixture <- inddepvar.mixture[-ii]
            id.X_mixture <- 1
	}else{
            inddepvar.mixture <- nom.mixture <- inddepvar.mixture.nom <- NULL
            id.X_mixture <- 0
	}
        ## random
	if(random[[2]] != "-1"){
            X_random <- model.matrix(random,data=newdata)	
            if(any(colnames(X_random)=="(Intercept)")){
                ii <- which(colnames(X_random)=="(Intercept)")
                colnames(X_random)[ii] <- "intercept"
                int.random <- 1
            }
            inddepvar.random <- inddepvar.random.nom <- colnames(X_random)
            if(int.random>0) inddepvar.random <- inddepvar.random[-ii]
            id.X_random <- 1
	}else{
            ## ad: add inddepvar.random.nom2 <- NULL 10/04/2012
            inddepvar.random <- inddepvar.random.nom <- NULL
            id.X_random <- 0
	}
        ## classmb
	if(classmb[[2]] != "-1"){ 
            if(attr(terms(classmb),"intercept")==0)
                {
                    classmb <- paste("~",classmb[2],"+1") 
                }
            X_classmb <- model.matrix(as.formula(classmb),data=newdata)
            if(any(colnames(X_classmb)=="(Intercept)")){
                ii <- which(colnames(X_classmb)=="(Intercept)")
                colnames(X_classmb)[ii] <- "intercept"
                int.classmb <- 1
            }
            id.X_classmb <- 1
            if(int.classmb>0) inddepvar.classmb <- colnames(X_classmb)[-ii]
            inddepvar.classmb.nom <- colnames(X_classmb)
	}
        else{
            inddepvar.classmb <- inddepvar.classmb.nom <- "intercept"
            id.X_classmb <- 0
	}	
        
        
                                        #7/05/2012
##############   COVARIATES       ##########################
                                        # intercept is always in inddepvar.classmb
        var.exp <- NULL
        var.exp <- c(var.exp,colnames(X_fixed))
        if(id.X_mixture == 1) var.exp <- c(var.exp,colnames(X_mixture))
        if(id.X_random == 1)var.exp <- c(var.exp,colnames(X_random))
        if(id.X_classmb == 1)var.exp <- c(var.exp,colnames(X_classmb))
        var.exp <- unique(var.exp)    
        if(ncor0>0) 
            { if(!(cor.var.time %in% var.exp)) 
                  {var.exp <- c(var.exp, cor.var.time)} #si la varaible de temps dans cor n'est dan sles variables expl, on l'ajoute
          }             
        
                                        #if(!(all(nom.mixture %in% nom.fixed))) stop("The covariates in mixture should be also included in the argument fixed")
                                        # controler si les variables de mixture sont toutes dans fixed : 
        z.fixed <- strsplit(nom.fixed,split=":",fixed=TRUE)
        z.fixed <- lapply(z.fixed,sort)
        
        if(id.X_mixture==1)
            {
                z.mixture <- strsplit(nom.mixture,split=":",fixed=TRUE)
                z.mixture <- lapply(z.mixture,sort)
            }
        else z.mixture <- list()

        if(!all(z.mixture %in% z.fixed))  stop("The covariates in mixture should also be included in the argument fixed")


        ## var dependante
        Y.name <- as.character(attributes(terms(fixed))$variables[2])
        Y0 <- newdata[,Y.name]

        ## var expli

        X0 <- X_fixed
        oldnames <- colnames(X0)

        z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
        z.X0 <- lapply(z.X0,sort)


        if(id.X_mixture == 1)
            {
                z.mixture <- strsplit(colnames(X_mixture),split=":",fixed=TRUE)
                z.mixture <- lapply(z.mixture,sort)
                for(i in 1:length(colnames(X_mixture)))
                    {
                                        #if((colnames(X_mixture)[i] %in% colnames(X0))==F){ 
                        if(!isTRUE(z.mixture[i] %in% z.X0))
                            {
                                X0 <- cbind(X0,X_mixture[,i])
                                colnames(X0) <- c(oldnames, colnames(X_mixture)[i])
                                oldnames <- colnames(X0)
                                
                                z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
                                z.X0 <- lapply(z.X0,sort)			
                            }
                    }
            }
        else
            {
                z.mixture <- list()
            }

        if(id.X_random == 1)
            {
                z.random <- strsplit(colnames(X_random),split=":",fixed=TRUE)
                z.random <- lapply(z.random,sort)
                for(i in 1:length(colnames(X_random)))
                    {
                                        #if((colnames(X_random)[i] %in% colnames(X0))==F){
                        if(!isTRUE(z.random[i] %in% z.X0))
                            {		
                                X0 <- cbind(X0,X_random[,i])
                                colnames(X0) <- c(oldnames, colnames(X_random)[i])
                                oldnames <- colnames(X0)
                                
                                z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
                                z.X0 <- lapply(z.X0,sort)				
                            }	                                                            
                    }
            }
        else
            {
                z.random <- list()
            }

        if(id.X_classmb == 1)
            {
                z.classmb <- strsplit(colnames(X_classmb),split=":",fixed=TRUE)
                z.classmb <- lapply(z.classmb,sort)
                for(i in 1:length(colnames(X_classmb)))
                    {
                                        #		if((colnames(X_classmb)[i] %in% colnames(X0))==F){
                        if(!isTRUE(z.classmb[i] %in% z.X0))
                            {			
                                X0 <- cbind(X0,X_classmb[,i])
                                colnames(X0) <- c(oldnames, colnames(X_classmb)[i])
                                oldnames <- colnames(X0)
                                
                                z.X0 <- strsplit(colnames(X0),split=":",fixed=TRUE)
                                z.X0 <- lapply(z.X0,sort)	         	 
                            }	
                    }
            }
        else
            {
                z.classmb <- list()
            }

        if(ncor0>0) 
            { 
                if(!(cor.var.time %in% colnames(X0)))   #cor.var.time jamais ne interaction donc ok (pas de z.cor)
                    {
                        X0 <- cbind(X0, newdata[,cor.var.time])
                        colnames(X0) <- c(oldnames, cor.var.time)
                    }
            }  

                                        #colnames(X0) <- var.exp # a remettre si on enleve les z.fixed etc

                                        #X0 <- X0[,-which(colnames(X0)=="intercept")]
### ad

        if((any(is.na(X0))==TRUE)|(any(is.na(Y0))==TRUE))stop("The data should not contain any missing value")
                                        # 
        n <- dim(data)[1]
                                        #if ((int.fixed+int.random)>0) X0<- cbind(intercept=rep(1,n),X0)
                                        #if (!((int.fixed+int.random)>0) & ("intercept" %in% colnames(X0))) X0 <- as.data.frame(X0[,-which(colnames(X0)=="intercept")])
        nom.X0 <- colnames(X0)
        nvar.exp <- length(nom.X0)

        IND <- newdata[,nom.subject]

        #IDnum <- as.numeric(IND)


#### INCLUSION PRIOR 
        #if(missing(prior)){ PRIOR <- seq(0,length=length(IDnum))} 
        if(missing(prior)){ PRIOR <- seq(0,length=length(IND))} 
        if(!missing(prior)){ 
            PRIOR <- newdata[,nom.prior]
            PRIOR[(is.na(PRIOR))] <- 0
        }
####

        ng0 <- ng
        idiag0 <- as.integer(idiag)
        nwg0 <- as.integer(nwg)

        idea0 <- rep(0,nvar.exp)
        idprob0 <- rep(0,nvar.exp)
        idg0 <- rep(0,nvar.exp)
        
        z.X0 <- strsplit(nom.X0,split=":",fixed=TRUE)
        z.X0 <- lapply(z.X0,sort)

        
        for (i in 1:nvar.exp){  
                                        #idea0[i] <- nom.X0[i]%in%inddepvar.random.nom
                                        #idprob0[i] <- nom.X0[i]%in%inddepvar.classmb.nom
                                        #if(nom.X0[i]%in%nom.fixed & !(nom.X0[i]%in%nom.mixture)) idg0[i] <- 1 
                                        #if(nom.X0[i]%in%nom.fixed & nom.X0[i]%in%nom.mixture) idg0[i] <- 2 

            idea0[i] <- z.X0[i] %in% z.random
            idprob0[i] <- z.X0[i] %in% z.classmb   
            if((z.X0[i] %in% z.fixed) & !(z.X0[i] %in% z.mixture)) idg0[i] <- 1 
            if((z.X0[i] %in% z.fixed) & (z.X0[i] %in% z.mixture)) idg0[i] <- 2 
        }

        idcor0 <- rep(0,length(z.X0))  
        if (ncor0!=0) idcor0 <- as.numeric(nom.X0 %in% cor.var.time) 

                                        #if((int.fixed+int.random)>0) idprob0[1] <- 0     
        if(("intercept" %in% nom.X0) | ("(Intercept)" %in% nom.X0))
            {
                ii <- which(nom.X0 %in% c("intercept","(Intercept)"))
                idprob0[ii] <- 0  
            }


                                        # on ordonne les donn es suivants la variable IND
        #matYX <- cbind(IDnum,IND,PRIOR,Y0,X0)
        #matYXord <- matYX[sort.list(matYX[,1]),]
        #Y0 <- matYXord[,4]  
        #X0 <- matYXord[,-c(1,2,3,4)]
        #IDnum <- matYXord[,1]
        #IND <- matYXord[,2]

        matYX <- cbind(IND,PRIOR,Y0,X0)
        matYXord <- matYX[sort.list(matYX[,1]),]
        Y0 <- as.numeric(matYXord[,3])
        X0 <- apply(matYXord[,-c(1,2,3),drop=FALSE],2,as.numeric)
        IND <- matYXord[,1]


#### INCLUSION PRIOR 
        PRIOR <- as.numeric(matYXord[,2])
        PRIOR <-as.integer(as.vector(PRIOR))
####

        X0<-as.numeric(as.matrix(X0))
        Y0<-as.numeric(as.matrix(Y0))
        #nmes0<-as.vector(table(IDnum))
        nmes0<-as.vector(table(IND))
        ns0<-length(nmes0)



##### INCLUSION PRIOR 
                                        # definition de prior a 0 pour l'analyse G=1
        prior2 <- as.integer(rep(0,ns0))
        prior0 <- prior2 
                                        # si prior pas missing alors mettre dedans la classe a priori. Attention tester q les valeurs sont dans 0, G
        if(!missing(prior)){ 
            prior0 <- PRIOR[cumsum(nmes0)]
        }
        INDuniq <- IND[cumsum(nmes0)]
        seqnG <- 0:ng0
        if (!(all(prior0  %in% seqnG))) stop ("The argument prior should contain integers between 0 and ng")
#####


        loglik <- as.double(0)
        ni <- 0
        istop <- 0
        gconv <-rep(0,3)
        ppi0 <- rep(0,ns0*ng0)
        nv0<-nvar.exp
        nobs0<-length(Y0)
        resid_m <- rep(0,nobs0)
        resid_ss <- rep(0,nobs0)
        pred_m_g <- rep(0,nobs0*ng0)
        pred_ss_g <- rep(0,nobs0*ng0)
        nea0 <- sum(idea0==1)
        predRE <- rep(0,nea0*ns0)

        ##---------------------------------------------------------------------------
        ##definition du vecteur de parametres + initialisation
        ##---------------------------------------------------------------------------

## gestion de B=random(mod)

        Brandom <- FALSE
        if(length(cl$B)==2)
            {
                if(class(eval(cl$B[[2]]))!="hlme") stop("The model specified in B should be of class hlme")
                if(as.character(cl$B[1])!="random") stop("Please use random() to specify random initial values")
                
                Brandom <- TRUE
                B <- eval(cl$B[[2]])

                if(length(posfix)) stop("Argument posfix is not compatible with random intial values")
            }

        

#####cas 1 : ng=1
        b <- NULL
        b1 <- NULL
        NPROB <- 0
        if(ng0==1| missing(B)){
            NEF<-sum(idg0!=0)  
            b1[1:NEF] <- 0
            if(int.fixed > 0)  b1[1] <- mean(Y0)

            if(idiag0==1){
                NVC<-sum(idea0==1)
                b1[(NEF+1):(NEF+NVC)] <- 1}

            if(idiag0==0){
                kk<-sum(idea0==1) 
                NVC<-(kk*(kk+1))/2
                indice<-cumsum(1:kk)
                bidiag<-rep(0,NVC)
                bidiag[indice] <- 1
                b1[(NEF+1):(NEF+NVC)] <- bidiag
            }
            if(ncor0==1)
                {b1[NEF+NVC+1] <- 1 }
            if(ncor0==2)
                {b1[(NEF+NVC+1):(NEF+NVC+ncor0)] <- c(0,1) }
            b1[NEF+NVC+ncor0+1] <- 1
            NPM <- length(b1)
            NW <- 0
            V <- rep(0,NPM*(NPM+1)/2) 
        }

#####cas 2 : ng>=2
        if(ng0>1){
            NPROB <- (sum(idprob0==1)+1)*(ng0-1)
            b[1:NPROB] <- 0
            NEF <- sum(idg0==1)+(sum(idg0==2))*ng0
            if(idiag0==1) NVC <- sum(idea0==1)
            if(idiag0==0){
                kk <- sum(idea0==1) 
                NVC <- (kk*(kk+1))/2}
            NW <- nwg0*(ng0-1)
            if(NW>0) b[(NPROB+NEF+NVC+1):(NPROB+NEF+NVC+NW)] <- 1   
            if(ncor0==1)
                {b[NEF+NVC+NW+1] <- 1 }
            if(ncor0==2)
                {b[(NPROB+NEF+NVC+NW+1):(NPROB+NEF+NVC+NW+ncor0)] <- c(0,1) }
            NPM <- NPROB+NEF+NVC+NW+ncor0+1
            V <- rep(0,NPM*(NPM+1)/2)
        }

        
        ## prm fixes
        fix0 <- rep(0,NPM)
        if(length(posfix))
            {
                if(any(!(posfix %in% 1:NPM))) stop("Indexes in posfix are not correct")
                
                fix0[posfix] <- 1
            }
        if(length(posfix)==NPM) stop("No parameter to estimate")
        
        
        if(missing(B)){

            if(ng0>1){
                idea2 <- idea0
                idprob2 <- rep(0,nv0)  
                idg2 <- rep(0,nv0) 
                idg2[idg0!=0] <- 1
                NEF2<-sum(idg2==1)
                NPM2<-NEF2+NVC+ncor0+1
                nwg2<-0
                ng2<-1
                ppi2<- rep(0,ns0)
                pred_m_g2 <- rep(0,nobs0)
                pred_ss_g2 <- rep(0,nobs0)
                maxiter2 <- min(75,maxiter)
                convB2 <- max(0.01,convB)
                convL2 <- max(0.01,convL)
                convG2 <- max(0.01,convG)

                V2 <- rep(0,NPM2*(NPM2+1)/2)
                best <- rep(0,NPM2)
                
                init <- .Fortran(C_hetmixlin,
                                 as.double(Y0),
                                 as.double(X0),
                                 as.integer(prior2),
                                 as.integer(idprob2),
                                 as.integer(idea2),
                                 as.integer(idg2),
                                 as.integer(idcor0),
                                 as.integer(ns0),
                                 as.integer(ng2),
                                 as.integer(nv0),
                                 as.integer(nobs0),
                                 as.integer(nea0),
                                 as.integer(nmes0),
                                 as.integer(idiag0),
                                 as.integer(nwg2),
                                 as.integer(ncor0),
                                 npm=as.integer(NPM2),
                                 best=as.double(b1),
                                 V=as.double(V2),
                                 loglik=as.double(loglik),
                                 niter=as.integer(ni),
                                 conv=as.integer(istop),
                                 gconv=as.double(gconv),
                                 ppi2=as.double(ppi2),
                                 resid_m=as.double(resid_m),
                                 resid_ss=as.double(resid_ss),
                                 pred_m_g=as.double(pred_m_g2),
                                 pred_ss_g=as.double(pred_ss_g2),
                                 predRE=as.double(predRE),
                                 as.double(convB2),
                                 as.double(convL2),
                                 as.double(convG2),
                                 as.integer(maxiter2),
                                 as.integer(fix0))

                k <- NPROB
                l <- 0
                t<- 0
                for (i in 1:nvar.exp)    {
                    if(idg0[i]==1){
                        l <- l+1
                        t <- t+1
                        b[k+t] <- init$best[l]
                    }
                    if(idg0[i]==2){
                        l <- l+1
                        for (g in 1:ng){
                            t <- t+1
                            if(init$conv==1) b[k+t] <- init$best[l]+(g-(ng+1)/2)*sqrt(init$V[l*(l+1)/2])
                            else b[k+t] <- init$best[l]+(g-(ng+1)/2)*init$best[l]
                        }
                    }
                }
                b[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- init$best[(NEF2+1):(NEF2+NVC)]
                if (ncor0>0) {b[(NPROB+NEF+NVC+NW+1):(NPROB+NEF+NVC+NW+ncor0)] <- init$best[(NPM2-ncor0):(NPM2-1)]}
                b[NPROB+NEF+NVC+NW+ncor0+1] <- init$best[NPM2]
            } 
            if(ng0==1 ){
                b <- b1
            }
        }
        else
            {
                if(is.vector(B))
                    {
                        if(length(B)!=NPM) stop(paste("Vector B should be of length",NPM))
                        else {b <-B}
                    }
                else
                    { 
                        if(class(B)!="hlme") stop("B should be either a vector or an object of class hlme")

                        ## B est le meme modele mais pr ng=1 :
                        if(ng>1 & B$ng==1)
                            {
                                NEF2 <- sum(idg0!=0)
                                NPM2 <- NEF2+NVC+ncor0+1
                                if(length(B$best)!=NPM2) stop("B is not correct")


                                if(Brandom==FALSE)
                                    {
                                        ## B deterministe
                                        l <- 0
                                        t <- 0
                                        for (i in 1:nvar.exp)
                                            {
                                                if(idg0[i]==1)
                                                    {
                                                        l <- l+1
                                                        t <- t+1
                                                        b[NPROB+t] <- B$best[l]
                                                    }
                                                if(idg0[i]==2)
                                                    {
                                                        l <- l+1
                                                        for (g in 1:ng)
                                                            {
                                                                t <- t+1
                                                                if(B$conv==1) b[NPROB+t] <- B$best[l]+(g-(ng+1)/2)*sqrt(B$V[l*(l+1)/2])
                                                                else b[NPROB+t] <- B$best[l]+(g-(ng+1)/2)*B$best[l]
                                                            }
                                                    }
                                            }
                                        if(NVC>0)
                                            {
                                                if(idiag==TRUE)
                                                    {
                                                        b[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- B$cholesky[(1:nea0)*(2:(nea0+1))/2]
                                                        
                                                    }
                                                else
                                                    {
                                                        b[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- B$cholesky
                                                    }
                                            }
                                        if (ncor0>0) {b[(NPROB+NEF+NVC+NW+1):(NPROB+NEF+NVC+NW+ncor0)] <- B$best[(NPM2-ncor0):(NPM2-1)]}
                                        b[NPROB+NEF+NVC+NW+ncor0+1] <- B$best[NPM2]
                                    }
                                else
                                    {
                                        ## B random
                                        bb <- rep(0,NPM-NPROB-NW)
                                        vbb <- matrix(0,NPM-NPROB-NW,NPM-NPROB-NW)
                                        
                                        VB <- matrix(0,NPM2,NPM2)
                                        VB[upper.tri(VB,diag=TRUE)] <- B$V
                                        VB <- t(VB)
                                        VB[upper.tri(VB,diag=TRUE)] <- B$V

                                        nbg <- idg0[which(idg0!=0)]
                                        nbg[which(nbg==2)] <- ng
                                        nbgnef <- unlist(sapply(nbg,function(k) if(k>1) rep(2,k) else k))
                                        
                                        vbb[which(nbgnef==1),setdiff(1:ncol(vbb),which(nbgnef!=1))] <- VB[which(nbg==1),setdiff(1:ncol(VB),which(nbg!=1))]
                                        vbb[(NEF+1):nrow(vbb),(NEF+1):ncol(vbb)] <- VB[(NEF2+1):nrow(VB),(NEF2+1):ncol(VB)]
                                        
                                        l <- 0
                                        t <- 0
                                        for (i in 1:nvar.exp)
                                            {
                                                if(idg0[i]==1)
                                                    {
                                                        l <- l+1
                                                        t <- t+1
                                                        bb[t] <- B$best[l]
                                                    }
                                                if(idg0[i]==2)
                                                    {
                                                        l <- l+1
                                                        for (g in 1:ng)
                                                            {
                                                                t <- t+1
                                                                bb[t] <- B$best[l]
                                                                vbb[t,t] <- VB[l,l]
                                                            }
                                                    }
                                            }
                                        
                                        if(NVC>0)
                                            {
                                                if(idiag==TRUE)
                                                    {
                                                        bb[NEF+1:NVC] <- B$cholesky[(1:nea0)*(2:(nea0+1))/2]
                                                    }
                                                else
                                                    {
                                                        bb[NEF+1:NVC] <- B$cholesky
                                                    }
                                            }
                            
                                
                                        if (ncor0>0)
                                            {
                                                bb[NEF+NVC+1:ncor0] <- B$best[(NPM2-ncor0):(NPM2-1)]
                                            }
                                
                                        bb[NEF+NVC+ncor0+1] <- B$best[NPM2]
                                                            
                                        up <- vbb[upper.tri(vbb,diag=TRUE)]
                                        vbb <- t(vbb)
                                        vbb[upper.tri(vbb,diag=TRUE)] <- up
                                        Chol <- chol(vbb)
                                        Chol <- t(Chol)
                                        
                                        b[c((NPROB+1):(NPROB+NEF+NVC),(NPROB+NEF+NVC+NW+1):NPM)] <- bb + Chol %*% rnorm(NPM-NPROB-NW)

                                        b[1:NPROB] <- 0
                                        if(NW>0) b[NPROB+NEF+NVC+1:NW] <- 1
                                   
                                    } 
                            }
                    }
            }
                
  
        ##------------------------------------------
        ##------nom au vecteur best
        ##--------------------------------------------


        if(ng0>=2){
            nom <-rep(c("intercept",nom.X0[idprob0==1]),each=ng0-1)
            nom1 <- paste(nom," class",c(1:(ng0-1)),sep="")
            names(b)[1:NPROB]<-nom1
        }

        if(ng0==1) names(b)[1:NEF] <- nom.X0[idg0!=0]

        if(ng0>1){
            nom1<- NULL
            for (i in 1:nvar.exp) {
                if(idg0[i]==2){ nom <- paste(nom.X0[i]," class",c(1:ng0),sep="")
                                nom1 <- cbind(nom1,t(nom))}
                if(idg0[i]==1) nom1 <- cbind(nom1,nom.X0[i])
            }
            names(b)[(NPROB+1):(NPROB+NEF)]<- nom1
        }
        if(NVC!=0)names(b)[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- paste("varcov",c(1:(NVC)))
        if(NW!=0)names(b)[(NPROB+NEF+NVC+1):(NPROB+NEF+NVC+NW)] <- paste("varprop class",c(1:(ng0-1)))
        names(b)[NPM] <- "stderr"
        if(ncor0!=0) {names(b)[(NPROB+NEF+NVC+NW+1):(NPROB+NEF+NVC+NW+ncor0)] <- paste("cor",1:ncor0,sep="") }

        N <- NULL
        N[1] <- NPROB
        N[2] <- NEF
        N[3] <- NVC
        N[4] <- NW
        N[5] <- ncor0

        idiag <- as.integer(idiag0)
        idea <- as.integer(idea0)
        nv <- as.integer(nv0)


################ Sortie ###########################

        out <- .Fortran(C_hetmixlin,
                        as.double(Y0),
                        as.double(X0),
                        as.integer(prior0),
                        as.integer(idprob0),
                        as.integer(idea0),
                        as.integer(idg0),
                        as.integer(idcor0),
                        as.integer(ns0),
                        as.integer(ng0),
                        as.integer(nv0),
                        as.integer(nobs0),
                        as.integer(nea0),
                        as.integer(nmes0),
                        as.integer(idiag0),
                        as.integer(nwg0),
                        as.integer(ncor0),
                        as.integer(NPM),
                        best=as.double(b),
                        V=as.double(V),
                        loglik=as.double(loglik),
                        niter=as.integer(ni),
                        conv=as.integer(istop),
                        gconv=as.double(gconv),
                        ppi2=as.double(ppi0),
                        resid_m=as.double(resid_m),
                        resid_ss=as.double(resid_ss),
                        pred_m_g=as.double(pred_m_g),
                        pred_ss_g=as.double(pred_ss_g),
                        predRE=as.double(predRE),
                        as.double(convB),
                        as.double(convL),
                        as.double(convG),
                        as.integer(maxiter),
                        as.integer(fix0))
        
    ### mettre 0 pr les prm fixes
    if(length(posfix))
        {
            mr <- NPM-length(posfix)
            Vr <- matrix(0,mr,mr)
            Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
            Vr <- t(Vr)
            Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
            V <- matrix(0,NPM,NPM)
            V[setdiff(1:NPM,posfix),setdiff(1:NPM,posfix)] <- Vr
            V <- V[upper.tri(V,diag=TRUE)]
        }
    else
        {
            V <- out$V
        }


            
### Creation du vecteur cholesky
        Cholesky <- rep(0,(nea0*(nea0+1)/2))
        if(idiag0==0 & NVC>0){
            Cholesky[1:NVC] <- out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)]
### Construction de la matrice U 
            U <- matrix(0,nrow=nea0,ncol=nea0)
            U[upper.tri(U,diag=TRUE)] <- Cholesky[1:NVC]
            z <- t(U) %*% U
            out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- z[upper.tri(z,diag=TRUE)]
        }
        if(idiag0==1 & NVC>0){
            id <- 1:nea0
            indice <- rep(id+id*(id-1)/2)
            Cholesky[indice] <- out$best[(NPROB+NEF+1):(NPROB+NEF+nea0)]
            out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- out$best[(NPROB+NEF+1):(NPROB+NEF+NVC)]**2 
        } 

####################################################
        if (nea0>0) {
            predRE <- matrix(out$predRE,ncol=nea0,byrow=T)
            predRE <- data.frame(INDuniq,predRE)
            colnames(predRE) <- c(nom.subject,nom.X0[idea0!=0])
        }


        if(ng0>1) {
            ppi<- matrix(out$ppi2,ncol=ng0,byrow=TRUE)
        }
        else {
            ppi <- matrix(rep(1,ns0),ncol=ng0)
        }


        classif<-apply(ppi,1,which.max)
        ppi<-data.frame(INDuniq,classif,ppi)
        temp<-paste("prob",1:ng0,sep="")
        colnames(ppi) <- c(nom.subject,"class",temp)
        rownames(ppi) <- 1:ns0

        pred_m_g <- matrix(out$pred_m_g,nrow=nobs0)
        pred_ss_g <- matrix(out$pred_ss_g,nrow=nobs0)
        pred_m <- Y0-out$resid_m
        pred_ss <- Y0-out$resid_ss
        pred <- data.frame(IND,pred_m,out$resid_m,pred_ss,out$resid_ss,Y0,pred_m_g,pred_ss_g)

        temp<-paste("pred_m",1:ng0,sep="")
        temp1<-paste("pred_ss",1:ng0,sep="")
        colnames(pred)<-c(nom.subject,"pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1) 

        names(out$best)<-names(b)
        btest <- out$best[1:length(inddepvar.fixed.nom)]
        names(btest) <-inddepvar.fixed.nom

### ad 2/04/2012
        if (!("intercept" %in% nom.X0)) X0.names2 <- X0.names2[-1]
### ad
        res <-list(ns=ns0,ng=ng0,idea0=idea0,idprob0=idprob0,idg0=idg0,idcor0=idcor0,loglik=out$loglik,best=out$best,V=V,gconv=out$gconv,conv=out$conv,call=cl,niter=out$niter,dataset=args$data,N=N,idiag=idiag0,pred=pred,pprob=ppi,predRE=predRE,Xnames=nom.X0,Xnames2=X0.names2,cholesky=Cholesky,na.action=na.action,AIC=2*(length(out$best)-length(posfix)-out$loglik),BIC=(length(out$best)-length(posfix))*log(ns0)-2*out$loglik)
        class(res) <-c("hlme") 
        cost<-proc.time()-ptm
        if(verbose==TRUE) cat("The program took", round(cost[3],2), "seconds \n")

        
        res
    }
