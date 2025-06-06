###################### derniere mise a jour : 2012/03/16 ############"
.Ordlcmm <-
    function(fixed,mixture,random,subject,classmb,ng,idiag,nwg,data,B,convB,convL,convG,prior,pprior=NULL,maxiter,zitr,ide,call,Ydiscrete,subset=subset,na.action,posfix,partialH,verbose,returndata,var.time,nproc,clustertype)
{
    
    ptm<-proc.time()
 
    cl <- match.call()

    args <- as.list(match.call(.Ordlcmm))[-1]

    nom.subject <- as.character(subject)

####
    if(!missing(mixture) & ng==1) stop("No mixture can be specified with ng=1")
    if(missing(mixture) & ng>1) stop("The argument mixture has to be specified for ng > 1")
    if(!missing(classmb) & ng==1) stop("No classmb can be specified with ng=1")
    if(missing(random)) random <- ~-1
    if(missing(fixed)) stop("The argument Fixed must be specified in any model")
    if(missing(classmb) & ng==1) classmb <- ~-1
    if(missing(classmb) & ng>1) classmb <- ~1
    if(missing(mixture)) mixture <- ~-1
    if(ng==1&nwg==TRUE) stop ("The argument nwg should be FALSE for ng=1")


    if(!inherits(fixed,"formula")) stop("The argument fixed must be a formula")
    if(!inherits(mixture,"formula")) stop("The argument mixture must be a formula")
    if(!inherits(random,"formula")) stop("The argument random must be a formula")
    if(!inherits(classmb,"formula")) stop("The argument classmb must be a formula")

    if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")} 
    if(missing(subject)){ stop("The argument subject must be specified in any model even without random-effects")} 
    if(!is.numeric(data[[subject]])) stop("The argument subject must be numeric")
### ad 2/04/2012
    X0.names2 <- c("intercept")
### ad
    int.fixed <- 0
    int.mixture <- 0
    int.random <- 0
    int.classmb <- 0

    ## garder data tel quel pour le renvoyer
    if(returndata==TRUE)
    {
        datareturn <- data
    }
    else
    {
        datareturn <- NULL
    }
    
    ## Table sans donnees manquante: newdata
    ##prendre le subset :
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

    ##enlever les NA
    if(!is.null(na.action)){
        newdata <- newdata[-na.action,]
    }

    attributes(newdata)$terms <- NULL

### names of covariate in intial fit
    X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(fixed)),data=newdata))[-1]))
    if(mixture[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(mixture)),data=newdata))))
    if(random[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(random)),data=newdata))))
    if(classmb[[2]] != "-1")X0.names2 <- unique(c(X0.names2,colnames(get_all_vars(formula(terms(classmb)),data=newdata))))


    ## Construction de nouvelles var eplicatives sur la nouvelle table
    ## fixed
    X_fixed <- model.matrix(fixed,data=newdata)
    if(colnames(X_fixed)[1]=="(Intercept)"){
        colnames(X_fixed)[1] <- "intercept"
        int.fixed <- 1
    }else{
        stop ("Only models with an intercept can be estimated using lcmm. This is required for identifiability purposes")
    }
    nom.fixed <- colnames(X_fixed)
    if(int.fixed>0)inddepvar.fixed <- inddepvar.fixed.nom <- nom.fixed[-1]

    ## mixture
    if(mixture[[2]] != "-1"){
        X_mixture <- model.matrix(mixture,data=newdata)	
        if(colnames(X_mixture)[1]=="(Intercept)"){
            colnames(X_mixture)[1] <- "intercept"
            int.mixture <- 1
        }
        nom.mixture <- inddepvar.mixture <- inddepvar.mixture.nom <- colnames(X_mixture)
        if(int.mixture>0)inddepvar.mixture <- inddepvar.mixture[-1]
        id.X_mixture <- 1
    }else{
        inddepvar.mixture <- nom.mixture <- inddepvar.mixture.nom <- NULL
        id.X_mixture <- 0
    }
    ## random
    if(random[[2]] != "-1"){
        X_random <- model.matrix(random,data=newdata)	
        if(colnames(X_random)[1]=="(Intercept)"){
            colnames(X_random)[1] <- "intercept"
            int.random <- 1
        }
        inddepvar.random <- inddepvar.random.nom <- colnames(X_random)
        if(int.random>0) inddepvar.random <- inddepvar.random[-1]
        id.X_random <- 1
    }else{
	## ad: add inddepvar.random.nom2 <- NULL 10/04/2012
        inddepvar.random <- inddepvar.random.nom <- NULL
        id.X_random <- 0
    }

    ## classmb
    if(classmb[[2]] != "-1"){ 
        X_classmb <- model.matrix(classmb,data=newdata)
        colnames(X_classmb)[1] <- "intercept"
        id.X_classmb <- 1
        inddepvar.classmb <- colnames(X_classmb)[-1]
        inddepvar.classmb.nom <- colnames(X_classmb)
    }else{
        inddepvar.classmb <- inddepvar.classmb.nom <- "intercept"
        id.X_classmb <- 0
    }	

                                        # intercept is always in inddepvar.classmb
    var.exp <- NULL
    var.exp <- c(var.exp,colnames(X_fixed))
    if(id.X_mixture == 1) var.exp <- c(var.exp,colnames(X_mixture))
    if(id.X_random == 1)var.exp <- c(var.exp,colnames(X_random))
    if(id.X_classmb == 1)var.exp <- c(var.exp,colnames(X_classmb))
    var.exp <- unique(var.exp)
    ## ad
    timeobs <- rep(0, nrow(newdata))
    if(!is.null(var.time))
    {
        timeobs <- newdata[,var.time]
        #if(any(is.na(timeobs))) stop(paste("Cannot use",var.time,"as time variable because it contains missing data"))
    }
    
    ## controler si les variables de mixture sont toutes dans fixed : 
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

    if((any(is.na(X0))==TRUE)|(any(is.na(Y0))==TRUE))stop("The data should not contain any missing value")
    
    n <- dim(data)[1]
                                       
    if (!((int.fixed+int.random)>0)) X0 <- as.data.frame(X0[,-which(colnames(X0)=="intercept")])
    nom.X0 <- colnames(X0)
    nvar.exp <- length(nom.X0)

    IND <- newdata[,nom.subject]
    #IDnum <- as.numeric(IND)


#### INCLUSION PRIOR 
    if(missing(prior)){ PRIOR <- seq(0,length=length(IND))} 
    if(!missing(prior)){ 
        PRIOR <- newdata[,prior]
        PRIOR[(is.na(PRIOR))] <- 0
    }
####

        ##pprior : proba d'appartenance a chaque classe
        if(is.null(pprior))
        {
            pprior0 <- matrix(1, length(IND), ng)            
        }
        else
        {
            if(ng==1) stop("pprior is only useful if ng>1")
            if(!is.character(pprior)) stop("pprior should be a character vector")
            if(length(pprior) != ng) stop("pprior should be of length ng")
            if(!all(pprior %in% colnames(newdata))) stop("pprior variables should be included in data")

            pprior0 <- newdata[,pprior]
        }
    
    ng0 <- ng
    idiag0 <- as.integer(idiag)
    nwg0 <- as.integer(nwg)

    idea0 <- rep(0,nvar.exp)
    idprob0 <- rep(0,nvar.exp)
    idg0 <- rep(0,nvar.exp)

    z.X0 <- strsplit(nom.X0,split=":",fixed=TRUE)
    z.X0 <- lapply(z.X0,sort)

    for (i in 1:nvar.exp)    
        {            
            idea0[i] <- z.X0[i] %in% z.random
            idprob0[i] <- z.X0[i] %in% z.classmb   
            if((z.X0[i] %in% z.fixed) & !(z.X0[i] %in% z.mixture)) idg0[i] <- 1 
            if((z.X0[i] %in% z.fixed) & (z.X0[i] %in% z.mixture)) idg0[i] <- 2     
        }

    #if((int.fixed+int.random)>0) idprob0[1] <- 0

    ## on ordonne les donnees suivant la variable IND
    matYX <- cbind(IND,PRIOR,pprior0,Y0,X0)
    matYXord <- matYX[sort.list(matYX[,1]),]
    Y0 <- as.numeric(matYXord[,3+ng] )
    X0 <- apply(matYXord[,-c(1:(3+ng)),drop=FALSE],2,as.numeric)
    #IDnum <- matYXord[,1]
    IND <- matYXord[,1]

#### INCLUSION PRIOR 
    PRIOR <- as.numeric(matYXord[,2])
    PRIOR <-as.integer(as.vector(PRIOR))
####

    X0<-as.numeric(as.matrix(X0))
    Y0<-as.numeric(as.matrix(Y0))
    nmes0<-as.vector(table(IND))
    ns0<-length(nmes0)


##### INCLUSION PRIOR 
    ## definition de prior a 0 pour l'analyse G=1
    prior2 <- as.integer(rep(0,ns0))
    prior0 <- prior2 
    ## si prior pas missing alors mettre dedans la classe a priori. Attention tester q les valeurs sont dans 0, G
    if(!missing(prior)){ 
        prior0 <- PRIOR[cumsum(nmes0)]
    }
    INDuniq <- IND[cumsum(nmes0)]
    seqnG <- 0:ng0
    if (!(all(prior0  %in% seqnG))) stop ("The argument prior should contain integers between 0 and ng")
#####

    ## pprior de dim ns
    if(!is.null(pprior))
    {
        pprior0 <- matYX[cumsum(nmes0),2+1:ng, drop=FALSE]
        if(any(is.na(pprior0))) stop("pprior contains missing values")
    }
    else
    {
        pprior0 <- matrix(1, ns0, ng0)
    }
    pprior0 <- t(pprior0)       
    

    loglik <- as.double(0)
    rlindiv <- rep(0,ns0)
    UACV <- as.double(0)
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
    predRE <- rep(0,nea0*ns0*(1+ng0))

    ntrtot0 <- sum(ide==1) 
    minY <- zitr[1]
    maxY <- zitr[2]



    ##-------------------------------------------------------------------------------
    ##definition du vecteur de parametre + initialisation
    ##-------------------------------------------------------------------------------
## pas de cor en ordinal
    ncor0 <- 0

#####cas 1 : ng=1
    b<-NULL
    b1 <- NULL
    NPROB <- 0
    if(ng0==1| missing(B)){
        NEF<-sum(idg0!=0)-1
        b1[1:NEF]<-0

        NVC <- 0
        if(idiag0==1&nea0>0){
            NVC<-sum(idea0==1)
            b1[(NEF+1):(NEF+NVC)]<-1}

        if(idiag0==0&nea0>0){
            kk<-sum(idea0==1) 
            NVC<-(kk*(kk+1))/2
            indice<-cumsum(1:kk)
            bidiag<-rep(0,NVC)
            bidiag[indice]<-1
            b1[(NEF+1):(NEF+NVC)]<-bidiag
        }

        ## valeurs initiales pour les transfos #
        if (ntrtot0==1) {
            b1[(NEF+NVC+1):(NEF+NVC+ntrtot0)] <- 0
        }else{
            b1[(NEF+NVC+1)] <- 2*qnorm(0.98)*(-median(Y0)+minY+1)/(ntrtot0-1)
            b1[(NEF+NVC+2):(NEF+NVC+ntrtot0)] <- sqrt(2*qnorm(0.98)/(ntrtot0-1))
        }

        NPM <- length(b1)
        NW <- 0
        V <- rep(0,NPM*(NPM+1)/2) 
    }

#####cas 2 : ng>=2
    if(ng0>1){
        NPROB <- (sum(idprob0==1))*(ng0-1)
        NEF <- sum(idg0==1)+(sum(idg0==2))*ng0-1
        NVC <- 0
        if(idiag0==1&nea0>0) NVC <- sum(idea0==1)
        if(idiag0==0&nea0>0){
            kk <- sum(idea0==1) 
            NVC <- (kk*(kk+1))/2}
        NW <- nwg0*(ng0-1)
        
        NPM <- NPROB+NEF+NVC+NW+ntrtot0
        
        b <- rep(0,NPM)
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
    
    ## pour H restreint
    Hr0 <- as.numeric(partialH)
    pbH0 <- rep(0,NPM)
    if(is.logical(partialH))
    {
        if(isTRUE(partialH))
        {
            if(partialH) pbH0 <- rep(1,NPM)
            #pbH0[posfix] <- 0
            if(sum(pbH0)==0 & Hr0==1) stop("No partial Hessian matrix can be defined")
        }
    }
    else
    {
        if(!all(Hr0 %in% 1:NPM)) stop("Indexes in partialH are not correct")
        pbH0[Hr0] <- 1
        #pbH0[posfix] <- 0
    }
    indexHr <- NULL
    if(sum(pbH0)>0)
    {
        if(length(posfix)) pbH1 <- pbH0[-posfix] else pbH1 <- pbH0
        indexHr <- which(pbH1==1)
    }

    if(missing(B)){

        if(ng0==1 ){
            b <- b1
        } else {
            stop("Please specify initial values with argument 'B'")
        }
    }
    else
        {
            if(is.vector(B))
                {
                    if(length(B)!=NPM){
                        stop("The length of the vector B is not correct")
                    }
                    else
                    {
                        b <-B
                        
                        if(NVC>0)
                        {
                            ## remplacer varcov des EA par les prm a estimer
                            
                            if(idiag0==1)
                            {
                                b[NPROB+NEF+1:NVC] <- sqrt(b[NPROB+NEF+1:NVC])
                            }
                            else
                            {
                                varcov <- matrix(0,nrow=nea0,ncol=nea0)
                                varcov[upper.tri(varcov,diag=TRUE)] <- b[NPROB+NEF+1:NVC]
                                varcov <- t(varcov)
                                varcov[upper.tri(varcov,diag=TRUE)] <- b[NPROB+NEF+1:NVC]
                                
                                ch <- chol(varcov)
                                
                                b[NPROB+NEF+1:NVC] <- ch[upper.tri(ch,diag=TRUE)]
                                
                            }
                        }                   
                        
                    }
                }
            else
                {
                    if(!inherits(B,"lcmm")) stop("B should be either a vector or an object of class lcmm")
                         
                    if(ng>1 & B$ng==1)
                        {
                            NEF2 <- sum(idg0!=0)-1
                            NPM2 <- NEF2+NVC+ntrtot0
                            if(length(B$best)!=NPM2) stop("B is not correct")
                                

                            if(!length(B$Brandom))
                                {
                                    ## B deterministe
                                    k <- NPROB
                                    l <- 0
                                    t<- 0
                                    for (i in 1:nvar.exp)    {
                                        if(idg0[i]==1 & i>1){
                                            l <- l+1
                                            t <- t+1
                                            b[k+t] <- B$best[l]
                                        }
                                        if(idg0[i]==2){
                                            if (i==1){
                                                for (g in 2:ng){
                                                    t <- t+1
                                                    b[k+t] <- - 0.5*(g-1)
                                                }
                                            }
                                            if (i>1){
                                                l <- l+1
                                                for (g in 1:ng){
                                                    t <- t+1
                                                    if(B$conv==1) b[k+t] <- B$best[l]+(g-(ng+1)/2)*sqrt(B$V[l*(l+1)/2])
                                                    else  b[k+t] <- B$best[l]+(g-(ng+1)/2)*B$best[l]
                                                }
                                            }
                                        }
                                    }

                                    
                                    ## pr nvc
                                    if(NVC>0)
                                        {
                                            if(idiag==TRUE)
                                                {
                                                    b[(NPROB+NEF+1):(NPROB+NEF+NVC)] <-B$cholesky[(1:nea0)*(2:(nea0+1))/2]    
                                                }
                                            else
                                                {
                                                    b[(NPROB+NEF+1):(NPROB+NEF+NVC)] <-B$cholesky
                                                }
                                        }

                                    ## cor et transfo
                                    b[(NPROB+NEF+NVC+NW+1):NPM] <- B$best[(NEF2+NVC+1):NPM2]
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
                                    nbgnef <- nbgnef[-1]
                                    nbg <- nbg[-1]

                                    vbb[which(nbgnef==1),setdiff(1:ncol(vbb),which(nbgnef!=1))] <- VB[which(nbg==1),setdiff(1:ncol(VB),which(nbg!=1))]
                                    vbb[(NEF+1):nrow(vbb),(NEF+1):ncol(vbb)] <- VB[(NEF2+1):nrow(VB),(NEF2+1):ncol(VB)]
                                    
                                    
                                    t <- 0
                                    l <- 0
                                    for (i in 1:nvar.exp)
                                        {
                                            if(idg0[i]==1)
                                                {
                                                    if(i==1) next
                                                    l <- l+1
                                                    t <- t+1
                                                    bb[t] <- B$best[l]
                                                }
                                            if(idg0[i]==2)
                                                {
                                                    if(i==1)
                                                        {
                                                            t <- t+ng-1
                                                            next
                                                        }
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

                                        #transfo
                                    bb[NEF+NVC+1:ntrtot0] <- B$best[NEF2+NVC+1:ntrtot0]
                                   
                                    
                                    
                                    if(idg0[1]>1)
                                        {
                                            bb <- bb[-(1:(ng-1))]
                                            vbb <- vbb[-(1:(ng-1)),-(1:(ng-1))]
                                        }
                        
                                    up <- vbb[upper.tri(vbb,diag=TRUE)]
                                    vbb <- t(vbb)
                                    vbb[upper.tri(vbb,diag=TRUE)] <- up
                                    ##Chol <- chol(vbb)
                                    ##Chol <- t(Chol)

                                    if(idg0[1]>1)
                                        {
                                            b[c((NPROB+ng):(NPROB+NEF+NVC),(NPROB+NEF+NVC+NW+1):NPM)] <- rmvnorm(n=1,mean=bb,sigma = vbb) #bb + Chol %*% rnorm(length(bb))
                                            b[NPROB+1:(ng-1)] <- 0
                                        } 
                                    else
                                        {
                                            b[c((NPROB+1):(NPROB+NEF+NVC),(NPROB+NEF+NVC+NW+1):NPM)] <- rmvnorm(n=1,mean=bb,sigma = vbb) #bb + Chol %*% rnorm(length(bb))
                                        }

                                    if(NPROB>0) b[1:NPROB] <- 0
                                    if(NW>0) b[NPROB+NEF+NVC+1:NW] <- 1
                                                                       
                                }
                        }
                }
        }


            ##------------------------------------------
            ##------nom au vecteur best
            ##--------------------------------------------

            

    if(ng0>=2 & NPROB>0){
        nom <-rep(nom.X0[idprob0!=0],each=ng0-1)
        nom1 <- paste(nom," class",c(1:(ng0-1)),sep="")
        names(b)[1:NPROB]<-nom1
    }



    if(ng0==1) names(b)[1:(NEF)] <-  nom.X0[-1][idg0[-1]!=0]

    if(ng0>1){
	nom1<- NULL
	for (i in 1:nvar.exp) {
            if(idg0[i]==2){ 
                if (i==1){
                    nom <- paste(nom.X0[i]," class",c(2:ng0),sep="")
                    nom1 <- cbind(nom1,t(nom))
                }
                if (i>1){
                    nom <- paste(nom.X0[i]," class",c(1:ng0),sep="")
                    nom1 <- cbind(nom1,t(nom))
                }
            }
            if(idg0[i]==1 & i>1) nom1 <- cbind(nom1,nom.X0[i])
	}
        names(b)[(NPROB+1):(NPROB+NEF)]<- nom1
    }

    names(b)[(NPM-ntrtot0+1):NPM]<- paste("thresh. parm",c(1:ntrtot0),sep="")



    if(NVC!=0)names(b)[(NPROB+NEF+1):(NPROB+NEF+NVC)] <- paste("varcov",c(1:(NVC)))
    if(NW!=0)names(b)[(NPROB+NEF+NVC+1):(NPROB+NEF+NVC+NW)] <- paste("varprop class",c(1:(ng0-1)))
    names_best <- names(b)

    N <- NULL
    N[1] <- NPROB
    N[2] <- NEF
    N[3] <- NVC
    N[4] <- NW
    N[5] <- nobs0
    N[6] <- 0 #=ncor

    idiag <- as.integer(idiag0)
    idea <- as.integer(idea0)
    nv <- as.integer(nv0)


################ Sortie ###########################


    marker <- rep(0,(maxY-minY+1)*2)
    transfY <- rep(0,(maxY-minY+1)*2)


    
    nfix <- sum(fix0)
    bfix <- 0
    if(nfix>0)
    {
        bfix <- b[which(fix0==1)]
        b <- b[which(fix0==0)]
        NPM <- NPM-nfix
    }

    nbzitr0 <- 0
    zitr0 <- 0
    epsY0 <- 0
    idlink0 <- 3
    idcor0 <- 0
    ncor0 <- 0
    nsim <- 2*(maxY-minY+1)
    
    if(maxiter==0)
    {
        vrais <- logliklcmm(b,Y0,X0,prior0,pprior0,idprob0,idea0,idg0,idcor0,
                            ns0,ng0,nv0,nobs0,nea0,nmes0,idiag0,nwg0,ncor0,
                            NPM,epsY0,idlink0,nbzitr0,zitr0,minY,maxY,ide,
                            fix0,nfix,bfix)
        
        out <- list(conv=2, V=rep(0, length(b)), best=b,
                    ppi2=rep(NA,ns0*ng0), predRE=rep(NA,ns0*nea0*(1+ng0)),
                    Yobs=rep(NA,nobs0),
                    resid_m=rep(NA,nobs0), resid_ss=rep(NA,nobs0),
                    marker=rep(NA,nsim), transfY=rep(NA,nsim),
                    pred_m_g=rep(NA,nobs0*ng0), pred_ss_g=rep(NA,nobs0*ng0),
                    vraisdiscret=NA, UACV=NA, rlindiv=rep(NA,ns0),
                    gconv=rep(NA,3), niter=0, loglik=vrais)
    }
    else
    {   
        res <- mla(b=b, m=length(b), fn=logliklcmm,
                   clustertype=clustertype,.packages=NULL,
                   epsa=convB,epsb=convL,epsd=convG,partialH=indexHr,
                   digits=8,print.info=verbose,blinding=FALSE,
                   multipleTry=25,file="",
                   nproc=nproc,maxiter=maxiter,minimize=FALSE,
                   Y0=Y0,X0=X0,prior0=prior0,pprior0=pprior0,idprob0=idprob0,idea0=idea0,
                   idg0=idg0,idcor0=idcor0,ns0=ns0,ng0=ng0,nv0=nv0,nobs=nobs0,
                   nea0=nea0,nmes0=nmes0,idiag=idiag0,nwg0=nwg0,ncor0=ncor0,
                   npm0=NPM,epsY0=epsY0,idlink0=idlink0,nbzitr0=nbzitr0,zitr0=zitr0,
                   minY0=minY,maxY0=maxY,ide0=ide,
                   fix0=fix0,nfix0=nfix,bfix0=bfix)
        
        out <- list(conv=res$istop, V=res$v, best=res$b,
                    ppi2=rep(NA,ns0*ng0), predRE=rep(NA,ns0*nea0*(1+ng0)), 
                    Yobs=rep(NA,nobs0),
                    resid_m=rep(NA,nobs0), resid_ss=rep(NA,nobs0),
                    marker=rep(NA,nsim), transfY=rep(NA,nsim),
                    pred_m_g=rep(NA,nobs0*ng0), pred_ss_g=rep(NA,nobs0*ng0),
                    vraisdiscret=NA, UACV=NA, rlindiv=rep(NA,ns0),
                    gconv=c(res$ca, res$cb, res$rdm), niter=res$ni,
                    loglik=res$fn.value)
    }
    
    if(out$conv %in% c(1,2,3))
    {
        estim0 <- 0
        ll <- 0
        ppi0 <- rep(0,ns0*ng0)
        resid_m <- rep(0,nobs0)
        resid_ss <- rep(0,nobs0)
        pred_m_g <- rep(0,nobs0*ng0)
        pred_ss_g <- rep(0,nobs0*ng0)
        predRE <- rep(0,ns0*nea0*(1+ng0))
        marker <- rep(0,nsim)
        transfY <- rep(0,nsim)
        Yobs <- rep(0,nobs0)
        Ydiscret <- 0
        vraisdiscret <- 0
        UACV <- 0
        rlindiv <- rep(0,ns0)
        post <- .Fortran(C_logliklcmmord,
                         as.double(Y0),
                         as.double(X0),
                         as.integer(prior0),
                         as.double(pprior0),
                         as.integer(idprob0),
                         as.integer(idea0),
                         as.integer(idg0),
                         as.integer(ns0),
                         as.integer(ng0),
                         as.integer(nv0),
                         as.integer(nobs0),
                         as.integer(nea0),
                         as.integer(nmes0),
                         as.integer(idiag0),
                         as.integer(nwg0),
                         as.integer(NPM),
                         as.double(out$best),
                         ppi=as.double(ppi0),
                         resid_m=as.double(resid_m),
                         resid_ss=as.double(resid_ss),
                         pred_m_g=as.double(pred_m_g),
                         pred_ss_g=as.double(pred_ss_g),
                         predRE=as.double(predRE),
                         as.integer(minY),
                         as.integer(maxY),
                         as.integer(ide),
                         marker=as.double(marker),
                         transfY=as.double(transfY),
                         UACV=as.double(UACV),
                         rlindiv=as.double(rlindiv),
                         as.double(out$V),
                         as.integer(fix0),
                         as.integer(nfix),
                         as.double(bfix),
                         as.integer(estim0),
                         loglik=as.double(ll))
        
        out$ppi2 <- post$ppi
        out$predRE <- post$predRE
        out$Yobs <- post$Yobs
        out$resid_m <- post$resid_m
        out$resid_ss <- post$resid_ss
        out$marker <- post$marker
        out$transfY <- post$transfY
        out$pred_m_g <- post$pred_m_g
        out$pred_ss_g <- post$pred_ss_g
        out$UACV <- post$UACV
        out$rlindiv <- post$rlindiv
    }
        
    
    
    ## creer best a partir de b et bfix
    best <- rep(NA,length(fix0))
    best[which(fix0==0)] <- out$best
    best[which(fix0==1)] <- bfix
    names(best) <- names_best
    out$best <- best
    NPM <- NPM+nfix
        
    
    

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


    if(ng0>1) {
        ppi<- matrix(out$ppi2,ncol=ng0,byrow=TRUE)
    }
    else {
        ppi <- matrix(rep(1,ns0),ncol=ng0)
    }

    chooseClass <- function(ppi)
    {
        res <- which.max(ppi)
        if(!length(res)) res <- NA
        return(res)
    }
    classif<-apply(ppi,1,chooseClass)
    ppi<-data.frame(INDuniq,classif,ppi)
    temp<-paste("prob",1:ng0,sep="")
    colnames(ppi) <- c(nom.subject,"class",temp)
    rownames(ppi) <- 1:ns0
                                 

                                       
    estimlink <- cbind(out$marker,out$transfY)
    colnames(estimlink) <- c("Y","transfY")

### ad 2/04/2012
    if (!("intercept" %in% nom.X0)) X0.names2 <- X0.names2[-1]
### ad

    ## levels = modalites des variables dans X0 (si facteurs)
    levelsdata <- vector("list", length(X0.names2))
    levelsfixed <- vector("list", length(X0.names2))
    levelsrandom <- vector("list", length(X0.names2))
    levelsmixture <- vector("list", length(X0.names2))
    levelsclassmb <- vector("list", length(X0.names2))
    names(levelsdata) <- X0.names2
    names(levelsfixed) <- X0.names2
    names(levelsrandom) <- X0.names2
    names(levelsmixture) <- X0.names2
    names(levelsclassmb) <- X0.names2
    for(v in X0.names2)
    {
        if(v == "intercept") next
        
        if(is.factor(data[,v]))
        {
            levelsdata[[v]] <- levels(data[,v])
        }
        if(length(grep(paste("factor\\(",v,"\\)",sep=""), fixed)))
        {
            levelsfixed[[v]] <- levels(as.factor(data[,v]))
        }
        if(length(grep(paste("factor\\(",v,"\\)",sep=""), random)))
        {
            levelsrandom[[v]] <- levels(as.factor(data[,v]))
        }
        if(length(grep(paste("factor\\(",v,"\\)",sep=""), mixture)))
        {
            levelsmixture[[v]] <- levels(as.factor(data[,v]))
        }
        if(length(grep(paste("factor\\(",v,"\\)",sep=""), classmb)))
        {
            levelsclassmb[[v]] <- levels(as.factor(data[,v]))
        }
        
    }
    levels <- list(levelsdata=levelsdata, levelsfixed=levelsfixed, levelsrandom=levelsrandom,
                   levelsmixture=levelsmixture, levelsclassmb=levelsclassmb)
    
    cost<-proc.time()-ptm
    
    res <-list(ns=ns0,ng=ng0,idea0=idea0,idprob0=idprob0,idg0=idg0,idcor0=rep(0,nvar.exp),loglik=out$loglik,best=out$best,V=V,gconv=out$gconv,conv=out$conv,call=call,niter=out$niter,N=N,idiag=idiag0,pprob=ppi,Xnames=nom.X0,Xnames2=X0.names2,cholesky=Cholesky,estimlink=estimlink,epsY=0,linktype=3,linknodes=zitr,ide=ide,Ydiscrete=Ydiscrete,discrete_loglik=out$loglik,UACV=out$UACV,IndivContrib=out$rlindiv,na.action=na.action,AIC=2*(length(out$best)-length(posfix)-out$loglik),BIC=(length(out$best)-length(posfix))*log(ns0)-2*out$loglik,data=datareturn, levels=levels, var.time=var.time, runtime=cost[3])
    class(res) <-c("lcmm")  

    if(verbose==TRUE) cat("The program took", round(cost[3],2), "seconds \n")
    res
}


