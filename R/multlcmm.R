
multlcmm <- function(fixed,mixture,random,subject,classmb,ng=1,idiag=FALSE,nwg=FALSE,randomY=FALSE,link="linear",intnodes=NULL,epsY=0.5,cor=NULL,data,B,convB=0.0001,convL=0.0001,convG=0.0001,maxiter=100,nsim=100,prior,range=NULL,subset=NULL,na.action=1,posfix=NULL,partialH=FALSE,verbose=TRUE)
{
    ptm<-proc.time()
    if(verbose==TRUE) cat("Be patient, multlcmm is running ... \n")

    cl <- match.call()

    nom.subject <- as.character(subject)

#### INCLUSION PRIOR
    nom.prior <- NULL
    if(!missing(prior)) nom.prior <- as.character(prior)
####

    if(!missing(mixture) & ng==1) stop("No mixture can be specified with ng=1")
    if(missing(mixture) & ng>1) stop("The argument mixture has to be specified for ng > 1")
    if(!missing(classmb) & ng==1) stop("No classmb can be specified with ng=1")
    if(missing(random)) stop("At least one random effect is required")
    if(random==~-1) stop("At least one random effect is required")
    if(missing(fixed)) stop("The argument Fixed must be specified in any model")
    if(missing(classmb) & ng==1) classmb <- ~-1
    if(missing(classmb) & ng>1) classmb <- ~1
    if(missing(mixture)) mixture <- ~-1
    if(ng==1&nwg==TRUE) stop ("The argument nwg should be FALSE for ng=1")


    if(class(fixed)!="formula") stop("The argument fixed must be a formula")
    if(class(mixture)!="formula") stop("The argument mixture must be a formula")
    if(class(random)!="formula") stop("The argument random must be a formula")
    if(class(classmb)!="formula") stop("The argument classmb must be a formula")
    if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
    if(nrow(data)==0) stop("Data should not be empty")
    if(missing(subject)){ stop("The argument subject must be specified")}
    if(any(link=="thresholds"))  stop("The link function thresholds is not available in multivariate case")
    if(all(link %in% c("linear","beta")) & !is.null(intnodes)) stop("Intnodes should only be specified with splines links")

    if(!(na.action%in%c(1,2)))stop("only 1 for 'na.omit' or 2 for 'na.fail' are required in na.action argument")

    if(length(posfix) & missing(B)) stop("A set of initial parameters must be specified if some parameters are not estimated")

### test de l'argument cor
    ncor0 <- 0
    cor.type <- cl$cor[1]
    cor.time <- cl$cor[2]
    cor <- paste(cor.type,cor.time,sep="-")
    if (!isTRUE(all.equal(cor,character(0))))
        {
            if (substr(cor,1,2)=="AR") { ncor0 <- 2 }
            else if (substr(cor,1,2)=="BM") { ncor0 <- 1  }
            else { stop("The argument cor must be of type AR or BM") }

            if(!(strsplit(cor,"-")[[1]][2] %in% colnames(data))) stop("Unable to find time variable from argument 'cor' in 'data'")
            else { cor.var.time <- strsplit(cor,"-")[[1]][2] }
        }
### fin test argument cor

    ##pour acces aux attributs des formules
    afixed <- terms(fixed, specials=c("factor","contrast"))
    if(attr(afixed,"intercept")==0) stop("An intercept should appear in fixed for identifiability purposes")
    amixture <- terms(mixture, specials=c("factor"))
                                        
    arandom <- terms(random, specials=c("factor"))
    aclassmb <- terms(classmb, specials=c("factor"))
    ##fixed sans contrast
    fixed2 <- gsub("contrast","",fixed)
    fixed2 <- formula(paste(fixed2[2],fixed2[3],sep="~"))   
    afixed2 <- terms(fixed2)
                                        
    ##verifier si totes les varialbes sont dans data
    variables <- c(attr(afixed,"variables"),attr(arandom,"variables"),attr(amixture,"variables"),attr(aclassmb,"variables"))
    variables <- unlist(lapply(variables,all.vars))  
    if(!all(variables %in% colnames(data))) stop(paste("Data should contain the variables",paste(unique(variables),collapse=" ")))


    ##contrast
    contr <- ~-1
    if(!is.null(attr(afixed,"specials")$contrast))
        {
            vcontr <- attr(afixed,"term.labels")[setdiff(attr(afixed,"specials")$contrast-1,untangle.specials(afixed,"contrast",2)$terms)]
            vcontr <- gsub("contrast","",vcontr)
            contr <- as.formula(paste("~-1+",paste(vcontr,collapse="+")))
        }
    acontr <- terms(contr)
    
    ##tjrs intercept dans classmb
    if(attr(aclassmb,"intercept")==0 & ng>1)
        {
            attr(aclassmb,"intercept") <- 1
            cat("The formula in classmb should always include an intercept. An intercept has been added.")
        }

###liste des outcomes
    nomsY <- as.character(attr(afixed,"variables")[2])
    nomsY <- strsplit(nomsY,split=" + ",fixed=TRUE)
    nomsY <- as.vector(nomsY[[1]])
    ny0 <- length(nomsY)

    ##pas de contrast ni randomY si un seul Y
    if(ny0<2 & length(attr(afixed,"specials")$contrast)) stop("No contrast can be included with less than two outcomes")
    if(ny0<2 & randomY==TRUE) stop("With less than 2 outcomes randomY should be FALSE")

###liste des variables utilisees  (sans les interactions et sans les Y)
    ttesLesVar <- colnames(get_all_vars(afixed,data=data[1,]))
    ttesLesVar <- c(ttesLesVar, colnames(get_all_vars(amixture,data=data[1,])))
    ttesLesVar <- c(ttesLesVar, colnames(get_all_vars(arandom,data=data[1,])))
    ttesLesVar <- c(ttesLesVar, colnames(get_all_vars(aclassmb,data=data[1,])))
    if (ncor0>0) ttesLesVar <- unique(c(ttesLesVar,cor.var.time))
    else ttesLesVar <- unique(ttesLesVar)
    ttesLesVar <- setdiff(ttesLesVar, nomsY)

### argument subset
    form1 <- paste(c(nom.subject,nomsY,ttesLesVar,nom.prior),collapse="+")
    if(!isTRUE(all.equal(as.character(cl$subset),character(0))))
        {
            cc <- cl
            cc <- cc[c(1,which(names(cl)=="subset"))]
            cc[[1]] <- as.name("model.frame")
            cc$formula <- formula(paste("~",form1))
            cc$data <- data
            cc$na.action <- na.pass
            data <- eval(cc)
        }

    attributes(data)$terms <- NULL

### si subject est un factor
    if(is.factor(data[,nom.subject]))
        {
            data[,nom.subject] <- as.numeric(data[,nom.subject])
        }

    
###subset de data avec les variables utilisees
    
    newdata <- data[,c(nom.subject,nomsY,ttesLesVar,nom.prior)]
    if(!is.null(nom.prior))
        {
            prior <- newdata[,nom.prior]
            newdata[which(is.na(prior)),nom.prior] <- 0
        }

###un data frame par outcome et creation Y0
    dataY <- paste("data",nomsY,sep=".")
    Y0 <- NULL
    IND <- NULL
    outcome <- NULL
    prior <- NULL
    data0 <- NULL
    nayk <- vector("list",ny0)
    for (k in 1:ny0)
        {
            dtemp <- newdata[,c(nom.subject,nomsY[k],ttesLesVar,nom.prior)]
            ##enlever les NA
            linesNA <- apply(dtemp,2,function(v) which(is.na(v)))
            linesNA <- unique(unlist(linesNA))
            if(length(linesNA)) nayk[[k]] <- linesNA
            if(na.action==1 & length(linesNA)>0) dtemp <- dtemp[-linesNA,]
            if(na.action==2 & length(linesNA)>0) stop("Data contains missing values")
            assign(dataY[k],dtemp)
            Y0 <- c(Y0, dtemp[,nomsY[k]])
            IND <- c(IND, dtemp[,nom.subject])
            outcome <- c(outcome,rep(nomsY[k],nrow(dtemp)))
            if(!is.null(nom.prior)) prior <- c(prior, dtemp[,nom.prior])
            data0 <- rbind(data0, dtemp[,setdiff(colnames(dtemp),nomsY[k]),drop=FALSE])   #dataset sans NA avec les covariables utilisees; obs ordonnees par outcome
        }

###prior=0 si pas specifie
    if(is.null(prior)) prior <- rep(0,length(Y0))

###creation de X0 (ttes les var + interactions)


    Xfixed <- model.matrix(fixed2[-2], data=data0)
    Xmixture <- model.matrix(mixture, data=data0)
    Xrandom <- model.matrix(random, data=data0)
    Xclassmb <- model.matrix(classmb, data=data0)
    Xcontr <- model.matrix(contr,data=data0)

    z.fixed <- strsplit(colnames(Xfixed),split=":",fixed=TRUE)
    z.fixed <- lapply(z.fixed,sort)

    z.random <- strsplit(colnames(Xrandom),split=":",fixed=TRUE)
    z.random <- lapply(z.random,sort)
    
    if(mixture != ~-1)
        {
            z.mixture <- strsplit(colnames(Xmixture),split=":",fixed=TRUE)
            z.mixture <- lapply(z.mixture,sort)
        }
    else
        {
            z.mixture <- list()
        }
    
    if(classmb != ~-1)
        {
            z.classmb <- strsplit(colnames(Xclassmb),split=":",fixed=TRUE)
            z.classmb <- lapply(z.classmb,sort)
        }
    else
        {
            z.classmb <- list()
        }
    
    if(contr != ~-1)
        {
            z.contr <- strsplit(colnames(Xcontr),split=":",fixed=TRUE)
            z.contr <- lapply(z.contr,sort)
        }
    else
        {
            z.contr <- list()
        }
    
    if(!all(z.mixture %in% z.fixed))  stop("The covariates in mixture should also be included in the argument fixed")
    if(!all(z.contr %in% z.fixed))  stop("The covariates in contrast should also appear in fixed")
    
    X0 <- cbind(Xfixed, Xrandom, Xclassmb)        
    nom.unique <- unique(colnames(X0))
    X0 <- X0[,nom.unique,drop=FALSE]
    
    if (ncor0>0)
        {
            if(!(cor.var.time %in% colnames(X0)))
                {
                    X0 <- cbind(X0, data0[,cor.var.time])
                    colnames(X0) <- c(nom.unique, cor.var.time)
                    nom.unique <- c(nom.unique,cor.var.time)
                }
        }

    X0 <- as.matrix(X0)
###X0 fini


###test de link
    if (length(link)!=1 & length(link)!=ny0) stop("One link per outcome should be specified")
    if(any(link %in% c("splines","Splines")))
        {
            link[which(link %in% c("splines","Splines"))] <- "5-quant-splines"
        }
    if(length(link)==1 & ny0>1)
        {
            link <- rep(link, ny0)
        }

    idlink0 <- rep(2,ny0)
    idlink0[which(link=="linear")] <- 0
    idlink0[which(link=="beta")] <- 1

    
    spl <- strsplit(link[which(idlink0==2)],"-")
    if(any(sapply(spl,length)!=3)) stop("Invalid argument 'link'")

    nySPL <- length(spl)
    nybeta <- sum(idlink0==1)
    ##remplir range si pas specifie
    if(!is.null(range) & length(range)!=2*(nySPL+nybeta)) stop("Length of vector range is not correct.")
    if((length(range)==2*(nySPL+nybeta)) & (nySPL+nybeta>0))
        {
            ind12 <- which(idlink0==1 | idlink0==2)
            for (i in 1:(nySPL+nybeta))
                {
                    rg <- range(get(dataY[ind12[i]])[,nomsY[ind12[i]]])
                    if(rg[1]<range[2*(i-1)+1] | rg[2]>range[2*(i-1)+2]) stop("The range specified do not cover the entire range of the data")
                }
        }
    if((is.null(range) & (nybeta+nySPL)>0) | length(range)!=2*(nySPL+nybeta))
        {
            range <- NULL
            for(k in which(idlink0!=0))
                {
                    min1 <- min(get(dataY[k])[,nomsY[k]])
                    min2 <- round(min1,3)
                    if(min1<min2) min2 <- min2-0.001

                    max1 <- max(get(dataY[k])[,nomsY[k]])
                    max2 <- round(max1,3)
                    if(max1>max2) max2 <- max2+0.001
                    
                    range <- c(range, min2, max2)
                }
        }


    ## epsY
    if (any(idlink0==1))
        {
            if (any(epsY<=0))
                {
                    stop("Argument 'epsY' should be positive.")
                }

            if(length(epsY)==1) epsY <- rep(epsY,nybeta)
            
            if(length(epsY)!=nybeta) stop(paste("Argument 'epsY' should be of length",nybeta))
            if(nybeta!=ny0)
                {
                    epsY2 <- rep(0,ny0)
                    epsY2[which(idlink0==1)] <- epsY
                    epsY <- epsY2
                }

        } 
    

    nbzitr0 <- rep(2,ny0) #nbzitr0 = nb de noeuds si splines, 2 sinon
    nbnodes <- NULL  #que pour les splines
    spltype <- NULL
    if(nySPL>0)
        {
            for (i in 1:nySPL)
                {
                    nbnodes <- c(nbnodes, spl[[i]][1])
                    spltype <- c(spltype, spl[[i]][2])
                    if(spl[[i]][3] != "splines") stop("Invalid argument link")
                }
        }
    nbnodes <- as.numeric(nbnodes)
    nbzitr0[which(idlink0==2)] <- nbnodes

    ##test splines
    if(!(all(spltype %in% c("equi","quant","manual")))) stop("The location of the nodes should be 'equi', 'quant' or 'manual'")

    ##tester longueur de intnodes
    if(!is.null(intnodes))
        {  
            if(length(intnodes) != sum(nbnodes[which(spltype=="manual")]-2)) stop(paste("Vector intnodes should be of length",sum(nbnodes[which(spltype=="manual")]-2)))
        }

    ##intnodes2 : contient tous les noeuds interieurs (pas seulement ceux de manual)
    intnodes2 <- rep(NA,sum(nbnodes-2))
    nb <- 0
    nbspl <- 0
    for (k in 1:ny0)
        {
            if (idlink0[k]!=2) next
            else
                {                                         
                    nbspl <- nbspl+1

                    if(spltype[nbspl]=="manual")
                        {
                            nodes <- intnodes[(nb+1):(nb+nbnodes[nbspl]-2)]
                            if(!length(nodes)) stop("The length of intnodes is not correct")
                            intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <-  nodes
                            nb <- nb+nbnodes[nbspl]-2

                            idrg <- length(which(idlink0[1:k] != 0))
                            if(any(nodes <= range[2*(idrg-1)+1]) | any(nodes >= range[2*idrg])) stop("Interior nodes must be in the range of the outcome")
                        }

                    if(spltype[nbspl]=="equi")
                        {
                            nodes <- seq(range[2*(nbspl-1)+1], range[2*nbspl], length.out=nbnodes[nbspl])
                            nodes <- nodes[-nbnodes[nbspl]]
                            nodes <- nodes[-1]
                            intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <- nodes
                        }

                    if(spltype[nbspl]=="quant")
                        {
                            nodes <- quantile(get(dataY[k])[,nomsY[k]], probs=seq(0,1,length.out=nbnodes[nbspl]))
                            if(length(unique(nodes)) != length(nodes)) stop(paste("Some nodes are equal for link number",k,"; Please try to reduce the number of nodes or use manual location."))
                            nodes <- nodes[-nbnodes[nbspl]]
                            nodes <- nodes[-1]
                            intnodes2[(sum(nbnodes[1:nbspl]-2)-(nbnodes[nbspl]-2)+1):sum(nbnodes[1:nbspl]-2)] <- as.vector(nodes)
                        }
                }
        }

    if(nb != length(intnodes)) stop(paste("The vector intnodes should be of length",nb))

    ##remplir zitr
    m <- 0
    if(nySPL>0) m <- max(nbnodes)
    zitr <- matrix(0,max(m,2),ny0)
    nb12 <- 0
    nbspl <- 0
    for (k in 1:ny0)
        {
            if(idlink0[k]==0) zitr[1:2,k] <- c(min(get(dataY[k])[,nomsY[k]]),max(get(dataY[k])[,nomsY[k]]))

            if(idlink0[k]==1)
                {
                    nb12 <- nb12 + 1
                    zitr[1:2,k] <- range[2*(nb12-1)+1:2]
                }

            if(idlink0[k]==2)
                {
                    nb12 <- nb12+1
                    nbspl <- nbspl+1
                    zitr[2:(nbzitr0[k]-1),k] <- intnodes2[ifelse(nbspl==1,0,1)*sum(nbnodes[1:(nbspl-1)]-2) + 1:(nbnodes[nbspl]-2)]
                    zitr[1,k] <- range[2*(nb12-1)+1]
                    zitr[nbnodes[nbspl],k]  <- range[2*nb12]
                    
                    ##verifier s'il y a des obs entre les noeuds
                    hcounts <- hist(get(dataY[k])[,nomsY[k]],breaks=zitr[1:nbnodes[nbspl],k],plot=FALSE,include.lowest=TRUE,right=TRUE)$counts
                    if(any(hcounts==0)) stop(paste("Link function number",k,"can not be estimated. Please try other nodes such that there are observations in each interval."))    
                }
        }

###uniqueY0 et indiceY0
    uniqueY0 <- NULL
    indiceY0 <- NULL
    nvalSPL0 <- NULL
    nb <- 0
    for (k in 1:ny0)
        {
            if(idlink0[k]!=2)
                {
                    indiceY0 <- c(indiceY0, rep(0,length(get(dataY[k])[,nomsY[k]])))
                    next
                }

            yk <- get(dataY[k])[,nomsY[k]]
            uniqueTemp <- sort(unique(yk))
            permut <- order(order(yk))  # sort(y)[order(order(y))] = y
            if(length(as.vector(table(yk)))==length(uniqueTemp))
                {
                    indice <- rep(1:length(uniqueTemp), as.vector(table(yk)))
                    indiceTemp <- nb + indice[permut]

                    nb <- nb + length(uniqueTemp)
                    uniqueY0 <- c(uniqueY0, uniqueTemp)
                    indiceY0 <- c(indiceY0, indiceTemp)
                    nvalSPL0 <- c(nvalSPL0, length(uniqueTemp))
                }
            else
                {
                    uniqueY0 <- yk
                    indiceY0 <- c(1:length(yk))
                    nvalSPL0 <- length(yk)
                }
        }


###ordonner les mesures par individu
    #IDnum <- as.numeric(IND)
    matYX <- cbind(IND,prior,Y0,indiceY0,outcome,X0)
    matYXord <- matYX[order(IND),]
    Y0 <- as.numeric(matYXord[,3])
    X0 <- apply(matYXord[,-c(1,2,3,4,5),drop=FALSE],2,as.numeric)
                                        #X0 <- as.matrix(X0)  a remettre si X0 <- as.data.frame(X0) remis l.211
    IND <- matYXord[,1]
    outcome <- matYXord[,5]
    indiceY0 <- as.numeric(matYXord[,4])
    prior0 <- as.numeric(matYXord[,2])


###parametres pour hetmixContMult
    ns0 <- length(unique(IND))
    ng0 <- ng
    nv0 <- dim(X0)[2]
    nobs0 <- length(Y0)
    idiag0 <- ifelse(idiag==TRUE,1,0)
    nwg0 <- ifelse(nwg==TRUE,1,0)
    nalea0 <- ifelse(randomY==TRUE,ny0,0)

    loglik <- 0
    ni <- 0
    istop <- 0
    gconv <- rep(0,3)
    ppi0 <- rep(0,ns0*ng0)
    resid_m <- rep(0,nobs0)
    resid_ss <- rep(0,nobs0)
    pred_m_g <- rep(0,nobs0*ng0)
    pred_ss_g <- rep(0,nobs0*ng0)
    Yobs <- rep(0,nobs0)
                                        
    predRE_Y <- rep(0,ns0*nalea0)
    rlindiv <- rep(0,ns0)
    marker <- rep(0,nsim*ny0)
    transfY <- rep(0,nsim*ny0)
    Ydiscrete <- 0
    UACV <- 0
    vraisdiscret <- 0


###nmes0
    nmes0 <- matrix(0,ns0,ny0)
    for (k in 1:ny0)
        {
            INDpresents <- which(unique(IND) %in% get(dataY[k])[,nom.subject])
            nmes0[INDpresents,k] <- as.vector(table(get(dataY[k])[,nom.subject]))
        }


###remplir idprob, etc
    z.X0 <- strsplit(nom.unique,split=":",fixed=TRUE)
    z.X0 <- lapply(z.X0,sort)
    
    idprob0 <- z.X0 %in% z.classmb + 0
    idea0 <- z.X0 %in% z.random + 0
    idg0 <- (z.X0 %in% z.fixed) + (z.X0 %in% z.mixture)
    idcontr0 <- z.X0 %in% z.contr + 0
    

    if (ncor0>0) idcor0 <- colnames(X0) %in% cor.var.time +0
    else idcor0 <- rep(0,nv0)


    nea0 <- sum(idea0)
    predRE <- rep(0,ns0*nea0)
    
    ##nombre total de parametres
    NPM <- (ng0-1)*sum(idprob0) + sum(idg0==1)-1 + ng0*sum(idg0==2) + ncor0 + (ny0-1)*sum(idcontr0) +
        ifelse(idiag0==1,nea0,nea0*(nea0+1)/2)-1 + (ng0-1)*nwg0 + nalea0 + ny0 +
            2*sum(idlink0==0) + 4*sum(idlink0==1) + sum(nbnodes+2)

    V <- rep(0, NPM*(NPM+1)/2)  #pr variance des parametres

    nef <- (ng0-1)*sum(idprob0) + sum(idg0==1)-1 + ng0*sum(idg0==2) + (ny0-1)*sum(idcontr0)
    ncontr <- (ny0-1)*sum(idcontr0)
    nvc <- ifelse(idiag0==1,nea0,nea0*(nea0+1)/2)-1
    nw <- (ng0-1)*nwg0
    ntrtot0 <- nbzitr0+2
    ntrtot0[which(idlink0==0)] <- 2
    nprob <- sum(idprob0)*(ng0-1)

    
## gestion de B=random(mod)

        Brandom <- FALSE
        if(length(cl$B)==2)
            {
                if(class(eval(cl$B[[2]]))!="multlcmm") stop("The model specified in B should be of class multlcmm")
                if(as.character(cl$B[1])!="random") stop("Please use random() to specify random initial values")
                
                Brandom <- TRUE
                B <- eval(cl$B[[2]])

                if(length(posfix)) stop("Argument posfix is not compatible with random intial values")
            }
    
###valeurs initiales
    if(!(missing(B)))
        {
            if(is.vector(B))
                {
                    if (length(B)==NPM) b <- B
                    else stop(paste("Vector B should be of length",NPM))
                }
            else
                {
                    if(class(B)!="multlcmm") stop("B should be either a vector or an object of class multlcmm")

                    if(ng==1 & B$ng==1)
                        {
                            if(length(B$best)!=NPM) stop("B is not correct")
                            b <- B$best
                        }

                    if(ng>1 & B$ng==1)
                            {
                                nef2 <- sum(idg0!=0)-1 + (ny0-1)*sum(idcontr0)
                                NPM2 <- nef2+ nvc+ncor0+nalea0+ny0+sum(ntrtot0)
                              
                                if(length(B$best)!=NPM2) stop("B is not correct")

                                if(Brandom==FALSE)
                                    {
                                        ### B deterministe
                                        b <- rep(0,NPM)

                                        ## calcul des valeurs initiales pour les effets fixes
                                        l <- 0
                                        t <- 0
                                        for (i in 1:nv0)
                                            {
                                                if(idg0[i]==1 & i>1)
                                                    {
                                                        l <- l+1
                                                        t <- t+1
                                                        b[nprob+t] <- B$best[l]
                                                    }
                                                if(idg0[i]==2)
                                                    {
                                                        if (i==1)
                                                            {
                                                                for (g in 2:ng0)
                                                                    {
                                                                        t <- t+1
                                                                        b[nprob+t] <- -0.5*(g-1)
                                                                    }
                                                            }
                                                        if (i>1)
                                                            {
                                                                l <- l+1
                                                                for (g in 1:ng0)
                                                                    {
                                                                        t <- t+1
                                                                        if(B$conv==1) b[nprob+t] <- B$best[l]+(g-(ng0+1)/2)*sqrt(B$V[l*(l+1)/2])
                                                                        else b[nprob+t] <- B$best[l]+(g-(ng0+1)/2)*B$best[l]
                                                                    }
                                                            }
                                                    }
                                            }

                                        ## remplacer varcov par cholesky pour les effets aleatoires
                                        if(nvc>0)
                                            {
                                                if(idiag==TRUE)
                                                    {
                                                        b[nef+1:nvc] <- B$cholesky[(1:nea0)*(2:(nea0+1))/2][-1]
                                                        
                                                    }
                                                else
                                                    {
                                                        b[nef+1:nvc] <- B$cholesky[-1]
                                                    }
                                            }

                                        ## les autres parametres sont inchanges
                                        if (ncor0>0) {b[nef+nvc+nw+1:ncor0] <- B$best[nef2+nvc+1:ncor0]}
                                        b[nef+nvc+nw+ncor0+1:ny0] <- B$best[nef2+nvc+ncor0+1:ny0]
                                        b[nef+nvc+nw+ncor0+ny0+1:nalea0] <- B$best[nef2+nvc+ncor0+ny0+1:nalea0]
                                        b[(nef+nvc+nw+ncor0+ny0+nalea0+1):NPM] <-B$best[(nef2+nvc+ncor0+ny0+nalea0+1):NPM2]

                                    }
                                else
                                    {
                                        ### B random

                                        ## initialiser le vecteur bb contenant les prm (avec repetition) du modele dans B et sa variance vbb
                                        bb <- rep(0,NPM-nprob-nw)
                                        vbb <- matrix(0,NPM-nprob-nw,NPM-nprob-nw)
                                        
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
                                        vbb[(nef-nprob+1):nrow(vbb),(nef-nprob+1):ncol(vbb)] <- VB[(nef2+1):nrow(VB),(nef2+1):ncol(VB)]

                                        

                                        ## remplir les effets fixes (avec repetition si effet specifique a la classe)
                                        l <- 0
                                        t <- 0
                                        for (i in 1:nv0)
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

                                        ## remplacer les varcov par la cholesky
                                        if(nvc>0)
                                            {
                                                if(idiag==TRUE)
                                                    {
                                                        bb[nef-nprob+1:nvc] <- B$cholesky[(1:nea0)*(2:(nea0+1))/2][-1]
                                                    }
                                                else
                                                    {
                                                        bb[nef-nprob+1:nvc] <- B$cholesky[-1]
                                                    }
                                            }
                            
                                        ##les autres parametres sont inchanges
                                        if (ncor0>0)
                                            {
                                                bb[nef-nprob+nvc+1:ncor0] <- B$best[(NPM2-ncor0):(NPM2-1)]
                                            }

                                                                      
                                        bb[nef-nprob+nvc+ncor0+1:ny0] <- B$best[nef2+nvc+ncor0+1:ny0]

                                        if(nalea0>0)
                                            {
                                              bb[nef-nprob+nvc+ncor0+ny0+1:nalea0] <- B$best[nef2+nvc+ncor0+ny0+1:nalea0]  
                                            }

                                        bb[nef-nprob+nvc+ncor0+ny0+nalea0+1:sum(ntrtot0)] <- B$best[nef2+nvc+ncor0+ny0+nalea0+1:sum(ntrtot0)]

                                        ## on enleve les intercepts car ils seront tous initialises a 0
                                        if(idg0[1]>1)
                                            {
                                                bb <- bb[-(1:(ng-1))]
                                                vbb <- vbb[-(1:(ng-1)),-(1:(ng-1))]
                                            }
                                                            
                                        up <- vbb[upper.tri(vbb,diag=TRUE)]
                                        vbb <- t(vbb)
                                        vbb[upper.tri(vbb,diag=TRUE)] <- up
                                        Chol <- chol(vbb)
                                        Chol <- t(Chol)

                                        ## vecteur final b cree a partir de bb et vbb
                                        b <- rep(0,NPM)
                                        
                                        if(idg0[1]>1)
                                            {
                                                b[c((nprob+ng):(nef+nvc),(nef+nvc+nw+1):NPM)] <- bb + Chol %*% rnorm(length(bb))
                                                b[nprob+1:(ng-1)] <- 0
                                            } 
                                        else
                                            {                                        
                                                b[c((nprob+1):(nef+nvc),(nef+nvc+nw+1):NPM)] <- bb + Chol %*% rnorm(NPM-nprob-nw)
                                            }

                                        ## les prm de classmb et nwg sont toujours initalises a 0 et 1
                                        b[1:nprob] <- 0
                                        if(nw>0) b[nef+nvc+1:nw] <- 1

                                        ## if(nvc>0)
                                        ##     {
                                        ##         cholRE <- matrix(0,nea0,nea0)
                                        ##         cholRE[upper.tri(cholRE,diag=TRUE)] <- c(1,b[nef+1:nvc])
                                        ##         varcovRE <- t(cholRE) %*% cholRE
                                        ##         b[nef+1:nvc] <- varcovRE[upper.tri(varcovRE,diag=TRUE)][-1]   
                                        ##     }
                                                                                
                                    }
                            }
                                
                }
        }
    else ## B missing
        {
            b <- rep(0,NPM)
            if (nvc>0)
                {
                    if(idiag==1) b[nef+1:nvc] <- rep(1,nvc)
                    if(idiag==0)
                        {
                            init.nvc <- diag(nea0)
                            init.nvc <- init.nvc[upper.tri(init.nvc, diag=TRUE)]
                            b[nef+1:nvc] <- init.nvc[-1]
                        }
                }
            if(nwg0>0) b[nef+nvc+1:nw] <- 1
            if(ncor0==1) b[nef+nvc+nw+1] <- 1
            if(ncor0==2) b[nef+nvc+nw+1:2] <- c(0,1)

            b[nef+nvc+nw+ncor0+1:ny0] <-  1

            if(nalea0>0) b[nef+nvc+nw+ncor0+ny0+1:nalea0] <- 1

            for(k in 1:ny0)
                {
                    if(idlink0[k]==0)
                        {
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-1] <- mean(get(dataY[k])[,nomsY[k]])
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])] <- 1
                        }
                    if(idlink0[k]==1)
                        {
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-3] <- 0
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-2] <- -log(2)
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-1] <- 0.7
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])] <- 0.1
                        }
                    if(idlink0[k]==2)
                        {
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-ntrtot0[k]+1] <- -2
                            b[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:k])-ntrtot0[k]+2:ntrtot0[k]] <- 0.1
                        }
                }
        }

    ##------------------------------------------
    ##------nom au vecteur best
    ##--------------------------------------------

    nom.X0 <- colnames(X0)
    nom.X0[nom.X0=="(Intercept)"] <- "intercept"
    if(ng0>=2)
        {
            nom <-rep(nom.X0[idprob0==1],each=ng0-1)
            nom1 <- paste(nom," class",c(1:(ng0-1)),sep="")
            names(b)[1:nprob]<-nom1
        }


    if(ng0==1) names(b)[1:(nef-ncontr)] <- nom.X0[-1][idg0[-1]!=0]
    if(ng0>1){
 	nom1<- NULL
 	for (i in 1:nv0) {
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
        names(b)[(nprob+1):(nef-ncontr)]<- nom1
    }

    if(idlink0[1]==0) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+1:ntrtot0[1]]<- c("Linear 1","Linear 2")
    if(idlink0[1]==1) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+1:ntrtot0[1]]<- paste("Beta",c(1:ntrtot0[1]),sep="")
    if(idlink0[1]==2) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+1:ntrtot0[1]]<- paste("I-splines",c(1:ntrtot0[1]),sep="")
    if(ny0>1)
        {
            for (yk in 2:ny0)
                {
                    if(idlink0[yk]==0) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:(yk-1)])+1:ntrtot0[yk]]<- c("Linear 1","Linear 2")
                    if(idlink0[yk]==1) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:(yk-1)])+1:ntrtot0[yk]]<- paste("Beta",c(1:ntrtot0[yk]),sep="")
                    if(idlink0[yk]==2) names(b)[nef+nvc+nw+ncor0+ny0+nalea0+sum(ntrtot0[1:(yk-1)])+1:ntrtot0[yk]]<- paste("I-splines",c(1:ntrtot0[yk]),sep="")
                }
        }
    if(nvc!=0)names(b)[nef+1:nvc] <- paste("varcov",c(1:nvc))
    if(nw!=0)names(b)[nef+nvc+1:nw] <- paste("varprop class",c(1:(ng0-1)))

    names(b)[nef+nvc+nw+ncor0+1:ny0] <- paste("std.err",1:ny0)

    if(ncor0>0) {names(b)[nef+nvc+nw+1:ncor0] <- paste("cor",1:ncor0,sep="")}
    if(nalea0!=0) names(b)[nef+nvc+nw+ncor0+ny0+1:nalea0] <- paste("std.randomY",1:ny0,sep="")
    if(ncontr!=0) names(b)[(nef-ncontr+1):nef] <- paste("contrast",paste(rep(1:sum(idcontr0),each=ny0-1),rep(1:(ny0-1),sum(idcontr0)),sep=""),sep="")


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
    pbH0[grep("I-splines",names(b))] <- 1
    pbH0[posfix] <- 0
    if(sum(pbH0)==0 & Hr0==1) stop("No partial Hessian matrix can be defined")

    
    ##initialisation pour ng>1
    if(missing(B) & ng0>1)
        {
            prior02 <- rep(0,nobs0)
            idprob02 <- rep(0,nv0)
            idg02 <- idg0
            idg02[idg02==2] <- 1
            idcontr02 <- rep(0,nv0)
            ng02 <- 1
            nwg02 <- 0
            nalea02 <- 0

            nef2 <- sum(idg0!=0)-1
            NPM2 <- nef2+ nvc+ncor0+ny0+sum(ntrtot0)
            b2 <- c(rep(0,nef2),b[(nef+1):NPM])
            ind1 <- which(substr(names(b2),1,7)=="varprop")
            if(length(ind1)) b2 <- b2[-ind1]
            ind2 <- which(substr(names(b2),1,11)=="std.randomY")
            if(length(ind2)) b2 <- b2[-ind2]

            V2 <- rep(0,NPM2*(NPM2+1)/2)
            loglik2 <- 0
            ppi02 <- rep(0,ns0)
            pred_m_g2 <- rep(0,nobs0)
            pred_ss_g2 <- rep(0,nobs0)
            predRE_Y2 <- rep(0,ns0)
            maxiter2 <- min(75,maxiter)
            convB2 <- max(0.01,convB)
            convL2 <- max(0.01,convL)
            convG2 <- max(0.01,convG)
            Hr02 <- 0

            init <- .Fortran(C_hetmixcontmult,
                             as.double(Y0),
                             as.double(X0),
                             as.integer(prior02),
                             as.integer(idprob02),
                             as.integer(idea0),
                             as.integer(idg02),
                             as.integer(idcor0),
                             as.integer(idcontr02),
                             as.integer(ny0),
                             as.integer(ns0),
                             as.integer(ng02),
                             as.integer(nv0),
                             as.integer(nobs0),
                             as.integer(nea0),
                             as.integer(nmes0),
                             as.integer(idiag0),
                             as.integer(nwg02),
                             as.integer(ncor0),
                             as.integer(nalea02),
                             as.integer(NPM2),
                             best=as.double(b2),
                             V=as.double(V2),
                             loglik=as.double(loglik2),
                             niter=as.integer(ni),
                             conv=as.integer(istop),
                             gconv=as.double(gconv),
                             ppi2=as.double(ppi02),
                             resid_m=as.double(resid_m),
                             resid_ss=as.double(resid_ss),
                             pred_m_g=as.double(pred_m_g2),
                             pred_ss_g=as.double(pred_ss_g2),
                             predRE=as.double(predRE),
                             predRE_Y=as.double(predRE_Y2),
                             as.double(convB2),
                             as.double(convL2),
                             as.double(convG2),
                             as.integer(maxiter2),
                             as.double(epsY),
                             as.integer(idlink0),
                             as.integer(nbzitr0),
                             as.double(zitr),
                             as.double(uniqueY0),
                             as.integer(indiceY0),
                             as.integer(nvalSPL0),
                             marker=as.double(marker),
                             transfY=as.double(transfY),
                             as.integer(nsim),
                             Yobs=as.double(Yobs),
                             as.integer(Ydiscrete),
                             vraisdiscret=as.double(vraisdiscret),
                             UACV=as.double(UACV),
                             rlindiv=as.double(rlindiv),
                             as.integer(pbH0),
                             as.integer(fix0))


            l <- 0
            t <- 0
            for (i in 1:nv0)
                {
                    if(idg0[i]==1 & i>1)
                        {
                            l <- l+1
                            t <- t+1
                            b[nprob+t] <- init$best[l]
                        }
                    if(idg0[i]==2)
                        {
                            if (i==1)
                                {
                                    for (g in 2:ng0)
                                        {
                                            t <- t+1
                                            b[nprob+t] <- -0.5*(g-1)
                                        }
                                }
                            if (i>1)
                                {
                                    l <- l+1
                                    for (g in 1:ng0)
                                        {
                                            t <- t+1
                                            if(init$conv==1) b[nprob+t] <- init$best[l]+(g-(ng0+1)/2)*sqrt(init$V[l*(l+1)/2])
                                            else b[nprob+t] <- init$best[l]+(g-(ng0+1)/2)*init$best[l]
                                        }
                                }
                        }
                }

            if(nvc>0) b[nef+1:nvc] <-init$best[nef2+1:nvc]
            if (ncor0>0) {b[nef+nvc+nw+1:ncor0] <- init$best[nef2+nvc+1:ncor0]}
            b[nef+nvc+nw+ncor0+1:ny0] <- init$best[nef2+nvc+ncor0+1:ny0]
            b[(nef+nvc+nw+ncor0+ny0+nalea0+1):NPM] <-init$best[(nef2+nvc+ncor0+ny0+1):NPM2]

        }

    
###estimation
    out <- .Fortran(C_hetmixcontmult,
                    as.double(Y0),
                    as.double(X0),
                    as.integer(prior0),
                    as.integer(idprob0),
                    as.integer(idea0),
                    as.integer(idg0),
                    as.integer(idcor0),
                    as.integer(idcontr0),
                    as.integer(ny0),
                    as.integer(ns0),
                    as.integer(ng0),
                    as.integer(nv0),
                    as.integer(nobs0),
                    as.integer(nea0),
                    as.integer(nmes0),
                    as.integer(idiag0),
                    as.integer(nwg0),
                    as.integer(ncor0),
                    as.integer(nalea0),
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
                    predRE_Y=as.double(predRE_Y),
                    as.double(convB),
                    as.double(convL),
                    as.double(convG),
                    as.integer(maxiter),
                    as.double(epsY),
                    as.integer(idlink0),
                    as.integer(nbzitr0),
                    as.double(zitr),
                    as.double(uniqueY0),
                    as.integer(indiceY0),
                    as.integer(nvalSPL0),
                    marker=as.double(marker),
                    transfY=as.double(transfY),
                    as.integer(nsim),
                    Yobs=as.double(Yobs),
                    as.integer(Ydiscrete),
                    vraisdiscret=as.double(vraisdiscret),
                    UACV=as.double(UACV),
                    rlindiv=as.double(rlindiv),
                    as.integer(pbH0),
                    as.integer(fix0))
    

### mettre NA pour les variances et covariances non calculees et  0 pr les prm fixes
    if(length(posfix))
        {
            if(out$conv==3)
                {
                    mr <- NPM-sum(pbH0)-length(posfix)
                    Vr <- matrix(0,mr,mr)
                    Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
                    Vr <- t(Vr)
                    Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
                    V <- matrix(NA,NPM,NPM)
                    V[setdiff(1:NPM,c(which(pbH0==1),posfix)),setdiff(1:NPM,c(which(pbH0==1),posfix))] <- Vr
                    V[,posfix] <- 0
                    V[posfix,] <- 0
                    V <- V[upper.tri(V,diag=TRUE)]
                }
            else
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
        }
    else
        {
            if(out$conv==3)
                {
                    mr <- NPM-sum(pbH0)
                    Vr <- matrix(0,mr,mr)
                    Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
                    Vr <- t(Vr)
                    Vr[upper.tri(Vr,diag=TRUE)] <- out$V[1:(mr*(mr+1)/2)]
                    V <- matrix(NA,NPM,NPM)
                    V[setdiff(1:NPM,which(pbH0==1)),setdiff(1:NPM,which(pbH0==1))] <- Vr
                    V <- V[upper.tri(V,diag=TRUE)]
                }
            else
                {
                    V <- out$V
                }
        }
    

### Creation du vecteur cholesky
    Cholesky <- rep(0,(nea0*(nea0+1)/2))
    if(idiag0==0 & nvc>0)
        {
            Cholesky[1:(nvc+1)] <- c(1,out$best[nef+1:nvc])
            ## Construction de la matrice U
            U <- matrix(0,nrow=nea0,ncol=nea0)
            U[upper.tri(U,diag=TRUE)] <- Cholesky[1:(nvc+1)]
            z <- t(U) %*% U
            out$best[nef+1:nvc] <- z[upper.tri(z,diag=TRUE)][-1]
        }
    if(idiag0==1 & nvc>0)
        {
            id <- 1:nea0
            indice <- rep(id+id*(id-1)/2)
            Cholesky[indice] <- c(1,out$best[nef+1:nvc])
            out$best[nef+1:nvc] <- out$best[nef+1:nvc]**2
        }

###predictions
    predRE <- matrix(out$predRE,ncol=nea0,byrow=T)
    predRE <- data.frame(unique(IND),predRE)
    colnames(predRE) <- c(nom.subject,nom.X0[idea0!=0])

    if (nalea0!=0)
        {
            predRE_Y <- matrix(out$predRE_Y,ncol=ny0,byrow=TRUE)
            predRE_Y <- data.frame(unique(IND),predRE_Y)
            colnames(predRE_Y)  <- c(nom.subject,nomsY)
        }
    else
        {
            predRE_Y <- rep(NA,nalea0*ns0)
        }

###ppi
    if(ng0>1) {
        ppi<- matrix(out$ppi2,ncol=ng0,byrow=TRUE)
    }
    else {
        ppi <- matrix(rep(1,ns0),ncol=ng0)
    }

    classif<-apply(ppi,1,which.max)
    ppi<-data.frame(unique(IND),classif,ppi)
    temp<-paste("prob",1:ng0,sep="")
    colnames(ppi) <- c(nom.subject,"class",temp)
    rownames(ppi) <- 1:ns0


###pred
    pred_m_g <- matrix(out$pred_m_g,nrow=nobs0)
    pred_ss_g <- matrix(out$pred_ss_g,nrow=nobs0)
    pred_m <- out$Yobs-out$resid_m
    pred_ss <- out$Yobs - out$resid_ss
    pred <- data.frame(IND,outcome,pred_m,out$resid_m,pred_ss,out$resid_ss,out$Yobs,pred_m_g,pred_ss_g)

    temp<-paste("pred_m",1:ng0,sep="")
    temp1<-paste("pred_ss",1:ng0,sep="")
    colnames(pred)<-c(nom.subject,"Yname","pred_m","resid_m","pred_ss","resid_ss","obs",temp,temp1)
    rownames(pred) <- NULL

###estimlink
    ysim <- matrix(out$marker,nsim,ny0)
    transfo <- matrix(out$transfY,nsim,ny0)
    estimlink <- as.vector(rbind(ysim,transfo))
    estimlink <- matrix(estimlink,nsim,2*ny0)
    colnames(estimlink) <- paste(c("","transf"),rep(nomsY, each=2),sep="")


    N <- NULL
    N[1] <- (ng0-1)*sum(idprob0)
    N[2] <- (ny0-1)*sum(idcontr0)
    N[3] <- (ng0-1)*sum(idprob0) + sum(idg0==1)-1 + ng0*sum(idg0==2) + (ny0-1)*sum(idcontr0)  #nef
    N[4] <- ifelse(idiag0==1,nea0,nea0*(nea0+1)/2)-1  #nvc
    N[5] <- (ng0-1)*nwg0
    N[6] <- nalea0
    N[7] <- ncor0
    N[8] <- ny0
    N[9] <- nobs0

    nom.X0[nom.X0=="(Intercept)"] <- "Intercept"

    res <-list(ns=ns0,ng=ng0,idea0=idea0,idprob0=idprob0,idg0=idg0,idcontr0=idcontr0,
               idcor0=idcor0,loglik=out$loglik,best=out$best,V=V,gconv=out$gconv,conv=out$conv,
               call=cl,niter=out$niter,N=N,idiag=idiag0,pred=pred,pprob=ppi,predRE=predRE,
               predRE_Y=predRE_Y,Ynames=nomsY,Xnames=nom.X0,Xnames2=ttesLesVar,cholesky=Cholesky,
               estimlink=estimlink,epsY=epsY,linktype=idlink0,linknodes=zitr,nbnodes=nbnodes,
               na.action=nayk,AIC=2*(length(out$best)-length(posfix)-out$loglik),BIC=(length(out$best)-length(posfix))*log(ns0)-2*out$loglik)
    
    names(res$best) <- names(b)
    class(res) <-c("multlcmm")

    cost<-proc.time()-ptm
    if(verbose==TRUE) cat("The program took", round(cost[3],2), "seconds \n")

    res
}



mlcmm <- multlcmm
