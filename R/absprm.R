absprm <- function(model, cholesky=TRUE)
{
    ## replace varcovRE with cholesky
    if(cholesky) prm <- estimates(model) else prm <- model$best
    ## attention : en idiag, la chol peut avoir des coef <0

    ## take the abs value of parameters :
    ##  - risq (if logscale=FALSE)
    ##  - wg
    ##  - sigma measurement error
    ##  - sigma cor AR/BM
    ##  - link (all prm but the first for splines and thresholds, only the last for other links)
    ##  - sigma randomY for multlcmm

    ################### hlme ###################
    if(inherits(model, "hlme"))
    {
        nprob <- model$N[1]
        nef <- model$N[2]
        nvc <- model$N[3]
        nw <- model$N[4]
        ncor <- model$N[5]

        if(nw>0) prm[nprob+nef+nvc+1:nw] <- abs(prm[nprob+nef+nvc+1:nw])
        
        if(ncor>0) prm[nprob+nef+nvc+nw+ncor] <- abs(prm[nprob+nef+nvc+nw+ncor])

        prm[nprob+nef+nvc+nw+ncor+1] <- abs(prm[nprob+nef+nvc+nw+ncor+1])
    }

    
    ################### lcmm ###################
    if(inherits(model, "lcmm"))
    {
        nprob <- model$N[1]
        nef <- model$N[2]
        nvc <- model$N[3]
        nw <- model$N[4]
        ncor <- model$N[6]
        ntrtot <- length(model$best) - (nprob+nef+nvc+nw+ncor)

        if(nw>0) prm[nprob+nef+nvc+1:nw] <- abs(prm[nprob+nef+nvc+1:nw])
        
        if(ncor>0) prm[nprob+nef+nvc+nw+ntrtot+ncor] <- abs(prm[nprob+nef+nvc+nw+ntrtot+ncor])

        if(model$linktype %in% c(2,3))
        {
            prm[nprob+nef+nvc+nw+2:ntrtot] <- abs(prm[nprob+nef+nvc+nw+2:ntrtot])
        }
        else
        {
            prm[nprob+nef+nvc+nw+ntrtot] <- abs(prm[nprob+nef+nvc+nw+ntrtot])
        }
    }

    
    ################### multlcmm ###################
    if(inherits(model, "multlcmm"))
    {
        nef <- model$N[3] #nprob and ncontr included in nef
        nvc <- model$N[4]
        nw <- model$N[5]
        nalea <- model$N[6]
        ncor <- model$N[7]
        ny <- model$N[8]
        ntrtot <- length(model$best) - (nef+nvc+nw+ncor+ny+nalea)

        if(nw>0) prm[nef+nvc+1:nw] <- abs(prm[nef+nvc+1:nw])
        
        if(ncor>0) prm[nef+nvc+nw+ncor] <- abs(prm[nef+nvc+nw+ncor])

        prm[nef+nvc+nw+ncor+1:ny] <- abs(prm[nef+nvc+nw+ncor+1:ny])
        
        if(nalea>0) prm[nef+nvc+nw+ncor+ny+1:nalea] <- abs(prm[nef+nvc+nw+ncor+ny+1:nalea])

        sumntr <- 0
        numspl <- 0
        for(k in 1:ny)
        {
            if(model$linktype[k] %in% c(2,3))
            {               
                if(model$linktype[k] == 2)
                {
                    numspl <- numspl + 1
                    ntr <- model$nbnodes[numspl]+2
                }            
                if(model$linktype[k] == 3)
                {
                    ntr <- model$nbmod[k]-1
                }
                                
                prm[nef+nvc+nw+ncor+ny+nalea+sumntr+2:ntr] <- abs(prm[nef+nvc+nw+ncor+ny+nalea+sumntr+2:ntr])
            }
            else
            {            
                if(model$linktype[k] == 0)
                {
                    ntr <- 2
                }  
                if(model$linktype[k] == 1)
                {
                    ntr <- 4
                }
                
                prm[nef+nvc+nw+ncor+ny+nalea+sumntr+ntr] <- abs(prm[nef+nvc+nw+ncor+ny+nalea+sumntr+ntr])
            }

            sumntr <- sumntr + ntr
        }
    }

    ################### Jointlcmm ###################
    if(inherits(model, "Jointlcmm"))
    {
        nprob <- model$N[1]
        nrisqtot <- model$N[2]

        ## survival part
        if(!(model$logspecif))
        {
            nbevt <- length(model$N) - 9
            typrisq <- model$hazard[[1]]
            hazardtype <- model$hazard[[2]]
            nz <- model$hazard[[4]]

            sumnrisq <- 0
            for(ke in 1:nbevt)
            {
                nprisq <- 2
                if(typrisq[ke] == 1) nprisq <- nz[ke] - 1
                if(typrisq[ke] == 3) nprisq <- nz[ke] + 2
                
                if(hazardtype[ke] == "Specific")
                {
                    nrisq <- model$ng*nprisq
                    prm[nprob+sumnrisq+1:nrisq] <- abs(prm[nprob+sumnrisq+1:nrisq])
                }
                else
                {
                    nrisq <- nprisq # "Common" risq
                    if(hazardtype[ke] == "PH") nrisq <- nprisq + model$ng - 1
                    prm[nprob+sumnrisq+1:nprisq] <- abs(prm[nprob+sumnrisq+1:nprisq])
                }

                sumnrisq <- sumnrisq + nrisq
            }
        }

        ## longitudinal part
        nvarxevt <- model$N[3]
        nef <- model$N[4]
        nvc <- model$N[5]
        nw <- model$N[6]
        ncor <- model$N[7]
        ntrtot <- model$N[8]
        
        if(nw > 0) prm[nprob+nrisqtot+nvarxevt+nef+nvc+1:nw] <- abs(prm[nprob+nrisqtot+nvarxevt+nef+nvc+1:nw])

        if(ncor > 0) prm[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor] <- abs(prm[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor])

        if(model$linktype < 2)
        {
            prm[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+ntrtot] <- abs(prm[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+ntrtot])
        }
        else
        {
            prm[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+2:ntrtot] <- abs(prm[nprob+nrisqtot+nvarxevt+nef+nvc+nw+ncor+2:ntrtot])
        }
    }

    ################### mpjlcmm ###################
    if(inherits(model, "mpjlcmm"))
    {
        nprob <- model$N[1]
        nrisqtot <- model$N[2]
        nbevt <- model$nbevt

        ## survival part
        if((nbevt>0) & !(model$logspecif))
        {
            typrisq <- model$typrisq
            hazardtype <- model$hazardtype
            nz <- model$nz

            sumnrisq <- 0
            for(ke in 1:nbevt)
            {
                nprisq <- 2
                if(typrisq[ke] == 1) nprisq <- nz[ke] - 1
                if(typrisq[ke] == 3) nprisq <- nz[ke] + 2

                if(hazardtype[ke] == "Specific")
                {
                    nrisq <- model$ng*nprisq
                    prm[nprob+sumnrisq+1:nrisq] <- abs(prm[nprob+sumnrisq+1:nrisq])
                }
                else
                {
                    nrisq <- nprisq # "Common" risq
                    if(hazardtype[ke] == "PH") nrisq <- nprisq + model$ng - 1
                    prm[nprob+sumnrisq+1:nprisq] <- abs(prm[nprob+sumnrisq+1:nprisq])
                }

                sumnrisq <- sumnrisq + nrisq
            }
        }

        ## longitudinal part
        K <- model$K
        ny <- model$ny
        nvarxevt <- model$N[3]
        if(nbevt == 0) nbevt <- 1 # for index in Nprm
        nef <- model$Nprm[2+nbevt+1:K]
        ncontr <- model$Nprm[2+nbevt+K+1:K]
        nvc <- model$Nprm[2+nbevt+2*K+1:K]
        nw <- model$Nprm[2+nbevt+3*K+1:K]
        ncor <- model$Nprm[2+nbevt+4*K+1:K]
        nerr <- model$Nprm[2+nbevt+5*K+1:K]
        nalea <- model$Nprm[2+nbevt+6*K+1:K]
        ntr <- model$Nprm[2+nbevt+7*K+1:sum(ny)]

        sumnpm <- nprob+nrisqtot+nvarxevt
        for(k in 1:K)
        {
            if(nw[k] > 0) prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+1:nw[k]] <- abs(prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+1:nw[k]])

            if(ncor[k] > 0) prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]] <- abs(prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]])

            if(nerr[k] > 0) prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+1:nerr[k]] <- abs(prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+1:nerr[k]])
            
            if(nalea[k] > 0) prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+1:nalea[k]] <- abs(prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+1:nalea[k]])

            sumntr <- 0
            for(m in 1:ny[k])
            {
                ym <- sum(ny[1:k])-ny[k]+m
                
                if(model$linktype[ym] %in% c(0,1)) prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+sumntr+ntr[ym]] <- abs(prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+sumntr+ntr[ym]])
                
                if(model$linktype[ym] == 2) prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+sumntr+2:ntr[ym]] <- abs(prm[sumnpm+nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+sumntr+2:ntr[ym]])
                
                sumntr <- sumntr + ntr[ym]
            }

            sumnpm <- sumnpm + nef[k]+ncontr[k]+nvc[k]+nw[k]+ncor[k]+nerr[k]+nalea[k]+sumntr
        }   
    }
    
    ## put in model object
    model$best <- prm

    return(model)
}
