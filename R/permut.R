#' Permutation of the latent classes
#'
#' This function allows a re-ordering of the latent classes of an estimated model.
#'
#' @param m an object inheriting from classes \code{hlme}, \code{lcmm}, \code{multlcmm} or \code{Jointlcmm}
#' @param order a vector (integer between 1 and ng) containing the new order of the latent classes
#' @param estim optional boolean specifying if the model should estimated with the reordered parameters as initial values. By default, the model is estimated. If FALSE, only the coefficients in \code{$best} are modified. All other outputs are not changed.
#' @return An object of the same class as m, with reordered classes, or the initial object with new coefficients if estim is FALSE.
#' @author Viviane Philipps and Cecile Proust-Lima
#' @examples
#'
#' ## Estimation of a hlme model with 2 classes
#' m2 <- hlme(Y~Time*X1,mixture=~Time,random=~Time,classmb=~X2+X3,subject='ID',
#'          ng=2,data=data_hlme,B=c(0.11,-0.74,-0.07,20.71,
#'                                  29.39,-1,0.13,2.45,-0.29,4.5,0.36,0.79,0.97))
#'
#' ## Exchange class 2 and class 1
#' m2b <- permut(m2,order=c(2,1))
#' 
#' @export
permut <- function(m,order,estim=TRUE)
    {
        ng <- m$ng

        if(!(class(m) %in% c("hlme","lcmm","multlcmm","Jointlcmm"))) stop("Please use this function only with hlme, lcmm, multlcmm or Jointlcmm models")
        if(ng==1) stop("Please use this function only with latent classes")
        if(length(order)!=ng) stop(paste("Argument 'order' should be of length",ng))
        if(!all(order %in% 1:ng)) stop(paste("Please specify classes between 1 and",ng,"in argument 'order'",collapse=" "))
        if(!all(c(1:ng) %in% order)) stop("Please specify all the classes in argument 'order'")

        bnew <- m$best

        ## coef de classmb reordonnes:
        nv <- m$N[1]/(ng-1) # nb de variables dans classmb (avec l'intercept)
        
        coefold <- rep(0,ng*nv) # coef 0 pour la classe G
        coefold[sapply(c(0:(nv-1)*ng),"+",1:(ng-1))] <- m$best[1:m$N[1]] # coef de best pour les autres classes

        coeford <- coefold[as.vector(sapply(c((0:(nv-1))*ng),"+",order))] #ordonner selon l'ordre donne dans l'argument order
        coefref <- coeford-rep(coeford[sapply(c(0:(nv-1)*ng),"+",ng)],each=ng) # soustraire le coef de la nouvelle classe de reference
        coefnew <- coefref[sapply(c(0:(nv-1)*ng),"+",1:(ng-1))] # enlever les coef de la nouvelle classe de reference
                
        bnew[1:m$N[1]] <- coefnew

        ## coef survie pour Jointlcmm
        if(class(m)=="Jointlcmm")
            {
                ## risque de base
                nbevt <- length(m$hazard[[1]])
                nprisq <- rep(NA,nbevt)
                nrisq <- rep(NA,nbevt)
                
                typrisq <- m$hazard[[1]]
                hazardtype <- m$hazard[[2]]
                nz <- m$hazard[[4]] 
                
                for(ke in 1:nbevt)
                    {
                        if(typrisq[ke]==1) nprisq[ke] <- nz[ke]-1
                        if(typrisq[ke]==2) nprisq[ke] <- 2
                        if(typrisq[ke]==3) nprisq[ke] <- nz[ke]+2
                        
                        if(hazardtype[ke]=="Common") nrisq[ke] <- nprisq[ke]
                        if(hazardtype[ke]=="PH") nrisq[ke] <- nprisq[ke]+m$ng-1
                        if(hazardtype[ke]=="Specific") nrisq[ke] <- nprisq[ke]*m$ng

                        if(hazardtype[ke]=="Specific")
                            {
                                ## echanger les coef
                                bnew[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+1:nrisq[ke]] <- m$best[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+sapply(order,function(x,np) {(x-1)*np+1:np},np=nprisq[ke])]
                            }
                        if(hazardtype=="PH")
                        {
                            if(order[ng]==ng)
                            {
                                wref <- 0
                            }
                            else
                            {
                                wref <- m$best[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+nprisq[ke]+order[ng]] 
                            }
                            
                            if(typrisq[ke]==2) #Weibull
                            {
                                if(m$logspecif==1)
                                {
                                    bnew[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+1] <- m$best[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+1] + wref # a devient a+wref 
                                }
                                else
                                {
                                    bnew[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+1] <- m$best[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+1] * exp(wref/(2*m$best[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+2]^2)) # a devient a*exp(wref/(2*b^2))
                                }
                                ## b ne change pas
                            }
                            else #piecewise ou splines
                            {
                                if(m$logspecif==1)
                                {
                                    bnew[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+1:nprisq[ke]] <- m$best[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+1:nprisq[ke]] + wref
                                }
                                else
                                {
                                    bnew[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+1:nprisq[ke]] <- m$best[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+1:nprisq[ke]] * exp(wref/2)
                                }
                            }
                            
                                ## PH : wg devient wg-wref
                                coefold <- c(m$best[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+nprisq[ke]+1:(ng-1)],0) # car wG=0
                                coeford <- coefold[order] # ordonner selon order
                                coefref <- coeford-coeford[ng] # soustraire le coef de la nouvelle classe de reference
                                coefnew <- coefref[1:(ng-1)] # enlever le dernier
                                bnew[m$N[1]+sum(nrisq[1:ke])-nrisq[ke]+nprisq[ke]+1:(ng-1)] <- coefnew
                            }
                    }


                ## variables explicatives
                if(m$N[3]>0)
                    {
                        idspecif <- matrix(m$idspecif,nrow=nbevt,byrow=TRUE)
                        avt <- 0
                        for(j in 1:ncol(idspecif))
                            {
                                if((m$idcom[j]==1) & (idspecif[1,j]==1)) # coef commun a tous les evts et pas en mixture
                                    {
                                        avt <- avt+1
                                    }
                                if((m$idcom[j]==1) & (idspecif[1,j]==2)) # coef commum et en mixture
                                    {
                                        bnew[m$N[1]+m$N[2]+avt+1:ng] <- m$best[m$N[1]+m$N[2]+avt+order] # echanger les coef
                                        avt <- avt+ng
                                    }
                                if(m$idcom[j]==0) # coef specifique a chaque evt 
                                    {
                                        for(ke in 1:nbevt)
                                            {
                                                if(idspecif[ke,j]==1) # pas en mixture
                                                    {
                                                        avt <- avt+1
                                                    }
                                                if(idspecif[ke,j]==2) # coef en mixture
                                                    {
                                                        bnew[m$N[1]+m$N[2]+avt+1:ng] <- m$best[m$N[1]+m$N[2]+avt+order] # echanger les coef
                                                        avt <- avt+ng
                                                    }
                                            }
                                    }
                            }
                       
                    }
                #if(avt!=m$N[3]) stop("pb varxevt")
            } # fin partie survie

        ## effets fixes modele mixte:
        if(class(m) %in% c("hlme","lcmm"))
            {
                avtnef <- m$N[1]
                nef <- m$N[2]
            }
        if(class(m)=="multlcmm")
            {
                avtnef <- m$N[1]
                nef <- m$N[3]-m$N[1]-m$N[2]
            }
        if(class(m)=="Jointlcmm")
            {
                avtnef <- m$N[1]+m$N[2]+m$N[3]
                nef <- m$N[4]
            }
        tmpnef <- 0
        if(m$idg[1]==2) # intercept en mixture
            {
                if((class(m)=="hlme") | ((class(m)=="Jointlcmm") & isTRUE(m$linktype==-1)))
                    {
                        bnew[avtnef+1:ng] <- m$best[avtnef+order] # echanger
                        tmpnef <- tmpnef+ng
                    }
                else
                    {
                        coefold <- c(0,m$best[avtnef+1:(ng-1)]) # car intercept fixe a 0 dans la premiere classe
                        coeford <- coefold[order] # ordonner
                        coefref <- coeford-coeford[1] # soustraire le nouveau coef de ref
                        coefnew <- coefref[2:ng] # enlever le premier coef
                        bnew[avtnef+1:(ng-1)] <- coefnew
                        tmpnef <- tmpnef+ng-1
                    }
            }
        if(nef>tmpnef) # autres coef
            {
                for(k in 2:length(m$idg))
                    {
                        if(m$idg[k]==1)
                            {
                                tmpnef <- tmpnef+1
                            }
                        if(m$idg[k]==2)
                            {
                                bnew[avtnef+tmpnef+1:ng] <- m$best[avtnef+tmpnef+order] #echanger
                                tmpnef <- tmpnef+ng
                            }
                    }
            }

        ## coef nw:
        if(class(m) %in% c("hlme","lcmm"))
        {
            avtnw <- m$N[1]+m$N[2]+m$N[3]
            nw <- m$N[4]
            nvc <- m$N[3]
        }
        if(class(m)=="multlcmm")
        {
            avtnw <- m$N[3]+m$N[4]
            nw <- m$N[5]
            nvc <- m$N[4]
        }
        if(class(m)=="Jointlcmm")
        {
            avtnw <- m$N[1]+m$N[2]+m$N[3]+m$N[4]+m$N[5]
            nw <- m$N[6]
            nvc <- m$N[5]
        }
        if(nw>0)
        {
            if(nvc>0)
            {
                reold <- m$best[(avtnw-nvc+1):avtnw]
                renew <- reold*m$best[avtnw+order[ng]]^2
                bnew[(avtnw-nvc+1):avtnw] <- renew
            }
            
            coefold <- c(m$best[avtnw+1:(ng-1)],1)
            coeford <- coefold[order]
            coefref <- coeford/coeford[ng]
            coefnew <- coefref[1:(ng-1)]
            bnew[avtnw+1:(ng-1)] <- coefnew
        }
        
        ## modele avec les nouveaux coef
        if(estim==TRUE)
            {
                z <- m$call
                z$B <- bnew
                mnew <- eval(z)
            }
        else
            {
                mnew <- m
                mnew$best <- bnew
            }
        
        return(mnew)
    }

