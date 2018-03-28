
WaldMult <- function(Mod,pos=NULL,contrasts=NULL,name=NULL,value=NULL)
{ 
    if (!(class(Mod) %in% c("hlme","lcmm","multlcmm","Jointlcmm"))) stop("applies to \"hlme\" or \"lcmm\" or \"multlcmm\" or \"Jointlcmm\" objects only")
    
    if(inherits(Mod,"hlme") | inherits(Mod,"lcmm"))
        {
            nea<-sum(Mod$idea) # Nombre d'effets aleatoires
            nef<-Mod$N[2]      # Nombre d'effets fixes
            nvc<-Mod$N[3]      # Nombre de parametres effets aleatoires
            nprob<-Mod$N[1]
            idiag <- ifelse(Mod$idiag==1,TRUE,FALSE)
        }
    if(inherits(Mod,"multlcmm"))
        {
            nea <- sum(Mod$idea0)
            nef <- Mod$N[3]
            nvc <- Mod$N[4]
            nprob <- 0 #nef contient deja nprob
            idiag <- ifelse(Mod$idiag==1,TRUE,FALSE)
        }
    
    if(inherits(Mod,"Jointlcmm"))
        {
            nea <- sum(Mod$idea)
            nef <- Mod$N[4]
            nvc <- Mod$N[5]
            nprob <- sum(Mod$N[1:3]) #nprob+nvarxevt+nristot
            idiag <- ifelse(Mod$idiag==1,TRUE,FALSE)
        }

    
    ## On remplace les varcov de $best par les parametres de cholesky
    if(nvc>0)
        {   
            debut <- nprob+nef+1
            fin <- nprob+nef+nvc
            cholesky <- Mod$cholesky
            if(inherits(Mod,"multlcmm")) cholesky[1] <- NA

            if(isTRUE(idiag)) cholesky[setdiff(1:(nea*(nea+1)/2),1:nea*(1:nea+1)/2)] <- NA
            
            Mod$best[debut:fin] <- na.omit(cholesky)
        }
    
    ## On rend symetrique la matrice des varcov des parametres estimes par le modele
    l <- length(Mod$best)
    V <- matrix(0,nrow=l,ncol=l)
    V[upper.tri(V,diag=TRUE)] <- Mod$V  
    V[lower.tri(V,diag=FALSE)] <- t(V)[lower.tri(V,diag=FALSE)]
    
    ## On teste si pos est un vecteur
    if (is.null(pos))
        {
            stop("pos must be specified")
        }
    else
        {
            if (!is.vector(pos)) stop("Error : pos must be a numeric vector")
            
            ## On cree la matrice qui recevra les varcov des parametres de pos
            Mat <- matrix(0,nrow=length(pos),ncol=length(pos))
            
            ## Remplissage de Mat sans boucles
            Mat <- V[pos,pos]
            
            ## Wald Multivarie, sans le vecteur contrasts
            Vect <- Mod$best[pos]
            
            
            if (is.null(contrasts))
                { 
                    
                    if (!is.null(value))
                        {
                            if (!is.vector(value)) stop("Error : value must be a numeric vector")
                            if (length(value)!=length(pos)) stop("value must have the same length as the vector pos")
                          
                            Vect <- Mod$best[pos]-value
                        }
                    
                    
                    Wald <- t(Vect)%*%solve(Mat)%*%Vect
                    
                    ## Nombre de degre de liberte
                    ddl <- length(pos)
                    ## Pvalue
                    p_value <- 1-pchisq(Wald,df=ddl)
                    ##cat("pvalue",p_value1,"\n")
                    
                    Results <- matrix(NA,nrow=1,ncol=2)
                    colnames(Results)<-c("Wald Test","p_value")
                    if (is.null(name)) 
                        {
                            if (!is.null(value))
                                {
                                    rownames(Results)<-paste(names(Mod$best[pos])," = ",value,collapse=" and ",sep="")
                                }
                            else
                                {
                                    rownames(Results)<-paste(paste(names(Mod$best[pos]),collapse=" = "),"= 0")
                                }
                        }
                    ## paste(names(Mod$best[pos])," = 0 ",collapse=" and ",sep="")} 
                    else
                        {
                            rownames(Results)<-name
                        }
                    
                    Results[,1] <- round(Wald,5)
                    Results[,2] <- round(p_value,5)
                }
            
            ## Wald Univarie avec le vecteur contrasts
            else
                {
                    ## Conditions d'application
                    if (length(contrasts)!=length(pos))
                        {
                            stop("contrasts must have the same length as the vector pos")
                        } 
                    if (sum(abs(contrasts))==0) 
                        {
                            stop("The absolute value of the sum of contratsts components must be different from 0")
                        }
                
                

                    Scalaire <- sum(Vect*contrasts)
                    
                    ## Utilisation de value
                    if (!is.null(value))
                        {
                            if (!is.vector(value)) stop("value must be a numeric vector")
                            if (length(value)!=1) stop("value must be a vector with a unique argument")
                        
                            Scalaire <- sum(Vect*contrasts)-value
                        }
                
                    ## Calcul de la variance de scalaire sans boucle
                    Var <- t(contrasts)%*%Mat%*%contrasts
                

                    Wald <- Scalaire/sqrt(Var)
                    p_value <- 2*(1-pnorm(abs(Wald)))
                    
                    Results <- matrix(NA,nrow=1,ncol=4)
                    colnames(Results)<-c("coef","Se","Wald Test","p_value")
                    if (is.null(name)) 
                    {
                        if(is.null(value)) value <- 0
                        rownames(Results)<-paste(paste(names(Mod$best[pos]),"*",contrasts,collapse=" + "),"= ",value)
                    } 
                    else
                        {
                            rownames(Results)<-name
                        }
                
                
                    Results[,1] <- round(sum(Vect*contrasts),5)
                    Results[,2] <- round(sqrt(Var),5)
                    Results[,3] <- round(Wald,5)
                    Results[,4] <- round(p_value,5)
                }
            
            return(Results)
            
        }  
}


