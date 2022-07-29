#' Multivariate Wald Test
#' 
#' This function provides multivariate and univariate Wald tests for
#' combinations of parameters from \code{hlme}, \code{lcmm}, \code{multlcmm},
#' \code{Jointlcmm} or \code{mpjlcmm} models.
#' 
#' 
#' @param Mod an object of class \code{hlme}, \code{lcmm}, \code{multlcmm},
#' \code{Jointlcmm} or \code{mpjlcmm}
#' @param pos a vector containing the indices in \code{Mod$best} of the parameters to
#' test
#' @param contrasts a numeric vector of same length as pos.  If NULL (the
#' default), a simultaneous test of the appropriate parameters is realised.  If
#' contrasts is specified, the quantity to test is the dot product of pos and
#' contrasts.
#' @param name a character containing the name the user wants to give to the
#' test. By default, the name's test is the null hypothesis.
#' @param value the value(s) to test against. By default, test against 0.
#' @return If contrasts is NULL, the function returns a matrix with 1 row and 2
#' columns containing the value of the Wald test's statistic and the associated
#' p-value.
#' 
#' If contrasts is not NULL, the function returns a matrix with 1 row and 4
#' columns containing the value of the coefficient (dot product of pos and
#' contrasts), his standard deviation, the value of the Wald test's statistic
#' and the associated p-value.
#' @author Cecile Proust-Lima, Lionelle Nkam and Viviane Philipps
#' 
#' @export
#' 
WaldMult <- function(Mod,pos=NULL,contrasts=NULL,name=NULL,value=NULL)
{ 
    if (!(class(Mod) %in% c("hlme","lcmm","multlcmm","Jointlcmm","mpjlcmm"))) stop("applies to \"hlme\" or \"lcmm\" or \"multlcmm\" or \"Jointlcmm\" or \"mpjlcmm\" objects only")


    ## prm estimes et leur variance
    Mod$best <- estimates(Mod)
    V <- VarCov(Mod)
    
    
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


