VarCovRE.multlcmm <- function(Mod)
{
  if(missing(Mod)) stop("The model should be specified")
  if(!inherits(Mod,"multlcmm")) stop("applies to \"multlcmm\" objects only")

   nea <- sum(Mod$idea0)
   nef <- Mod$N[3]
   nvc <- Mod$N[4]
   nprob <- 0 #nef contient deja nprob
   idiag <- ifelse(Mod$idiag==1,TRUE,FALSE)

  if(nvc==0) return(NA)

  # Recuperation des parametres de Cholesky dans un vecteur
  cholesky<-Mod$cholesky

  # Transformation de Cholesky telle que UU'=B
  tU <- matrix(0,nrow=nea,ncol=nea)
  tU[upper.tri(tU,diag=TRUE)] <- cholesky  # La partie superieure de la matrice U y compris la diagonale recoit le vecteur cholesky
  U<-t(tU)

  B<-U%*%tU   # Matrice de variance-covariance des effets aleatoires

  # On end symetrique la matrice de variance-covariance V des parametres estimes par le modele
  l<-length(Mod$best)
  V<- matrix(0,nrow=l,ncol=l)
  V[upper.tri(V,diag=TRUE)] <- Mod$V
  V[lower.tri(V,diag=FALSE)] <- t(V)[lower.tri(V,diag=FALSE)]


  # On recupere la sous matrice de var-cov des parametres de cholesky dans V
  debut<-nprob+nef+1
  fin<-nprob+nef+nvc
  Vc<-V[debut:fin,debut:fin]

  # Calcul de la derivee de B en fonction des parametres de cholesky : D=f'(c)/c
  # D est different selon la valeur de idiag
  if (idiag == FALSE)
  {
    # calcul de derivee
    B[lower.tri(B,diag=FALSE)]<-0
    tB<-B
    Sup<-B[upper.tri(B,diag=TRUE)]
    D<-matrix(0,nrow=length(Sup),ncol=length(Sup))
    # Vaut mieux eviter les boucles et faire le moins de if possible

    for (j in 1:ncol(tU))
    {
     for (i in 1:j)
     {
      for (l in 1:ncol(tB))
      {
       for (k in 1:l)
       {
        m<-k+(l*(l-1))/2
        n<-i+(j*(j-1))/2
        if (l==j && k==l)
        {
         D[m,n]<-2*tU[i,j]
        }

        if (l!=j && k==j)
        {
         D[m,n]<-tU[i,l]
        }

        if (k>=i && l==j && k!=l)
        {
         D[m,n]<-tU[i,k]
        }

       }
      }
     }
    }
    
   D <- D[-1,-1]  #car premier parametre pas estime
  }

  if (idiag==TRUE)
  {
   D<- matrix(0,nrow=nvc,ncol=nvc)
   diag(D)<- 2*diag(B)[-1]
  }

  # Delta-Method
  VFC<- D%*%Vc%*%t(D)


  # Tableau presentant les resultats
  #Results<- matrix(nrow=length(cholesky),ncol=4)
  Results<- matrix(nrow=nvc,ncol=4)
  colnames(Results)<- c("coef","Se","Wald test","p_value")

  # Noms des variables presents dans Xnames[Mod$idea==1]
   nom<-Mod$Xnames[Mod$idea==1]

  n<-c(0:nrow(Results))
  if(!isTRUE(idiag))
  {
   for(i in 1:length(nom))
   {
    for (j in 1:i)
    {if (i!=j)
     {
      n[i*(i-1)/2+j]<-paste(" Cov(",nom[j],"," ,nom[i], ")" ,sep="")
     }
     else
     {
      n[i*(i-1)/2+j]<-paste(" Var(",nom[i],")",sep="")
     }
    }
   }
  }
  else
  {
   n <- paste(" Var(",nom,")",sep="")
  }
  
  rownames(Results) <- n[-1]

  Results[,1]<- round(Mod$best[debut:fin],5)
  Results[,2]<- round(sqrt(diag(VFC)),5)
  # Le test de Wald est a effectuer uniquement sur les covariances
  # On commence par le calculer pour tout le monde, ensuite on remplace les resultats des variances par NA
  # Le meme principe est applique pour le calcul des p-values
  Results[,3]<- round(Results[,1]/Results[,2],5)

  # Le calcul de la p-value differe selon le signe de la stat de test de Wald
  # On va traiter un seul cas en utilisant la valeur absolue de la stat du test de Wald
  Results[,4]<-round(2*(1-pnorm(abs(Results[,3]))),5)
  if(!isTRUE(idiag))
  {
   for (i in 1:nea)
   {
     pos<-i*(i+1)/2
       {
       Results[pos-1,3]<-NA    # suit une N(0,1) sous H0
       Results[pos-1,4]<-NA
       }
   }
  }
  else #cas diagonal, on a que des variances
  {
   Results[,3:4] <- NA
  }
  print(Results,na.print="")
  return(invisible(Results))

}

VarCovRE <- function(Mod) UseMethod("VarCovRE")
