VarCovRE.Jointlcmm <- function(Mod)
{
  if(missing(Mod)) stop("The model should be specified")
  if(!inherits(Mod,"Jointlcmm")) stop("applies to \"Jointlcmm\" objects only")

   nea <- sum(Mod$idea) # Nombre d'effets aleatoires
   nef <- Mod$N[4]      # Nombre d'effets fixes
   nvc <- Mod$N[5]      # Nombre de parametres  effets aleatoires
   nprob <- Mod$N[1]


  if(nvc==0) return(NA)

  # Recuperation des parametres de Cholesky dans un vecteur
  cholesky<-Mod$cholesky

  # Transformation de Cholesky telle que UU'=B
  tU <- matrix(0,nrow=nea,ncol=nea)
  if(Mod$idiag==0)
      {
          tU[upper.tri(tU,diag=TRUE)] <- cholesky
      }
  if(Mod$idiag==1)
      {
          diag(tU) <- cholesky
      }

  U <- t(tU)

  B <- U%*%tU   # Matrice de variance-covariance des effets aleatoires

  # On rend symetrique la matrice de variance-covariance V des parametres
  l <- length(Mod$best)
  V <- matrix(0,nrow=l,ncol=l)
  V[upper.tri(V,diag=TRUE)] <- Mod$V
  V[lower.tri(V,diag=FALSE)] <- t(V)[lower.tri(V,diag=FALSE)]


  # On recupere la sous matrice de var-cov des parametres de cholesky dans V
  debut <- sum(Mod$N[1:4])+1
  fin <- sum(Mod$N[1:4])+nvc
  Vc <- V[debut:fin,debut:fin]

  # Calcul de la derivee de B en fonction des parametres de cholesky : D=f'(c)/c
  # D est different selon la valeur de idiag
  if (Mod$idiag==0)
      {
          # calcul de derivee
          B[lower.tri(B,diag=FALSE)] <- 0
          tB <- B
          Sup <- B[upper.tri(B,diag=TRUE)]
          D <- matrix(0,nrow=length(Sup),ncol=length(Sup))
          
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
      }
    
  if (Mod$idiag==1)
      {
          D<- matrix(0,nrow=nvc,ncol=nvc)
          diag(D)<- 2*diag(B)
      }
  


  # Delta-Method
  VFC<- D%*%Vc%*%t(D)


  # Tableau presentant les resultats
  #Results<- matrix(nrow=length(cholesky),ncol=4)
  Results <- matrix(nrow=nvc,ncol=4)
  colnames(Results) <- c("coef","Se","Wald test","p_value")

  # Noms des variables presents dans Xnames[Mod$idea==1]
   nom <- Mod$Names$Xnames[Mod$idea==1]

  n<-c(1:nrow(Results))
  if(!isTRUE(Mod$idiag))
      {
          for(i in 1:length(nom))
              {
                  for (j in 1:i)
                      if (i!=j)
                          {
                              n[i*(i-1)/2+j]<-paste(" Cov(",nom[j],"," ,nom[i], ")" ,sep="")
                          }
                      else
                          {
                              n[i*(i-1)/2+j]<-paste(" Var(",nom[i],")",sep="")
                          }
              }
      }
  else
      {
          n <- paste(" Var(",nom,")",sep="")
      }
   
  rownames(Results) <- n

  Results[,1]<- round(Mod$best[debut:fin],5)
  Results[,2]<- round(sqrt(diag(VFC)),5)
  # Le test de Wald est a effectuer uniquement sur les covariances
  # On commence par le calculer pour tout le monde, ensuite on remplace les resultats des variances par NA
  # Le meme principe est applique pour le calcul des p-values
  Results[,3]<- round(Results[,1]/Results[,2],5)

  # Le calcul de la p-value differe selon le signe de la stat de test de Wald
  # On va traiter un seul cas en utilisant la valeur absolue de la stat du test de Wald
  Results[,4]<-round(2*(1-pnorm(abs(Results[,3]))),5)
  if(!isTRUE(Mod$idiag))
  {
      for (i in 1:nea)
          {
              pos<-i*(i+1)/2
              {
                  Results[pos,3]<-NA    # suit une N(0,1) sous H0
                  Results[pos,4]<-NA
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
