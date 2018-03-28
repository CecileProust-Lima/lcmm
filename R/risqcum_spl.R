risqcum_spl <-
function(t,z,b)
    {
        ## t vecteur de temps
        ## z vecteur de noeuds (sans repetition des noeuds externes)
        ## b vectuer des coef beta


        nz <- length(z)
        z <- sort(z)
        if(any(t<z[1]) | any(t>z[nz])) stop("t en dehors de z")

        ## repeter les noeuds externes
        zz <- rep(0,nz+6)
        zz[1:3] <- rep(z[1],3)
        zz[3+1:nz] <- z
        zz[3+nz+1:3] <- z[nz]

        ## declaration du resultat
        res <- matrix(0,nrow=length(t),ncol=nz+2)
        
        ## boucle sur les temps
        for(j in 1:length(t))
            {
                ## encadrer t par les noeuds
                tz <- c(t[j],z)
                itz <- c(1,rep(0,nz))
                l <- 3+which(itz[order(tz)]==1)-1

                ## si t=min, on laisse le resultat a zero
                if(t[j]==z[1]) next


                ## calcul des Tmm(t) non nuls

                if(t[j]<z[nz])
                    {
                        Tim3 <- ((t[j]-zz[l-3])*(zz[l+1]-t[j])**3)/((zz[l+1]-zz[l])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l-2])*(zz[l+1]-zz[l-3])) +
                            ((t[j]-zz[l-2])*(zz[l+1]-t[j])**2)/((zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])*(zz[l+1]-zz[l-2])) -
                                ((t[j]-zz[l+2])*(t[j]-zz[l-1])*(zz[l+1]-t[j]))/((zz[l+2]-zz[l-1])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])) +
                                    ((t[j]-zz[l+2])**2*(t[j]-zz[l]))/((zz[l+2]-zz[l])*(zz[l+1]-zz[l])*(zz[l+2]-zz[l-1])) +                        
                                        ((t[j]-zz[l-1])**2*(zz[l+1]-t[j]))/((zz[l+2]-zz[l-1])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])) -
                                            ((t[j]-zz[l-1])*(t[j]-zz[l])*(t[j]-zz[l+2]))/((zz[l+2]-zz[l])*(zz[l+1]-zz[l])*(zz[l+2]-zz[l-1])) +
                                                ((zz[l+3]-t[j])*(t[j]-zz[l])**2)/((zz[l+3]-zz[l])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l])) +
                                                    ((t[j]-zz[l])**3)/((zz[l+3]-zz[l])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l]))

                        Tim2 <- ((t[j]-zz[l-2])**2*(zz[l+1]-t[j])**2)/((zz[l+2]-zz[l-2])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])*(zz[l+1]-zz[l-2])) -
                            ((t[j]-zz[l-2])*(t[j]-zz[l+2])*(t[j]-zz[l-1])*(zz[l+1]-t[j]))/((zz[l+2]-zz[l-2])*(zz[l+2]-zz[l-1])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])) +
                                ((t[j]-zz[l-2])*(t[j]-zz[l+2])**2*(t[j]-zz[l]))/((zz[l+2]-zz[l-2])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l])*(zz[l+2]-zz[l-1])) +
                                    ((t[j]-zz[l-1])**2*(zz[l+1]-t[j]))/((zz[l+2]-zz[l-1])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])) -
                                        ((t[j]-zz[l-1])*(t[j]-zz[l])*(t[j]-zz[l+2]))/((zz[l+2]-zz[l])*(zz[l+1]-zz[l])*(zz[l+2]-zz[l-1])) +
                                            ((zz[l+3]-t[j])*(t[j]-zz[l])**2)/((zz[l+3]-zz[l])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l])) +
                                                ((t[j]-zz[l])**3)/((zz[l+3]-zz[l])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l]))


                        Tim1 <- ((t[j]-zz[l-1])**3*(zz[l+1]-t[j]))/((zz[l+3]-zz[l-1])*(zz[l+2]-zz[l-1])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])) -
                            ((t[j]-zz[l-1])**2*(t[j]-zz[l])*(t[j]-zz[l+2]))/((zz[l+3]-zz[l-1])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l])*(zz[l+2]-zz[l-1])) +
                                ((t[j]-zz[l-1])*(zz[l+3]-t[j])*(t[j]-zz[l])**2)/((zz[l+3]-zz[l-1])*(zz[l+3]-zz[l])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l])) +
                                    ((t[j]-zz[l])**3)/((zz[l+3]-zz[l])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l]))


                        Tim <- ((t[j]-zz[l])**4)/((zz[l+4]-zz[l])*(zz[l+3]-zz[l])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l]))
                        
                                                                              

                        if(t[j]>z[2]) res[j,1:(l-4)] <- 1
                        res[j,(l-3):l] <- c(Tim3,Tim2,Tim1,Tim)     
                    }
            

                if(t[j]==z[nz])
                    {
                        Tim3 <- (zz[l+4]-zz[l])/(z[nz]-z[nz-1])
                        Tim2 <- (zz[l+4]-zz[l])/(z[nz]-z[nz-1])
                        Tim1 <- (zz[l+4]-zz[l])/(z[nz]-z[nz-1])
                        Tim <- (t[j]-zz[l])/(z[nz]-z[nz-1])

                        res[j,1:(l-4)] <- 1
                        res[j,(l-3):l] <- c(Tim3,Tim2,Tim1,Tim)
                    }
                
            }

        return(res%*%b)
    }
