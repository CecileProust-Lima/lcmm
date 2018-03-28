risq_spl <-
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
                        Tmm3 <- (4*(zz[l+1]-t[j])**3)/((zz[l+1]-zz[l])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l-2])*(zz[l+1]-zz[l-3]))

                        Tmm2 <- (4*(t[j]-zz[l-2])*(zz[l+1]-t[j])**2)/((zz[l+2]-zz[l-2])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])*(zz[l+1]-zz[l-2])) -
                            (4*(t[j]-zz[l+2])*(t[j]-zz[l-1])*(zz[l+1]-t[j]))/((zz[l+2]-zz[l-2])*(zz[l+2]-zz[l-1])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])) +
                                (4*(t[j]-zz[l+2])**2*(t[j]-zz[l]))/((zz[l+2]-zz[l-2])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l])*(zz[l+2]-zz[l-1]))


                        Tmm1 <- (4*(t[j]-zz[l-1])**2*(zz[l+1]-t[j]))/((zz[l+3]-zz[l-1])*(zz[l+2]-zz[l-1])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])) -
                            (4*(t[j]-zz[l-1])*(t[j]-zz[l])*(t[j]-zz[l+2]))/((zz[l+3]-zz[l-1])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l])*(zz[l+2]-zz[l-1])) +
                                (4*(zz[l+3]-t[j])*(t[j]-zz[l])**2)/((zz[l+3]-zz[l-1])*(zz[l+3]-zz[l])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l]))
                        
                        Tmm <- (4*(t[j]-zz[l])**3)/((zz[l+4]-zz[l])*(zz[l+3]-zz[l])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l]))
                    }

                if(t[j]==z[nz])
                    {
                        Tmm3 <- 0
                        Tmm2 <- 0
                        Tmm1 <- 0
                        Tmm <- (4*(t[j]-zz[l])**3)/((zz[l+4]-zz[l])*(zz[l+3]-zz[l])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l]))
                    }


                res[j,(l-3):l] <- c(Tmm3,Tmm2,Tmm1,Tmm)
                
            }

        return(res%*%b)
    }
