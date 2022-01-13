transfo_spl <- function(y,z,b)
    {
        nz <- length(z)
        z <- sort(z)
        if(any(na.omit(y)<z[1]) | any(na.omit(y)>z[nz])) stop("y en dehors de z")

        b <- c(b[1],b[-1]**2)

        ## repeter les noeuds externes
        zz <- rep(0,nz+4)
        zz[1:2] <- rep(z[1],2)
        zz[2+1:nz] <- z
        zz[2+nz+1:2] <- z[nz]

        ## declaration du resultat
        res <- matrix(0,nrow=length(y),ncol=nz+1)

        ## boucle sur les y
        for(j in 1:length(y))
            {
                if(!is.na(y[j]))
                    {
                        ## encadrer y par les noeuds
                        yz <- c(y[j],z)
                        iyz <- c(1,rep(0,nz))
                        l <- 2+which(iyz[order(yz)]==1)-1
                        
                        ## si y=min, on laisse le resultat a zero
                        if(y[j]==z[1]) next


                        ## calcul des im(y) non nuls
                        
                        if(y[j]<z[nz])
                            {
                                im2 <- ((y[j]-zz[l-2])*(zz[l+1]-y[j])**2)/((zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])*(zz[l+1]-zz[l-2])) +
                                    ((y[j]-zz[l-1])*(zz[l+1]-y[j]))/((zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])) +
                                        ((y[j]-zz[l])*(zz[l+2]-y[j]))/((zz[l+2]-zz[l])*(zz[l+1]-zz[l])) +
                                            ((y[j]-zz[l])**2)/((zz[l+2]-zz[l])*(zz[l+1]-zz[l]))
                                
                                im1 <- ((y[j]-zz[l-1])**2*(zz[l+1]-y[j]))/((zz[l+2]-zz[l-1])*(zz[l+1]-zz[l-1])*(zz[l+1]-zz[l])) +
                                    ((y[j]-zz[l-1])*(y[j]-zz[l])*(zz[l+2]-y[j]))/((zz[l+2]-zz[l])*(zz[l+1]-zz[l])*(zz[l+2]-zz[l-1])) +
                                        ((y[j]-zz[l])**2)/((zz[l+2]-zz[l])*(zz[l+1]-zz[l]))
                                
                                im <- ((y[j]-zz[l])**3)/((zz[l+3]-zz[l])*(zz[l+2]-zz[l])*(zz[l+1]-zz[l]))

                                if(y[j]>z[2]) res[j,1:(l-3)] <- 1
                                res[j,(l-2):(l)] <- c(im2,im1,im)
                            }
                        
                        if(y[j]==z[nz])
                            {
                                res[j,] <- 1
                            }
                    }
                else
                    {
                        res[j,] <- NA
                    }
            }

        return(b[1]+res%*%b[-1])
    }
