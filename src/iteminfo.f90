
subroutine iteminfo(X0,condRE_Y,nalea,ny,maxmes,npm,b1,debut,nbzitr,idlink,&
    modalite,nbmod,nsim,ninfo,info)  

!  use optim
  IMPLICIT NONE

  ! in input
  integer,intent(in)::maxmes,npm,nsim,nalea,ny,debut,condRE_Y,ninfo
  double precision,dimension(maxmes*ny),intent(in) ::X0
  integer,dimension(ny),intent(in)::idlink,nbmod,nbzitr
  double precision,dimension(npm),intent(in)::b1
  integer, dimension(sum(nbmod)),intent(in)::modalite

  ! for output
  double precision,dimension(maxmes*ninfo),intent(out) ::info

  ! for computation
  integer::j,m,l,jj,yk,sumntrtot
  integer::sumnbmod
  double precision,dimension(:),allocatable::ysim
  double precision::x22,ai,asim
  double precision::ytemp,SX
  double precision::aa1,aa
  double precision,dimension(:),allocatable::grandPhi,petitphi
  integer,dimension(ny)::ntrtot
  double precision::alnorm

  m=nbmod(1)-1
  do yk=2,ny
     if(nbmod(yk)+1 .gt. m) m = nbmod(yk)-1
  end do
  
  allocate(ysim(maxmes*ny), grandPhi(maxmes*m), petitphi(maxmes*m))

  
  do yk=1,ny 
     if (idlink(yk).eq.0) ntrtot(yk)=2
     if (idlink(yk).eq.1) ntrtot(yk)=4
     if (idlink(yk).eq.2) ntrtot(yk)=nbzitr(yk)+2
     if (idlink(yk).eq.3) ntrtot(yk)=nbmod(yk)-1
  end do

  ysim=0.d0

  info=0.d0  

  do l=1,nsim

     sumntrtot=0
     sumnbmod = 0
     do yk=1,ny

        grandPhi=0.d0
        petitphi=0.d0

        if (idlink(yk).eq.3) then ! ordinal

           if(condRE_Y.eq.1 .and.l.eq.1) then
              !! sans randomY, pas besoin de MC

              !! calculer Phi((seuil-Lambda)/sigma) et phi((seuil-Lambda)/sigma)
              aa = b1(debut+ny+nalea+sumntrtot+1)
              do jj=1,nbmod(yk)-1
                 if(jj.gt.1) aa = aa + b1(debut+ny+nalea+sumntrtot+jj)**2
                 do j=1,maxmes
                    grandPhi(maxmes*(jj-1)+j) = alnorm((aa-X0(maxmes*(yk-1)+j))/abs(b1(debut+yk)),.false.)
                    petitphi(maxmes*(jj-1)+j) = exp(-((aa-X0(maxmes*(yk-1)+j))/abs(b1(debut+yk)))**2&
                         / 2.d0) / sqrt(dble(2*3.14159265))
                 end do
              end do
                            
              !! calculer P_kl et I_kl*P_kl
              !! pour modalite min
              aa = b1(debut+ny+nalea+sumntrtot+1)
              do j=1,maxmes
                 info(maxmes*(2*sumnbmod+yk-1)+j) = GrandPhi(j) ! P_kl
                 info(maxmes*(2*sumnbmod+yk-1)+maxmes*nbmod(yk)+j) = &
                      (petitphi(j)**2) / (abs(b1(debut+yk))**2 * GrandPhi(j)) + &
                      (aa-X0(maxmes*(yk-1)+j))/abs(b1(debut+yk))**3 * petitphi(j) ! I_klP_kl
                 info(maxmes*(2*sumnbmod+yk-1)+maxmes*2*nbmod(yk)+j) = &
                      info(maxmes*(2*sumnbmod+yk-1)+maxmes*nbmod(yk)+j) ! I_k
              end do

              !! pour modalite intermediaire
              if(nbmod(yk).gt.3) then
                 do jj=2,nbmod(yk)-1
                    aa1 = aa ! seuil precedent
                    aa = aa + b1(debut+ny+nalea+sumntrtot+jj)**2 ! seuil courant
                    do j=1,maxmes
                       info(maxmes*(2*sumnbmod+yk-1)+maxmes*(jj-1)+j) =  &
                            GrandPhi(maxmes*(jj-1)+j) - GrandPhi(maxmes*(jj-2)+j) ! P_kl
                       info(maxmes*(2*sumnbmod+yk-1)+maxmes*(nbmod(yk)+jj-1)+j) = &
                            ((petitphi(maxmes*(jj-1)+j) - petitphi(maxmes*(jj-2)+j))**2)/&
                            (abs(b1(debut+yk))**2 * info(maxmes*(2*sumnbmod+yk-1)+maxmes*(jj-1)+j)) + &
                            (aa-X0(maxmes*(yk-1)+j))/abs(b1(debut+yk))**3 * petitphi(maxmes*(jj-1)+j) - &
                            (aa1-X0(maxmes*(yk-1)+j))/abs(b1(debut+yk))**3 * petitphi(maxmes*(jj-2)+j)! I_klP_kl
                       info(maxmes*(2*sumnbmod+yk-1)+maxmes*2*nbmod(yk)+j) = &
                            info(maxmes*(2*sumnbmod+yk-1)+maxmes*2*nbmod(yk)+j) + &
                            info(maxmes*(2*sumnbmod+yk-1)+maxmes*(nbmod(yk)+jj-1)+j) ! I_k
                    end do
                 end do
              end if
              
              !! pour modalite max
              do j=1,maxmes
                 info(maxmes*(2*sumnbmod+yk-1)+maxmes*(nbmod(yk)-1)+j) = 1 - GrandPhi(maxmes*(nbmod(yk)-2)+j) ! P_kl
                 info(maxmes*(2*sumnbmod+yk-1)+maxmes*(2*nbmod(yk)-1)+j) = &
                      (petitphi(maxmes*(nbmod(yk)-2)+j)**2) / &
                      (abs(b1(debut+yk))**2 * info(maxmes*(2*sumnbmod+yk-1)+maxmes*(nbmod(yk)-1)+j)) - &
                      (aa-X0(maxmes*(yk-1)+j))/abs(b1(debut+yk))**3 * petitphi(maxmes*(nbmod(yk)-2)+j)! I_klP_kl
                 info(maxmes*(2*sumnbmod+yk-1)+maxmes*2*nbmod(yk)+j) = &
                      info(maxmes*(2*sumnbmod+yk-1)+maxmes*2*nbmod(yk)+j) + &
                      info(maxmes*(2*sumnbmod+yk-1)+maxmes*(2*nbmod(yk)-1)+j) ! I_k
              end do
              

              
              
           else if(condRE_Y.eq.0) then ! integrer sur randomY
              
              !! simuler l'EA specifique au test
              ai=0.d0
              call bgos(SX,0,asim,x22,0.d0)
              ai = b1(debut+ny+yk)*asim
              do j=1,maxmes
                 ysim(maxmes*(yk-1)+j) = X0(maxmes*(yk-1)+j) + ai
              end do
           
              aa = b1(debut+ny+nalea+sumntrtot+1)
              do j=1,maxmes
                 ytemp = modalite(sumnbmod + nbmod(yk)) - &
                      alnorm((aa-ysim(maxmes*(yk-1)+j))/abs(b1(debut+yk)),.false.)
                 info(maxmes*(yk-1)+j) = info(maxmes*(yk-1)+j) + ytemp / dble(nsim)
              end do
              
              if(nbmod(yk).gt.2) then
                 do jj=2,nbmod(yk)-1
                    aa = aa + b1(debut+ny+nalea+sumntrtot+jj)**2
                    do j=1,maxmes
                       ytemp = (modalite(sumnbmod + jj) - modalite(sumnbmod + jj+1))* &
                            alnorm((aa-ysim(maxmes*(yk-1)+j))/abs(b1(debut+yk)),.false.)
                       info(maxmes*(yk-1)+j) = info(maxmes*(yk-1)+j) + ytemp / dble(nsim)                            
                    end do
                 end do
              end if
           end if
           
           sumnbmod = sumnbmod + nbmod(yk)

        end if
        sumntrtot = sumntrtot+ntrtot(yk) 
     end do
  end do

  deallocate(ysim, grandPhi, petitphi)


  return

end subroutine iteminfo
