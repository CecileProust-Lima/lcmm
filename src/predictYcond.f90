
subroutine predictcondmult(X0,nalea,ny,nerr,maxmes,npm,b1,debut,epsY,idlink &
     ,nbzitr,zitr0,nsim,Ymarg)  

  use optim
  IMPLICIT NONE

  ! in input
  integer,intent(in)::maxmes,npm,nsim,nalea,ny,nerr,debut
  double precision,dimension(maxmes*ny),intent(in) ::X0
  integer,dimension(ny),intent(in)::idlink,nbzitr
  double precision,dimension(npm),intent(in)::b1
  double precision,dimension(ny),intent(in) ::epsY
  double precision,dimension(maxval(nbzitr),ny),intent(in)::zitr0

  ! for output
  double precision,dimension(maxmes*ny),intent(out) ::Ymarg

  ! for computation
  integer ::j,k,kk,m,l,jj,j1,j2,yk,sumntrtot
  integer ::ier,niter
  double precision,dimension(:,:),allocatable::VC
  double precision,dimension(:),allocatable::Vi,ysim,usim
  double precision :: eps,x22
  double precision ::ytemp,diff,SX,beta
  double precision::aa1,bb1,dd1,aa,bb,cc1
  double precision,dimension(:),allocatable::zitr,splaa
  double precision::xinbta,INV_ISPLINES,beta_ln
  double precision,dimension(ny)::minY,maxY
  integer,dimension(ny)::ntrtot

  allocate(ysim(maxmes*ny),usim(maxmes*ny),Vi(maxmes*ny*(maxmes*ny+1)/2), &
       VC(maxmes*ny,maxmes*ny))


  do yk=1,ny 

     minY(yk)=zitr0(1,yk)
     maxY(yk)=zitr0(nbzitr(yk),yk)

     !  if(verbose==1) print*,"min=",minY(yk),"max=",maxY(yk)

     if (idlink(yk).eq.0) ntrtot(yk)=2
     if (idlink(yk).eq.1) ntrtot(yk)=4
     if (idlink(yk).eq.2) ntrtot(yk)=nbzitr(yk)+2
  end do

  !  if(verbose==1) print*,"ntrtot=",ntrtot

  if(all(idlink/=2)) then
     allocate(zitr(1),splaa(1))
  else
     allocate(zitr(-1:maxval(ntrtot)),splaa(-1:(maxval(ntrtot)-3)))
  end if

  !  print*,"apres allocate2"
  zitr=0.d0
  splaa=0.d0
  ysim=0.d0
  usim=0.d0
  VC=0.d0

  eps=1.d-20

  ymarg=0.d0  


  do yk=1,ny
     do j1=1,maxmes
        do j2=1,maxmes  
           if(nalea.eq.ny) VC(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2) = VC(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2) +&
                b1(debut+ny+yk)**2 !variance de alpha_k

           if(j1.eq.j2) then
              if(nerr.gt.0) then
                 VC(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2) = VC(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2)+&
                      b1(debut+yk) **2 ! variance de l'erreur
              else
                 VC(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2) = VC(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2)+1 ! variance de l'erreur fixee a 1
              end if
           end if
        end do
     end do
  end do

  jj=0
  Vi=0.d0
  do j=1,maxmes*ny
     do k=j,maxmes*ny
        jj=j+k*(k-1)/2
        Vi(jj)=VC(j,k)
     end do
  end do

  CALL DMFSD(Vi,maxmes*ny,EPS,IER)
  if (ier.eq.-1) then
     ymarg=9999.d0
     goto 654
  end if

  VC=0.d0
  do j=1,maxmes*ny
     do k=1,j
        VC(j,k)=Vi(k+j*(j-1)/2)
     end do
  end do



  do l=1,nsim

     if(any(idlink.ne.0)) then
        usim=0.d0
        ysim=0.d0
        do m=1,maxmes*ny
           SX=1.d0
           call bgos(SX,0,usim(m),x22,0.d0)
        end do
        ysim=X0+MATMUL(VC,usim)
     end if

     sumntrtot=0
     do yk=1,ny
        !          print*,"boucle sur idlink, yk=",yk, "nsim=",l
        if (idlink(yk).eq.0 .and. l.eq.1) then  ! Linear link

           aa = b1(debut+ny+nalea+sumntrtot+1)
           bb = b1(debut+ny+nalea+sumntrtot+2)
           do j=1,maxmes
              Ymarg(maxmes*(yk-1)+j) = X0(maxmes*(yk-1)+j)*bb+aa
           end do

           !sumntrtot = sumntrtot +ntrtot(yk) 

        else if (idlink(yk).eq.1) then  ! Beta link

           aa1=exp(b1(debut+ny+nalea+sumntrtot+1))/ &
                (1+exp(b1(debut+ny+nalea+sumntrtot+1)))
           bb1=exp(b1(debut+ny+nalea+sumntrtot+2))/ &
                (1+exp(b1(debut+ny+nalea+sumntrtot+2)))
           bb1=aa1*(1.d0-aa1)*bb1

           cc1=abs(b1(debut+ny+nalea+sumntrtot+3))

           dd1=abs(b1(debut+ny+nalea+sumntrtot+4))

           aa=aa1*aa1*(1-aa1)/bb1-aa1
           bb=aa*(1-aa1)/aa1
           beta=beta_ln(aa,bb)


           do j=1,maxmes
              ytemp=ysim(maxmes*(yk-1)+j)*dd1+cc1
              if (ytemp.lt.0) then
                 ytemp=0.d0
              end if
              if (ytemp.gt.1) then
                 ytemp=1.d0
              end if
              ier=0
              ymarg(maxmes*(yk-1)+j)=ymarg(maxmes*(yk-1)+j)+xinbta(aa,bb,beta,ytemp,ier)/dble(nsim)
              if (ier.ne.0.or.ymarg(maxmes*(yk-1)+j).eq.9999.d0) then
                 ymarg(maxmes*(yk-1)+j)=9999.d0
                 !                    if (verbose==1) print*,"pb beta"
              end if
           end do

           !sumntrtot = sumntrtot +ntrtot(yk) 

        else if (idlink(yk).eq.2) then ! Splines link

           zitr=0.d0
           zitr(1:nbzitr(yk))=zitr0(1:nbzitr(yk),yk)
           zitr(-1)=zitr(1)
           zitr(0)=zitr(1)
           zitr(ntrtot(yk)-1)=zitr(ntrtot(yk)-2)
           zitr(ntrtot(yk))=zitr(ntrtot(yk)-1)

           !       if(verbose==1 .and. l==1) print*,"yk=",yk,"zitr=",zitr

           bb=b1(debut+ny+nalea+sumntrtot+1)
           do kk=2,ntrtot(yk)
              splaa(kk-3)=b1(debut+ny+nalea+sumntrtot+kk)**2
           end do

           do j=1,maxmes
              niter=0
              diff=0.d0
              ier=0
              ytemp=INV_ISPLINES(ysim(maxmes*(yk-1)+j),splaa,bb,nbzitr(yk),zitr,ier,niter,diff)
              if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3).or.ymarg(maxmes*(yk-1)+j).eq.9999.d0) then
                 ymarg(maxmes*(yk-1)+j)=9999.d0
                 !                       if (verbose==1) print*,"pb inversion splines"
              else
                 ymarg(maxmes*(yk-1)+j)=ymarg(maxmes*(yk-1)+j)+ytemp/dble(nsim)
              end if

           end do

           !sumntrtot = sumntrtot+ntrtot(yk)     

        end if
        sumntrtot = sumntrtot+ntrtot(yk) 
     end do
  end do


  do yk=1,ny
     ! pour les betas !
     if (idlink(yk).eq.1) then
        do j=1,maxmes
           if (ymarg(maxmes*(yk-1)+j).ne.9999.d0) then
              ymarg(maxmes*(yk-1)+j)=ymarg(maxmes*(yk-1)+j)*(maxY(yk)-minY(yk)+2*epsY(yk))+minY(yk)-epsY(yk)
              if (ymarg(maxmes*(yk-1)+j).lt.minY(yk)) ymarg(maxmes*(yk-1)+j)=minY(yk)
              if (ymarg(maxmes*(yk-1)+j).gt.maxY(yk)) ymarg(maxmes*(yk-1)+j)=maxY(yk)
           end if
        end do
     end if

  end do


654 continue
  !        print*,"avant delocate"
  deallocate(zitr,splaa)
  deallocate(ysim,usim,VC,Vi)


  return

end subroutine predictcondmult
