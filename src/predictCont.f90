






!====================================================================

! PROGRAMME PRINCIPAL

!-----------------------------------------------------------
!           PREDICTION by PROFILE of COVARIATES
!------------------------------------------------------------



!===================================================================


subroutine predictcont(X0,idprob,idea,idg,idcor &
     ,ng,ncor,nv,maxmes,idiag,nwg,npm,b1,epsY,idlink &
     ,nbzitr,zitr0,nsim,methInteg,Ydiscret,Ymarg)


!  use optim
  IMPLICIT NONE

  ! in input
  integer,intent(in)::ng,ncor,nv,maxmes,idiag,nwg,npm,idlink,nbzitr,nsim,Ydiscret,methInteg
  double precision,dimension(maxmes*nv),intent(in) ::X0
  integer,dimension(nv),intent(in)::idprob,idea,idg,idcor
  double precision,dimension(npm),intent(in)::b1
  double precision,intent(in) ::epsY
  double precision,dimension(nbzitr),intent(in)::zitr0

  ! for computation
  integer ::j,k,l,m,g,l2,m2,jj,ntrtot,npm2,j1,j2
  integer ::ier,nmoins,kk,nrisq,nvarxevt,niter,nvarprob,ncg,ncssg,nea,nef,nvc,nprob
!  double precision,dimension(nv) ::Xprob,  bprob
!  double precision::temp
!  double precision,dimension(ng) :: pi
  double precision,dimension(nv,nv) ::Ut,Ut1

  double precision,dimension(:,:),allocatable::VC,Z,P,R,X00,X2
  double precision,dimension(:),allocatable::Vi,ysim,usim,mu,tcor
  double precision,dimension(nv) :: b0,b2
  double precision :: eps,ytemp2,x22
  double precision ::ytemp,diff,SX,beta
  double precision::aa1,bb1,dd1,aa,bb,cc1,rangeY,minY,maxY
  double precision,dimension(:),allocatable::zitr,splaa
  double precision::xinbta,INV_ISPLINES,beta_ln
  double precision,dimension(2,51)::gauss

  ! for output
  double precision,dimension(maxmes*ng),intent(out) ::Ymarg

  
  !! call GetRNGstate
  call getrand()


  !============= recup des places de parametres

  allocate(ysim(maxmes),usim(maxmes),mu(maxmes),tcor(maxmes),Vi(maxmes*(maxmes+1)/2), &
       VC(maxmes,maxmes),Z(maxmes,nv),P(maxmes,nv),R(maxmes,maxmes),X00(maxmes,nv),X2(maxmes,nv))

  ! en prevision de l'extension au conjoint
  nrisq=0
  nvarxevt=0
  ! fin en prevision

  minY=zitr0(1)
  maxY=zitr0(nbzitr)

  rangeY=0
  if (Ydiscret.eq.1) rangeY=maxY-minY

  if (idlink.eq.0) ntrtot=2
  if (idlink.eq.1) ntrtot=4
  if (idlink.eq.3) ntrtot=INT(zitr0(nbzitr))-INT(zitr0(1))
  if (idlink.eq.2) then
     ntrtot=nbzitr+2
     allocate(zitr(-1:(ntrtot)),splaa(-1:(ntrtot-3)))
     zitr(1:nbzitr)=zitr0(1:nbzitr)
     zitr(-1)=zitr(1)
     zitr(0)=zitr(1)
     zitr(ntrtot-1)=zitr(ntrtot-2)
     zitr(ntrtot)=zitr(ntrtot-1)
  else
     allocate(zitr(1),splaa(1))
  end if


  eps=1.d-20


!  write(*,*)'entree'
!  write(*,*)idprob
!  write(*,*)idea
!  write(*,*)idg
!  write(*,*)ng,nv,maxmes,idiag,nwg,npm
!  write(*,*)b1
!  write(*,*)'epsY',epsY,idlink,nbzitr
!  write(*,*)'zitr',zitr0
!  write(*,*)nsim,methInteg,Ydiscret
!  write(*,*)'revu'




  ! creation des parametres

  nea=0
  ncg=0
  ncssg=0
  nprob=0 !ng-1
  nvarprob=min(ng-1,1)
  do k=1,nv
     if (idg(k).eq.1) then
        ncssg=ncssg+1      ! nb var. sans melange
     else if (idg(k).eq.2) then
        ncg=ncg+1      ! nb var. sans melange
     end if
     nea=nea+idea(k)
     nprob=nprob+(idprob(k))*(ng-1)
     nvarprob=nvarprob+idprob(k)
  end do

  if((ng.eq.1.and.ncg.gt.0).or.(ng.eq.1.and.nprob.gt.0)) then
     ymarg=9999.d0
     go to 654
  end if


  !  nb effets fixes = nb effets fixes sans melange
  !                  + ng fois le nb de var dans melange

  nvc=nea
  if(idiag.eq.0) then
     nvc=nea*(nea+1)/2
  end if

  nef=nprob+ncssg+ncg*ng+nvarxevt+nrisq-1
  npm2=nef+nvc+nwg+ntrtot+ncor



  if (npm.ne.npm2) then
     ymarg=9999.d0
     goto 654
  end if


!========== debut des calculs =======================
!   Variance of random-effects


  Ut=0.d0
  If (idiag.eq.1) then
     do j=1,nea
        do k=1,nea
           if (j.eq.k) then
              Ut(j,k)=b1(nef+j)
           else
              Ut(j,k)=0.d0
           end if
        end do
     end do
  end if

  If (idiag.eq.0) then
     do j=1,nea
        do k=1,j
           Ut(j,k)=b1(nef+k+j*(+j-1)/2)
        end do
     end do
  end if
  Ut1=Ut



  ymarg=0.d0





! -------- creation de Vi = ZiGZi'+se*seIni ----------
! creation de Zi

  Z=0.d0
  l=0
  do k=1,nv
     if (idea(k).eq.1) then
        l=l+1
        do j=1,maxmes
           Z(j,l)=dble(X0(maxmes*(k-1)+j))
        end do
     end if
  end do
  
  ! creation de Ri
        R=0.d0
        tcor=0.d0
        if (ncor.gt.0) then
           do k=1,nv
              if (idcor(k).eq.1) then
                 do j=1,maxmes
                    tcor(j) = X0(maxmes*(k-1)+j)
                 end do
              end if
           end do
         end if
         do j1=1,maxmes
            do j2=1,maxmes
               if (ncor.eq.1) then 
                  R(j1,j2) = b1(npm)*b1(npm)*min(tcor(j1),tcor(j2))
               else if (ncor.eq.2) then
                  R(j1,j2) = b1(npm)*b1(npm)*exp(-b1(npm-1)*abs(tcor(j1)-tcor(j2)))
               end if
            end do
         end do     




! creation de P=Zi*Ut et V=P*P' que si non spec aux classes

  if (nwg.eq.0.OR.NG.EQ.1) then
     P=0.d0
     P=MATMUL(Z,Ut)
     VC=0.d0
     VC=MATMUL(P,transpose(P))+R
     ! j'ajoute l'erreur a 1 apres


     ! Vi en vecteur

     jj=0
     Vi=0.d0
     do j=1,maxmes
        VC(j,j)=VC(j,j)+1.d0
        do k=j,maxmes
           jj=j+k*(k-1)/2
           Vi(jj)=VC(j,k)
           ! ajouter 1 si diagonale
        end do
     end do

     if (methInteg.eq.1) then 
        CALL DMFSD(Vi,maxmes,EPS,IER)
        if (ier.eq.-1) then
           ymarg=9999.d0
           goto 654
        end if
        
        VC=0.d0
        do j=1,maxmes
           do k=1,j
              VC(j,k)=Vi(k+j*(j-1)/2)
           end do
        end do
     end if
  end if
     

  ! contribution individuelle a la vraisemblance

  ! cas 1 : ng=1

  if (ng.eq.1) then

     b0=0.d0
     l=0
     X00=0.d0
     do k=1,nv
        if (idg(k).ne.0) then
           l=l+1
           do j=1,maxmes
              X00(j,l)=dble(X0(maxmes*(k-1)+j))
           end do
           ! idg ne 0 pour l'iintercept forcement donc on met le parm a 0
           if (k.eq.1) then
              b0(l)=0.d0
           else
              b0(l)=b1(nprob+l-1)
           end if
        end if
     end do

     mu=0.d0
     mu=matmul(X00,b0)



     if (idlink.eq.0) then  ! Linear link
        ymarg=9999.d0
        goto 654
        ! ne devrait pas venir ici

     else if (idlink.eq.1) then  ! Beta link

        aa1=exp(b1(nef+nvc+nwg+1))/ &
             (1+exp(b1(nef+nvc+nwg+1)))
        bb1=exp(b1(nef+nvc+nwg+2))/ &
             (1+exp(b1(nef+nvc+nwg+2)))
        bb1=aa1*(1.d0-aa1)*bb1

        cc1=abs(b1(nef+nvc+nwg+3))

        dd1=abs(b1(npm-ncor))

        aa=aa1*aa1*(1-aa1)/bb1-aa1
        bb=aa*(1-aa1)/aa1
        beta=beta_ln(aa,bb)



        if (methInteg.eq.1) then
           ! i.e. methode de MC a faire
           do l=1,nsim
              usim=0.d0
              ysim=0.d0
              do m=1,maxmes
                 SX=1.d0
                 call bgos(SX,0,usim(m),x22,0.d0)
              end do

              ysim=mu+MATMUL(VC,usim)
              do j=1,maxmes
                 ytemp=ysim(j)*dd1+cc1
                 if (ytemp.lt.0) then
                    ytemp=0.d0
                 end if
                 if (ytemp.gt.1) then
                    ytemp=1.d0
                 end if
                 ymarg(j)=ymarg(j)+xinbta(aa,bb,beta,ytemp,ier)/dble(nsim)
                 if (ier.ne.0.or.ymarg(j).eq.9999.d0) then
                    ymarg(j)=9999.d0
                 end if
              end do
           end do

        else

           call gausshermite(gauss,nsim)
           do j=1,maxmes
              do l=1,nsim
                 ytemp=mu(j)+sqrt(VC(j,j))*gauss(1,l)

!                 ytemp=X0(maxmes+j)

                 ytemp=ytemp*dd1+cc1

                 if (ytemp.lt.0) then
                    ytemp=0.d0
                 end if
                 if (ytemp.gt.1) then
                    ytemp=1.d0
                 end if
                 ymarg(j)=ymarg(j)+xinbta(aa,bb,beta,ytemp,ier)*gauss(2,l)
                if (ier.ne.0.or.ymarg(j).eq.9999.d0) then
                    ymarg(j)=9999.d0
                 end if
              end do
           end do

        end if

     else if (idlink.eq.2) then ! Splines link


        bb=b1(nef+nvc+nwg+1)
        do kk=2,ntrtot
           splaa(kk-3)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
        end do

        if (methInteg.eq.1) then
           ! i.e. methode de MC a faire
           do l=1,nsim
              usim=0.d0
              ysim=0.d0
              do m=1,maxmes
                 SX=1.d0
                 call bgos(SX,0,usim(m),x22,0.d0)
              end do

              ysim=mu+MATMUL(VC,usim)
              do j=1,maxmes
                 niter=0
                 diff=0.d0
                 ier=0
                 ytemp=INV_ISPLINES(ysim(j),splaa,bb,nbzitr,zitr,ier,niter,diff)
                    if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3).or.ymarg(j).eq.9999.d0) then
                       ymarg(j)=9999.d0
                    else
                        ymarg(j)=ymarg(j)+ytemp/dble(nsim)
                    end if

              end do
           end do


        else

           call gausshermite(gauss,nsim)


           do j=1,maxmes
              do l=1,nsim
                 niter=0
                 diff=0.d0
                 ier=0
                 ytemp=mu(j)+sqrt(VC(j,j))*gauss(1,l)

!                ytemp=X0(maxmes+j)

                 ytemp2=INV_ISPLINES(ytemp,splaa,bb,nbzitr,zitr,ier,niter,diff)
                   if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3).or.ymarg(j).eq.9999.d0) then
                       ymarg(j)=9999.d0
                    else
                        ymarg(j)=ymarg(j)+ytemp2*gauss(2,l)
                    end if

              end do
           end do
        end if
     end if


! cas 2 :  ng>1  composantes
  else


     ! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
!     Xprob=0.d0
!     Xprob(1)=1
!     l=0
!     do k=1,nv
!        if (idprob(k).eq.1) then
!           l=l+1
!           Xprob(1+l)=X0(maxmes*(k-1)+1)
!        end if
!     end do
!     pi=0.d0
!     temp=0.d0
!     Do g=1,ng-1
!        bprob=0.d0
!        do k=1,nvarprob
!           bprob(k)=b1((k-1)*(ng-1)+g)
!        end do

!        temp=temp+exp(DOT_PRODUCT(bprob,Xprob))

!        pi(g)=exp(DOT_PRODUCT(bprob,Xprob))

!     end do

!     pi(ng)=1/(1+temp)

!     do g=1,ng-1
!        pi(g)=pi(g)*pi(ng)
!     end do

     ! creation des vecteurs de variables explicatives
     l=0
     m=0
     X00=0.d0
     X2=0.d0
     do k=1,nv
        if (idg(k).eq.2) then
           l=l+1
           do j=1,maxmes
              X2(j,l)=dble(X0(maxmes*(k-1)+j))
           end do
        else if (idg(k).eq.1) then
           m=m+1
           do j=1,maxmes
              X00(j,m)=dble(X0(maxmes*(k-1)+j))
           end do
        end if
     end do

     b2=0.d0
     b0=0.d0
     do g=1,ng
        nmoins=0
        l2=0
        m2=0
        do k=1,nv
           if (idg(k).eq.1) then
              ! parametre a 0 pour l'intercept
              if (k.eq.1) then
                 m2=m2+1
                 b0(m2)=0.d0
              else
                 m2=m2+1
                 b0(m2)=b1(nprob+nmoins+1)
                 nmoins=nmoins+1
              end if
           else if (idg(k).eq.2) then
              ! parametre a 0 pour l'intercept de la premiere classe
              if (k.eq.1) then
                 if (g.eq.1) then
                    l2=l2+1
                    b2(l2)=0.d0
                    nmoins=nmoins+ng-1
                 else
                    l2=l2+1
                    b2(l2)=b1(nprob+nmoins+g-1)
                    nmoins=nmoins+ng-1
                 end if
              else
                 l2=l2+1
                 b2(l2)=b1(nprob+nmoins+g)
                 nmoins=nmoins+ng
              end if
           end if
        end do

!        write(*,*)'b2',b2
!        write(*,*)'b0',b0

        ! variance covariance si spec aux classes :

        Ut1=Ut
        if (nwg.ne.0) then
           Ut1=0.d0
           if (g.eq.ng) then
              Ut1=Ut
           else
              Ut1=Ut*abs(b1(nef+nvc+g))
           end if


           P=0.d0
           P=MATMUL(Z,Ut1)
           VC=0.d0
           VC=MATMUL(P,transpose(P))+R

           ! Vi en vecteur
           Vi=0.d0
           jj=0
           do j=1,maxmes
              VC(j,j)=VC(j,j)+1.d0
              do k=j,maxmes
                 jj=j+k*(k-1)/2
                 Vi(jj)=VC(j,k)
              end do
           end do

           if (methInteg.eq.1) then 
              CALL DMFSD(Vi,maxmes,EPS,IER)
              if (ier.eq.-1) then
                 ymarg=9999.d0
                 goto 654
              end if
              
              VC=0.d0
              do j=1,maxmes
                 do k=1,j
                    VC(j,k)=Vi(k+j*(j-1)/2)
                 end do
              end do
           end if
        end if

        mu=0.d0
        mu=matmul(X00,b0)+matmul(X2,b2)
                
!        write(*,*)'mu'

        if (idlink.eq.0) then  ! Linear link
           ymarg=9999.d0
           goto 654
           ! ne devrait pas venir ici

        else if (idlink.eq.1) then  ! Beta link

        aa1=exp(b1(nef+nvc+nwg+1))/ &
             (1+exp(b1(nef+nvc+nwg+1)))
        bb1=exp(b1(nef+nvc+nwg+2))/ &
             (1+exp(b1(nef+nvc+nwg+2)))
        bb1=aa1*(1.d0-aa1)*bb1

        cc1=abs(b1(nef+nvc+nwg+3))

        dd1=abs(b1(npm-ncor))

        aa=aa1*aa1*(1-aa1)/bb1-aa1
        bb=aa*(1-aa1)/aa1
        beta=beta_ln(aa,bb)

!                 write(*,*)'b1 beta',b1((nef+nvc+nwg+1):(nef+nvc+nwg+4)),b1(npm)

           if (methInteg.eq.1) then
              ! i.e. methode de MC a faire
              do l=1,nsim
                 usim=0.d0
                 ysim=0.d0
                 do m=1,maxmes
                    SX=1.d0
                    call bgos(SX,0,usim(m),x22,0.d0)
                 end do

                 ysim=mu+MATMUL(VC,usim)

                 do j=1,maxmes
                    ytemp=ysim(j)*dd1+cc1
                    if (ytemp.lt.0) then
                       ytemp=0.d0
                    end if
                    if (ytemp.gt.1) then
                       ytemp=1.d0
                    end if
                    ier=0
                    ymarg(maxmes*(g-1)+j)=ymarg(maxmes*(g-1)+j)+xinbta(aa,bb,beta,ytemp,ier)/dble(nsim)
!                   ymarg(maxmes*(g-1)+j)=ymarg(maxmes*(g-1)+j)+ytemp/dble(nsim) 
                    if (ier.ne.0.or.ymarg(maxmes*(g-1)+j).eq.9999.d0) then
                       ymarg(maxmes*(g-1)+j)=9999.d0
                    end if
                 end do
              end do


           else

              call gausshermite(gauss,nsim)
              do j=1,maxmes
                 do l=1,nsim
                    ytemp=mu(j)+sqrt(VC(j,j))*gauss(1,l)
                    ytemp=ytemp*dd1+cc1
                    if (ytemp.lt.0) then
                       ytemp=0.d0
                    end if
                    if (ytemp.gt.1) then
                       ytemp=1.d0
                    end if
                    ier=0
                    ymarg(maxmes*(g-1)+j)=ymarg(maxmes*(g-1)+j)+xinbta(aa,bb,beta,ytemp,ier)*gauss(2,l)
!                    ymarg(maxmes*(g-1)+j)=ymarg(maxmes*(g-1)+j)+ytemp*gauss(2,l)
                 if (ier.ne.0.or.ymarg(maxmes*(g-1)+j).eq.9999.d0) then
                    ymarg(maxmes*(g-1)+j)=9999.d0
                 end if
                 end do
              end do

           end if

        else if (idlink.eq.2) then ! Splines link


           bb=b1(nef+nvc+nwg+1)
           do kk=2,ntrtot
              splaa(kk-3)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
           end do

!write(*,*)'bb',bb,b1((nef+nvc+nwg+2):(nef+nvc+nwg+ntrtot))
                   
                   
           if (methInteg.eq.1) then

              ! i.e. methode de MC a faire
              do l=1,nsim
                 usim=0.d0
                 ysim=0.d0
                 do m=1,maxmes
                    SX=1.d0
                    call bgos(SX,0,usim(m),x22,0.d0)
                 end do

                 ysim=mu+MATMUL(VC,usim)
                 do j=1,maxmes
                    niter=0
                    diff=0.d0
                    ier=0
                    ytemp=INV_ISPLINES(ysim(j),splaa,bb,nbzitr,zitr,ier,niter,diff)
                    if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3).or.ymarg(maxmes*(g-1)+j).eq.9999.d0) then
                       ymarg(maxmes*(g-1)+j)=9999.d0
                    else
                        ymarg(maxmes*(g-1)+j)=ymarg(maxmes*(g-1)+j)+ytemp/dble(nsim)
                    end if
                 end do
              end do


           else

              call gausshermite(gauss,nsim)


              do j=1,maxmes
                 do l=1,nsim
                    niter=0
                    diff=0.d0
                    ier=0
                    ytemp=mu(j)+sqrt(VC(j,j))*gauss(1,l)
                     !   write(*,*)'avant',ll,j,mu(ll+j),sqrt(VC(ll+j,ll+j)),gauss(1,l),VC(ll+j,ll+j)

                    ytemp=INV_ISPLINES(ytemp,splaa,bb,nbzitr,zitr,ier,niter,diff)
                    if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3).or.ymarg(maxmes*(g-1)+j).eq.9999.d0) then
                       ymarg(maxmes*(g-1)+j)=9999.d0
                    else
                        ymarg(maxmes*(g-1)+j)=ymarg(maxmes*(g-1)+j)+ytemp*gauss(2,l)
                    end if

                 end do
              end do
           end if

        end if

     end do



  end if

  if (idlink.eq.1) then
     do j=1,maxmes*ng
        if (ymarg(j).ne.9999.d0) then
        ymarg(j)=ymarg(j)*(maxY-minY+2*epsY)+minY-epsY
        if (ymarg(j).lt.minY) ymarg(j)=minY
        if (ymarg(j).gt.maxY) ymarg(j)=maxY
        end if
     end do
  elseif (idlink.eq.2) then
     do j=1,maxmes*ng
        if (ymarg(j).ne.9999.d0) then
        ymarg(j)=ymarg(j)
        end if
     end do
  end if



654 continue

  deallocate(zitr,splaa)
  deallocate(Z,P,R,X00,X2,ysim,usim,mu,tcor,VC,Vi)

  !! call PutRNGstate
  call putrand()


  return

end subroutine predictcont
















