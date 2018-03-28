






!====================================================================

! PROGRAMME PRINCIPAL

!-----------------------------------------------------------
!           PREDICTION by PROFILE of COVARIATES
!------------------------------------------------------------



!===================================================================


subroutine predictcont(X0,idprob,idea,idg,idcor &
     ,ng,ncor,nv,maxmes,idiag,nwg,npm,b1,epsY,idlink &
     ,nbzitr,zitr0,nsim,methInteg,Ydiscret,Ymarg)


  use optim
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
  nprob=ng-1
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



  return

end subroutine predictcont


















!===========================================================
!
      !SUBROUTINES POUR BETA
!
!===========================================================





!C********************************************************************
!C            calcul d'une inverse de Beta incomplète
!C********************************************************************




      double precision function xinbta(p,q,beta,alpha,ifault)



      implicit double precision (a-h,o-z)


      double precision :: beta,alpha
      integer :: ifault
      double precision ::iex


!c     algorithm as 109 appl. statist. (1977), vol.26, no.1
!c     (replacing algorithm as 64  appl. statist. (1973),
!c     vol.22, no.3)
!c
!c     Remark AS R83 and the correction in vol40(1) p.236 have been
!c     incorporated in this version.
!c
!c     Computes inverse of the incomplete beta function
!c     ratio for given positive values of the arguments
!c     p and q, alpha between zero and one.
!c     log of complete beta function, beta, is assumed to be known.
!c
!c     Auxiliary function required: BETAIN = algorithm AS63
!c
      logical indx
!c
!c     Define accuracy and initialise.
!c     SAE below is the most negative decimal exponent which does not
!c     cause an underflow; a value of -308 or thereabouts will often be
!c     OK in double precision.
!c  variable SAE in XINBTA changed from -37D.0 to -308D.0 to avoid
!c  infinite loop (only a problem on digital unix).

!c
!C      data acu/1.0d-14/
!C      data SAE/-308.D0/
!C      data zero/0.0d0/, one/1.0d0/, two/2.0d0/
!C      data three/3.0d0/, four/4.0d0/, five/5.0d0/, six/6.0d0/



      double precision ::SAE=-308.D0,zero=0.0d0,one=1.0d0,two=2.0d0     &
   ,three=3.0d0,four=4.0d0,five=5.0d0,six=6.0d0

      double precision ::a,pp,qq,p,q,y,r,t,s,h,w,yprev,sq,prev,ACU
      double precision ::betain,xin,g,adj,tx,fpu

      fpu = 10.d0 ** sae
      xinbta = alpha
!c
!c     test for admissibility of parameters
!c
      ifault = 1
      if (p.le.zero .or. q.le.zero) return
      ifault = 2
      if (alpha.lt.zero .or. alpha.gt.one) return
      ifault = 0
      if (alpha.eq.zero .or. alpha.eq.one) return
!c
!c     change tail if necessary
!c

      if (alpha.le.0.5d0) goto 1
      a = one-alpha
      pp = q
      qq = p
      indx = .true.
      goto 2
    1 a = alpha
      pp = p
      qq = q
      indx = .false.
!c
!c     calculate the initial approximation
!c
    2 r = dsqrt(-dlog(a*a))
      y = r-(2.30753d0+0.27061d0*r)/(one+(0.99229d0+0.04481d0*r)*r)
      if(pp.gt.one .and. qq.gt.one) goto 5
      r = qq+qq
      t = one/(9.0d0*qq)
      t = r*(one-t+y*dsqrt(t))**3
      if(t.le.zero) goto 3
      t = (four*pp+r-two)/t
      if(t.le.one) goto 4
      xinbta = one-two/(t+one)
      goto 6
    3 xinbta = one-dexp((dlog((one-a)*qq)+beta)/qq)
      goto 6
    4 xinbta = dexp((dlog(a*pp)+beta)/pp)
      goto 6
    5 r = (y*y-three)/six
      s = one/(pp+pp-one)
      t = one/(qq+qq-one)
      h = two/(s+t)
      w = y*dsqrt(h+r)/h-(t-s)*(r+five/six-two/(three*h))
      xinbta = pp/(pp+qq*dexp(w+w))

!c
!c     solve for x by a modified newton-raphson method,
!c     using the function betain
!c
    6 r = one-pp
      t = one-qq
      yprev = zero
      sq = one
      prev = one
      if(xinbta.lt.0.0001d0) xinbta = 0.0001d0
      if(xinbta.gt.0.9999d0) xinbta = 0.9999d0
      IEX = MAX(-5.D0/PP**2 - 1.D0/A**.2 - 13.D0, SAE)
      ACU = 10.D0 ** IEX

      ACU=1.0D-30


    7 y = betain(xinbta,pp,qq,beta,ifault)

      if(ifault.eq.0) goto 8
      ifault = 3
      return
    8 continue
      xin = xinbta
      y = (y-a)*exp(beta+r*log(xin)+t*log(one-xin))
      if(y*yprev.le.zero) prev = max(sq, fpu)
      g = one
    9 adj = g*y
      sq = adj*adj
      if(sq.ge.prev) goto 10
      tx = xinbta-adj
      if(tx.ge.zero .and. tx.le.one) goto 11
   10 g = g/three
      goto 9
   11 if(prev.le.acu) goto 12
      if(y*y.le.acu) goto 12
      if(tx.eq.zero .or. tx.eq.one) goto 10
      if(tx.eq.xinbta) goto 12
      xinbta = tx
      yprev = y
      goto 7
   12 if (indx) xinbta = one-xinbta
      return
      end function xinbta



!C****************************************************************



      double precision function betain(x, p, q, beta, ifault)
      implicit double precision (a-h, o-z)
!c
!c     algorithm as 63  appl. statist. (1973), vol.22, no.3
!c
!c     computes incomplete beta function ratio for arguments
!c     x between zero and one, p and q positive.
!c     log of complete beta function, beta, is assumed to be known
!c
      logical indx
      integer :: ifault,ns
      double precision :: beta
!c
!c     define accuracy and initialise
!c

      double precision ::zero=0.0d0,one=1.0d0,acu=0.1d-14

      double precision :: x,p,q,psq,cx,xx,pp,qq,term,ai,rx,temp




      betain=x

!c     test for admissibility of arguments
!c
      ifault=1
      if(p.le.zero .or. q.le.zero) return
      ifault=2
      if(x.lt.zero .or. x.gt.one) return
      ifault=0
      if(x.eq.zero .or. x.eq. one) return
!c
!c     change tail if necessary and determine s
!c
      psq=p+q
      cx=one-x
      if(p.ge.psq*x) goto 1
      xx=cx
      cx=x
      pp=q
      qq=p
      indx=.true.
      goto 2
    1 xx=x
      pp=p
      qq=q
      indx=.false.
    2 term=one
      ai=one
      betain=one
      ns=INT(qq+cx*psq)
!c
!c     user soper's reduction formulae.
!c
      rx=xx/cx
    3 temp=qq-ai
      if(ns.eq.0) rx=xx
    4 term=term*temp*rx/(pp+ai)
      betain=betain+term
      temp=abs(term)
      if(temp.le.acu .and. temp.le.acu*betain) goto 5
      ai=ai+one
      ns=ns-1
      if(ns.ge.0) goto 3
      temp=psq
      psq=psq+one
      goto 4
!cc     calculate result
!c
    5 betain=betain*exp(pp*log(xx)+(qq-one)*log(cx)-beta)/pp
      if(indx) betain=one-betain
      return
      end function betain


!C****************************************************************

      DOUBLE PRECISION function beta_ln(z,w)

      implicit none
      double precision :: z,w
      double precision :: gammln2

      beta_ln=gammln2(z)+gammln2(w)-gammln2(z+w)
      return
      end function beta_ln

      double precision Function gammln2(xx)

      implicit none
!C retourne la valeur ln(gamma(xx)) pour xx>0
      double precision ::xx
      integer j
      double precision ser,stp,tmp,x,y,cof(6)
      save cof,stp

      data cof,stp/76.18009172947146d0,-86.50532032941677d0,     &
      24.01409824083091d0, -1.231739572450155d0,.1208650973866179d-2,     &
      -.5395239384953d-5,2.5066282746310005d0/

      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
         y=y+1.d0
         ser=ser+cof(j)/y
      end do
      gammln2=tmp+log(stp*ser/x)

      return
      end function gammln2


!==================================================================
!
      !  SUBROUTINES INVERSION ISPLINES
!
!==================================================================




!C ------------------------------------------------------------
!C
!C     INVERSE de I-splines
!C
!C     BY USING NEWTON-RAPHSON METHOD
!C ------------------------------------------------------------


!C ------------------------------------------------------------
!C
!C     INVERSE de I-splines
!C
!C     BY USING NEWTON-RAPHSON METHOD
!C ------------------------------------------------------------

      double precision FUNCTION INV_ISPLINES(X00,splaa,bb,nztr,zi_eval,istop,iter,eps)


      implicit none
      integer::iter,istop,nztr
      double precision,dimension(-1:nztr+2)::zi_eval
      double precision::X0,X00,X1,fx0,f1x0,eps,bb,bb1
      double precision,dimension(-1:nztr-1)::splaa


!       write(*,*)'X00',X00
      eps=1.d-5
      iter=1
      X0=1.d10
      call eval_splines(X0,fx0,f1x0,splaa,bb,nztr,zi_eval)
      if (fx0.eq.1.d9.and.f1x0.eq.1.d9) then
         INV_ISPLINES=fx0
         istop=3
         goto 1234
      end if

      if (X00.ge.fx0) then
         INV_ISPLINES=zi_eval(nztr)
         istop=1
         goto 1234
      end if
      X0=-1.d10
      call eval_splines(X0,fx0,f1x0,splaa,bb,nztr,zi_eval)
      if (fx0.eq.1.d9.and.f1x0.eq.1.d9) then
         INV_ISPLINES=fx0
         istop=3
         goto 1234
      end if
      if (X00.le.fx0) then
         INV_ISPLINES=zi_eval(1)
         istop=1
         goto 1234
      end if
      bb1=bb-X00
      X0=0
      call eval_splines(X0,fx0,f1x0,splaa,bb1,nztr,zi_eval)
      if (fx0.eq.1.d9.and.f1x0.eq.1.d9) then
         INV_ISPLINES=fx0
         istop=3
         goto 1234
      end if
      X1=X0-fx0/f1x0
      do while (ABS((X1-X0)/X0).GT.EPS.and.iter.lt.500)
         iter=iter+1
         X0=X1
         call eval_splines(X0,fx0,f1x0,splaa,bb1,nztr,zi_eval)
        if (fx0.eq.1.d9.and.f1x0.eq.1.d9) then
            INV_ISPLINES=fx0
            istop=3
            goto 1234
        end if
        X1=X0-fx0/f1x0
      end do
      INV_ISPLINES=zi_eval(1)+(zi_eval(nztr)-zi_eval(1))*exp(X1)/(1.d0+exp(X1))

      if (ABS((X1-X0)/X0).le.EPS) then
         istop=1
      else if (iter.ge.500) then
         istop=2
      else
         istop=3
      end if

      eps=ABS((X1-X0)/X0)

1234  continue

      return
      end function INV_ISPLINES



!C ------------------------------------------------------------
!C
!C     EVALUATION OF I-SPLINES and M-splines
!C ------------------------------------------------------------




      SUBROUTINE eval_splines(X00,Ispl,Mspl,splaa,bb,nztr,zi_eval)

      implicit none
      integer ::k,l,i,nztr
      double precision::X00,X0,Ispl,Mspl,som,bb
      double precision,dimension(-1:nztr+2)::zi_eval
      double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n &
      ,hn,hht,mmeval,mm1eval,mm2eval,imeval,im1eval,im2eval

      double precision,dimension(-1:nztr-1)::splaa


!C ou se trouve la valeur de X


      mmeval=0.d0
      mm1eval=0.d0
      mm2eval=0.d0
      imeval=0.d0
      im1eval=0.d0
      im2eval=0.d0

      X0=zi_eval(1)+(zi_eval(nztr)-zi_eval(1))*(1.d0-1.d0/(1.d0+exp(X00)))
      l=0
      do k = 2,nztr
         if ((X0.ge.zi_eval(k-1)).and.(X0.lt.zi_eval(k))) then
            l=k-1
         endif
      end do

      !if (X0.eq.zi_eval(nztr)) then
      if (abs(X0-zi_eval(nztr)).lt.1.d-10) then
         l=nztr-1
      end if

      ht2 = zi_eval(l+1)-X0
      htm= X0-zi_eval(l-1)
      ht = X0-zi_eval(l)
      ht3 = zi_eval(l+2)-X0
      hht = X0-zi_eval(l-2)
      h = zi_eval(l+1)-zi_eval(l)
      hh= zi_eval(l+1)-zi_eval(l-1)
      hn= zi_eval(l+1)-zi_eval(l-2)
      h2n=zi_eval(l+2)-zi_eval(l-1)
      h2= zi_eval(l+2)-zi_eval(l)
      h3= zi_eval(l+3)-zi_eval(l)


      if (h.eq.0.or.hh.eq.0.or.hn.eq.0.or.h2n.eq.0.or.h2.eq.0.or.h3.eq.0)  then
         Mspl=1.d9
         Ispl=1.d9
         go to 587
      end if

      !if (X0.ne.zi_eval(nztr)) then 
      if (abs(X0-zi_eval(nztr)).ge.1.d-10) then
         mm2eval = (3.d0*ht2*ht2)/(hh*h*hn)
         mm1eval = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
         mmeval  = (3.d0*ht*ht)/(h3*h2*h)
      end if
      !if (X0.eq.zi_eval(nztr)) then
      if (abs(X0-zi_eval(nztr)).lt.1.d-10) then


         mm2eval = 0.d0
         mm1eval = 0.d0
         mmeval  = 3.d0/h

      end if

      if (mm2eval.lt.0.or.mm1eval.lt.0.or.mmeval.lt.0) then
         Mspl=1.d9
         Ispl=1.d9
         go to 587
      end if

      im2eval=hht*mm2eval/(3.d0)+ h2n*mm1eval/(3.d0)+h3*mmeval/(3.d0)
      im1eval=htm*mm1eval/(3.d0)+h3*mmeval/(3.d0)
      imeval=ht*mmeval/(3.d0)

      som=0.d0
      if (l.gt.1) then
         do i=2,l
            som=som+splaa(i-3)
         end do
      end if

      Ispl=bb+ som +splaa(l-2)*im2eval+splaa(l-1)*im1eval+splaa(l)*imeval

      Mspl= (splaa(l-2)*mm2eval+splaa(l-1)*mm1eval+splaa(l)*mmeval)*      &
            (1.d0-1.d0/((1.d0+exp(X00))**2))*(zi_eval(nztr)-zi_eval(1))



 587   continue


      end subroutine eval_splines





!=============================================================
!
      !SUBROUTINES simulation
!
!=============================================================




!C ******************** BGOS ********************************


      SUBROUTINE BGOS(SX,ID,X1,X2,RO)


!C     ID=1:U(0,SX); ID DIFF DE 1 :N(0,SX)

      implicit none
      double precision ::RO,SX
      integer ::ID
      double precision ::F,V1,V2,S,DLS,RO2
      double precision ::X1,X2,UNIRAN
!C     write(*,*)'dans bgos'


 5    CONTINUE

!C     write(*,*)'avant rand :'

!C     X1=RAND()
!C     X2=RAND()

      X1=UNIRAN()
      X2=UNIRAN()

      IF(ID.NE.1) GO TO 10
      F=2.*SQRT(3.)
      X1=(X1-0.5)*F
      X2=(X2-0.5)*F
      GO TO 20
10    CONTINUE
      V1=2.*X1-1
      V2=2.*X2-1
      S=V1*V1+V2*V2
      IF(S.GE.1.) GO TO 5
      DLS=SQRT(-2.*LOG(S)/S)
      X1=V1*DLS
      X2=V2*DLS
20    CONTINUE
      RO2=RO*RO
      IF(ABS(RO).GT.1.E-10) X2=(X1+X2*SQRT(1./RO2-1.))*RO
      X1=X1*SX
      X2=X2*SX

!C      write(*,*) 'X1 ',X1,' X2 ',X2
!C OK, X1 et X2 sont créés

!C      write(*,*)'fin bgos'

      RETURN
      END subroutine bgos


!C ------------------- FIN SUBROUTINE BGOS -----------------


!C ------------------------------------------------------

      DOUBLE PRECISION FUNCTION UNIRAN()
!C
!C     Random number generator(RCARRY), adapted from F. James
!C     "A Review of Random Number Generators"
!C      Comp. Phys. Comm. 60(1990), pp. 329-344.
!C
      implicit none
      DOUBLE PRECISION SEEDS(24), TWOM24, CARRY, ONE
      PARAMETER ( ONE = 1, TWOM24 = ONE/16777216 )
      INTEGER I, J
      SAVE I, J, CARRY, SEEDS
      DATA I, J, CARRY / 24, 10, 0.0 /
      DATA SEEDS /      &
    0.8804418, 0.2694365, 0.0367681, 0.4068699, 0.4554052, 0.2880635,      &
    0.1463408, 0.2390333, 0.6407298, 0.1755283, 0.7132940, 0.4913043, &
    0.2979918, 0.1396858, 0.3589528, 0.5254809, 0.9857749, 0.4612127, &
    0.2196441, 0.7848351, 0.4096100, 0.9807353, 0.2689915, 0.5140357/
      UNIRAN = SEEDS(I) - SEEDS(J) - CARRY
      IF ( UNIRAN .LT. 0 ) THEN
         UNIRAN = UNIRAN + 1
         CARRY = TWOM24
      ELSE
         CARRY = 0
      ENDIF
      SEEDS(I) = UNIRAN
      I = 24 - MOD( 25-I, 24 )
      J = 24 - MOD( 25-J, 24 )

      return
      END function uniran










! ===================== GAUSSHERMITE NODES =======================


      subroutine gausshermite(gauss,npg)

!
!C
!C     Gauss-Hermite pour f(x)*exp(-x*x/2)/(rac(2*pi))
!C       = somme (f(xg)w(xg))
!C
!C
!C
      implicit none

      double precision, dimension(2,51)::Gauss
      double precision, dimension(51,51)::Wg,Tg
      integer :: npg,i


      DATA ( Wg(I, 5), Tg(I, 5), I = 1, 3) / &
       0.1125741132772071D-01, 0.2856970013872805D+01, &
       0.2220759220056126D+00, 0.1355626179974265D+01, &
       0.5333333333333342D+00, 0.9386691848789097D-16/
      DATA ( Wg(I, 7), Tg(I, 7), I = 1, 4) / &
       0.5482688559722184D-03, 0.3750439717725742D+01, &
       0.3075712396758645D-01, 0.2366759410734542D+01, &
       0.2401231786050126D+00, 0.1154405394739968D+01, &
       0.4571428571428575D+00, 0.2669848554723344D-16/
      DATA ( Wg(I, 9), Tg(I, 9), I = 1, 5) / &
      0.2234584400774664D-04, 0.4512745863399781D+01, &
       0.2789141321231769D-02, 0.3205429002856470D+01, &
       0.4991640676521780D-01, 0.2076847978677829D+01, &
       0.2440975028949394D+00, 0.1023255663789133D+01, &
       0.4063492063492066D+00, 0.0000000000000000D+00/
     DATA ( Wg(I,15), Tg(I,15), I = 1, 8) / &
       0.8589649899633300D-09, 0.6363947888829836D+01, &
       0.5975419597920602D-06, 0.5190093591304780D+01, &
       0.5642146405189029D-04, 0.4196207711269018D+01, &
       0.1567357503549958D-02, 0.3289082424398766D+01, &
       0.1736577449213763D-01, 0.2432436827009758D+01, &
       0.8941779539984458D-01, 0.1606710069028730D+01, &
       0.2324622936097322D+00, 0.7991290683245483D+00, &
       0.3182595182595181D+00, 0.0000000000000000D+00/
      DATA ( Wg(I,20), Tg(I,20), I = 1,10) / &
       0.1257800672437914D-12, 0.7619048541679760D+01, &
       0.2482062362315163D-09, 0.6510590157013660D+01, &
       0.6127490259983006D-07, 0.5578738805893195D+01, &
       0.4402121090230841D-05, 0.4734581334046057D+01, &
       0.1288262799619300D-03, 0.3943967350657311D+01, &
       0.1830103131080496D-02, 0.3189014816553389D+01, &
       0.1399783744710099D-01, 0.2458663611172367D+01, &
       0.6150637206397690D-01, 0.1745247320814126D+01, &
       0.1617393339840001D+00, 0.1042945348802752D+01, &
       0.2607930634495551D+00, 0.3469641570813557D+00/
      DATA ( Wg(I,30), Tg(I,30), I = 1,15) / &
      0.1640807008117853D-20, 0.9706235997359524D+01, &
       0.1585560944966296D-16, 0.8680837722732207D+01, &
       0.1624080129972436D-13, 0.7825051744352813D+01, &
       0.4573425871326147D-11, 0.7055396866960296D+01, &
       0.5178459467189710D-09, 0.6339997686869597D+01, &
       0.2882175154047618D-07, 0.5662381850082873D+01, &
       0.8909088868621158D-06, 0.5012600596486518D+01, &
       0.1657998163067346D-04, 0.4384020365898051D+01, &
       0.1965129439848249D-03, 0.3771894423159236D+01, &
       0.1544707339866097D-02, 0.3172634639420402D+01, &
       0.8295747557723240D-02, 0.2583402100229274D+01, &
       0.3111177018350134D-01, 0.2001858612956431D+01, &
       0.8278683671562172D-01, 0.1426005658374115D+01, &
       0.1580469532090208D+00, 0.8540733517109733D+00, &
       0.2179999718155776D+00, 0.2844387607362094D+00/
      DATA ( Wg(I,40), Tg(I,40), I = 1,20) / &
       0.1461839873869467D-28, 0.1145337784154873D+02, &
       0.4820467940200524D-24, 0.1048156053467427D+02, &
       0.1448609431551587D-20, 0.9673556366934033D+01, &
       0.1122275206827074D-17, 0.8949504543855559D+01, &
       0.3389853443248306D-15, 0.8278940623659475D+01, &
       0.4968088529197761D-13, 0.7646163764541459D+01, &
       0.4037638581695192D-11, 0.7041738406453829D+01, &
       0.1989118526027766D-09, 0.6459423377583766D+01, &
       0.6325897188548972D-08, 0.5894805675372016D+01, &
       0.1360342421574886D-06, 0.5344605445720084D+01, &
       0.2048897436081474D-05, 0.4806287192093873D+01, &
       0.2221177143247582D-04, 0.4277826156362752D+01, &
       0.1770729287992397D-03, 0.3757559776168985D+01, &
       0.1055879016901825D-02, 0.3244088732999869D+01, &
       0.4773544881823334D-02, 0.2736208340465433D+01, &
       0.1653784414256937D-01, 0.2232859218634873D+01, &
       0.4427455520227679D-01, 0.1733090590631720D+01, &
       0.9217657917006089D-01, 0.1236032004799159D+01, &
       0.1499211117635710D+00, 0.7408707252859313D+00, &
       0.1910590096619904D+00, 0.2468328960227240D+00/
      DATA ( Wg(I,50), Tg(I,50), I = 1,25) / &
       0.1034607500576990D-36, 0.1298588445541555D+02, &
       0.9443414659584510D-32, 0.1205301838092448D+02, &
       0.6856280758924735D-28, 0.1127923332148262D+02, &
       0.1206044550761014D-24, 0.1058738174919177D+02, &
       0.7995094477915292D-22, 0.9948035709637500D+01, &
       0.2522482807168144D-19, 0.9346039593575728D+01, &
       0.4368171816201588D-17, 0.8772299579514598D+01, &
       0.4566698246800344D-15, 0.8220815907982127D+01, &
       0.3083828687005300D-13, 0.7687362406712500D+01, &
       0.1414228936126661D-11, 0.7168814837853899D+01, &
       0.4576636712310442D-10, 0.6662775399018720D+01, &
       0.1077060789389039D-08, 0.6167347388659921D+01, &
       0.1888225976835208D-07, 0.5680992291033284D+01, &
       0.2514609880838772D-06, 0.5202434993399912D+01, &
       0.2584937658949391D-05, 0.4730598550228594D+01, &
       0.2078485175734569D-04, 0.4264557843038109D+01, &
       0.1321726328668984D-03, 0.3803505741742012D+01, &
       0.6708280619787080D-03, 0.3346727774732429D+01, &
       0.2738160896935348D-02, 0.2893582727707738D+01, &
       0.9045054154849623D-02, 0.2443487452654017D+01, &
       0.2430481286424306D-01, 0.1995904709795124D+01, &
       0.5334352453170102D-01, 0.1550333214338771D+01, &
       0.9593054035810168D-01, 0.1106299289397183D+01, &
       0.1416854132499443D+00, 0.6633496795082918D+00, &
       0.1721258519924433D+00, 0.2210451816445435D+00/



!      if(npg.ne.5.and.npg.ne.7.and.npg.ne.9.and.npg.ne.15.   &
!          and.npg.ne.20.and.npg.ne.30.and.npg.ne.40.and.npg.ne.50) then
!         write(*,*)'nb pts GH = 5,7,9,15,20,30,40, ou 50'
!         stop
!      end if


!ccccccccccccccccccccccccccccccccccccccccccccccc
      DO I = 1, NPG/2
         GAUSS(1,I) = -Tg(I,NPG)
         GAUSS(2,I) =  Wg(I,NPG)
         GAUSS(1,NPG-I+1) = Tg(I,NPG)
         GAUSS(2,NPG-I+1) = Wg(I,NPG)
      END DO
      IF ( MOD( NPG, 2 ) .EQ. 1 ) THEN
         GAUSS(1, NPG/2 + 1 ) = 0.D0
         GAUSS(2, NPG/2 + 1 ) = Wg( NPG/2 + 1, NPG )
      END IF

!ccccccccccccccccccccccccccccccccccccccccccccccc

      end subroutine gausshermite

