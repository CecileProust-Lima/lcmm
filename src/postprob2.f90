


!-------------------------------------------------------------
!
!          Computation of the individual posterior probabilty of 
!          surviving without the event in each component 
!          for longitudinal data truncated in s
!-------------------------------------------------------------
! -> if no measures at times < s , then ppi=-5




subroutine postprob2(Y0,Xdata0,ns0,nmes0,nobs0,ng0,nv0,idiag0,nwg0, &
     ncor0,logspecif0,zi0,idea0,idg0,idprob0,idcor0,idcom0,idspecif0,idtdv0, &
     risqcom0,nvdepsurv0,typrisq0,nz0,nprisq0,nbzitr0, &
     zitr0,idlink0,indiceY0,uniqueY0,nvalSPL0,epsY0,nbevt0, &
     tsurv00,idtrunc0,b,npm0,timemes,landmark,nbland,ppi)


  use commun_comp
  use optim

  implicit none

  !Declaration des variables en entree
  integer,intent(in)::nobs0,ns0,nv0,ng0,npm0,nbland,nbevt0,nbzitr0
  integer,intent(in)::nvalSPL0,idlink0,epsY0
  double precision,dimension(nobs0),intent(in):: Y0
  integer,dimension(nobs0),intent(in):: indiceY0
  double precision,dimension(nobs0*nv0),intent(in)::Xdata0
  integer,dimension(ns0),intent(in)::nmes0
  double precision,dimension(ns0),intent(in)::tsurv00
  integer,intent(in)::idiag0,nwg0,ncor0,logspecif0,nvdepsurv0,idtrunc0
  integer,dimension(nbevt0),intent(in)::risqcom0,typrisq0,nprisq0,nz0
  integer,dimension(nv0),intent(in)::idcom0,idea0,idg0,idprob0,idcor0,idtdv0
  integer,dimension(nv0*nbevt0),intent(in)::idspecif0
  double precision,dimension(maxval(nz0),nbevt0),intent(in)::zi0
  double precision,dimension(npm0),intent(in)::b
  double precision,dimension(nobs0),intent(in)::timemes
  double precision,dimension(nbland),intent(in)::landmark
  double precision,dimension(nbzitr0+2),intent(in)::zitr0
  double precision,dimension(nvalSPL0),intent(in)::uniqueY0

  integer ::i,j,k,l,m,g,l2,m2,jj,it,ier,nmoins,kk,j1,j2,ll
  integer::dejaspl,jtemp,ks,ktemp,sumnmes,nxevtcurr
  double precision,dimension(maxval(nmes0),nv0) ::Z,P,X0,X2
  double precision,dimension(nv0) ::Xprob
  double precision,dimension(nv0,nv0) ::Ut,Ut1
  double precision,dimension(:,:),allocatable ::VC,Corr
  double precision,dimension(npm0) :: b1
  double precision,dimension(:),allocatable ::Vi
  double precision,dimension(nv0) :: b0,b2,bprob
  double precision :: eps,det,eta0,ytemp
  double precision,dimension(-1:nbzitr0-1)::splaa
  double precision ::temp,aa,bb,aa1,bb1,cc1,dd1,som,betai
  double precision ::Y4,fevt
  double precision,dimension(ng0) ::fi,pi,fi1
  double precision,dimension(ns0)::tsurvint0

  double precision,dimension(maxval(nmes0)) :: mu,Y1,Y2,Y3,tcor

  double precision,dimension(nv0*nbevt0)::bevt,Xevt
  double precision,dimension(nbevt0)::bevtint
  double precision::surv_glob,varexpsurv
  integer::sumnrisq,ii,ke
  double precision,dimension(maxval(nprisq0))::brisq
  double precision,dimension(ng0,nbevt0)::risq,surv,surv0,survint


  !sortie
  double precision,dimension(ns0*nbland,ng0),intent(out)::ppi

  eps=1.D-20
  ppi=-5.D0
  b1=0.d0
  do k=1,npm0
     b1(k)=b(k)
  end do


  ! alloc pour partie modele mixte
  minY=zitr0(1)
  maxY=zitr0(nbzitr0)

  rangeY=0
  !if (Ydiscret.eq.1) rangeY=maxY-minY


  epsY=epsY0
  idlink=idlink0
  nvalSPL=nvalSPL0
  if (idlink.eq.0) ntrtot=2  ! linear
  if (idlink.eq.1) ntrtot=4  ! beta
  if (idlink.eq.2) then      ! splines nbzitr0 noeuds
     ntrtot=nbzitr0+2
     allocate(zitr(-1:(ntrtot)))

     allocate(mm(nvalSPL),mm1(nvalSPL),mm2(nvalSPL))
     allocate(im(nvalSPL),im1(nvalSPL),im2(nvalSPL))

     allocate(indiceY(nobs0),uniqueY(nvalSPL0))

     zitr(1:nbzitr0)=zitr0(1:nbzitr0)
     zitr(-1)=zitr(1)
     zitr(0)=zitr(1)
     zitr(ntrtot-1)=zitr(ntrtot-2)
     zitr(ntrtot)=zitr(ntrtot-1)

     indiceY(1:nobs0)=indiceY0(1:nobs0)
     uniqueY(1:nvalSPL)=uniqueY0(1:nvalSPL0)
  else
     allocate(zitr(1))
     allocate(mm(1),mm1(1),mm2(1),im(1),im1(1),im2(1))
     allocate(indiceY(1),uniqueY(1))

     zitr(1)=0.d0
     mm(1)=0.d0
     mm1(1)=0.d0
     mm2(1)=0.d0
     im(1)=0.d0
     im1(1)=0.d0
     im2(1)=0.d0
     indiceY(1)=0
     uniqueY(1)=0.d0
  end if

  !print*,"ntrtot",ntrtot
  allocate(Y(nobs0),idprob(nv0),X(nobs0,nv0) &
       ,idea(nv0),idg(nv0),idcor(nv0),nmes(ns0) &
       ,idcom(nv0),idspecif(nv0*nbevt0),idtdv(nv0))



  !  alloc pour partie survie

  allocate(Tsurv0(ns0),Tsurv(ns0),Tsurvint(ns0),ind_survint(ns0),Devt(ns0))
  allocate(risqcom(nbevt0),typrisq(nbevt0),nz(nbevt0),nprisq(nbevt0),nrisq(nbevt0))
  allocate(nevtparx(nv0),nxevtspec(nbevt0),nxcurr(nv0))

  nbevt=nbevt0    
  typrisq=typrisq0
  risqcom=risqcom0
  idtrunc=idtrunc0 
  Tsurv0=tsurv00      
  devt=0    
  ind_survint=0
  logspecif=logspecif0  
  nvdepsurv=nvdepsurv0


  ! zi : contient noeuds pour hazard (ou min et max si Weibull)
  if(any(typrisq.eq.3)) then
     allocate(zi(-2:maxval(nz0)+3,nbevt))
  else
     allocate(zi(maxval(nz0),nbevt))
  end if

  ! recopier le reste des parametres dans module commun_comp
  ns=ns0
  ng=ng0
  ncor=ncor0
  nv=nv0
  nobs=nobs0
  if (nwg0.eq.0) then
     nwg=0
  else
     nwg=ng-1
  end if
  idiag=idiag0


  Y=0.d0
  X=0.d0
  idprob=0 !classmb
  idea=0   !random
  idg=0    !fixed
  idcor=0  !cor
  idcom=0 !survival
  idspecif(1:(nv*nbevt))=idspecif0(1:(nv*nbevt))
  tsurvint0=0.d0

  ktemp=0
  do k=1,nv
     idprob(k)=idprob0(k)
     idea(k)=idea0(k)
     idg(k)=idg0(k)
     idcor(k)=idcor0(k)
     idcom(k)=idcom0(k)
     idtdv(k)=idtdv0(k)

     jtemp=0
     it=0
     DO i=1,ns
        if (k.eq.1) then
           do j=1,nmes0(i)
              jtemp=jtemp+1
              Y(jtemp)=Y0(jtemp)
           end do
        end if

        do j=1,nmes0(i)
           ktemp=ktemp+1
           it=it+1
           X(it,k)=Xdata0(ktemp)
        end do

        if(idtdv(k).eq.1) then
           tsurvint0(i)=X(it,k)
        end if

     end do
  end do

  nprisq=0
  nrisq=0
  nrisqtot=0

  do ke=1,nbevt

     nz(ke)=nz0(ke) ! nb de noeuds pour hazard (ou 2 pr Weibull)

     if (typrisq(ke).eq.1) then
        nprisq(ke)=nz(ke)-1
     end if
     if (typrisq(ke).eq.2) then
        nprisq(ke)=2
     end if
     if (typrisq(ke).eq.3) then
        nprisq(ke)=nz(ke)+2
     end if

     if (risqcom(ke).eq.1) then
        nrisq(ke) = nprisq(ke)
     end if
     if (risqcom(ke).eq.2) then
        nrisq(ke) = nprisq(ke)+ng-1
     end if
     if (risqcom(ke).eq.0) then
        nrisq(ke) = nprisq(ke)*ng
     end if

     nrisqtot = nrisqtot+nrisq(ke)  ! nb total de prm pour hazards
     zi(1:nz(ke),ke)=zi0(1:nz(ke),ke)
  end do


  ! nvarxevt etc
  nxevt=0
  nxevtspec=0
  nevtparx=0
  do j=1,nv
     if (idtdv(j).eq.1) then 
        if(idcom(j).eq.1.and.idspecif(j).eq.1) then
           nevtparx(j)=1
        end if
        if(idcom(j).eq.1.and.idspecif(j).eq.2) then 
           nevtparx(j)=ng
        end if
        if(idcom(j).eq.0) then
           do ke=1,nbevt
              if(idspecif((ke-1)*nv+j).eq.1) then
                 nevtparx(j)=nevtparx(j)+1
              end if
              if(idspecif((ke-1)*nv+j).eq.2) then
                 nevtparx(j)=nevtparx(j)+ng
              end if
           end do
        end if

     else


        if(idcom(j).eq.1.and.idspecif(j).eq.1) then
           nevtparx(j)=1
           do ke=1,nbevt
              nxevtspec(ke)=nxevtspec(ke)+1
           end do
        end if
        if(idcom(j).eq.1.and.idspecif(j).eq.2) then
           nevtparx(j)=ng
           do ke=1,nbevt
              nxevtspec(ke)=nxevtspec(ke)+1
           end do
        end if
        if(idcom(j).eq.0) then
           do ke=1,nbevt
              if(idspecif((ke-1)*nv+j).eq.1) then
                 nevtparx(j)=nevtparx(j)+1
                 nxevtspec(ke)=nxevtspec(ke)+1
              end if
              if(idspecif((ke-1)*nv+j).eq.2) then
                 nevtparx(j)=nevtparx(j)+ng
                 nxevtspec(ke)=nxevtspec(ke)+1
              end if
           end do
        end if
     end if
     nvarxevt = sum(nevtparx)
     nxevt=sum(nxevtspec)

  end do


  nea=0
  ncg=0
  ncssg=0
  nprob=0 
  nvarprob=0
  do k=1,nv
     if (idg(k).eq.1) then
        ncssg=ncssg+1      ! nb var. sans melange
     else if (idg(k).eq.2) then
        ncg=ncg+1      ! nb var. dans melange
     end if
     nea=nea+idea(k)
     nprob=nprob+(idprob(k))*(ng-1)
     nvarprob=nvarprob+idprob(k)
  end do

  if (idiag.eq.1) then
     nvc=nea
  else if(idiag.eq.0) then
     nvc=nea*(nea+1)/2
  end if


  nef=ncssg+ncg*ng
  if(idlink.ne.-1) nef=nef-1


  ! creer base de splines si au moins un hazard splines
  if(any(typrisq.eq.3)) then
     allocate(Tmm(ns),Tmm1(ns),Tmm2(ns),Tmm3(ns),Tim(ns)         &
          ,Tim1(ns),Tim2(ns),Tim3(ns),Tmm0(ns),Tmm01(ns),Tmm02(ns)  &
          ,Tmm03(ns),Tim0(ns),Tim01(ns),Tim02(ns),Tim03(ns),          &
          Tmmt(ns),Tmmt1(ns),Tmmt2(ns),Tmmt3(ns),Timt(ns),Timt1(ns) &
          ,Timt2(ns),Timt3(ns))
  end if

  ! creer base de splines si transfo splines
  if (idlink.eq.2) then 
     call design_splines_comp(ier)
     if (ier.eq.-1) then
        go to 123
     end if
  end if


  ! nmbre max de mesures pour 1 sujet
  maxmes=maxval(nmes0)

  ! modules ok (sauf tsurv,nmes,Tmm,tsurvint,ind_survint  a faire pour chaque landmark)



  allocate(VC(maxmes,maxmes),Corr(maxmes,maxmes),Vi(maxmes*(maxmes+1)/2))


  do ks=1,nbland

     VC=0.d0
     Vi=0.d0
     Corr=0.d0


     ! print*,"ks=",ks
     tsurv(1:ns0)=landmark(ks)
     tsurvint(1:ns0)=landmark(ks)
     ind_survint(1:ns0)=0

     nmes=0
     sumnmes=0
     do i=1,ns
        do j=1,nmes0(i)
           !print*,"timemes",timemes(sumnmes+j),i,j
           if(timemes(sumnmes+j).le.landmark(ks)) then
              nmes(i) = nmes(i) + 1 
           end if
        end do
        ! print*,"nmes=",nmes

        if(tsurvint0(i)<landmark(ks)) then
           tsurvint(i) = tsurvint0(i)
           ind_survint(i) = 1
        end if
        sumnmes = sumnmes+nmes0(i)
     end do


     dejaspl=0
     if (any(typrisq.eq.3)) then
        do ke=1,nbevt
           if(typrisq(ke).eq.3 .and. dejaspl.eq.0) then
              call splines(ke)
              dejaspl=1
           end if
        end do
        !print*,"Tmm",Tmm
     end if

!!! debut postprob !!!

     ! calcul de brisq en chaque composante et risq, surv et surv0 pour chaque evt et tous les g

     Ut=0.d0
     if (nea>0) then

        If (idiag.eq.1) then
           do j=1,nea
              do k=1,nea
                 if (j.eq.k) then
                    Ut(j,k)=b1(nprob+nrisqtot+nvarxevt+nef+j)
                 else
                    Ut(j,k)=0.d0
                 end if
              end do
           end do
        end if

        If (idiag.eq.0) then
           do j=1,nea
              do k=1,j
                 Ut(j,k)=b1(nprob+nrisqtot+nvarxevt+nef+k+j*(j-1)/2)
              end do
           end do
        end if

     end if



     ! boucle sur les sujets


     it=0
     do i=1,ns

        if(nmes(i)>0) then

           ! calcul de brisq en chaque composante et risq, surv et surv0 pour chaque evt et tous les g

           risq=0.d0
           surv=0.d0
           surv0=0.d0
           survint=0.d0
           sumnrisq=0  
           do ke=1,nbevt

              do g=1,ng

                 brisq=0.d0
                 if (logspecif.eq.1) then
                    if (risqcom(ke).eq.0) then
                       do k=1,nprisq(ke)
                          brisq(k)=exp(b1(nprob+sumnrisq+nprisq(ke)*(g-1)+k))
                       end do
                    elseif(risqcom(ke).eq.1) then
                       do k=1,nprisq(ke)
                          brisq(k)=exp(b1(nprob+sumnrisq+k))
                       end do
                    elseif (risqcom(ke).eq.2) then
                       do k=1,nprisq(ke)
                          brisq(k)=exp(b1(nprob+sumnrisq+k))
                       end do
                    end if

                 else
                    if (risqcom(ke).eq.0) then
                       do k=1,nprisq(ke)
                          brisq(k)=b1(nprob+sumnrisq+nprisq(ke)*(g-1)+k) &
                               *b1(nprob+sumnrisq+nprisq(ke)*(g-1)+k)
                       end do
                    elseif(risqcom(ke).eq.1) then
                       do k=1,nprisq(ke)
                          brisq(k)=b1(nprob+sumnrisq+k)*b1(nprob+sumnrisq+k)
                       end do
                    elseif (risqcom(ke).eq.2) then
                       do k=1,nprisq(ke)
                          brisq(k)=b1(nprob+sumnrisq+k)*b1(nprob+sumnrisq+k)
                       end do
                    end if
                 end if


                 call fct_risq_comp_i(i,ke,brisq,g,risq,surv,surv0,survint)


                 if (risqcom(ke).eq.2.and.ng.gt.1.and.g.lt.ng) then
                    risq(g,ke)=risq(g,ke)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
                    surv(g,ke)=surv(g,ke)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
                    survint(g,ke)=survint(g,ke)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
                 end if

              end do ! fin boucle classe
              sumnrisq= sumnrisq+nrisq(ke)
           end do


           ! matrice de variance pour le sujet i : VC=ZGZ+Cor+sigma2I

           Z=0.d0
           l=0
           do k=1,nv
              if (idea(k).eq.1) then
                 l=l+1
                 do j=1,nmes(i)
                    Z(j,l)=dble(X(it+j,k))
                 end do
              end if
           end do

           Corr=0.d0
           tcor=0.d0
           if (ncor.gt.0) then
              do k=1,nv
                 if (idcor(k).eq.1) then
                    do j=1,nmes(i)
                       tcor(j) = X(it+j,k)
                    end do
                 end if
              end do
           end if

           do j1=1,nmes(i)
              do j2=1,nmes(i)
                 if (j1.eq.j2) then
                    if(idlink.ne.-1) Corr(j1,j2) = 1
                    if(idlink.eq.-1) Corr(j1,j2) = b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)**2
                 end if
                 if (ncor.eq.1) then 
                    Corr(j1,j2) = Corr(j1,j2)+ &
                         b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor)**2 * &
                         min(tcor(j1),tcor(j2))
                 else if (ncor.eq.2) then
                    Corr(j1,j2) = Corr(j1,j2)+ &
                         b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor)**2 * &
                         exp(-b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+1)* &
                         abs(tcor(j1)-tcor(j2)))
                 end if
              end do
           end do

           ! creation de P=Zi*Ut et V=P*P' que si non spec aux classes

           if (nwg.eq.0) then
              P=0.d0
              P=MATMUL(Z,Ut)
              VC=0.d0
              VC=MATMUL(P,transpose(P))+Corr

              ! Vi en vecteur
              jj=0
              do j=1,nmes(i)
                 do k=j,nmes(i)
                    jj=j+k*(k-1)/2
                    Vi(jj)=VC(j,k)
                 end do
              end do

              CALL DSINV(Vi,nmes(i),eps,ier,det)
              if (ier.eq.-1) then
                 ppi=-1.d0
                 go to 456
              end if

              ! retransformation du vecteur Vi en matrice :

              do j=1,nmes(i)
                 do k=1,nmes(i)
                    if (k.ge.j) then
                       VC(j,k)=Vi(j+k*(k-1)/2)
                    else
                       VC(j,k)=Vi(k+j*(j-1)/2)
                    end if
                 end do
              end do
           end if
           ! VC ok si non spec aux classes


           ! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
           Xprob=0.d0
           l=0
           do k=1,nv
              if (idprob(k).eq.1) then
                 l=l+1
                 Xprob(l)=X(it+1,k)
              end if
           end do
           pi=0.d0
           temp=0.d0
           Do g=1,ng-1
              bprob=0.d0
              do k=1,nvarprob
                 bprob(k)=b1((k-1)*(ng-1)+g)
              end do

              temp=temp+exp(DOT_PRODUCT(bprob,Xprob))

              pi(g)=exp(DOT_PRODUCT(bprob,Xprob))

           end do

           pi(ng)=1/(1+temp)

           do g=1,ng-1
              pi(g)=pi(g)*pi(ng)
           end do
           ! pig ok


           l=0
           m=0
           X0=0.d0
           X2=0.d0
           do k=1,nv
              if (idg(k).eq.2) then
                 l=l+1
                 do j=1,nmes(i)
                    X2(j,l)=dble(X(it+j,k))
                 end do
              else if (idg(k).eq.1) then
                 m=m+1
                 do j=1,nmes(i)
                    X0(j,m)=dble(X(it+j,k))
                 end do
              end if
           end do
           ! X0 et X2 ok

           ! creation des donnees transformees

           Y1=0.d0
           splaa=0.d0

           if (idlink.eq.0) then  ! Linear link

              do j=1,nmes(i)
                 Y1(j)=(dble(Y(it+j))- &
                      b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1))/ &
                      abs(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2))
              end do

           else if (idlink.eq.1) then  ! Beta link


              aa1=exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1))/ &
                   (1+exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)))
              bb1=exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2))/ &
                   (1+exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2)))
              bb1=aa1*(1.d0-aa1)*bb1

              cc1=abs(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+3))

              dd1=abs(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+4))

              aa=aa1*aa1*(1-aa1)/bb1-aa1
              bb=aa*(1-aa1)/aa1

              do j=1,nmes(i)

                 ytemp=(dble(Y(it+j))-minY+epsY)/(maxY-minY+2*epsY)
                 Y1(j)=(betai(aa,bb,ytemp)-cc1)/dd1


                 if (Y1(j).eq.999.d0) then
                    ppi=-1.d0
                    goto 456
                 end if

              end do

           else if (idlink.eq.2) then ! Splines link

              splaa=0.d0
              eta0=0.d0
              eta0=b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)

              do kk=2,ntrtot
                 splaa(kk-3)=b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+kk)**2
              end do

              do j=1,nmes(i)
                 ll=0

                 if (Y(it+j).eq.zitr(ntrtot-2)) then
                    ll=ntrtot-3
                 end if

                 som=0.d0
                 do kk = 2,ntrtot-2
                    if ((Y(it+j).ge.zitr(kk-1)).and. &
                         (Y(it+j).lt.zitr(kk))) then
                       ll=kk-1
                    end if
                 end do

                 if (ll.lt.1.or.ll.gt.ntrtot-3) then
                    ppi=-1.d0
                    goto 456
                 end if
                 if (ll.gt.1) then
                    do ii=2,ll
                       som=som+splaa(ii-3)
                    end do
                 end if


                 Y1(j)=eta0+som +splaa(ll-2)*im2(indiceY(it+j)) &
                      +splaa(ll-1)*im1(indiceY(it+j))&
                      + splaa(ll)*im(indiceY(it+j))

              end do
           else if (idlink.eq.-1) then
              do j=1,nmes(i)
                 Y1(j)=dble(Y(it+j))
              end do
           end if
           ! Y1 ok

           ! calcul des proba a posteriori de chaque classe
           !sumig=0
           fi=0.d0
           fi1=0.d0


           do g=1,ng

              ! creer tous les beta
              bevt=0.d0
              bevtint=0.d0
              Xevt=0.d0
              if (nxevt.ne.0) then
                 l=0
                 ! for each covariate, number of previous effects before ke
                 nxcurr=0
                 do ke=1,nbevt     
                    ! number of previous effects in b1 before k    
                    nxevtcurr=0
                    do k=1,nv                  
                       if (idtdv(k).ne.1) then
                          if (idcom(k).eq.1) then  
                             if (idspecif(k).eq.2) then 
                                l=l+1
                                bevt(l)=b1(nprob+nrisqtot+nxevtcurr+g)
                                Xevt(l)=X(it+1,k)
                             else 
                                l=l+1
                                bevt(l)=b1(nprob+nrisqtot+nxevtcurr+1)
                                Xevt(l)=X(it+1,k)
                             end if

                          else 

                             if (idspecif((ke-1)*nv+k).eq.1) then   
                                l=l+1
                                bevt(l)=b1(nprob+nrisqtot+nxevtcurr+nxcurr(k)+ 1)
                                Xevt(l)=X(it+1,k)
                                nxcurr(k)=nxcurr(k)+1
                             end if
                             if (idspecif((ke-1)*nv+k).eq.2) then  
                                l=l+1
                                bevt(l)=b1(nprob+nrisqtot+nxevtcurr+nxcurr(k)+g)
                                Xevt(l)=X(it+1,k)
                                nxcurr(k)=nxcurr(k)+ng
                             end if
                          end if

                       else

                          if (idcom(k).eq.1) then  
                             if(idspecif(k).eq.1) then
                                bevtint(ke)=b1(nprob+nrisqtot+nxevtcurr+1)
                             end if
                             if(idspecif(k).eq.2) then
                                bevtint(ke)=b1(nprob+nrisqtot+nxevtcurr+g)
                             end if
                          else 
                             if (idspecif((ke-1)*nv+k).eq.1) then 
                                bevtint(ke)=b1(nprob+nrisqtot+nxevtcurr+nxcurr(k)+1)
                                nxcurr(k)=nxcurr(k)+1
                             end if
                             if (idspecif((ke-1)*nv+k).eq.2) then
                                bevtint(ke)=b1(nprob+nrisqtot+nxevtcurr+nxcurr(k)+g)
                                nxcurr(k)=nxcurr(k)+ng
                             end if
                          end if
                       end if

                       nxevtcurr=nxevtcurr+nevtparx(k)
                    end do
                 end Do

              end if


              b0=0.d0
              b2=0.d0
              nmoins=0
              l2=0
              m2=0
              do k=1,nv
                 if (idg(k).eq.1) then
                    ! parametre a 0 pour l'intercept
                    if (k.eq.1.and.idlink.ne.-1) then
                       m2=m2+1
                       b0(m2)=0.d0
                    else
                       m2=m2+1
                       b0(m2)=b1(nprob+nrisqtot+nvarxevt+nmoins+1)
                       nmoins=nmoins+1
                    end if
                 else if (idg(k).eq.2) then
                    ! parametre a 0 pour l'intercept de la premiere classe
                    if (k.eq.1.and.idlink.ne.-1) then
                       if (g.eq.1) then
                          l2=l2+1
                          b2(l2)=0.d0
                          nmoins=nmoins+ng-1
                       else
                          l2=l2+1
                          b2(l2)=b1(nprob+nrisqtot+nvarxevt+nmoins+g-1)
                          nmoins=nmoins+ng-1
                       end if
                    else
                       l2=l2+1
                       b2(l2)=b1(nprob+nrisqtot+nvarxevt+nmoins+g)
                       nmoins=nmoins+ng
                    end if
                 end if
              end do



              !     variance covariance si spec aux classes :
              if (nwg.ne.0) then
                 Ut1=0.d0
                 if (g.eq.ng) then
                    Ut1=Ut
                 else
                    Ut1=Ut*abs(b1(nprob+nrisqtot+nvarxevt+nef+nvc+g))
                 end if

                 P=0.d0
                 P=MATMUL(Z,Ut1)
                 VC=0.d0
                 VC=MATMUL(P,transpose(P))+Corr

                 !     Vi en vecteur
                 jj=0
                 do j=1,nmes(i)
                    do k=j,nmes(i)
                       jj=j+k*(k-1)/2
                       Vi(jj)=VC(j,k)
                    end do
                 end do

                 CALL DSINV(Vi,nmes(i),eps,ier,det)
                 if (ier.eq.-1) then
                    ppi=-1.d0
                    goto 456
                 end if

                 !     retransformation du vecteur Vi en matrice :

                 do j=1,nmes(i)
                    do k=1,nmes(i)
                       if (k.ge.j) then
                          VC(j,k)=Vi(j+k*(k-1)/2)
                       else
                          VC(j,k)=Vi(k+j*(j-1)/2)
                       end if
                    end do
                 end do

              end if


              mu=0.d0
              mu=matmul(X0,b0)+matmul(X2,b2)

              Y2=Y1-mu
              Y3=MATMUL(VC,Y2)
              Y4=DOT_PRODUCT(Y2,Y3)

              fi(g)=fi(g) - nmes(i)*log(dble(2*3.14159265))
              fi(g)=fi(g) - det
              fi(g)=fi(g) - Y4
              fi(g)=fi(g)/(2.d0)

              surv_glob=0.d0
              nxevtcurr=0
              fevt=0.d0
              do ke=1,nbevt
                 varexpsurv=0.d0
                 if (nxevt.ne.0) then 
                    varexpsurv=DOT_PRODUCT(Xevt((nxevtcurr+1):(nxevtcurr+nxevtspec(ke)))&
                         ,bevt((nxevtcurr+1):(nxevtcurr+nxevtspec(ke))))
                 end if

                 Surv_glob=surv_glob+ survint(g,ke)*exp(varexpsurv)+ &
                      exp(bevtint(ke)+varexpsurv)*(surv(g,ke)-survint(g,ke))


                 if (Devt(i).eq.ke) then
                    fevt=log(risq(g,ke))+varexpsurv              
                    if (ind_survint(i).eq.1) then
                       fevt=fevt+ bevtint(ke)
                    end if
                 end IF
                 nxevtcurr=nxevtcurr+nxevtspec(ke)
              end do

              fi1(g)=exp(fi(g))

              fi(g)=fi(g)+fevt-surv_glob            
              fi(g)=exp(fi(g))

              !ppi(ns*(ks-1)+i,g) = pi(g)*fi1(g)

           end do ! fin boucle classe

           !sumig=DOT_PRODUCT(pi,fi)

           do g=1,ng
              ppi(ns*(ks-1)+i,g)=pi(g)*fi1(g)/DOT_PRODUCT(pi,fi)
           end do

        end if

        it = it + nmes0(i)


     end do ! fin boucle sujet


  end do ! fin boucle landmark

456 continue



  deallocate(VC,Vi,Corr)

123 continue

  if(any(typrisq.eq.3)) then

     deallocate(Tmm,Tmm1,Tmm2,Tmm3,Tim         &
          ,Tim1,Tim2,Tim3,Tmm0,Tmm01,Tmm02  &
          ,Tmm03,Tim0,Tim01,Tim02,Tim03,          &
          Tmmt,Tmmt1,Tmmt2,Tmmt3,Timt,Timt1 &
          ,Timt2,Timt3)

  end if

  deallocate(tsurv0,tsurv,tsurvint,ind_survint,Devt,risqcom,typrisq, &
       nz,nprisq,nrisq,Y,idprob,X,idea,idg,idcor,nmes,idcom,idspecif,idtdv, &
       zitr,mm,mm1,mm2,im,im1,im2,indiceY,uniqueY,zi,nevtparx,nxevtspec,nxcurr)

  return

end subroutine postprob2
