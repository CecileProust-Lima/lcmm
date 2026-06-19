subroutine postREmult(g0,Y0,X0,idprob0,idea0,idg0,idcor0,idcontr0 &
     ,ny0,ns0,ng0,nv0,nobs0,nea0,nmes0,idiag0,nwg0,ncor0,nalea0&
     ,npm0,b1,RE,ncontr0,nvc0,ntrtot0 &
     ,epsY0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0 &
     ,nvalSPLORD0,nMC0,dimMC0,seqMC0,chol0,vraispost)

  use communmo

  IMPLICIT NONE

  !Declaration des variables en entree
  integer,intent(in)::g0,nv0,ny0,chol0,ncontr0,nvc0,nMC0,dimMC0
  integer, intent(in)::ns0,ng0,nobs0,idiag0,nwg0,nea0,ncor0,nalea0,npm0
  double precision,dimension(ny0),intent(in)::epsY0
  integer, dimension(ny0),intent(in)::idlink0,nbzitr0,nvalSPLORD0,ntrtot0
  double precision,dimension(maxval(nbzitr0),ny0),intent(in)::zitr0
  integer,dimension(nobs0),intent(in)::indiceY0
  double precision,dimension(sum(nvalSPLORD0(:))),intent(in)::uniqueY0
  integer, dimension(nv0),intent(in)::idea0,idg0,idprob0,idcor0,idcontr0
  integer,dimension(ny0)::nmes0   
  double precision,dimension(nobs0),intent(in)::Y0
  double precision,dimension(nobs0*nv0),intent(in)::X0
  double precision,dimension(dimMC0*nMC0),intent(in)::seqMC0
  double precision, dimension(npm0), intent(in) :: b1

  double precision, dimension(nea0), intent(in)::RE
  double precision, intent(out)::vraispost

  !Variables locales
  integer ::i,j,k,l,m,g,l2,m2,jj,ll,ii,numSPL,ykord,nmescur
  integer ::ier,nmoins,kk,j1,j2,sumMesYk,yk,sumntrtot,sumntr
  integer::jtemp,ktemp,k1,k2
  double precision,dimension(nobs0,nv0) ::X00,X2
  double precision,dimension(nobs0,nea0) ::Z
  double precision,dimension(nobs0,(ncontr0+sum(idcontr0)))::X01
  double precision,dimension(ncontr0+sum(idcontr0))::b01
  double precision,dimension(nea0,nea0) ::Ut,varB
  double precision,dimension(nobs0,nobs0) ::VC,Corr
  double precision,dimension(nobs0*(nobs0+1)/2) ::Vi
  double precision,dimension(nv0) :: b0,b2
  double precision,dimension(nvc0+1)::mvc
  double precision,dimension(nea0*(nea0+1)/2)::vectVarB
  double precision::det,eta0,jacobien,beta_densite,ytemp,Y4,YRE2
  double precision,dimension(nobs0) :: mu,Y1,Y2,Y3,tcor,wi,wsim
  double precision :: asim, ai
  double precision,dimension(-1:maxval(ntrtot0)-3)::splaa
  double precision::aa1,bb1,dd1,aa,bb,betai,cc1
  double precision,dimension(nea0)::ui,YRE
  double precision::binf,bsup,som,eps
  double precision::div,vrais_l,vrais_u
  double precision,external::alnorm



  !      write(*,*)'indice entres',indiceY0

  !          print*,"Y0=",Y0(1:10)
  !          print*,"X0=",X0(nobs0-1:nobs+12)
  !          print*,"zitr0=",zitr0

  ! sorties initialisees




  allocate(rangeY(ny0),minY(ny0),maxY(ny0),idlink(ny0),ntrtot(ny0),epsY(ny0))



  nySPL=0
  nyORD=0
  rangeY=0
  epsY=epsY0
  do k=1,ny0
     idlink(k)=idlink0(k)
     minY(k)=zitr0(1,k)
     maxY(k)=zitr0(nbzitr0(k),k)
     if (idlink(k).eq.2) then
        nySPL=nySPL+1
     end if
     if (idlink(k).eq.3) then
        nyORD=nyORD+1
     end if
  end do

  !        print*,"min,max",minY,maxY
  !      if (verbose==1) write(*,*)'nySPL',nySPL

  if(nySPL.gt.0) then 
     allocate(nvalSPL(nySPL))
     nvalSPL=0
  else
     allocate(nvalSPL(1))
     nvalSPL(1) = 0
  end if

  if(nyORD.gt.0) then
     allocate(nvalORD(nyORD))
     nvalORD=0
  else
     allocate(nvalORD(1))
     nvalORD(1) = 0
  end if

  k1=0
  k2=0
  do k=1,ny0
     if(idlink(k).eq.2) then
        k1=k1+1
        nvalSPL(k1)=nvalSPLORD0(k)
     else if (idlink(k).eq.3) then
        k2=k2+1
        nvalORD(k2)=nvalSPLORD0(k)
     end if
  end do
  ntotvalSPL=sum(nvalSPL(:))
  ntotvalORD=sum(nvalORD(:))


  if(all(idlink.ne.2)) then
     allocate(zitr(1,1))
     allocate(mm(1),mm1(1),mm2(1),im(1),im1(1),im2(1))
     mm(1)=0.d0
     mm1(1)=0.d0
     mm2(1)=0.d0
     im(1)=0.d0
     im1(1)=0.d0
     im2(1)=0.d0
  else
     allocate(zitr(-1:(maxval(nbzitr0)+2),nySPL))
     allocate(mm(ntotvalSPL),mm1(ntotvalSPL),mm2(ntotvalSPL),im(ntotvalSPL),im1(ntotvalSPL),im2(ntotvalSPL))
  end if


  zitr=0.d0  
  k1=0
  k2=0
  do k=1,ny0
     if (idlink(k).eq.0) ntrtot(k)=2     
     if (idlink(k).eq.1) ntrtot(k)=4
     if (idlink(k).eq.2) then
        k1=k1+1
        ntrtot(k)=nbzitr0(k)+2

        zitr(1:nbzitr0(k),k1)=zitr0(1:nbzitr0(k),k)
        zitr(-1,k1)=zitr(1,k1)
        zitr(0,k1)=zitr(1,k1)
        zitr(ntrtot(k)-1,k1)=zitr(ntrtot(k)-2,k1)
        zitr(ntrtot(k),k1)=zitr(ntrtot(k)-1,k1)
     end if
     if (idlink(k).eq.3) then
        k2 = k2+1
        ntrtot(k) = nvalORD(k2)-1
     end if
  end do

  !if (verbose==1)       write(*,*)'zitr',zitr



  allocate(Y(nobs0),idprob(nv0),X(nobs0,nv0),uniqueY(ntotvalSPL+ntotvalORD) &
       ,idea(nv0),idg(nv0),idcor(nv0),idcontr(nv0),nmes(ns0,ny0),indiceY(nobs0))

  ! enregistrement pour les modules
  ny=ny0
  ns=ns0
  ng=ng0
  nv=nv0
  nobs=nobs0
  ncor=ncor0
  nalea=nalea0
  if (nwg0.eq.0) then
     nwg=0
  else
     nwg=ng-1
  end if

  idiag=idiag0

  !     if (verbose==1) write(*,*)'ntotvalSPL',ntotvalSPL

  if (ntotvalSPL+ntotvalORD.gt.0) uniqueY(1:ntotvalSPL+ntotvalORD)=uniqueY0(1:ntotvalSPL+ntotvalORD)

  nmes=0
  Y=0.d0
  X=0.d0
  idprob=0
  idea=0
  idg=0
  idcor=0
  idcontr=0
  ktemp=0



  do k=1,nv
     idprob(k)=idprob0(k)
     idea(k)=idea0(k)
     idg(k)=idg0(k)
     idcor(k)=idcor0(k)
     idcontr(k)=idcontr0(k)

     jtemp=0
     DO i=1,ns
        do yk=1,ny            
           if (k.eq.1) then
              nmes(1,yk)=nmes0(yk)   !dim(nmes)=ny    
              do j=1,nmes(1,yk)
                 jtemp=jtemp+1
                 Y(jtemp)=Y0(jtemp)
                 indiceY(jtemp)=indiceY0(jtemp)
                 ktemp=ktemp+1
                 X(jtemp,k)=X0(ktemp)
              end do
           else
              do j=1,nmes(1,yk)
                 ktemp=ktemp+1
                 jtemp=jtemp+1
                 X(jtemp,k)=X0(ktemp)
              end do
           end if
        end do
     end do
  end do
  !         write(*,*)'X k:',X(1:50,k)


  ! creation des parametres

  nea=0
  ncg=0
  ncssg=0
  nprob=0
  ncontr=0
  do k=1,nv
     if (idg(k).eq.1) then
        ncssg=ncssg+1      ! nb var. sans melange
     else if (idg(k).eq.2) then
        ncg=ncg+1      ! nb var. dans melange
     end if
     nea=nea+idea(k)
     nprob=nprob+(idprob(k))*(ng-1)
     ncontr=ncontr+idcontr(k)*(ny-1) 
  end do



  !  nb effets fixes = nb effets fixes sans melange
  !                  + ng fois le nb de var dans melange


  if (idiag.eq.1) then
     nvc=nea-1
  else if(idiag.eq.0) then
     nvc=(nea+1)*nea/2-1
  end if

  nef=ncssg+ncg*ng-1
  npmtot=nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sum(ntrtot(:))


  chol=chol0

!  print*,"npmtot=", npmtot

  ! points qmc
  nMC = nMC0
  allocate(seqMC(dimMC0*nMC))
  seqMC = seqMC0(1:dimMC0*nMC) 
  

  ! base de splines
  if (any(idlink.eq.2)) then 
     call design_splines_multo(ier)
     if (ier.eq.-1) then
        go to 1589
     end if
  end if

  if((ng.eq.1.and.ncg.gt.0).or.(ng.eq.1.and.nprob.gt.0)) then
     !      if(verbose==1) write(*,*)"ng",ng,"ncg",ncg,"nprob",nprob
     go to 1589
  end if


!!!!! calcul de f(Y/b)f(b)
  !! on suppose que ns=1
  nmescur=0

  !! on est sur la classe g
  g = g0

  Ut=0.d0
  Ut(1,1)=1.d0
  if (nea>1) then 

     If (idiag.eq.1) then
        do j=2,nea
           do k=2,nea
              if (j.eq.k) then
                 Ut(j,k)=b1(nprob+nef+ncontr+j-1)
              else
                 Ut(j,k)=0.d0
              end if
           end do
        end do
     end if

     If (idiag.eq.0) then
        do j=2,nea
           do k=1,j
              ! Ut triangulaire inferieure
              Ut(j,k)=b1(nprob+nef+ncontr+k-1+j*(j-1)/2)
           end do
        end do
     end if
        
  end if

!print*,"Ut=",Ut
  ! creation de Zi
  Z=0.d0
  l=0
  do k=1,nv
     if (idea(k).eq.1) then
        l=l+1
        do j=1,sum(nmes(1,:))
           Z(j,l)=dble(X(nmescur+j,k))
        end do
     end if
  end do


  !matrice Corr variance de BM/AR
  Corr=0.d0
  tcor=0.d0
  if (ncor.gt.0) then
     do k=1,nv
        if (idcor(k).eq.1) then
           do j=1,sum(nmes(1,:))
              tcor(j) = X(nmescur+j,k)
           end do
        end if
     end do
     do j1=1,sum(nmes(1,:))
        do j2=1,sum(nmes(1,:))
           if (ncor.eq.1) then 
              Corr(j1,j2) = Corr(j1,j2)+b1(nprob+nef+ncontr+nvc+nwg+ncor)* &
                   b1(nprob+nef+ncontr+nvc+nwg+ncor)*min(tcor(j1),tcor(j2))
           else if (ncor.eq.2) then
              Corr(j1,j2) = Corr(j1,j2)+b1(nprob+nef+ncontr+nvc+nwg+ncor)* &
                   b1(nprob+nef+ncontr+nvc+nwg+ncor)* &
                   exp(-b1(nprob+nef+ncontr+nvc+nwg+1)*abs(tcor(j1)-tcor(j2)))
           end if
        end do
     end do
     
     ! passer en cholesky si on a de l ordinal
     if(any(idlink.eq.3) .or. nMC.ne.0) then
        jj=0
        Vi=0.d0
        do j=1,sum(nmes(1,:))
           do k=j,sum(nmes(1,:))
              jj=j+k*(k-1)/2
              Vi(jj)=Corr(j,k)
           end do
        end do

        CALL DMFSD(Vi,sum(nmes(1,:)),EPS,IER)

        Corr=0.d0
        do j=1,sum(nmes(1,:))
           do k=1,j
              Corr(j,k)=Vi(k+j*(j-1)/2)
           end do
        end do
     end if
  end if

  ! creation de Y1
  Y1=0.d0
  splaa=0.d0

  sumMesYk = 0
  sumntrtot=0
  numSPL=0
  do yk=1,ny

     if (idlink(yk).eq.0) then  ! Linear link

        do j=1,nmes(1,yk)
           Y1(sumMesYk+j)=(dble(Y(nmescur+sumMesYk+j))-b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1)) &
                /abs(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+2))

           jacobien = jacobien - log(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+2))
        end do

     else if (idlink(yk).eq.1) then  ! Beta link


        aa1=exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1))/ &
             (1+exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1)))
        bb1=exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+2))/ &
             (1+exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+2)))
        bb1=aa1*(1.d0-aa1)*bb1

        cc1=abs(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+3))

        dd1=abs(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+4))

        aa=aa1*aa1*(1-aa1)/bb1-aa1
        bb=aa*(1-aa1)/aa1

        do j=1,nmes(1,yk)

           ytemp=(dble(Y(nmescur+sumMesYk+j))-minY(yk)+epsY(yk))/(maxY(yk)-minY(yk)+2*epsY(yk))
           Y1(sumMesYk+j)=(betai(aa,bb,ytemp)-cc1)/dd1


           if (abs(Y1(sumMesYk+j) - 999.d0).lt.1.d-8) then
              vraispost=-1.d9
              !print*,"-1.d9 Y1=999"
              goto 1589
           end if

           jacobien = jacobien + log(abs(beta_densite(ytemp,aa,bb))/dd1)
           jacobien=jacobien-log(abs(maxY(yk)-minY(yk)+2*epsY(yk)))
        end do

     else if (idlink(yk).eq.2) then ! Splines link
        numSPL=numSPL+1

        splaa=0.d0
        eta0=0.d0
        eta0=b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1)

        do kk=2,ntrtot(yk)
           splaa(kk-3)=b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+kk)&
                *b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+kk)
        end do
        do j=1,nmes(1,yk)
           ll=0
           if (abs(Y(nmescur+sumMesYk+j) - zitr(ntrtot(yk)-2,numSPL)).lt.1.d-8) then
              ll=ntrtot(yk)-3
           end if

           som=0.d0
           do kk = 2,ntrtot(yk)-2
              if ((Y(nmescur+sumMesYk+j).ge.zitr(kk-1,numSPL)).and. &
                   (Y(nmescur+sumMesYk+j).lt.zitr(kk,numSPL))) then
                 ll=kk-1
              end if
           end do

           if (ll.lt.1.or.ll.gt.ntrtot(yk)-3) then          
              vraispost=-1.d9
              !print*,"-1.d9 ll<1 ou ll>ntrtot-3",ll!," ntrtot=",ntrtot(yk)," numSPL=",numSPL," y=",Y(nmescur+sumMesYk+j)
              goto 1589
           end if
           if (ll.gt.1) then
              do ii=2,ll
                 som=som+splaa(ii-3)
              end do
           end if



           Y1(sumMesYk+j)=eta0+som +splaa(ll-2)*im2(indiceY(nmescur+sumMesYk+j)) &
                +splaa(ll-1)*im1(indiceY(nmescur+sumMesYk+j))&
                + splaa(ll)*im(indiceY(nmescur+sumMesYk+j))

           jacobien = jacobien + log(splaa(ll-2)*mm2(indiceY(nmescur+sumMesYk+j)) &
                +splaa(ll-1)*mm1(indiceY(nmescur+sumMesYk+j))&
                +splaa(ll)*mm(indiceY(nmescur+sumMesYk+j)))

        end do
     else if (idlink(yk).eq.3) then
        do j=1,nmes(1,yk)
           Y1(sumMesYk+j)=Y(nmescur+sumMesYk+j)
        end do
     end if
     sumMesYk=sumMesYk+nmes(1,yk)
     sumntrtot=sumntrtot+ntrtot(yk)
  end do !fin boucle yk

!  print*,"Y1=",Y1
!  print*,"indiceY=",indiceY


  if(chol.eq.0) then
     ! parametrisation sd et cor
     varB=0.d0
     do j1=1,nea
        do j2=1,j1
           if(j1.eq.j2) then
              varB(j1,j2) = Ut(j1,j2)*Ut(j1,j2)
           else
              varB(j1,j2) = (exp(Ut(j1,j2))-1)/(exp(Ut(j1,j2))+1)
              varB(j2,j1) = varB(j1,j2)
           end if
        end do
     end do
     do j1=1,nea
        do j2=1,nea
           if(j1.ne.j2) then
              varB(j1,j2) = varB(j1,j2)*sqrt(varB(j1,j1)*varB(j2,j2))
           end if
        end do
     end do

     ! calculer la cholesky pour les points MC
     mvc=0.d0
     j=0
     do j1=1,nea
        do j2=1,j1
           j = j+1
           mvc(j)=varB(j1,j2)
        end do
     end do
     CALL dmfsd(mvc,nea,EPS,IER)
     Ut=0.d0
     j=0
     do j1=1,nea
        do j2=1,j1
           j = j+1
           Ut(j1,j2) = mvc(j)
        end do
     end do

  else
     ! parametrisation cholesky
     varB = matmul(Ut,transpose(Ut))
  end if


  
  
  
  vraispost=0.d0


     ui=0.d0
     if(nea.gt.0) then
        do k=1,nea
           ui(k) = RE(k)
        end do
     end if
  !   print*,"ui=", ui
     g = g0
     
     b0=0.d0
     b01=0.d0
     l=0
     m=0
     X00=0.d0
     X01=0.d0
     X2=0.d0
     nmoins = 0
     l2=0
     m2=0
     do k=1,nv
        ! mixture
        if (idg(k).eq.2) then 
           l=l+1
           do j=1,sum(nmes(1,:))
              X2(j,l)=dble(X(nmescur+j,k))
           end do
           
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
           
           ! fixed 
        else if (idg(k).eq.1) then 
           m=m+1
           do j=1,sum(nmes(1,:))
              X00(j,m)=dble(X(nmescur+j,k))
           end do
           
           ! parametre a 0 pour l'intercept
           if (k.eq.1) then
              m2=m2+1
              b0(m2)=0.d0
           else
              m2=m2+1
              b0(m2)=b1(nprob+nmoins+1)
              nmoins=nmoins+1
           end if
        end if

        !contrast : 
        if (idcontr(k).ne.0) then
           m=m+1
           sumMesYk=0
           do yk=1,ny
              ! creation matrice design des contrastes: X01
              do j=1,nmes(1,yk)
                 X01(sumMesYk+j,(m-1)*ny+yk) = dble(X(nmescur+sumMesYk+j,k))
              end do
              sumMesYk=sumMesYk+nmes(1,yk)
              ! creation vecteur parms des contrastes: b01
              if (yk<ny) then
                 b01((m-1)*ny+yk)=b1(nprob+nef+(m-1)*(ny-1)+yk)
              else
                 b01((m-1)*ny+ny) =-sum(b1(nprob+nef+(m-1)*(ny-1)+1 &
                      :nprob+nef+(m-1)*(ny-1)+ny-1))
              end if
           end do
        end if
     end do
     

     ! esperance conditionnelle
     mu = matmul(X00,b0)+matmul(X2,b2)+matmul(X01,b01)+matmul(Z,ui)
     !print*,"mu = ", mu
     !print*,"nMC = ", nMC
     som = 0.d0
     do l = 1,nMC
        
        ! simuler le BM ou AR
        if(ncor.gt.0) then
           
           wsim=0.d0
           do j=1,sum(nmes(1,:))
              wsim(j)=seqMC(nMC*(j-1)+l)
           end do
           wi=0.d0
           wi=matmul(Corr,wsim)

           mu = mu + wi
           
        end if

        
        
     !if(ncor.gt.0) mu = mu+wi
!     print*,"mu=",mu
!     print*,"b1=",b1
!     print*,"u=",ui
!     print*,"a=",ai
     sumMesYk=0
     sumntr=0
     ykord=0
     vrais_l=1.d0
     do yk =1,ny

        !! simuler l'EA specifique au test
        ai=0.d0
        if(nalea.gt.0) then
           asim = seqMC(nMC * (sum(nmes(1,:)) + yk - 1) + l)
           ai = b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)*asim
        end if
        !print*,"ai=",ai
        if(idlink(yk).eq.3) then
           !! yk est ordinal
           ykord = ykord + 1

           do j=1,nmes(1,yk)

              !! trouver binf et bsup tq binf < lambda + epsilon < bsup

              !! initialiser au premier seuil
              binf = b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntr+1)
              bsup = binf

              !! si Y>minY ajouter b1(..)^2
              if(indiceY(nmescur+sumMesYk+j).gt.1) then
                 do ll=2,min(indiceY(nmescur+sumMesYk+j),ntrtot(yk))
                    bsup = bsup + b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntr+ll)**2
                    if(ll.lt.indiceY(nmescur+sumMesYk+j)) then
                       binf = binf + b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntr+ll)**2
                    end if
                 end do
              end if

              !! centrer et standardiser
              binf = (binf - mu(sumMesYk+j) - ai) / b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)
              bsup = (bsup - mu(sumMesYk+j) - ai) / b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)

              if(indiceY(nmescur+sumMesYk+j).eq.1) then
                 !! si Y=minY
                 vrais_l = vrais_l * alnorm(binf,.false.)

              else if(indiceY(nmescur+sumMesYk+j).eq.nvalORD(ykord)) then
                 !! si Y=maxY
                 vrais_l = vrais_l * (1.d0-alnorm(bsup,.false.))

              else
                 !! minY < Y < maxY
                 vrais_l = vrais_l * (alnorm(bsup,.false.)-alnorm(binf,.false.))

              end if
!print*,"vrais_l=",vrais_l
           end do


        else
           !! yk est continu

           if(nmes(1,yk).gt.0) then
              
              !! inverse de variance de Y|ui,wi
              VC=0.d0
              det=0.d0
              do j1=1,nmes(1,yk)
                 VC(j1,j1) = 1.d0 / b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)**2 !variance de l'erreur yk
                 det = det + dlog(b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)**2)
              end do

!print*,"det =", det
              ! calcul de la vrais
              Y2=0.d0
              Y3=0.d0
              Y4=0.d0
              do j=1,nmes(1,yk)
                 Y2(j) = Y1(sumMesYk+j)-mu(sumMesYk+j)-ai
              end do
              Y3=matmul(VC,Y2)
              Y4=DOT_PRODUCT(Y2,Y3)

              div = (dble(2*3.14159265)**(dble(nmes(1,yk)/2.d0)))*sqrt(exp(det))
!print*, "Y4=",Y4, "  div=",div
vrais_l = vrais_l * exp(-Y4/2.d0)/div
!print*,"vrais_l=", vrais_l
           end if

        end if

        sumMesYk = sumMesYk+nmes(1,yk)
        sumntr = sumntr + ntrtot(yk)
     end do ! fin boucle yk

     som = som + vrais_l
     
  end do ! fin boucle MC
  
  vraispost = log(som)- dlog(dble(nMC))

  !! calcul de f(u)
  if(nea.gt.0) then
     
     ! VarB en vecteur
     jj=0
     VectVarB=0.d0
     do j1=1,nea
        do j2=j1,nea
           jj=j1+j2*(j2-1)/2
           if(nwg.eq.0.or.g.eq.ng) then
              VectVarB(jj) = VarB(j1,j2)
           else
              VectVarB(jj) = b1(nprob+nef+ncontr+nvc+g)*b1(nprob+nef+ncontr+nvc+g)*VarB(j1,j2)
           end if
        end do
     end do

     ! inversion
     CALL dsinv(VectVarB,nea,eps,ier,det)
     if (ier.eq.-1) then
        vraispost=-1.d9
        !print*,"-1.d9 dsinv continu MC"
        goto 1589
     end if
     !print*, "inv B" , VectVarB
     ! retransformation du vecteur en matrice :
     VarB=0.d0
     do j1=1,nea
        do j2=1,nea
           if (j2.ge.j1) then
              VarB(j1,j2)=VectVarB(j1+j2*(j2-1)/2)
           else
              VarB(j1,j2)=VectVarB(j2+j1*(j1-1)/2)
           end if
        end do
     end do
     
     
     YRE = MATMUL(VarB,ui)
     YRE2 = sum(ui*YRE)
     
     vrais_u = -nea*dlog(dble(2*3.14159265))-det
     vrais_u = vrais_u - YRE2
     vrais_u = vrais_u/2.d0
     !print*,"vrais_u=", vrais_u
     vraispost = vraispost + vrais_u
      !      print*,"vraispost=",vraispost
  end if
  
  
1589 continue
  
  deallocate(Y,X,idprob,idea,idg,idcor,idcontr,nmes,uniqueY,indiceY,ntrtot,seqMC)
  
  
  deallocate(zitr,mm,mm1,mm2,im,im1,im2,minY,maxY,rangeY,idlink,nvalSPL,nvalORD,epsY)
  
  
  !write(*,*)'fin'
  return
end subroutine postREmult
