
module communmo

  implicit none
  integer,save ::ny,ns,ng,nv,idiag,ncssg,nvc,nea,ncg,nwg,ncor &
       ,nprob,nvarprob,maxmes,nobs,nef,ncontr,nalea &
       ,nySPL,ntotvalSPL,nyORD,ntotvalORD,npmtot &
       ,nMC,methInteg,chol
  double precision,dimension(:),allocatable,save::Y,uniqueY,minY,maxY,rangeY
  double precision,dimension(:,:),allocatable,save ::X
  integer,dimension(:),allocatable,save ::idea,idg,idprob,idcor,idcontr,indiceY
  integer,dimension(:),allocatable,save ::prior
  integer,dimension(:),allocatable,save ::idlink,ntrtot
  integer,dimension(:,:),allocatable,save ::nmes
  integer,dimension(:),allocatable,save::nvalSPL,nvalORD
  double precision,dimension(:),allocatable,save :: seqMC
  double precision,dimension(:),allocatable,save :: epsY
  double precision,dimension(:,:),allocatable,save::zitr
  double precision,dimension(:),allocatable,save::mm,mm1,mm2,im,im1,im2
  integer,dimension(:),allocatable,save::fix
  double precision,dimension(:),allocatable,save::bfix
end module communmo


module donnees_indivmo

  implicit none
  double precision,dimension(:,:),allocatable::Ut1
  double precision,dimension(:),allocatable::mu
  double precision,dimension(:,:),allocatable::Z
  integer::numpat,nmescur
  integer,parameter ::nf=1
  double precision,dimension(:),allocatable,save::seuils
end module donnees_indivmo



!-----------------------------------------------------------
!                       VRAIS_MULTO_i
!------------------------------------------------------------


double precision function vrais_multo_i(b,npm,id,thi,jd,thj,i) 

!  use parameters
  use communmo
!  use optim
  use donnees_indivmo,only:nmescur

  IMPLICIT NONE
  integer ::i,j,k,l,m,g,l2,m2,id,jd,jj,npm,ll,ii,numSPL,ykord
  integer ::ier,nmoins,kk,j1,j2,q,sumMesYk,yk,sumntrtot,sumntr
  double precision,dimension(maxmes,nv) ::X00,X2
  double precision,dimension(maxmes,nea) ::Z,P
  double precision,dimension(maxmes,(ncontr+sum(idcontr)))::X01
  double precision,dimension(ncontr+sum(idcontr))::b01
  double precision,dimension(nprob) ::Xprob,bprob
  double precision,dimension(nea,nea) ::Ut,Ut1,varB
  double precision,dimension(maxmes,maxmes) ::VC,Corr
  double precision,dimension(npm) :: b
  double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
  double precision,dimension(nv) :: b0,b2
  double precision,dimension(npmtot)::b1
  double precision,dimension(nvc+1)::mvc

  double precision :: vrais,eps,det,som,thi,thj,temp,eta0,vrais_Y
  double precision ::Y4,jacobien,beta_densite,ytemp
  double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
  double precision,dimension(ng) :: pi
  double precision,dimension(-1:maxval(ntrtot)-3)::splaa
  double precision::aa1,bb1,dd1,aa,bb,betai,cc1
  double precision,dimension(nea)::ui,usim
  double precision,dimension(maxmes)::wi,wsim
  double precision,dimension(ny)::asim
  double precision::ai,binf,bsup
  double precision::SX,x22,div,vrais_l
  double precision,external::alnorm


  ! definir le nombre total de mesures pour le sujet i : nmestot (valable que pour cette fonction)

  ! if (verbose==1) write(*,*)'i',i 

  !if(i==1 .and. id==0 .and. jd==0) then
  !print*, "b=",b   
  !end if
  b1=0.d0
  eps=1.D-20
  l=0
  m=0
  do k=1,npmtot
     if(fix(k).eq.0) then
        l=l+1
        b1(k)=b(l)
        if(id.eq.l) b1(k)=b1(k)+thi
        if(jd.eq.l) b1(k)=b1(k)+thj
     end if
     if(fix(k).eq.1) then
        m=m+1
        b1(k)=bfix(m)
     end if
  end do


  !----------- rappel des parametres utilises ---------



  !        write(*,*)'i',i,nmescur,nmescur + nmes(i)


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

  vrais_Y=0.d0
  jacobien=0.d0
  ! -------- creation de Vi = ZiGZi'+se*seIni ----------
  ! creation de Zi

  Z=0.d0
  l=0
  do k=1,nv
     if (idea(k).eq.1) then
        l=l+1
        do j=1,sum(nmes(i,:))
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
           do j=1,sum(nmes(i,:))
              tcor(j) = X(nmescur+j,k)
           end do
        end if
     end do
     do j1=1,sum(nmes(i,:))
        do j2=1,sum(nmes(i,:))
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
        do j=1,sum(nmes(i,:))
           do k=j,sum(nmes(i,:))
              jj=j+k*(k-1)/2
              Vi(jj)=Corr(j,k)
           end do
        end do

        CALL DMFSD(Vi,sum(nmes(i,:)),EPS,IER)

        Corr=0.d0
        do j=1,sum(nmes(i,:))
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

     ! si que du continu, ajouter alpha et sigma dans corr
     if(all(idlink.ne.3) .and. nMC.eq.0) then
        do j1=1,nmes(i,yk)
           Corr(sumMesYk+j1,sumMesYk+j1) =  Corr(sumMesYk+j1,sumMesYk+j1)+b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)**2 !variance de l'erreur k
           if (nalea.eq.ny) then
              do j2=1,nmes(i,yk)
                 Corr(sumMesYk+j1,sumMesYk+j2) = Corr(sumMesYk+j1,sumMesYk+j2)+b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)**2 ! effet aleatoire specfique au test k
              end do
           end if
        end do
     end if

     if (idlink(yk).eq.0) then  ! Linear link

        do j=1,nmes(i,yk)
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

        do j=1,nmes(i,yk)

           ytemp=(dble(Y(nmescur+sumMesYk+j))-minY(yk)+epsY(yk))/(maxY(yk)-minY(yk)+2*epsY(yk))
           Y1(sumMesYk+j)=(betai(aa,bb,ytemp)-cc1)/dd1


           if (Y1(sumMesYk+j).eq.999.d0) then
              vrais_multo_i=-1.d9
              !print*,"-1.d9 Y1=999"
              goto 654
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
        do j=1,nmes(i,yk)
           ll=0
           if (Y(nmescur+sumMesYk+j).eq.zitr(ntrtot(yk)-2,numSPL)) then
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
              vrais_multo_i=-1.d9
              !print*,"-1.d9 ll<1 ou ll>ntrtot-3",ll!," ntrtot=",ntrtot(yk)," numSPL=",numSPL," y=",Y(nmescur+sumMesYk+j)
              goto 654
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
        do j=1,nmes(i,yk)
           Y1(sumMesYk+j)=Y(nmescur+sumMesYk+j)
        end do
     end if
     sumMesYk=sumMesYk+nmes(i,yk)
     sumntrtot=sumntrtot+ntrtot(yk)
  end do !fin boucle yk



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
     varB = matmul(transpose(Ut),Ut)
  end if
  


  ! contribution individuelle a la vraisemblance
  vrais_multo_i=0.d0

  if (ng.eq.1) then
     ! sans classes latentes : ng=1
     
     b0=0.d0
     b01=0.d0
     l=0
     m=0
     X00=0.d0
     X01=0.d0
     do k=1,nv
        if (idg(k).ne.0) then
           l=l+1
           do j=1,sum(nmes(i,:))
              X00(j,l)=dble(X(nmescur+j,k))  
           end do
           ! idg ne 0 pour l'intercept forcement donc on met le parm a 0
           if (k.eq.1) then
              b0(l)=0.d0
           else
              b0(l)=b1(nprob+l-1)
           end if
        end if

        !contrast : 
        if (idcontr(k).ne.0) then
           m=m+1
           sumMesYk=0
           do yk=1,ny
              ! creation matrice design des contrastes: X01
              do j=1,nmes(i,yk)
                 X01(sumMesYk+j,(m-1)*ny+yk) = dble(X(nmescur+sumMesYk+j,k))
              end do
              sumMesYk=sumMesYk+nmes(i,yk)
              ! creation vecteur parms des contrastes: b01
              if (yk<ny) THEN
                 b01((m-1)*ny+yk)=b1(nprob+nef+(m-1)*(ny-1)+yk)
              else
                 b01((m-1)*ny+ny) =-sum(b1(nprob+nef+(m-1)*(ny-1)+1 &
                      :nprob+nef+(m-1)*(ny-1)+ny-1))
              end if
           end do
        end if
     end do

     if(all(idlink.ne.3) .and. nMC.eq.0) then
!!!!!!!!!!! que du continu !!!!!!!!!!!!!!!
        P=0.d0
        P=MATMUL(Z,Ut)
        VC=0.d0
        VC=MATMUL(P,transpose(P))+Corr

        if(chol.eq.0) then
           P=0.d0
           P=MATMUL(Z,varB)
           VC=0.d0
           VC=MATMUL(P,transpose(Z))+Corr
        end if
        
        ! Vi en vecteur
        jj=0
        Vi=0.d0
        do j=1,sum(nmes(i,:))
           do k=j,sum(nmes(i,:))
              jj=j+k*(k-1)/2
              Vi(jj)=VC(j,k)
           end do
        end do

        CALL dsinv(Vi,sum(nmes(i,:)),eps,ier,det)
        if (ier.eq.-1) then
           vrais_multo_i=-1.d9
           !print*,"-1.d9 dsinv que continu"
           goto 654
        end if

        !     retransformation du vecteur Vi en matrice :
        VC=0.d0
        do j=1,sum(nmes(i,:))
           do k=1,sum(nmes(i,:))
              if (k.ge.j) then
                 VC(j,k)=Vi(j+k*(k-1)/2)
              else
                 VC(j,k)=Vi(k+j*(j-1)/2)
              end if
           end do
        end do

        mu=0.d0
        y2=0.d0
        y3=0.d0
        y4=0.d0
        mu=matmul(X00,b0)+matmul(X01,b01)
        Y2=Y1-mu
        Y3=matmul(VC,Y2)
        Y4=DOT_PRODUCT(Y2,Y3)

        vrais_Y=-sum(nmes(i,:))*dlog(dble(2*3.14159265))-det
        vrais_Y = vrais_Y-Y4
        vrais_Y = vrais_Y/2.d0

     else
!!!!!!!!!! au moins 1 ordinal !!!!!!!!!!!!!

        som=0.d0
        do l=1,nMC
           
           vrais_l=1.d0

!!!!!!!!! MC pour EA et BM/AR !!!!!!!!!

           if(methInteg.eq.1) then 
              ! !!!!!!!!!!!!! MCO !!!!!!!!!!!!!

              ! simuler les effets aleatoires
              if(nea.gt.0) then
                 x22=0.d0
                 SX=1.d0
                 do j=1,nea
                    call bgos(SX,0,usim(j),x22,0.d0)
                 end do
                 ui=0.d0
                 ui=matmul(Ut,usim)
              end if

              ! simuler le BM ou AR
              if(ncor.gt.0) then
                 x22=0.d0
                 SX=1.d0
                 do j=1,sum(nmes(i,:))
                    call bgos(SX,0,wsim(j),x22,0.d0)
                 end do
                 wi=0.d0
                 wi=matmul(Corr,wsim)
              end if

           else if(methInteg.eq.2) then 
              ! !!!!!!!!!!!!! MCA !!!!!!!!!!!!!


              if(mod(l,2).eq.0) then
                 ! si l est pair on prend l'oppose des precedents
                 usim = -usim
                 wsim = -wsim
              else
                 ! sinon on simule des nouveaux

                 ! simuler les effets aleatoires
                 if(nea.gt.0) then
                    x22=0.d0
                    SX=1.d0
                    do j=1,nea
                       call bgos(SX,0,usim(j),x22,0.d0)
                    end do
                    ui=0.d0
                    ui=matmul(Ut,usim)
                 end if

                 ! simuler le BM ou AR
                 if(ncor.gt.0) then
                    x22=0.d0
                    SX=1.d0
                    do j=1,sum(nmes(i,:))
                       call bgos(SX,0,wsim(j),x22,0.d0)
                    end do
                    wi=0.d0
                    wi=matmul(Corr,wsim)
                 end if

              end if

           else 
              ! !!!!!!!!!!!!! QMC !!!!!!!!!!!!!

              ! simuler les effets aleatoires
              if(nea.gt.0) then
                 usim=0.d0
                 do j=1,nea
                    usim(j)=seqMC(nMC*(j-1)+l)
                 end do
                 ui=0.d0
                 ui=matmul(Ut,usim)
              end if
                 
              ! simuler le BM ou AR
              if(ncor.gt.0) then
                 wsim=0.d0
                 do j=1,sum(nmes(i,:))
                    wsim(j)=seqMC(nMC*(nea+j-1)+l)
                 end do
                 wi=0.d0
                 wi=matmul(Corr,wsim)
              end if
              

           end if ! fin if methInteg


           ! esperance conditionnelle
           mu = matmul(X00,b0)+matmul(X01,b01)+matmul(Z,ui)
           if(ncor.gt.0) mu = mu+wi

           sumMesYk=0
           sumntr=0
           ykord=0
           do yk =1,ny

              if(idlink(yk).eq.3) then
                 !! yk est ordinal
                 ykord = ykord + 1

                 !! MC pour simuler l'EA specifique au test
                 ai=0.d0
                 if(nalea.gt.0) then
                    if(methInteg.eq.1) then
                       !! MCO
                       call bgos(SX,0,asim(yk),x22,0.d0)
                       ai = b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)*asim(yk)
                    else if(methInteg.eq.2) then
                       !! MCA
                       if(mod(l,2).eq.0) then
                          ! si l est pair on prend l'oppose du precedent
                          ai = -ai
                       else
                          call bgos(SX,0,asim(yk),x22,0.d0)
                          ai = b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)*asim(yk)
                       end if
                    else
                       !! QMC
                       asim(yk) = seqMC(nMC*(nea+sum(nmes(i,:)))+l)
                       ai = b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)*asim(yk)
                    end if
                 end if

                 do j=1,nmes(i,yk)

                    !! on ajoute ai a mu
                    if(nalea.gt.0) mu(sumMesYk+j) = mu(sumMesYk+j)+ai

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
                    binf = (binf - mu(sumMesYk+j))/b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)
                    bsup = (bsup - mu(sumMesYk+j))/b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)

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
                    
                 end do


              else
                 !! yk est continu

                 if(nmes(i,yk).gt.0) then
                 !! variance de Y|ui,wi
                 VC=0.d0
                 do j1=1,nmes(i,yk)
                    VC(j1,j1) = b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)**2 !variance de l'erreur yk
                    if (nalea.eq.ny) then ! intercept aleatoire de yk
                       do j2=1,nmes(i,yk)
                          VC(j1,j2) = VC(j1,j2) + &
                               b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)**2
                       end do
                    end if
                 end do

                 ! Vi en vecteur
                 jj=0
                 Vi=0.d0
                 do j1=1,nmes(i,yk)
                    do j2=j1,nmes(i,yk)
                       jj=j1+j2*(j2-1)/2
                       Vi(jj)=VC(j1,j2)
                    end do
                 end do

                 ! inversion
                 CALL dsinv(Vi,nmes(i,yk),eps,ier,det)
                 if (ier.eq.-1) then
                    vrais_multo_i=-1.d9
                    !print*,"-1.d9 dsinv continu MC"
                    goto 654
                 end if

                 ! retransformation du vecteur Vi en matrice :
                 VC=0.d0
                 do j1=1,nmes(i,yk)
                    do j2=1,nmes(i,yk)
                       if (j2.ge.j1) then
                          VC(j1,j2)=Vi(j1+j2*(j2-1)/2)
                       else
                          VC(j1,j2)=Vi(j2+j1*(j1-1)/2)
                       end if
                    end do
                 end do

                 ! calcul de la vrais
                 Y2=0.d0
                 Y3=0.d0
                 Y4=0.d0
                 do j=1,nmes(i,yk)
                    Y2(j) = Y1(sumMesYk+j)-mu(sumMesYk+j)
                 end do
                 Y3=matmul(VC,Y2)
                 Y4=DOT_PRODUCT(Y2,Y3)
                 
                 div = (dble(2*3.14159265)**(dble(nmes(i,yk)/2.d0)))*sqrt(exp(det))

                 vrais_l = vrais_l * exp(-Y4/2.d0)/div
              end if

           end if

              sumMesYk = sumMesYk+nmes(i,yk)
              sumntr = sumntr + ntrtot(yk)
           end do ! fin boucle yk

           som = som + vrais_l
           
        end do ! fin boucle MC

        vrais_Y = log(som) - dlog(dble(nMC))

     end if
     ! avec classes latentes:  ng>1  
  else

     if (prior(i).ne.0) then
        pi=0.d0
        pi(prior(i))=1.d0
     else


        ! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
        Xprob=0.d0
        l=0
        do k=1,nv
           if (idprob(k).eq.1) then
              l=l+1
              Xprob(l)=X(nmescur+1,k)
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

     end if

     ! creation des vecteurs de variables explicatives
     l=0
     m=0
     q=0
     X00=0.d0
     X2=0.d0
     X01=0.d0
     b01=0.d0
     do k=1,nv
        if (idg(k).eq.2) then
           l=l+1
           do j=1,sum(nmes(i,:))
              X2(j,l)=dble(X(nmescur+j,k))
           end do
        else if (idg(k).eq.1) then
           m=m+1
           do j=1,sum(nmes(i,:))
              X00(j,m)=dble(X(nmescur+j,k))
           end do
        end if

        !contrast : 
        if (idcontr(k).ne.0) then
           q=q+1
           sumMesYk=0
           do yk=1,ny
              ! creation matrice design des contrastes: X01
              do j=1,nmes(i,yk)
                 X01(sumMesYk+j,(q-1)*ny+yk) = dble(X(nmescur+sumMesYk+j,k))
              end do
              sumMesYk=sumMesYk+nmes(i,yk)
              ! creation vecteur parms des contrastes: b01
              if (yk<ny) THEN
                 b01((q-1)*ny+yk)=b1(nprob+nef+(q-1)*(ny-1)+yk)
              else
                 b01((q-1)*ny+ny) =-sum(b1((nprob+nef+(q-1)*(ny-1)+1) &
                      :(nprob+nef+(q-1)*(ny-1)+ny-1)))
              end if
           end do
        end if
     end do
     

     
     if(all(idlink.ne.3).and.nMC.eq.0) then
        b2=0.d0
        b0=0.d0
        vrais=0.d0
        !! boucle sur les classes
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


           ! variance covariance si spec aux classes :
           Ut1=Ut
           if (nwg.ne.0) then
              Ut1=0.d0
              if (g.eq.ng) then
                 Ut1=Ut
              else
                 Ut1=Ut*abs(b1(nprob+nef+ncontr+nvc+g))
              end if
           end if


           P=0.d0
           P=MATMUL(Z,Ut1)
           VC=0.d0
           VC=MATMUL(P,transpose(P))+Corr

           ! Vi en vecteur
           Vi=0.d0
           jj=0
           do j=1,sum(nmes(i,:))
              do k=j,sum(nmes(i,:))
                 jj=j+k*(k-1)/2
                 Vi(jj)=VC(j,k)
              end do
           end do


           CALL dsinv(Vi,sum(nmes(i,:)),eps,ier,det)
           if (ier.eq.-1) then
              vrais_multo_i=-1.d9
              !print*,"-1.d9 dsinv 2"
              goto 654
           end if

           !     retransformation du vecteur Vi en matrice :
           VC=0.d0
           do j=1,sum(nmes(i,:))
              do k=1,sum(nmes(i,:))
                 if (k.ge.j) then
                    VC(j,k)=Vi(j+k*(k-1)/2)
                 else
                    VC(j,k)=Vi(k+j*(j-1)/2)
                 end if
              end do
           end do


           mu=0.d0
           y2=0.d0
           y3=0.d0
           y4=0.d0
           mu=matmul(X00,b0)+matmul(X2,b2)+matmul(X01,b01)
           Y2=Y1-mu
           Y3=Matmul(VC,Y2)
           Y4=0.d0
           Y4=DOT_PRODUCT(Y2,Y3)
           vrais_Y = vrais_Y + pi(g)*exp((-det-Y4)/2.d0)
        end do

        vrais_Y = -sum(nmes(i,:))*dlog(dble(2*3.14159265)) + log(vrais_Y)

     else
        !! au moins 1 ordinal
        !ppi(i,:) = 0.d0
        som = 0.d0   
        do l=1,nMC
           
           
           if(methInteg.eq.1) then 
              ! !!!!!!!!!!!!! MCO !!!!!!!!!!!!!
              
              ! simuler les effets aleatoires
              if(nea.gt.0) then
                 usim=0.d0
                 x22=0.d0
                 SX=1.d0
                 do j=1,nea
                    call bgos(SX,0,usim(j),x22,0.d0)
                    !print*,"usim=",usim(j)
                 end do
                 ui=0.d0
                 !ui=matmul(Ut,usim)
                 !print*," ui=",ui, " Ut=",Ut, "  nea=",nea(k)
              end if

              ! simuler le BM ou AR
              if(ncor.gt.0) then
                 wsim=0.d0
                 x22=0.d0
                 SX=1.d0
                 do j=1,sum(nmes(i,:))
                    call bgos(SX,0,wsim(j),x22,0.d0)
                 end do
                 wi=0.d0
                 wi=matmul(Corr,wsim)
              end if

              ! simuler EA specif
              if(nalea.gt.0) then
                 asim=0.d0
                 x22=0.d0
                 SX=1.d0
                 do j=1,ny
                    call bgos(SX,0,asim(j),x22,0.d0)
                 end do
              end if

           else if(methInteg.eq.2) then 
              ! !!!!!!!!!!!!! MCA !!!!!!!!!!!!!


              if(mod(l,2).eq.0) then
                 ! si l est pair on prend l'oppose des precedents
                 usim = -usim
                 wsim = -wsim
                 asim = -asim
              else
                 ! sinon on simule des nouveaux

                 ! simuler les effets aleatoires
                 if(nea.gt.0) then
                    usim=0.d0
                    x22=0.d0
                    SX=1.d0
                    do j=1,nea
                       call bgos(SX,0,usim(j),x22,0.d0)
                    end do
                    ui=0.d0
                    !ui=matmul(Ut,usim)
                 end if

                 ! simuler le BM ou AR
                 if(ncor.gt.0) then
                    wsim=0.d0
                    x22=0.d0
                    SX=1.d0
                    do j=1,sum(nmes(i,:))
                       call bgos(SX,0,wsim(j),x22,0.d0)
                    end do
                    wi=0.d0
                    wi=matmul(Corr,wsim)
                 end if

                 ! simuler EA specif
                 if(nalea.gt.0) then
                    asim=0.d0
                    x22=0.d0
                    SX=1.d0
                    do j=1,ny
                       call bgos(SX,0,asim(j),x22,0.d0)
                    end do
                 end if

              end if

           else 
              ! !!!!!!!!!!!!! QMC !!!!!!!!!!!!!

              ! simuler les effets aleatoires
              if(nea.gt.0) then
                 usim=0.d0
                 do j=1,nea
                    usim(j)=seqMC(nMC*(j-1)+l)
                 end do
                 ui=0.d0
                 !ui=matmul(Ut,usim)
              end if

              ! simuler le BM ou AR
              if(ncor.gt.0) then
                 wsim=0.d0
                 do j=1,sum(nmes(i,:))
                    wsim(j)=seqMC(nMC*(nea+j-1)+l)
                 end do
                 wi=0.d0
                 wi=matmul(Corr,wsim)
              end if

              ! simuler EA specif
              if(nalea.gt.0) then
                 asim=0.d0
                 do j=1,ny
                    asim(j) = seqMC(nMC*(nea+sum(nmes(i,:)))+l)
                 end do
              end if
           endif

           b2=0.d0
           b0=0.d0
           vrais=0.d0
           !! boucle sur les classes
           do g=1,ng
              vrais_l = 1.d0
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


              ! variance covariance si spec aux classes :
              Ut1=Ut
              if (nwg.ne.0) then
                 Ut1=0.d0
                 if (g.eq.ng) then
                    Ut1=Ut
                 else
                    Ut1=Ut*abs(b1(nprob+nef+ncontr+nvc+g))
                 end if
              end if

              !! EA avec variance specif a la classe
              ui=matmul(Ut1,usim)


              ! esperance conditionnelle
              mu = matmul(X00,b0)+matmul(X01,b01)+matmul(X2,b2)
              mu = mu+matmul(Z,ui)
              if(ncor.gt.0) mu = mu + wi

              sumMesYk=0
              sumntr=0
              ykord=0
              do yk =1,ny

                 if(idlink(yk).eq.3) then
                    !! yk est ordinal
                    ykord=ykord+1

                    !! EA specifique au test
                    ai=0.d0
                    if(nalea.gt.0) then                       
                       ai = b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)*asim(yk)
                    end if

                    do j=1,nmes(i,yk)

                       !! on ajoute ai a mu
                       if(nalea.gt.0) mu(sumMesYk+j) = mu(sumMesYk+j)+ai

                       !! trouver binf et bsup tq binf < lambda + epsilon < bsup :

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
                       binf = (binf - mu(sumMesYk+j))/b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)
                       bsup = (bsup - mu(sumMesYk+j))/b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)

                       if(indiceY(nmescur+sumMesYk+j).eq.1) then
                          !! si Y=minY
                          vrais_l = vrais_l*alnorm(binf,.false.)

                       else if(indiceY(nmescur+sumMesYk+j).eq.nvalORD(ykord)) then
                          !! si Y=maxY
                          vrais_l = vrais_l*(1.d0-alnorm(bsup,.false.))

                       else
                          !! minY < Y < maxY
                          vrais_l = vrais_l*(alnorm(bsup,.false.)-alnorm(binf,.false.))

                       end if

                    end do


                 else
                    !! yk est continu

                    if(nmes(i,yk).gt.0) then
                    !! variance de Y|ui,wi
                    VC=0.d0
                    do j1=1,nmes(i,yk)
                       VC(j1,j1) = b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)**2 !variance de l'erreur yk
                       if (nalea.eq.ny) then ! intercept aleatoire de yk
                          do j2=1,nmes(i,yk)
                             VC(j1,j2) = VC(j1,j2) + &
                                  b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)**2
                          end do
                       end if
                    end do

                    ! Vi en vecteur
                    jj=0
                    Vi=0.d0
                    do j1=1,nmes(i,yk)
                       do j2=j1,nmes(i,yk)
                          jj=j1+j2*(j2-1)/2
                          Vi(jj)=VC(j1,j2)
                       end do
                    end do

                    ! inversion
                    CALL dsinv(Vi,nmes(i,yk),eps,ier,det)
                    if (ier.eq.-1) then
                       vrais_multo_i=-1.d9
                       goto 654
                    end if

                    ! retransformation du vecteur Vi en matrice :
                    VC=0.d0
                    do j1=1,nmes(i,yk)
                       do j2=1,nmes(i,yk)
                          if (j2.ge.j1) then
                             VC(j1,j2)=Vi(j1+j2*(j2-1)/2)
                          else
                             VC(j1,j2)=Vi(j2+j1*(j1-1)/2)
                          end if
                       end do
                    end do

                    ! calcul de la vrais
                    Y2=0.d0
                    Y3=0.d0
                    Y4=0.d0
                    do j=1,nmes(i,yk)
                       Y2(j) = Y1(sumMesYk+j)-mu(sumMesYk+j)
                    end do
                    Y3=matmul(VC,Y2)
                    Y4=DOT_PRODUCT(Y2,Y3) 

                    !! densite gaussienne de yk
                    vrais_l = vrais_l * exp((-det-Y4)/2.d0)

                    if(l.eq.1) vrais_multo_i = vrais_multo_i -dble(nmes(i,yk)/2.d0)*dlog(dble(2*3.14159265))
                 end if

              end if

                 sumMesYk = sumMesYk+nmes(i,yk)
                 sumntr = sumntr + ntrtot(yk)
              end do ! fin boucle yk

              vrais = vrais + pi(g) * vrais_l
              
           end do ! fin boucle g
           
           som = som + vrais
              
        end do ! fin boucle MC
        
        vrais_Y = log(som) - dlog(dble(nMC))
     
     end if
     
  end if !fin if ng
  !print*,"i=",i," vrais_Y=",vrais_Y,"jac=",jacobien," vrais_multo_i=",vrais_multo_i
  vrais_multo_i = vrais_multo_i + vrais_Y + jacobien
!   print*,"vrais_multo_i=",vrais_multo_i
654 continue
  
  return
  
end function vrais_multo_i





double precision function vrais_multo(b,m,id,thi,jd,thj)


  use communmo,only:ns,nmes
  use donnees_indivmo,only:nmescur
  !use parameters, only:verbose

  implicit none

  integer::m,i,id,jd
  double precision::thi,thj,vrais_multo_i,temp
  double precision,dimension(m)::b

  nmescur=0
  vrais_multo=0.d0
  do i=1,ns
     temp=vrais_multo_i(b,m,id,thi,jd,thj,i)          
     vrais_multo = vrais_multo + temp
     if (temp.eq.-1.d9 .or. temp/temp.ne.1) then 
        !if (temp/temp.ne.1) write(*,*)"i=",i,"vrais= ",temp
        !if (temp.eq.-1.d9) then 
        vrais_multo = -1.d9
     ! print*,"i=",i," vrais_multo=",vrais_multo," m=",m," b=",b
        ! if(verbose==1) write(*,*)"i=",i,"vrais= ",temp
        goto 541
     end iF
     nmescur = nmescur + sum(nmes(i,:))
  end do
541 continue
  return

end function vrais_multo




! =============================================
! subroutine de creation de design matrix
! =============================================




      subroutine design_splines_multo (ier)

      use communmo

      implicit none

      integer ::jj,l,k,ier,yk,q,sumnval
      double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n,hn,hht

      ier=0
      jj=0
      l=0
      q=0
      sumnval=0
      do yk=1,ny
         if (idlink(yk).eq.2) then 
            q=q+1
            do jj=1,nvalSPL(q)      !     ou se trouve la valeur de zi

               do k = 2,ntrtot(yk)-2
                  if ((uniqueY(sumnval+jj).ge.zitr(k-1,q)).and.(uniqueY(sumnval+jj).lt.zitr(k,q))) then
                     l=k-1
                  end if
                End do  


            if (uniqueY(sumnval+jj).eq.zitr(ntrtot(yk)-2,q)) then
               l=ntrtot(yk)-3
            end if

            ht2 = zitr(l+1,q)-uniqueY(sumnval+jj)
            htm= uniqueY(sumnval+jj)-zitr(l-1,q)
            ht = uniqueY(sumnval+jj)-zitr(l,q)
            ht3 = zitr(l+2,q)-uniqueY(sumnval+jj)
            hht = uniqueY(sumnval+jj)-zitr(l-2,q)
            h = zitr(l+1,q)-zitr(l,q)
            hh= zitr(l+1,q)-zitr(l-1,q)
            hn= zitr(l+1,q)-zitr(l-2,q)
            h2n=zitr(l+2,q)-zitr(l-1,q)
            h2= zitr(l+2,q)-zitr(l,q)
            h3= zitr(l+3,q)-zitr(l,q)

            if (uniqueY(sumnval+jj).ne.zitr(ntrtot(yk)-2,q)) then
               mm2(sumnval+jj) = (3.d0*ht2*ht2)/(hh*h*hn)
               mm1(sumnval+jj) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
               mm(sumnval+jj)  = (3.d0*ht*ht)/(h3*h2*h)

            end if
            if (uniqueY(sumnval+jj).eq.zitr(ntrtot(yk)-2,q)) then
               mm2(sumnval+jj) = 0.d0
               mm1(sumnval+jj) = 0.d0
               mm(sumnval+jj)  = 3.d0/h
            end if

            if (mm2(sumnval+jj).lt.0.or.mm1(sumnval+jj).lt.0.or.mm(sumnval+jj).lt.0) then
                ier=-1
                goto 765
            end if

            im2(sumnval+jj)=hht*mm2(sumnval+jj)/(3.d0)+ h2n*mm1(sumnval+jj)/(3.d0) &
             +h3*mm(sumnval+jj)/(3.d0)
            im1(sumnval+jj)=htm*mm1(sumnval+jj)/(3.d0)+h3*mm(sumnval+jj)/(3.d0)
            im(sumnval+jj)=ht*mm(sumnval+jj)/(3.d0)

        end do
        sumnval = sumnval + nvalSPL(q)
        
   end if

end do


765     continue

      end subroutine design_splines_multo
!fin design_splines


      
      subroutine postprobmo(b,npm,PPI)
        use communmo
!        use optim
        implicit none

        integer ::i,j,k,l,m,g,l2,m2,jj,npm,ier,nmoins,kk,ii,ll,j1,j2,nmescur
        integer::q,numSPL,sumMesYk,yk,sumntrtot,sumntr,ykord
        double precision,dimension(maxmes,nea) ::Z,P
        double precision,dimension(maxmes,nv) ::X00,X2
        double precision,dimension(nprob) ::Xprob,bprob
        double precision,dimension(nea,nea) ::Ut,Ut1,varB
        double precision,dimension(maxmes,(ncontr+sum(idcontr)))::X01
        double precision,dimension(ncontr+sum(idcontr))::b01
        double precision,dimension(maxmes,maxmes) ::VC,Corr
        double precision,dimension(npm) :: b,b1
        double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
        double precision,dimension(nvc+1)::mvc
        double precision,dimension(nv) :: b0,b2
        double precision :: eps,det,som,temp,Y4,f,ytemp
        double precision,dimension(ng) ::fi,pi
        double precision,dimension(ns,ng) ::PPI
        double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
        double precision,dimension(-1:(maxval(ntrtot)-3))::splaa
        double precision::aa1,bb1,dd1,aa,bb,betai,cc1,eta0
        double precision,dimension(nea)::ui,usim
        double precision,dimension(maxmes)::wi,wsim
        double precision,dimension(ny)::asim
        double precision::ai,binf,bsup
        double precision::SX,x22,vrais_l
        double precision,external::alnorm
        

        eps=1.D-20

        PPI=0.D0

        do k=1,npm
           b1(k)=b(k)
        end do

        ! cholesky de la variance des EA
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

        ! ----------- boucle sur les individus -------------
        nmescur=0
        do i=1,ns

           ! creation de Zi
           Z=0.d0
           l=0
           do k=1,nv
              if (idea(k).eq.1) then
                 l=l+1
                 do j=1,sum(nmes(i,:))
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
                    do j=1,sum(nmes(i,:))
                       tcor(j) = X(nmescur+j,k)
                    end do
                 end if
              end do
              do j1=1,sum(nmes(i,:))
                 do j2=1,sum(nmes(i,:))
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
              if(any(idlink.eq.3)) then
                 jj=0
                 Vi=0.d0
                 do j=1,sum(nmes(i,:))
                    do k=j,sum(nmes(i,:))
                       jj=j+k*(k-1)/2
                       Vi(jj)=Corr(j,k)
                    end do
                 end do

                 CALL DMFSD(Vi,sum(nmes(i,:)),EPS,IER)

                 Corr=0.d0
                 do j=1,sum(nmes(i,:))
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

              ! si que du continu, ajouter alpha et sigma dans corr
              if(all(idlink.ne.3)) then
                 do j1=1,nmes(i,yk)
                    Corr(sumMesYk+j1,sumMesYk+j1) =  Corr(sumMesYk+j1,sumMesYk+j1)+b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)**2 !variance de l'erreur k
                    if (nalea.eq.ny) then
                       do j2=1,nmes(i,yk)
                          Corr(sumMesYk+j1,sumMesYk+j2) = Corr(sumMesYk+j1,sumMesYk+j2)+b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)**2 ! effet aleatoire specfique au test k
                       end do
                    end if
                 end do
              end if


              ! transformations
              if (idlink(yk).eq.0) then  ! Linear link

                 do j=1,nmes(i,yk)
                    Y1(sumMesYk+j)=(dble(Y(nmescur+sumMesYk+j))-b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1)) &
                         /abs(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+2))

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

                 do j=1,nmes(i,yk)

                    ytemp=(dble(Y(nmescur+sumMesYk+j))-minY(yk)+epsY(yk))/(maxY(yk)-minY(yk)+2*epsY(yk))
                    Y1(sumMesYk+j)=(betai(aa,bb,ytemp)-cc1)/dd1


                    if (Y1(sumMesYk+j).eq.999.d0) then
                       ppi=-1.d0
                       goto 654
                    end if

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

                 do j=1,nmes(i,yk)
                    ll=0
                    if (Y(nmescur+sumMesYk+j).eq.zitr(ntrtot(yk)-2,numSPL)) then
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
                       ppi=-1.d0
                       goto 654
                    end if
                    if (ll.gt.1) then
                       do ii=2,ll
                          som=som+splaa(ii-3)
                       end do
                    end if

                    Y1(sumMesYk+j)=eta0+som +splaa(ll-2)*im2(indiceY(nmescur+sumMesYk+j)) &
                         +splaa(ll-1)*im1(indiceY(nmescur+sumMesYk+j))&
                         + splaa(ll)*im(indiceY(nmescur+sumMesYk+j))

                 end do
              else if (idlink(yk).eq.3) then
                 do j=1,nmes(i,yk)
                    Y1(sumMesYk+j)=Y(nmescur+sumMesYk+j)
                 end do
              end if
              sumMesYk=sumMesYk+nmes(i,yk)
              sumntrtot=sumntrtot+ntrtot(yk)
           end do !fin boucle yk


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
              varB = matmul(transpose(Ut),Ut)
           end if


           !! proba de chaque classe 
           if (prior(i).ne.0) then
              pi=0.d0
              pi(prior(i))=1.d0
           else

              ! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
              Xprob=0.d0
              l=0
              do k=1,nv
                 if (idprob(k).eq.1) then
                    l=l+1
                    Xprob(l)=X(nmescur+1,k)
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

           end if

           ! creation des vecteurs de variables explicatives
           l=0
           m=0
           q=0
           X00=0.d0
           X2=0.d0
           X01=0.d0
           b01=0.d0
           do k=1,nv
              if (idg(k).eq.2) then ! en mixture
                 l=l+1
                 do j=1,sum(nmes(i,:))
                    X2(j,l)=dble(X(nmescur+j,k))
                 end do
              else if (idg(k).eq.1) then ! en fixed commun
                 m=m+1
                 do j=1,sum(nmes(i,:))
                    X00(j,m)=dble(X(nmescur+j,k))
                 end do
              end if

              !contrast : 
              if (idcontr(k).ne.0) then
                 q=q+1
                 sumMesYk=0
                 do yk=1,ny
                    ! creation matrice design des contrastes: X01
                    do j=1,nmes(i,yk)
                       X01(sumMesYk+j,(q-1)*ny+yk) = dble(X(nmescur+sumMesYk+j,k))
                    end do
                    sumMesYk=sumMesYk+nmes(i,yk)
                    ! creation vecteur parms des contrastes: b01
                    if (yk<ny) THEN
                       b01((q-1)*ny+yk)=b1(nprob+nef+(q-1)*(ny-1)+yk)
                    else
                       b01((q-1)*ny+ny) =-sum(b1((nprob+nef+(q-1)*(ny-1)+1) &
                            :(nprob+nef+(q-1)*(ny-1)+ny-1)))
                    end if
                 end do
              end if
           end do


           !! cas tous les y sont continus
           if(all(idlink.ne.3).and.nMC.eq.0) then
              b2=0.d0
              b0=0.d0
              f=0.d0
              fi=0.d0
              !! boucle sur les classes
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


                 ! variance covariance si spec aux classes :
                 Ut1=Ut
                 if (nwg.ne.0) then
                    Ut1=0.d0
                    if (g.eq.ng) then
                       Ut1=Ut
                    else
                       Ut1=Ut*abs(b1(nprob+nef+ncontr+nvc+g))
                    end if
                 end if


                 P=0.d0
                 P=MATMUL(Z,Ut1)
                 VC=0.d0
                 VC=MATMUL(P,transpose(P))+Corr

                 ! Vi en vecteur
                 Vi=0.d0
                 jj=0
                 do j=1,sum(nmes(i,:))
                    do k=j,sum(nmes(i,:))
                       jj=j+k*(k-1)/2
                       Vi(jj)=VC(j,k)
                    end do
                 end do


                 CALL dsinv(Vi,sum(nmes(i,:)),eps,ier,det)
                 if (ier.eq.-1) then
                    ppi=-1.d0
                    goto 654
                 end if

                 !     retransformation du vecteur Vi en matrice :
                 VC=0.d0
                 do j=1,sum(nmes(i,:))
                    do k=1,sum(nmes(i,:))
                       if (k.ge.j) then
                          VC(j,k)=Vi(j+k*(k-1)/2)
                       else
                          VC(j,k)=Vi(k+j*(j-1)/2)
                       end if
                    end do
                 end do


                 mu=0.d0
                 y2=0.d0
                 y3=0.d0
                 y4=0.d0
                 mu=matmul(X00,b0)+matmul(X2,b2)+matmul(X01,b01)
                 Y2=Y1-mu
                 Y3=Matmul(VC,Y2)
                 Y4=0.d0
                 Y4=DOT_PRODUCT(Y2,Y3)
                 
                 fi(g)=fi(g)- sum(nmes(i,:))*log(dble(2*3.14159265))
                 fi(g)=fi(g) -det
                 fi(g)=fi(g) - Y4
                 fi(g)=fi(g)/(2.d0)
                 fi(g)=exp(fi(g))
                 
              end do

              

           else
              !! cas au moins 1 y est ordinal
              
              som = 0.d0
              f=0.d0
              fi=0.d0
              do l=1,nMC


                 if(methInteg.eq.1) then 
                    ! !!!!!!!!!!!!! MCO !!!!!!!!!!!!!

                    ! simuler les effets aleatoires
                    if(nea.gt.0) then
                       usim=0.d0
                       x22=0.d0
                       SX=1.d0
                       do j=1,nea
                          call bgos(SX,0,usim(j),x22,0.d0)
                       end do
                       ui=0.d0
                    end if

                    ! simuler le BM ou AR
                    if(ncor.gt.0) then
                       wsim=0.d0
                       x22=0.d0
                       SX=1.d0
                       do j=1,sum(nmes(i,:))
                          call bgos(SX,0,wsim(j),x22,0.d0)
                       end do
                       wi=0.d0
                       wi=matmul(Corr,wsim)
                    end if

                    ! simuler EA specif
                    if(nalea.gt.0) then
                       asim=0.d0
                       x22=0.d0
                       SX=1.d0
                       do j=1,ny
                          call bgos(SX,0,asim(j),x22,0.d0)
                       end do
                    end if

                 else if(methInteg.eq.2) then 
                    ! !!!!!!!!!!!!! MCA !!!!!!!!!!!!!


                    if(mod(l,2).eq.0) then
                       ! si l est pair on prend l'oppose des precedents
                       usim = -usim
                       wsim = -wsim
                       asim = -asim
                    else
                       ! sinon on simule des nouveaux

                       ! simuler les effets aleatoires
                       if(nea.gt.0) then
                          usim=0.d0
                          x22=0.d0
                          SX=1.d0
                          do j=1,nea
                             call bgos(SX,0,usim(j),x22,0.d0)
                          end do
                          ui=0.d0
                       end if

                       ! simuler le BM ou AR
                       if(ncor.gt.0) then
                          wsim=0.d0
                          x22=0.d0
                          SX=1.d0
                          do j=1,sum(nmes(i,:))
                             call bgos(SX,0,wsim(j),x22,0.d0)
                          end do
                          wi=0.d0
                          wi=matmul(Corr,wsim)
                       end if

                       ! simuler EA specif
                       if(nalea.gt.0) then
                          asim=0.d0
                          x22=0.d0
                          SX=1.d0
                          do j=1,ny
                             call bgos(SX,0,asim(j),x22,0.d0)
                          end do
                       end if

                    end if

                 else 
                    ! !!!!!!!!!!!!! QMC !!!!!!!!!!!!!

                    ! simuler les effets aleatoires
                    if(nea.gt.0) then
                       usim=0.d0
                       do j=1,nea
                          usim(j)=seqMC(nMC*(j-1)+l)
                       end do
                       ui=0.d0
                    end if

                    ! simuler le BM ou AR
                    if(ncor.gt.0) then
                       wsim=0.d0
                       do j=1,sum(nmes(i,:))
                          wsim(j)=seqMC(nMC*(nea+j-1)+l)
                       end do
                       wi=0.d0
                       wi=matmul(Corr,wsim)
                    end if

                    ! simuler EA specif
                    if(nalea.gt.0) then
                       asim=0.d0
                       do j=1,ny
                          asim(j) = seqMC(nMC*(nea+sum(nmes(i,:)))+l)
                       end do
                    end if
                 endif

                 b2=0.d0
                 b0=0.d0
                 !! boucle sur les classes
                 do g=1,ng
                    vrais_l = 1.d0
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


                    ! variance covariance si spec aux classes :
                    Ut1=Ut
                    if (nwg.ne.0) then
                       Ut1=0.d0
                       if (g.eq.ng) then
                          Ut1=Ut
                       else
                          Ut1=Ut*abs(b1(nprob+nef+ncontr+nvc+g))
                       end if
                    end if

                    !! EA avec variance specif a la classe
                    ui=matmul(Ut1,usim)


                    ! esperance conditionnelle
                    mu = matmul(X00,b0)+matmul(X01,b01)+matmul(X2,b2)
                    mu = mu+matmul(Z,ui)
                    if(ncor.gt.0) mu = mu + wi

                    sumMesYk=0
                    sumntr=0
                    ykord=0
                    do yk =1,ny

                       if(idlink(yk).eq.3) then
                          !! yk est ordinal
                          ykord=ykord+1

                          !! EA specifique au test
                          ai=0.d0
                          if(nalea.gt.0) then                       
                             ai = b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)*asim(yk)
                          end if

                          do j=1,nmes(i,yk)

                             !! on ajoute ai a mu
                             if(nalea.gt.0) mu(sumMesYk+j) = mu(sumMesYk+j)+ai

                             !! trouver binf et bsup tq binf < lambda + epsilon < bsup :

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
                             binf = (binf - mu(sumMesYk+j))/b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)
                             bsup = (bsup - mu(sumMesYk+j))/b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)

                             if(indiceY(nmescur+sumMesYk+j).eq.1) then
                                !! si Y=minY
                                vrais_l = vrais_l*alnorm(binf,.false.)

                             else if(indiceY(nmescur+sumMesYk+j).eq.nvalORD(ykord)) then
                                !! si Y=maxY
                                vrais_l = vrais_l*(1.d0-alnorm(bsup,.false.))

                             else
                                !! minY < Y < maxY
                                vrais_l = vrais_l*(alnorm(bsup,.false.)-alnorm(binf,.false.))

                             end if

                          end do


                       else
                          !! yk est continu

                          !! variance de Y|ui,wi
                          VC=0.d0
                          do j1=1,nmes(i,yk)
                             VC(j1,j1) = b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)**2 !variance de l'erreur yk
                             if (nalea.eq.ny) then ! intercept aleatoire de yk
                                do j2=1,nmes(i,yk)
                                   VC(j1,j2) = VC(j1,j2) + &
                                        b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)**2
                                end do
                             end if
                          end do

                          ! Vi en vecteur
                          jj=0
                          Vi=0.d0
                          do j1=1,nmes(i,yk)
                             do j2=j1,nmes(i,yk)
                                jj=j1+j2*(j2-1)/2
                                Vi(jj)=VC(j1,j2)
                             end do
                          end do

                          ! inversion
                          CALL dsinv(Vi,nmes(i,yk),eps,ier,det)
                          if (ier.eq.-1) then
                             ppi=-1.d0
                             goto 654
                          end if

                          ! retransformation du vecteur Vi en matrice :
                          VC=0.d0
                          do j1=1,nmes(i,yk)
                             do j2=1,nmes(i,yk)
                                if (j2.ge.j1) then
                                   VC(j1,j2)=Vi(j1+j2*(j2-1)/2)
                                else
                                   VC(j1,j2)=Vi(j2+j1*(j1-1)/2)
                                end if
                             end do
                          end do

                          ! calcul de la vrais
                          Y2=0.d0
                          Y3=0.d0
                          Y4=0.d0
                          do j=1,nmes(i,yk)
                             Y2(j) = Y1(sumMesYk+j)-mu(sumMesYk+j)
                          end do
                          Y3=matmul(VC,Y2)
                          Y4=DOT_PRODUCT(Y2,Y3) 

                          !! densite gaussienne
                          vrais_l = vrais_l * dble(2*3.14159265)**(dble(nmes(i,yk)/2.d0))*exp((-det-Y4)/2.d0)
                      
                       end if

                       sumMesYk = sumMesYk+nmes(i,yk)
                       sumntr = sumntr + ntrtot(yk)
                    end do ! fin boucle yk

                    fi(g) = fi(g) + vrais_l / dble(nMC)
                 end do ! fin boucle g

              end do ! fin boucle MC

           end if


           f=DOT_PRODUCT(pi,fi)
           do g=1,ng
              PPI(i,g)=pi(g)*fi(g)/f
           end do
           
           nmescur = nmescur + sum(nmes(i,:))

        end do

654     continue
        return

      end subroutine postprobmo




          subroutine transfos_estimees_2(b,npm,nsim,marker,transfY)

      use communmo

      implicit none

      integer::kk,nsim,npm,j,k,yk,sumntrtot,numSPL,l
      double precision,dimension(nsim*ny)::marker,transfY
      double precision,dimension(maxval(ntrtot))::splaa
      double precision,dimension(maxval(ntrtot))::Xspl
      double precision,dimension(nsim)::mmm,mmm1,mmm2,iim,iim1,iim2
      double precision::aa1,bb1,dd1,aa,bb,betai,eps,pas,ytemp,cc1
      double precision, dimension(npm)::b,b1
      double precision ::ht,htm,ht2,ht3,hht,h,hh,h2,h3,h2n,hn


      b1=0.d0
      eps=1.D-20
      do k=1,npm
         b1(k)=b(k)
      end do

 !     print*,"dans trasnfos"
 !      write(*,*)'infos',minY,maxY,nsim,npm
 !      write(*,*)'b',(b1(j),j=1,npm)



       marker=0.d0
       transfY=0.d0
       
       sumntrtot = 0
       numSPL =0
       do yk=1,ny

       pas=(maxY(yk)-minY(yk))/dble(nsim-1)
       j=1
       marker((yk-1)*nsim+1)=minY(yk)
       do while(j.lt.nsim)
           j=j+1
           marker((yk-1)*nsim+j)=marker((yk-1)*nsim+j-1)+pas
       end do
       marker(yk*nsim)=maxY(yk)


       if (idlink(yk).eq.2) then
       numSPL = numSPL+1
        splaa=0.d0

        splaa(1)=b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1)
        do kk=2,ntrtot(yk)
           splaa(kk)=b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+kk)*b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+kk)
        end do

        !calcul de H(y) pour remplacer call estim_splines..
        do j=1,nsim
! ou se trouve la valeur
         l=0

         do k = 2,ntrtot(yk)-2
               if ((marker((yk-1)*nsim+j).ge.zitr(k-1,numSPL)).and.(marker((yk-1)*nsim+j).lt.zitr(k,numSPL))) then
                  l=k-1
               end if
        end do

           if (marker((yk-1)*nsim+j).eq.zitr(ntrtot(yk)-2,numSPL)) then
               l=ntrtot(yk)-3
            end if

!         if (l.lt.1.or.l.gt.ntrtot-1) then
!            write(*,*)'probleme estim splines',l
!            write(*,*)'j=',j,'test(j)',test(j)
!            stop
!         end if


               ht2 = zitr(l+1,numSPL)-marker((yk-1)*nsim+j)
               htm= marker((yk-1)*nsim+j)-zitr(l-1,numSPL)
               ht = marker((yk-1)*nsim+j)-zitr(l,numSPL)
               ht3 = zitr(l+2,numSPL)-marker((yk-1)*nsim+j)
               hht = marker((yk-1)*nsim+j)-zitr(l-2,numSPL)
               h = zitr(l+1,numSPL)-zitr(l,numSPL)
               hh= zitr(l+1,numSPL)-zitr(l-1,numSPL)
               hn= zitr(l+1,numSPL)-zitr(l-2,numSPL)
               h2n=zitr(l+2,numSPL)-zitr(l-1,numSPL)
               h2= zitr(l+2,numSPL)-zitr(l,numSPL)
               h3= zitr(l+3,numSPL)-zitr(l,numSPL)

               if (marker((yk-1)*nsim+j).ne.zitr(ntrtot(yk)-2,numSPL)) then
                  mmm2(j) = (3.d0*ht2*ht2)/(hh*h*hn)
                  mmm1(j) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
                  mmm(j)  = (3.d0*ht*ht)/(h3*h2*h)
               end if
               if (marker((yk-1)*nsim+j).eq.zitr(ntrtot(yk)-2,numSPL)) then
                  mmm2(j) = 0.d0
                  mmm1(j) = 0.d0
                  mmm(j)  = 3.d0/h
               end if

               iim2(j)=hht*mmm2(j)/(3.d0)+ h2n*mmm1(j)/(3.d0) &
      +h3*mmm(j)/(3.d0)
               iim1(j)=htm*mmm1(j)/(3.d0)+h3*mmm(j)/(3.d0)

               iim(j)=ht*mmm(j)/(3.d0)

!-------- transformation et IC de la transformation :

            Xspl=0.d0
            Xspl(1)=1
            do k=2,l
               Xspl(k)=1
            end do
            Xspl(l+1)=iim2(j)
            Xspl(l+2)=iim1(j)
            Xspl(l+3)=iim(j)
            transfY((yk-1)*nsim+j)= dot_product(Xspl,splaa)
      end do
        !fin H(y)


            !call   mult_estim_splines_ssstd(nsim,splaa(1:ntrtot(yk)),marker,transfY,yk) 
            
            
            

        else if (idlink(yk).eq.1) then

            aa1=exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1))/ &
             (1+exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1)))
            bb1=exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+2))/ &
             (1+exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+2)))
            bb1=aa1*(1.d0-aa1)*bb1
            cc1=b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+3)
            dd1=abs(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+4))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1
            do j=1,nsim
               ytemp=(marker((yk-1)*nsim+j)-minY(yk)+epsY(yk))/(maxY(yk)-minY(yk)+2*epsY(yk))
               transfY((yk-1)*nsim+j)=(betai(aa,bb,ytemp)-cc1)/dd1
               if (transfY((yk-1)*nsim+j).eq.999.d0) then
                  !                    write(*,*)'problem'
               end if
               
            end do


        else if (idlink(yk).eq.0) then

                 do j=1,nsim
                    transfY((yk-1)*nsim+j)=(marker((yk-1)*nsim+j)-b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1)) &
                    /abs(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+2))
                 end do
        end if
        sumntrtot = sumntrtot + ntrtot(yk)
        end do
 !     write(*,*)(marker(j),j=1,ny*nsim)
 !     write(*,*)(transfY(j),j=1,ny*nsim)          
        end subroutine transfos_estimees_2
        !fin transfos_estimeees




        
      subroutine residualsmo(b1,npm,ppi,resid_m,pred_m_g,resid_ss,pred_ss_g,pred_RE,pred_RE_Y,Yobs)

      use communmo

!      use optim

      implicit none
      integer ::i,j,k,l,m,g,l2,m2,jj,npm,ll,ii,j1,j2,sumMesYk,q,yk
      integer ::ier,nmoins,nmes_cur,n2,nmoins2,kk,numSPL,sumntrtot
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(maxmes,nv) ::X0,X2
      double precision,dimension(maxmes,(ncontr+sum(idcontr)))::X01
      double precision,dimension(nprob+1) ::Xprob,bprob
      double precision,dimension(nea) ::err2
      double precision,dimension(nea,nea) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Corr,VC1,Valea,SigmaE,CovDev
      double precision,dimension(npm) ::b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,b3
      double precision,dimension(ncontr+sum(idcontr))::b01
      double precision :: eps,det,temp
      double precision,dimension(nea,maxmes)::covUY
      double precision,dimension(maxmes) :: mu,Y1,Y2,pred1,err1,tcor
      double precision,dimension(nalea,maxmes)::err3
      double precision,dimension(nalea)::err4
      double precision,dimension(ng) :: pi
      double precision,dimension(nobs)::resid_m,resid_ss,Yobs
      double precision,dimension(nobs*ng)::pred_m_g,pred_ss_g
      double precision,dimension(ns,ng) ::PPI
      double precision,dimension(ns*nea)::pred_RE
      double precision,dimension(ns*nalea)::pred_RE_Y
      double precision,dimension(-1:(maxval(ntrtot)-3))::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,ytemp,som,cc1,eta0

      eps=1.D-20
      

!----------- rappel des parametres utilises ---------

! creation de Ut, decomposition de cholesky pour G
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
                 Ut(j,k)=b1(nprob+nef+ncontr+k-1+j*(j-1)/2)
            end do
         end do
      end if
      
      end if
      

! ----------- boucle sur les individus -------------

      resid_m=0.d0
      pred_m_g=0.d0
      resid_ss=0.d0
      pred_ss_g=0.d0
      Yobs=0.d0
      nmes_cur=0
      kk=0
      pred_RE=0.d0
      pred_RE_Y=0.d0

      do i=1,ns

!     -------- creation de Vi = ZiGZi'+Ri+se*seIni ----------
! creation de Zi

         Z=0.d0
         l=0
         do k=1,nv
            if (idea(k).eq.1) then
               l=l+1
               do j=1,sum(nmes(i,:))
                  Z(j,l)=dble(X(nmes_cur+j,k))
               end do
            end if

         end do
                 
!matrice Corri=Ri 
        Corr=0.d0
        tcor=0.d0
        if (ncor.gt.0) then
           do k=1,nv
              if (idcor(k).eq.1) then
                 do j=1,sum(nmes(i,:))
                    tcor(j) = X(nmes_cur+j,k)
                 end do
              end if
           end do

         do j1=1,sum(nmes(i,:))
            do j2=1,sum(nmes(i,:))
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
         end if  
                 
! creation de Valea et Y1
         Valea=0.d0
         Y1=0.d0

        SigmaE=0.d0
      sumMesYk=0
      sumntrtot=0
      numSPL=0
      do yk=1,ny
         do j1=1,nmes(i,yk)
            SigmaE(sumMesYk+j1,sumMesYk+j1) = b1(nprob+nef+ncontr+nvc+nwg+ncor+yk)**2 !erreur du test yk
         
            if (nalea.eq.ny) then
               do j2=1,nmes(i,yk)
                  Valea(sumMesYk+j1,sumMesYk+j2) = b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)**2
              end do
         end if
      end do

         if (idlink(yk).eq.0) then  ! Linear link

            do j=1,nmes(i,yk)
               Y1(sumMesYk+j)=dble(Y(nmes_cur+sumMesYk+j)-b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1))&
               /abs(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+2))
               Yobs(nmes_cur+sumMesYk+j)=Y1(sumMesYk+j)
            end do

         elseif (idlink(yk).eq.1) then  ! Beta link


            aa1=exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1))/ &
             (1+exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1)))
            bb1=exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+2))/ &
             (1+exp(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+2)))
            bb1=aa1*(1.d0-aa1)*bb1
            cc1=b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+3)
            dd1=abs(b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+4))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1


            do j=1,nmes(i,yk)

               !nmescur=nmescur+1

               ytemp=dble(Y(nmes_cur+sumMesYk+j)-minY(yk)+epsY(yk))/(maxY(yk)-minY(yk)+2*epsY(yk))
               Y1(sumMesYk+j)=(betai(aa,bb,ytemp)-cc1)/dd1

               if (Y1(sumMesYk+j).eq.999.d0) then
                  do k=1,nmes(i,yk)
                     resid_m(nmes_cur+sumMesYk+k)=9999.d0
                     pred_m_g(nmes_cur+sumMesYk+k)=9999.d0
                     Yobs(nmes_cur+sumMesYk+k)=9999.d0
                     resid_ss(nmes_cur+sumMesYk+k)=9999.d0
                     pred_ss_g(nmes_cur+sumMesYk+k)=9999.d0
                  end do
                  do k=1,nea
                     pred_RE((i-1)*nea+k)=9999.d0
                  end do
                  do k=1,nalea
                  pred_RE_Y((i-1)*nalea+k)=9999.d0
                  end do
                  goto 654
               else
                  Yobs(nmes_cur+sumMesYk+j)=Y1(sumMesYk+j)
               end if

            end do

         elseif (idlink(yk).eq.2) then ! Splines link
         
           numSPL=numSPL+1
           
           splaa=0.d0
           eta0=0.d0
           eta0=b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+1)

            do kk=2,ntrtot(yk)
               splaa(kk-3)=b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+kk)* &
                    b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+nalea+sumntrtot+kk)
            end do

            do j=1,nmes(i,yk)

               ll=0
               if (Y(nmes_cur+sumMesYk+j).eq.zitr(ntrtot(yk)-2,numSPL)) then
                  ll=ntrtot(yk)-3
               end if
               som=0.d0
               do kk = 2,ntrtot(yk)-2
                  if ((Y(nmes_cur+sumMesYk+j).ge.zitr(kk-1,numSPL)).and. &
                  Y(nmes_cur+sumMesYk+j).lt.zitr(kk,numSPL)) then
                     ll=kk-1
                  end if
               end do

               if (ll.lt.1.or.ll.gt.ntrtot(yk)-3) then
               do k=1,nmes(i,yk)
                  resid_m(nmes_cur+sumMesYk+k)=9999.d0
                  pred_m_g(nmes_cur+sumMesYk+k)=9999.d0
                  Yobs(nmes_cur+sumMesYk+k)=9999.d0
                  resid_ss(nmes_cur+sumMesYk+k)=9999.d0
                  pred_ss_g(nmes_cur+sumMesYk+k)=9999.d0
               end do
               do k=1,nea
                    pred_RE((i-1)*nea+k)=9999.d0
                end do
                do k=1,nalea
                    pred_RE_Y((i-1)*nalea+k)=9999.d0
                end do
               goto 654
               end if
               if (ll.ge.1) then
                  do ii=2,ll
                     som=som+splaa(ii-3)
                  end do
               end if

               Y1(sumMesYk+j)=eta0+som +splaa(ll-2)*im2(indiceY(nmes_cur+sumMesYk+j)) &
                 +splaa(ll-1)*im1(indiceY(nmes_cur+sumMesYk+j))&
                 +splaa(ll)*im(indiceY(nmes_cur+sumMesYk+j))

               Yobs(nmes_cur+sumMesYk+j)=Y1(sumMesYk+j)


            end do

         end if
         !fin if idlink
         sumMesYk = sumMesYk+nmes(i,yk)
         sumntrtot=sumntrtot+ntrtot(yk)
      end do 
      !fin boucle yk
      
      
!     debut du calcul de la vraisemblance

!     cas 1 : ng=1

         if (ng.eq.1) then


            covUY=0.d0
            VC=0.d0
            P=0.d0
            CovDev=0.d0


            P=MATMUL(Z,Ut)
            covUY=MATMUL(Ut,transpose(P))
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Corr+Valea+SigmaE

! covDev covariance matrix for individual deviation

            CovDev=MATMUL(P,transpose(P))+Corr+Valea

!     Vi en vecteur

            jj=0
            do j=1,sum(nmes(i,:))
               do k=j,sum(nmes(i,:))
                  jj=j+k*(k-1)/2
                  Vi(jj)=VC(j,k)
               end do
            end do

            CALL dsinv(Vi,sum(nmes(i,:)),eps,ier,det)
            if (ier.eq.-1) then
               do j=1,sum(nmes(i,:))
                  resid_m(nmes_cur+j)=9999.d0
                  pred_m_g(nmes_cur+j)=9999.d0
                  resid_ss(nmes_cur+j)=9999.d0
                  pred_ss_g(nmes_cur+j)=9999.d0
               end do
               do k=1,nea
                  pred_RE((i-1)*nea+k)=9999.d0
               end do
               do k=1,nalea
                  pred_RE_Y((i-1)*nalea+k)=9999.d0
               end do          
               goto 654
            end if

!     retransformation du vecteur Vi en matrice :
            VC1=0.d0
            do j=1,sum(nmes(i,:))
               do k=1,sum(nmes(i,:))
                  if (k.ge.j) then
                     VC1(j,k)=Vi(j+k*(k-1)/2)
                  else
                     VC1(j,k)=Vi(k+j*(j-1)/2)
                  end if
               end do
            end do


            b0=0.d0
            l=0
            m=0
            X0=0.d0
            X01=0.d0
            b01=0.d0
            do k=1,nv
               if (idg(k).ne.0) then
                  l=l+1
                  do j=1,sum(nmes(i,:))
                     X0(j,l)=dble(X(nmes_cur+j,k))
                  end do
                  if (k.gt.1) then
                     b0(l)=b1(nprob+l-1)
                  end if
               end if
               
               !contrast : 
               if (idcontr(k).ne.0) then
                  m=m+1
                  sumMesYk=0
                  do yk=1,ny
! creation matrice design des contrastes: X01
                     do j=1,nmes(i,yk)
                        X01(sumMesYk+j,(m-1)*ny+yk) = dble(X(nmes_cur+sumMesYk+j,k))
                     end do
                     sumMesYk=sumMesYk+nmes(i,yk)
! creation vecteur parms des contrastes: b01
                     if (yk<ny) THEN
                        b01((m-1)*ny+yk)=b1(nprob+nef+(m-1)*(ny-1)+yk)
                     else
                        b01((m-1)*ny+ny) =-sum(b1(nprob+nef+(m-1)*(ny-1)+1 &
                             :nprob+nef+(m-1)*(ny-1)+ny-1))
                     end if
                  end do                  
               end if
            end do
            !end do
            

            mu=matmul(X0,b0)+matmul(X01,b01)
            Y2=Y1-mu

            err1=0.d0
            err1=MATMUL(VC1,Y2)
            err2=0.d0
            err2=MATMUL(covUY,err1)

            pred1=0.d0
            pred1=MATMUL(CovDev,err1)



            do j=1,sum(nmes(i,:))
               resid_m(nmes_cur+j)=Y2(j)
               pred_m_g(nmes_cur+j)=mu(j)
               pred_ss_g(nmes_cur+j)=mu(j)+pred1(j)
               resid_ss(nmes_cur+j)=Y2(j)-pred1(j)
            end do


            do k=1,nea
               pred_RE((i-1)*nea+k)=err2(k)
            end do

           !pred_RE_Y 
           err3=0.d0
           if (nalea.eq.ny) then
           do yk=1,ny
            err3(yk,:) = Valea(sum(nmes(i,1:yk)),:)
           end do
           end if
           err4=0.d0
           err4= matmul(err3,err1)
           do k=1,nalea
              pred_RE_Y((i-1)*nalea+k) = err4(k)
           end do

!     cas 2 :  ng>1  composantes
         else

            if (prior(i).ne.0) then
                pi=0.d0
                pi(prior(i))=1.d0
            else

!     transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
            Xprob=0.d0
            Xprob(1)=1.d0
            l=0
            do k=1,nv
               if (idprob(k).eq.1) then
                  l=l+1
                  Xprob(1+l)=X(nmes_cur+1,k)
               end if
            end do
!     write(*,*)'l apres Xprob',l,(Xprob(j),j=1,10)
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
            pi(ng)=1.d0/(1.d0+temp)
            do g=1,ng-1
               pi(g)=pi(g)*pi(ng)
            end do


           end if
           !fin if prior
           
!     creation des vecteurs de variables explicatives
            l=0
            m=0
            q=0
            X0=0.d0
            X2=0.d0
            X01=0.d0
            b01=0.d0
            do k=1,nv
               if (idg(k).eq.2) then
                  l=l+1
                  do j=1,sum(nmes(i,:))
                     X2(j,l)=dble(X(nmes_cur+j,k))
                  end do
               else if (idg(k).eq.1) then
                  m=m+1
                  do j=1,sum(nmes(i,:))
                     X0(j,m)=dble(X(nmes_cur+j,k))
                  end do
               end if
               !contrast : 
               if (idcontr(k).ne.0) then
                  q=q+1
                  sumMesYk=0
                  do yk=1,ny
! creation matrice design des contrastes: X01
                     do j=1,nmes(i,yk)
                        X01(sumMesYk+j,(q-1)*ny+yk) = dble(X(nmes_cur+sumMesYk+j,k))
                     end do
                     sumMesYk=sumMesYk+nmes(i,yk)
! creation vecteur parms des contrastes: b01
                     if (yk<ny) THEN
                        b01((q-1)*ny+yk)=b1(nprob+nef+(q-1)*(ny-1)+yk)
                     else
                        b01((q-1)*ny+ny) =-sum(b1(nprob+nef+(q-1)*(ny-1)+1 &
                             :nef-ncontr+(q-1)*(ny-1)+ny-1))
                     end if
                  end do                  
               end if
            end do


            do g=1,ng
               nmoins=0
               l2=0
               m2=0
               n2=0
               b0=0.d0
               b2=0.d0
               b3=0.d0
               nmoins2=0

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
                  IF (IDEA(k).EQ.1.and.idg(k).eq.1) THEN
                  ! parametre a 0 pour l'intercept
                     if (k.eq.1) then
                        n2=n2+1
                        b3(n2)=0.d0
                     else
                        n2=n2+1
                        b3(n2)=b1(nprob+nmoins2+1)
                        nmoins2=nmoins2+1
                     end if
                  else if(IDEA(k).EQ.1.and.idg(k).eq.2) THEN
                    if (k.eq.1) then
                        if (g.eq.1) then
                            n2=n2+1
                            b3(n2)=0.d0
                            nmoins2=nmoins2+ng-1
                        else
                            n2=n2+1
                            b3(n2)=b1(nprob+nmoins2+g-1)
                            nmoins2=nmoins2+ng-1
                        end if
                     else
                        n2=n2+1
                        b3(n2)=b1(nprob+nmoins2+g)
                        nmoins2=nmoins2+ng
                     end if
                  end if
               end do

               mu=0.d0
               mu=matmul(X0,b0)+matmul(X01,b01)+matmul(X2,b2)

               VC=0.d0
               P=0.d0
               Ut1=Ut
               if (nwg.ne.0) then
                  if (g.eq.ng) then
                     Ut1=Ut
                  else
                     Ut1=Ut*b1(nprob+nef+ncontr+nvc+g)
                  end if
               end if
               P=0.d0
               covUY=0.d0
               VC=0.d0
               P=MATMUL(Z,Ut1)
               covUY=MATMUL(Ut1,transpose(P))
               VC=MATMUL(P,transpose(P))+Corr+Valea+SigmaE

! covDev covariance matrix for individual deviation

            CovDev=MATMUL(P,transpose(P))+Corr+Valea

!     Vi en vecteur
               jj=0
               do j=1,sum(nmes(i,:))
                  do k=j,sum(nmes(i,:))
                     jj=j+k*(k-1)/2
                     Vi(jj)=VC(j,k)
                  end do
               end do
               CALL dsinv(Vi,sum(nmes(i,:)),eps,ier,det)
               if (ier.eq.-1) then
                  do j=1,sum(nmes(i,:))
                     resid_m(nmes_cur+j)=9999.d0
                     resid_ss(nmes_cur+j)=9999.d0
                     do l=1,ng
                        pred_ss_g((l-1)*nobs+nmes_cur+j)=9999.d0
                        pred_m_g((l-1)*nobs+nmes_cur+j)=9999.d0
                     end do
                  end do
               do k=1,nea
                  pred_RE((i-1)*nea+k)=9999.d0
               end do
               do k=1,nalea
                  pred_RE_Y((i-1)*nalea+k)=9999.d0
               end do
                  goto 654
               end if

!     retransformation du vecteur Vi en matrice :
               VC1=0.d0
               do j=1,sum(nmes(i,:))
                  do k=1,sum(nmes(i,:))
                     if (k.ge.j) then
                        VC1(j,k)=Vi(j+k*(k-1)/2)
                     else
                        VC1(j,k)=Vi(k+j*(j-1)/2)
                     end if
                  end do
               end do


               Y2=Y1-mu

               err1=0.d0
               err1=MATMUL(VC1,Y2)
               err2=0.d0
               err2=MATMUL(covUY,err1)
               pred1=0.d0
               pred1=MATMUL(CovDev,err1)


               do j=1,sum(nmes(i,:))
                  pred_m_g((g-1)*nobs+nmes_cur+j)=mu(j)
                  pred_ss_g((g-1)*nobs+nmes_cur+j)=mu(j)+pred1(j)

                  resid_ss(nmes_cur+j)=resid_ss(nmes_cur+j) &
                      +ppi(i,g)*(Y1(j)-mu(j)-pred1(j))
                  resid_m(nmes_cur+j)=resid_m(nmes_cur+j)+pi(g)*(Y2(j))
               end do
               do k=1,nea
                  pred_RE((i-1)*nea+k)=pred_RE((i-1)*nea+k)+ppi(i,g)*err2(k)
               end do
               
           !pred_RE_Y
           err3=0.d0
           do yk=1,nalea
            err3(yk,:) = Valea(sum(nmes(i,1:yk)),:)
           end do
           err4=0.d0
           err4= matmul(err3,err1)         
           do k=1,nalea
              pred_RE_Y((i-1)*nalea+k) = pred_RE_Y((i-1)*nalea+k)+ppi(i,g)*err4(k)
           end do

            end do
            !fin do g

         end if
         !fin if ng

 654  continue

         nmes_cur=nmes_cur+sum(nmes(i,:))
      end do

! FIN BOUCLE SUJET

!write(*,*)"fin boucle sujets","i = ",i

      end subroutine residualsmo



      
      subroutine loglikmultlcmm(Y0,X0,Prior0,idprob0,idea0,idg0,idcor0,idcontr0 &
           ,ny0,ns0,ng0,nv0,nobs0,nea0,nmes0,idiag0,nwg0,ncor0,nalea0&
           ,npm0,b0,ppi0,resid_m &
           ,resid_ss,pred_m_g,pred_ss_g,pred_RE,pred_RE_Y &
           ,epsY0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0 &
           ,nvalSPLORD0,marker,transfY,nsim0,Yobs,Ydiscret,vraisdiscret,UACV,rlindiv&
           ,fix0,nfix0,bfix0,methInteg0,nMC0,dimMC0,seqMC0,chol0,estim0,vrais)

        !use parameters
        use communmo
        !use optim

        IMPLICIT NONE

        !Declaration des variables en entree
        integer,intent(in)::nv0,Ydiscret,ny0,nMC0,methInteg0,dimMC0,chol0,nfix0
        integer, intent(in)::ns0,ng0,nobs0,idiag0,nwg0,nea0,nsim0,ncor0,nalea0,npm0
        double precision,dimension(ny0),intent(in)::epsY0
        integer, dimension(ny0),intent(in)::idlink0,nbzitr0,nvalSPLORD0
        double precision,dimension(maxval(nbzitr0),ny0),intent(in)::zitr0
        integer,dimension(nobs0),intent(in)::indiceY0
        double precision,dimension(sum(nvalSPLORD0(:))),intent(in)::uniqueY0
        integer, dimension(nv0),intent(in)::idea0,idg0,idprob0,idcor0,idcontr0
        integer, dimension(ns0),intent(in)::prior0
        integer,dimension(ns0,ny0)::nmes0   
        double precision,dimension(nobs0),intent(in)::Y0
        double precision,dimension(nobs0*nv0),intent(in)::X0
        integer,dimension(npm0+nfix0),intent(in)::fix0
        double precision,dimension(nfix0),intent(in)::bfix0
        double precision,dimension(dimMC0*nMC0),intent(in)::seqMC0
        integer,intent(in)::estim0 ! estim0=1 pour estimer ; estim0=0 pour faire le postfit

        !Declaration des variable en entree et sortie
        double precision, dimension(npm0), intent(inout) :: b0

        !Declaration des variables en sortie
        double precision,intent(out)::vrais,vraisdiscret,UACV
        double precision,dimension(ns0*ng0),intent(out)::ppi0
        double precision,dimension(ns0),intent(out)::rlindiv
        double precision,dimension(nobs0),intent(out)::resid_m,resid_ss,Yobs
        double precision,dimension(nobs0*ng0),intent(out)::pred_m_g
        double precision,dimension(nobs0*ng0),intent(out)::pred_ss_g
        double precision,dimension(ns0*nea0),intent(out)::pred_RE ! commun sur proc latent
        double precision,dimension(ns0*nalea0),intent(out)::pred_RE_Y ! effets specifiques
        double precision,dimension(nsim0*ny0),intent(out)::marker,transfY 

        !Variables locales
        integer::jtemp,i,g,j,ij,ier,k,ktemp,ig,id,yk,k1,k2,mi,nbfix,jd,npmtot0
        double precision::thi,thj
        double precision,dimension(npm0+nfix0)::btot
        double precision,dimension(ns0,ng0)::PPI
        double precision,external::vrais_multo


        !      write(*,*)'indice entres',indiceY0

        !          print*,"Y0=",Y0(1:10)
        !          print*,"X0=",X0(nobs0-1:nobs+12)
        !          print*,"zitr0=",zitr0

        ! sorties initialisees

        ppi0=0.d0
        pred_ss_g=0.d0
        pred_m_g=0.d0
        pred_RE=0.d0
        pred_RE_Y=0.d0
        marker=0.d0
        transfY=0.d0
        resid_m=0.d0
        resid_ss=0.d0
        vrais=0.d0
        Yobs=0.d0
        vraisdiscret=0.d0
        UACV=0.d0
        rlindiv=0.d0



        ! fin en prevision

        maxmes=0
        do i=1,ns0
           mi=sum(nmes0(i,:))
           !         mi=0 
           !         do k=1,ny0
           ! mi = mi+nmes0(i,k)
           !         end do
           if (mi.gt.maxmes) then
              maxmes=mi
           end if
        end do


        allocate(rangeY(ny0),minY(ny0),maxY(ny0),idlink(ny0),ntrtot(ny0),epsY(ny0))



        nySPL=0
        nyORD=0
        rangeY=0
        epsY=epsY0
        do k=1,ny0
           idlink(k)=idlink0(k)
           minY(k)=zitr0(1,k)
           maxY(k)=zitr0(nbzitr0(k),k)
           if (Ydiscret.eq.1) rangeY(k)=maxY(k)-minY(k)
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

        methInteg = methInteg0
        nMC = nMC0

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
             ,idea(nv0),idg(nv0),idcor(nv0),idcontr(nv0),nmes(ns0,ny0),prior(ns0) &
             ,indiceY(nobs0))

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
                    prior(i)=prior0(i)
                    nmes(i,yk)=nmes0(i,yk)   !dim(nmes)=ns*ny    
                    do j=1,nmes(i,yk)
                       jtemp=jtemp+1
                       Y(jtemp)=Y0(jtemp)
                       indiceY(jtemp)=indiceY0(jtemp)
                       ktemp=ktemp+1
                       X(jtemp,k)=X0(ktemp)
                    end do
                 else
                    do j=1,nmes(i,yk)
                       ktemp=ktemp+1
                       jtemp=jtemp+1
                       X(jtemp,k)=X0(ktemp)
                    end do
                 end if
              end do
           end do
        end do
        !         write(*,*)'X k:',X(1:50,k)


        ! prm fixes
        npmtot0=npm0+nfix0
        allocate(fix(npmtot0))
        fix=0
        fix(1:npmtot0)=fix0(1:npmtot0)
        nbfix=sum(fix)
        if(nbfix.eq.0) then
           allocate(bfix(1))
        else
           allocate(bfix(nbfix))
        end if
        bfix(1:nbfix)=bfix0(1:nbfix)



        ! creation des parametres

        nea=0
        ncg=0
        ncssg=0
        !nprob=ng-1 (a remettre si on remet idprob[1]=0)
        nprob=0
        !nvarprob=min(ng-1,1)  (a remettre si on remet idprob[1]=0)
        nvarprob=0
        ncontr=0
        do k=1,nv
           if (idg(k).eq.1) then
              ncssg=ncssg+1      ! nb var. sans melange
           else if (idg(k).eq.2) then
              ncg=ncg+1      ! nb var. dans melange
           end if
           nea=nea+idea(k)
           nprob=nprob+(idprob(k))*(ng-1)
           nvarprob=nvarprob+idprob(k)
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
        
        !if (nwg.gt.0) then
        !   do i=1,nwg
        !      btot(nprob+nef+ncontr+nvc+i)=abs(btot(nprob+nef+ncontr+nvc+i))
        !   end do
        !end if !! le faire dans R !**


        ! points qmc
        if(methInteg.ne.3) then 
           allocate(seqMC(1))
        else
           allocate(seqMC(dimMC0*nMC))
           seqMC = seqMC0(1:dimMC0*nMC) 
        end if


        ! base de splines
        if (any(idlink.eq.2)) then 
           call design_splines_multo(ier)
           if (ier.eq.-1) then
              go to 1589
           end if
        end if

        !print*,"npmtot=",npmtot, " npm=",npm
        !print*,"b=",b, " bfix=",bfix
        ! lancement de l'optimisation

        if((ng.eq.1.and.ncg.gt.0).or.(ng.eq.1.and.nprob.gt.0)) then
           !      if(verbose==1) write(*,*)"ng",ng,"ncg",ncg,"nprob",nprob
           go to 1589
        end if


        if (estim0.eq.1) then

           IF (npm0.eq.1) then
              go to 1589
           else

              id=0
              jd=0
              thi=0.d0
              thj=0.d0
              
              vrais=vrais_multo(b0,npm0,id,thi,jd,thj)
              
           end if

        else
              
           !  injecter le b estime dans btot
           btot=0.d0
           k=0
           k2=0
           do j=1,npmtot
              if(fix0(j).eq.0) then
                 k=k+1
                 btot(j)=b0(k)
              else
                 k2=k2+1
                 btot(j)=bfix0(k2)
              end if
           end do
           
              !if (verbose==1) write(*,*)'avant transfo'
              call transfos_estimees_2(btot,npmtot,nsim0,marker,transfY)
              ! end if

              !if (istop.eq.1) then
              if (ng.gt.1) then
                 !                if(verbose==1)  write(*,*)'avant postprob'
                 call postprobmo(btot,npmtot,PPI)
              end if



              !            if(verbose==1) write(*,*)'avant residuals'

              call residualsmo(btot,npmtot,ppi,resid_m,pred_m_g,resid_ss &
                   ,pred_ss_g,pred_RE,pred_RE_Y,Yobs)


              ig=0
              ij=0
              do i=1,ns
                 do g=1,ng0
                    ig=ig+1
                    ppi0(ig)=PPI(i,g)
                 end do
              end do
              
           end if



        !write(*,*)'avant deallocate'


1589    continue

        deallocate(Y,X,idprob,idea,idg,idcor,idcontr,nmes,prior,uniqueY,indiceY,ntrtot)


        deallocate(zitr,mm,mm1,mm2,im,im1,im2,minY,maxY,rangeY,idlink,nvalSPL,nvalORD,epsY)

        deallocate(fix,bfix,seqMC)

        !write(*,*)'fin'
        return
      end subroutine loglikmultlcmm






