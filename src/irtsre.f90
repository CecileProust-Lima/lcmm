
module modirtsre

  implicit none
  integer,save ::ny,ns,nv,idiag,nvc,nea,ncor &
       ,maxmes,nobs,nef,ncontr,ntrtot,nalea &
       ,nySPL,ntotvalSPL,nyORD,ntotvalORD,npmtot &
       ,nMC,methInteg,nmescur &
       ,nvarxevt,nbevt,logspecif,idtrunc,nvdepsurv,nrisqtot,nxevt,nasso &
       ,expectancy
  integer,dimension(:),allocatable,save::typrisq,nz,nprisq,nevtparx,nxcurr
  double precision,dimension(:),allocatable,save::Y,uniqueY,minY,maxY,rangeY
  double precision,dimension(:,:),allocatable,save ::X,zi
  double precision,dimension(:),allocatable,save::Tsurv0,Tsurv,Tsurvint
  integer,dimension(:),allocatable,save::Devt,ind_survint
  integer,dimension(:),allocatable,save::idtdv,idsurv
  integer,dimension(:),allocatable,save::idea,idg,idcor,idcontr,indiceY
  integer,dimension(:),allocatable,save::idlink,ntr
  integer,dimension(:,:),allocatable,save::nmes
  integer,dimension(:),allocatable,save::nvalSPL,nvalORD
  double precision,dimension(:),allocatable,save :: seqMC
  double precision,dimension(:),allocatable,save :: epsY
  double precision,dimension(:,:),allocatable,save::zitr
  double precision,dimension(:),allocatable,save::mm,mm1,mm2,im,im1,im2
  double precision,dimension(:),allocatable,save::Tmm,Tmm1,&
       Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,Tmm01,Tmm02,Tmm03,Tim0, &
       Tim01,Tim02,Tim03,Tmmt,Tmmt1,Tmmt2,Tmmt3,Timt,Timt1,&
       Timt2,Timt3
!  double precision,dimension(:),allocatable,save::Tmm_est,Tmm1_est &
!       ,Tmm2_est,Tmm3_est,Tim_est,Tim1_est,Tim2_est,Tim3_est
  double precision,dimension(:),allocatable,save::brisq_est
  !  double precision,save::vrais_surv
  integer,dimension(:),allocatable,save::fix
  double precision,dimension(:),allocatable,save::bfix

end module modirtsre




subroutine irtsre(Y0,X0,Tentr0,Tevt0,Devt0,ind_survint0 &
     ,idea0,idg0,idcor0,idcontr0,idsurv0,idtdv0 &
     ,typrisq0,nz0,zi0,nbevt0,idtrunc0,logspecif0 &
     ,ny0,ns0,nv0,nobs0,nea0,nmes0,idiag0,ncor0,nalea0&
     ,npmtot0,btot,Vopt,vrais,ni,istop,gconv,resid_m &
     ,resid_ss,pred_RE,pred_RE_Y,convB,convL,convG &
     ,maxiter0,epsY0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0 &
     ,nvalSPLORD0,time,risq_est,risqcum_est,marker,transfY,nsim0,Yobs &
     ,rlindiv,pbH0,fix0,methInteg0,nMC0,dimMC0,seqMC0)

  use parameters
  use modirtsre
  use optim

  IMPLICIT NONE

  !Declaration des variables en entree
  integer,intent(in)::nv0,maxiter0,ny0,nMC0,methInteg0,dimMC0
  integer, intent(in)::ns0,nobs0,idiag0,npmtot0,nea0,nsim0,ncor0,nalea0
  integer,intent(in)::idtrunc0,logspecif0,nbevt0
  double precision, dimension(ns0),intent(in)::Tentr0,Tevt0
  integer, dimension(ns0),intent(in)::ind_survint0,Devt0
  integer, dimension(nv0),intent(in)::idtdv0,idsurv0
  integer,dimension(nbevt0),intent(in)::typrisq0,nz0
  double precision,dimension(maxval(nz0),nbevt0),intent(in)::zi0    
  double precision,dimension(ny0),intent(in)::epsY0
  integer, dimension(ny0),intent(in)::idlink0,nbzitr0,nvalSPLORD0
  double precision,dimension(maxval(nbzitr0),ny0),intent(in)::zitr0
  integer,dimension(nobs0),intent(in)::indiceY0
  double precision,dimension(sum(nvalSPLORD0(:))),intent(in)::uniqueY0
  integer, dimension(nv0),intent(in)::idea0,idg0,idcor0,idcontr0
  integer,dimension(ns0,ny0)::nmes0   
  double precision,dimension(nobs0),intent(in)::Y0
  double precision,dimension(nobs0*nv0),intent(in)::X0
  double precision,intent(in)::convB,convL,convG
  integer,dimension(npmtot0),intent(in)::pbH0,fix0
  double precision,dimension(dimMC0*nMC0),intent(in)::seqMC0

  !Declaration des variable en entree et sortie
  double precision, dimension(npmtot0), intent(inout) :: btot

  !Declaration des variables en sortie
  double precision,intent(out)::vrais
  double precision,dimension(3),intent(out)::gconv
  double precision,dimension(ns0),intent(out)::rlindiv
  double precision,dimension(nobs0),intent(out)::resid_m,resid_ss,Yobs
  double precision,dimension(ns0*nea0),intent(out)::pred_RE ! commun sur proc latent
  double precision,dimension(ns0*nalea0),intent(out)::pred_RE_Y ! effets specifiques
  double precision,dimension(nsim0*ny0),intent(out)::marker,transfY
  double precision,dimension(nsim0),intent(out)::time
  double precision,dimension(nsim0,nbevt0),intent(out)::risq_est,risqcum_est
  double precision,dimension(npmtot0*(npmtot0+1)/2),intent(out)::Vopt
  integer, intent(out)::ni,istop

  !Variables locales
  integer::jtemp,i,j,npm,ier,k,ktemp,yk,k1,k2,mi,nbfix
  integer::ke,sumnrisq,it
  double precision::eps,ca,cb,dd
  double precision,dimension(npmtot0)::mvc,b
  double precision,dimension(npmtot0*(npmtot0+3)/2)::V
  double precision,external::vrais_irtsre

  !      write(*,*)'indice entres',indiceY0

  !          print*,"Y0=",Y0(1:10)
  !          print*,"X0=",X0(nobs0-1:nobs+12)
  !          print*,"zitr0=",zitr0

  ! sorties initialisees

  Vopt=0.d0
  gconv=0.d0
  pred_RE=0.d0
  pred_RE_Y=0.d0
  marker=0.d0
  transfY=0.d0
  resid_m=0.d0
  resid_ss=0.d0
  vrais=0.d0
  ni=0
  Yobs=0.d0
  rlindiv=0.d0


  maxmes=0
  do i=1,ns0
     mi=sum(nmes0(i,:))
     if (mi.gt.maxmes) then
        maxmes=mi
     end if
  end do

  epsa=convB
  epsb=convL
  epsd=convG
  maxiter=maxiter0

  allocate(rangeY(ny0),minY(ny0),maxY(ny0),idlink(ny0),ntr(ny0),epsY(ny0))

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

  methInteg = methInteg0
  nMC = nMC0

!print*,"ok2"," nvalsplord=",nvalsplord0, " nvalord=",nvalord
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
     if (idlink(k).eq.0) ntr(k)=2     
     if (idlink(k).eq.1) ntr(k)=4
     if (idlink(k).eq.2) then
        k1=k1+1
        ntr(k)=nbzitr0(k)+2

        zitr(1:nbzitr0(k),k1)=zitr0(1:nbzitr0(k),k)
        zitr(-1,k1)=zitr(1,k1)
        zitr(0,k1)=zitr(1,k1)
        zitr(ntr(k)-1,k1)=zitr(ntr(k)-2,k1)
        zitr(ntr(k),k1)=zitr(ntr(k)-1,k1)
     end if
     if (idlink(k).eq.3) then
        k2 = k2+1
        ntr(k) = nvalORD(k2)-1
     end if
  end do

  ntrtot = sum(ntr)

  !if (verbose==1)       write(*,*)'zitr',zitr


  allocate(Y(nobs0),X(nobs0,nv0),uniqueY(ntotvalSPL+ntotvalORD) &
       ,idea(nv0),idg(nv0),idcor(nv0),idcontr(nv0),nmes(ns0,ny0) &
       ,indiceY(nobs0))

  allocate(Tsurv0(ns0),Tsurv(ns0),Tsurvint(ns0),ind_survint(ns0),Devt(ns0))
  allocate(typrisq(nbevt0),nz(nbevt0),nprisq(nbevt0),nevtparx(nv0),nxcurr(nv0))
  allocate(idsurv(nv0),idtdv(nv0))

  ! zi : contient noeuds pour hazard (ou min et max si Weibull)
  if(any(typrisq0.eq.3)) then
     allocate(zi(-2:maxval(nz0)+3,nbevt0))
  else
     allocate(zi(maxval(nz0),nbevt0))
  end if


  eps=1.d-20

  ! enregistrement pour les modules
  nbevt=nbevt0    
  typrisq=typrisq0
  idtrunc=idtrunc0
  Tsurv0=Tentr0   
  Tsurv=Tevt0    
  devt=devt0    
  ind_survint=ind_survint0
  logspecif=logspecif0 
  Tsurvint=Tsurv
  ny=ny0
  ns=ns0
  nv=nv0
  nobs=nobs0
  ncor=ncor0
  nalea=nalea0
  idiag=idiag0

  !     if (verbose==1) write(*,*)'ntotvalSPL',ntotvalSPL

  if (ntotvalSPL+ntotvalORD.gt.0) uniqueY(1:ntotvalSPL+ntotvalORD)=uniqueY0(1:ntotvalSPL+ntotvalORD)

  nmes=0
  Y=0.d0
  X=0.d0
  idea=0
  idg=0
  idcor=0
  idcontr=0
  idsurv=0
  idtdv=0
  ktemp=0

  do k=1,nv
     idsurv(k)=idsurv0(k)
     idtdv(k)=idtdv0(k)
     idea(k)=idea0(k)
     idg(k)=idg0(k)
     idcor(k)=idcor0(k)
     idcontr(k)=idcontr0(k)

     jtemp=0
     DO i=1,ns
        do yk=1,ny            
           if (k.eq.1) then
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
          ! write(*,*)'X k:',X(1:50,2)

  ! definition Tsurvint 
  nvdepsurv=0
  if(sum(ind_survint).gt.0) then
     nvdepsurv=1
     do k=1,nv
        if (idtdv(k).eq.1) then
           it=0
           do i=1,ns
              Tsurvint(i)=X(it+1,k)
              it=it+maxval(nmes(i,:))
           end do
        end if
     end do
  end if


  ! prm fixes
  allocate(fix(npmtot0))
  fix=0
  fix(1:npmtot0)=fix0(1:npmtot0)
  nbfix=sum(fix)
  if(nbfix.eq.0) then
     allocate(bfix(1))
  else
     allocate(bfix(nbfix))
  end if
  bfix=0.d0



  ! creation des parametres

  nprisq=0
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

     nrisqtot = nrisqtot+nprisq(ke)  ! nb total de prm pour hazards
     zi(1:nz(ke),ke)=zi0(1:nz(ke),ke)
  end do

  ! nvarxevt = nombre total de coef pour survie (sans prm hazard)
  nxevt=0
  nevtparx=0
  do j=1,nv

     if(idtdv(j).ne.1) then

        if(idsurv(j).eq.1) then
           nevtparx(j) = 1
           nxevt = nxevt + 1
        end if
        if(idsurv(j).eq.2) then 
           nevtparx(j) = nbevt
           nxevt = nxevt + 1
        end if
     end if

  end do

  nvarxevt = sum(nevtparx) + nvdepsurv

  nea=0
  nef=0
  ncontr=0
  do k=1,nv
     if (idg(k).eq.1) then
        nef = nef + 1
     end if
     nea=nea+idea(k)
     ncontr=ncontr+idcontr(k)*(ny-1) 
  end do
  nef = nef - 1 !intercept pas estime

  nasso = nbevt*nea


  if (idiag.eq.1) then
     nvc=nea-1
  else if(idiag.eq.0) then
     nvc=(nea+1)*nea/2-1
  end if

  npmtot = nrisqtot+nvarxevt+nef+ncontr+nvc+ncor+ntrtot+nalea+ny+nasso
  !print*,nrisqtot,nvarxevt,nef,ncontr,nvc,ncor,ntrtot,nalea,ny,nasso
  !print*,"npmtot=",npmtot
  !print*,"btot=",btot
  ! if (verbose==1 ) write(*,*)"npm0=",npm0,"npm=",npm

  !  write(*,*)'idlink',idlink
  !  write(*,*)'idea',idea
  !  write(*,*)'NVC',nvc


  if (idiag.eq.1) then
     DO j=1,nvc
        btot(nrisqtot+nvarxevt+nasso+nef+ncontr+j)=dsqrt(abs(btot(nrisqtot+nvarxevt+nasso+nef+ncontr+j)))
     END DO
  end if

  ! si idiag=0, on met dans le vecteur des parms, les parms
  ! de la transformee de Cholesky

  if (idiag.eq.0) then

     mvc(1)=1.d0
     DO j=1,nvc
        mvc(1+j)=btot(nrisqtot+nvarxevt+nasso+nef+ncontr+j)
     END DO

     CALL dmfsd(mvc,nea,EPS,IER)
     DO j=1,nvc
        btot(nrisqtot+nvarxevt+nasso+nef+ncontr+j)=mvc(1+j)
     END DO
  end if


  ! points qmc
  if(methInteg.ne.3) then 
     allocate(seqMC(1))
  else
     allocate(seqMC(dimMC0*nMC))
     seqMC = seqMC0(1:dimMC0*nMC) 
  end if

  ! creation du vecteur b avec slt les prm a estimer
  b=0.d0
  npm=0
  k=0
  do j=1,npmtot
     if(fix0(j).eq.0) then
        npm=npm+1
        b(npm)=btot(j)
     end if
     if(fix0(j).eq.1) then
        k=k+1
        bfix(k)=btot(j)
     end if
  end do

  allocate(pbH(npm))
  k=0
  do j=1,npmtot
     if(fix0(j).eq.0) then
        k=k+1
        pbH(k)=pbH0(j)
     end if
  end do



  ! creer base de splines si au moins un hazard splines
  if(any(typrisq.eq.3)) then
     allocate(Tmm(ns*nbevt),Tmm1(ns*nbevt),Tmm2(ns*nbevt),Tmm3(ns*nbevt),Tim(ns*nbevt)         &
          ,Tim1(ns*nbevt),Tim2(ns*nbevt),Tim3(ns*nbevt),Tmm0(ns*nbevt),Tmm01(ns*nbevt),Tmm02(ns*nbevt)  &
          ,Tmm03(ns*nbevt),Tim0(ns*nbevt),Tim01(ns*nbevt),Tim02(ns*nbevt),Tim03(ns*nbevt),          &
          Tmmt(ns*nbevt),Tmmt1(ns*nbevt),Tmmt2(ns*nbevt),Tmmt3(ns*nbevt),Timt(ns*nbevt),Timt1(ns*nbevt) &
          ,Timt2(ns*nbevt),Timt3(ns*nbevt))

     do ke=1,nbevt
        if(typrisq(ke).eq.3) then
           call splines_irtsre(ke)
        end if
     end do

  end if



  ! base de splines transfos
  if (any(idlink.eq.2)) then 
     call design_splines_irtsre(ier)
     if (ier.eq.-1) then
        istop=9
        go to 1589
     end if
  end if



  ! indicateur vraisemblance ou especrance de vie
  expectancy = 0 ! on va calculer une vraisemblance 

  !print*,"npmtot=",npmtot, " npm=",npm
  !print*,"b=",b, " bfix=",bfix
  ! lancement de l'optimisation


  IF (npm.eq.1) then
     istop=10
     go to 1589
  else
     ca=0.d0
     cb=0.d0
     dd=0.d0
     !      write(*,*)"before optimisation",npm,b
     !         if (verbose==1) write(*,*)"before optimisation",npm,b
     call marq98(b,npm,ni,V,vrais,ier,istop,ca,cb,dd,vrais_irtsre)

     !         if (verbose==1) write(*,*)"after optimisation",npm,b
     !         write(*,*)
     !      write(*,*)'    FIN OPTIMISATION  ..... '
     ! if(verbose==1) write(*,*)'istop',istop,'vrais',vrais

     gconv=0.d0
     gconv(1)=ca
     gconv(2)=cb
     gconv(3)=dd
     vopt(1:(npm*(npm+1)/2))=V(1:(npm*(npm+1)/2))

     !  injecter le b estime dans btot
     k=0
     do j=1,npmtot
        if(fix0(j).eq.0) then
           k=k+1
           btot(j)=b(k)
        end if
     end do


     if (istop.eq.1.or.istop.eq.2.or.istop.eq.3) then  
        !postfit

        ! grille de valeur pour les risques de base
        allocate(brisq_est(maxval(nprisq)))
        brisq_est=0.d0

        sumnrisq=0
        do ke=1,nbevt
           
           brisq_est=0.d0

           if (logspecif.eq.1) then
              do k=1,nprisq(ke)
                 brisq_est(k)=exp(btot(sumnrisq+k))
              end do
           else
              do k=1,nprisq(ke)
                 brisq_est(k)=btot(sumnrisq+k)*btot(sumnrisq+k)
              end do
           end if
              
           call fct_risq_estime_irtsre(ke,brisq_est,time,nsim0,risq_est,risqcum_est)

           sumnrisq = sumnrisq + nprisq(ke)
        end do

        deallocate(brisq_est)


        ! grille de valeurs pour la transfo
        call transfos_estimees_irtsre(btot,npmtot,nsim0,marker,transfY)

     end if

  end if


  !write(*,*)'avant deallocate'


1589 continue

  if (any(typrisq.eq.3)) then
     deallocate(Tmm,Tmm1,Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,    &
          Tmm01,Tmm02,Tmm03,Tim0,Tim01,Tim02,Tim03,Tmmt,Tmmt1,     &
          Tmmt2,Tmmt3,Timt,Timt1,Timt2,Timt3)
  endif

  deallocate(Tsurv0,Tsurv,Tsurvint &
       ,ind_survint,zi,devt,typrisq,nz,nprisq,idsurv,idtdv &
       ,nevtparx,nxcurr)

  deallocate(Y,X,idea,idg,idcor,idcontr,nmes,uniqueY,indiceY,ntr)


  deallocate(zitr,mm,mm1,mm2,im,im1,im2,minY,maxY,rangeY,idlink,nvalSPL,nvalORD,epsY)

  deallocate(fix,bfix,pbH,seqMC)


  !write(*,*)'fin'
  return
end subroutine irtsre






!-----------------------------------------------------------
!                       VRAIS_IRTSRE_i
!------------------------------------------------------------


double precision function vrais_irtsre_i(b,npm,id,thi,jd,thj,i) 

  use parameters
  use modirtsre
  use optim

  IMPLICIT NONE
  integer ::i,j,k,l,m,id,jd,jj,npm,ll,ii,numSPL,ykord
  integer ::ier,kk,j1,j2,sumMesYk,yk,sumntr,ke,sumnrisq
  integer::nevtxcurr,nxevtcurr
  double precision,dimension(maxmes,nv) ::X00
  double precision,dimension(maxmes,nea) ::Z
  double precision,dimension(maxmes,(ncontr+sum(idcontr)))::X01
  double precision,dimension(ncontr+sum(idcontr))::b01
  double precision,dimension(nea,nea) ::Ut
  double precision,dimension(maxmes,maxmes) ::VC,Corr
  double precision,dimension(npm) :: b
  double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
  double precision,dimension(nv) :: b0
  double precision,dimension(npmtot)::b1
  double precision,dimension(nxevt)::Xevt,bevt
  double precision,dimension(nbevt)::bevtint
  double precision,dimension(maxval(nprisq))::brisq

  double precision :: eps,det,som,thi,thj,eta0
  double precision ::Y4,jacobien,beta_densite,ytemp
  double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
  double precision,dimension(-1:maxval(ntr)-3)::splaa
  double precision::aa1,bb1,dd1,aa,bb,betai,cc1
  double precision,dimension(nea)::ui,usim
  double precision,dimension(maxmes)::wi,wsim
  double precision,dimension(ny)::asim
  double precision::ai,binf,bsup
  double precision,dimension(nbevt)::risq,surv,surv0,survint
  double precision::SX,x22,div,vrais_Y,vrais_surv,varexpsurv
  double precision::surv0_glob,surv_glob,fevt,easurv,trunc
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
                 Ut(j,k)=b1(nrisqtot+nvarxevt+nasso+nef+ncontr+j-1)
              else
                 Ut(j,k)=0.d0
              end if
           end do
        end do
     end if

     If (idiag.eq.0) then
        do j=2,nea
           do k=1,j
              Ut(j,k)=b1(nrisqtot+nvarxevt+nasso+nef+ncontr+k-1+j*(j-1)/2)
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
              Corr(j1,j2) = Corr(j1,j2)+b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)* &
                   b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)*min(tcor(j1),tcor(j2))
           else if (ncor.eq.2) then
              Corr(j1,j2) = Corr(j1,j2)+b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)* &
                   b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor)* &
                   exp(-b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+1)*abs(tcor(j1)-tcor(j2)))
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
  sumntr=0
  numSPL=0
  do yk=1,ny

     if (idlink(yk).eq.0) then  ! Linear link

        do j=1,nmes(i,yk)
           Y1(sumMesYk+j)=(dble(Y(nmescur+sumMesYk+j))-b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)) &
                /abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2))

           jacobien = jacobien - log(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2))
        end do

     else if (idlink(yk).eq.1) then  ! Beta link


        aa1=exp(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1))/ &
             (1+exp(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)))
        bb1=exp(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2))/ &
             (1+exp(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+2)))
        bb1=aa1*(1.d0-aa1)*bb1

        cc1=abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+3))

        dd1=abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+4))

        aa=aa1*aa1*(1-aa1)/bb1-aa1
        bb=aa*(1-aa1)/aa1

        do j=1,nmes(i,yk)

           ytemp=(dble(Y(nmescur+sumMesYk+j))-minY(yk)+epsY(yk))/(maxY(yk)-minY(yk)+2*epsY(yk))
           Y1(sumMesYk+j)=(betai(aa,bb,ytemp)-cc1)/dd1


           if (Y1(sumMesYk+j).eq.999.d0) then
              vrais_irtsre_i=-1.d9
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
        eta0=b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)

        do kk=2,ntr(yk)
           splaa(kk-3)=b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+kk)&
                *b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+kk)
        end do
        !if(i==1 .and. id==0 .and. jd==0) print*,"eta0=",eta0,"splaa=",sqrt(splaa)
        do j=1,nmes(i,yk)
           ll=0
           !if(i==1 .and. id==0 .and. jd==0) print*,"Y=",Y(nmescur+sumMesYk+j)
           if (Y(nmescur+sumMesYk+j).eq.zitr(ntr(yk)-2,numSPL)) then
              ll=ntr(yk)-3
           end if

           som=0.d0
           do kk = 2,ntr(yk)-2
              if ((Y(nmescur+sumMesYk+j).ge.zitr(kk-1,numSPL)).and. &
                   (Y(nmescur+sumMesYk+j).lt.zitr(kk,numSPL))) then
                 ll=kk-1
              end if
           end do

           if (ll.lt.1.or.ll.gt.ntr(yk)-3) then          
              vrais_irtsre_i=-1.d9
              print*,"-1.d9 ll<1 ou ll>ntr-3",ll!," ntr=",ntr(yk)," numSPL=",numSPL," y=",Y(nmescur+sumMesYk+j)
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

           !print*,"jac=",jacobien
           !print*,"ll =",ll
           !print*,"splaa =",splaa(ll-2)
           !print*,"nmescur+sumMesYk+j =",nmescur+sumMesYk+j
           !print*,"indiceY= ",indiceY(nmescur+sumMesYk+j)
           !print*,"mm =",mm(nmescur+sumMesYk+j)
           !print*,"mm1 =",mm1(nmescur+sumMesYk+j)
           !print*,"mm2 =",mm2(nmescur+sumMesYk+j)
           !print*,"Y =",Y(nmescur+sumMesYk+j)

           !                write(*,*)'Y',Y1(sumMesYk+j),sumMesYk,yk,j,jacobien
        end do
     else if (idlink(yk).eq.3) then
        do j=1,nmes(i,yk)
           Y1(sumMesYk+j)=Y(nmescur+sumMesYk+j)
           !if(nmescur.lt.15) then
           !   print*,"nmescur=",nmescur," sumMesYk=",sumMesYk," j=",j
           !   print*,"Y=",Y(nmescur+sumMesYk+j), "  Y1=",Y1(sumMesYk+j)
           !end if
        end do
     end if
     sumMesYk=sumMesYk+nmes(i,yk)
     sumntr=sumntr+ntr(yk)
  end do !fin boucle yk


  !         if (i.lt.3)then
  !            write(*,*)'nmes',nmes(i),b1((nef+ncontr+nvc+1):npm),nef+ncontr
  !            write(*,*)'Y1',Y1
  !         end if




  ! contribution individuelle a la vraisemblance
  ! print*,"i=",i," -ni*log(2pi)=",-sum(nmes(i,:))*dlog(dble(2*3.14159265)), " log(det)=",det
  ! print*,"Vi=",VC
  ! sans classes latentes : ng=1
  vrais_irtsre_i=0.d0


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
           b0(l)=b1(nrisqtot+nvarxevt+nasso+l-1)
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
              b01((m-1)*ny+yk)=b1(nrisqtot+nvarxevt+nasso+nef+(m-1)*(ny-1)+yk)
           else
              b01((m-1)*ny+ny) =-sum(b1(nrisqtot+nvarxevt+nasso+nef+(m-1)*(ny-1)+1 &
                   :nef+(m-1)*(ny-1)+ny-1))
           end if
        end do
     end if
  end do


  trunc=0
  som=0.d0
  do l=1,nMC

     vrais_Y=1.d0
     mu=0.d0

!!!!!!!!! MC pour EA et BM/AR !!!!!!!!!

     if(methInteg.eq.1) then 
        ! !!!!!!!!!!!!! MCO !!!!!!!!!!!!!

        ! simuler les effets aleatoires
        if(nea.gt.0) then
           x22=0.d0
           SX=1.d0
           do j=1,nea
              call bgos(SX,0,usim(j),x22,0.d0)
              !print*,"usim=",usim(j)
           end do
           ui=0.d0
           ui=matmul(Ut,usim)
           !print*,"usim=",usim(j)," ui=",ui, " Ut=",Ut, "  nea=",nea
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
           ui = -ui
           wi = -wi
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
     !if(i.lt.4) then
     !print*,"i=",i," nmes=",nmes(i,1)," nmescur=",nmescur
     !print*,"Xb=",matmul(X00,b0)+matmul(X01,b01)
     !print*,"mu=",mu
     !print*,"ui=",ui, " wi=",wi," xb=",matmul(X00,b0)+matmul(X01,b01)
     !end if
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
                 ai = b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk)*asim(yk)
              else if(methInteg.eq.2) then
                 !! MCA
                 if(mod(l,2).eq.0) then
                    ! si l est pair on prend l'oppose du precedent
                    ai = -b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk)*asim(yk)
                 else
                    call bgos(SX,0,asim(yk),x22,0.d0)
                    ai = b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk)*asim(yk)
                 end if
              else
                 !! QMC
                 asim(yk) = seqMC(nMC*(nea+sum(nmes(i,:)))+l)
                 ai = b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk)*asim(yk)
              end if
           end if
!if(i.lt.4) print*,"i=",i," avant do j, vrais_Y=",vrais_Y
           do j=1,nmes(i,yk)

              !! on ajoute ai a mu
              if(nalea.gt.0) mu(sumMesYk+j) = mu(sumMesYk+j)+ai

              !! trouver binf et bsup tq binf < lambda + epsilon < bsup

              !! initialiser au premier seuil
              binf = b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+1)
              bsup = binf

              !! si Y>minY ajouter b1(..)^2
              if(indiceY(nmescur+sumMesYk+j).gt.1) then
                 do ll=2,min(indiceY(nmescur+sumMesYk+j),ntr(yk))
                    bsup = bsup + b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+ll)**2
                    if(ll.lt.indiceY(nmescur+sumMesYk+j)) then
                       binf = binf + b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+sumntr+ll)**2
                    end if
                 end do
              end if
              !if(i.lt.4)print*,"y=",Y1(sumMesYk+j)," indiceY=",indiceY(nmescur+sumMesYk+j), " Y=",Y(nmescur+sumMesYk+j)," mu=",mu(sumMesYk+j)
              !print*," binf=",binf," bsup=",bsup
              !! centrer et standardiser
              binf = (binf - mu(sumMesYk+j))/abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+yk))
              bsup = (bsup - mu(sumMesYk+j))/abs(b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+yk))

              !if(i.lt.4) print*,"nmescur=", nmescur, "sumMesYk=", sumMesYk,"j=",j
              !if(i.lt.4) print*," nvalORD=",nvalORD(ykord)," binf=",binf," bsup=",bsup
              
              if(indiceY(nmescur+sumMesYk+j).eq.1) then
                 !! si Y=minY
                 vrais_Y = vrais_Y * alnorm(binf,.false.)
                 !if(expectancy.eq.1) then
                 !   print*,"Y=minY, P(Y=y)=",alnorm(binf,.false.)
                 !   print*,"t=",X00(sumMesYk+j,2), " mu=", mu(sumMesYk+j)
                 !   print*,"binf=", binf
                 !end if
                 ! P(y=ymin) = P(y<=ymin) donc pareil si vrais ou expect
              else if(indiceY(nmescur+sumMesYk+j).eq.nvalORD(ykord)) then
                 !! si Y=maxY
                 if(expectancy.eq.0) then
                    vrais_Y = vrais_Y * (1.d0-alnorm(bsup,.false.))
                 end if
                 ! si expectancy, P(y<=ymax)=1 donc on n a rien a calculer 

              else
                 !! minY < Y < maxY
                 if(expectancy.eq.0) then
                    vrais_Y = vrais_Y * (alnorm(bsup,.false.)-alnorm(binf,.false.))
!                    print*,"<Y<"
                 else
                    vrais_Y = vrais_Y * alnorm(bsup,.false.)
                 end if

              end if
              !if(expectancy.eq.1) then
              !   print*, "indiceY=",indiceY(nmescur+sumMesYk+j)
              !   print*,"bsup=", bsup, " binf=", binf
              !   print*,"vrais_Y=",vrais_Y
              !end if
           end do


        else
           !! yk est continu
           !print*,"MC continu"
           if(nmes(i,yk).gt.0) then
              !! variance de Y|ui,wi
              VC=0.d0
              do j1=1,nmes(i,yk)
                 VC(j1,j1) = b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+yk)**2 !variance de l'erreur yk
                 if (nalea.eq.ny) then ! intercept aleatoire de yk
                    do j2=1,nmes(i,yk)
                       VC(j1,j2) = VC(j1,j2) + &
                            b1(nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+yk)**2
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
                 vrais_irtsre_i=-1.d9
                 print*,"-1.d9 dsinv continu MC"
                 !print*,"b=",b
                 !print*,"bfix=",bfix
                 !print*,"id=",id,"thi=",thi
                 !print*,"jd=",jd,"thj=",thj
                 !print*,"fix=",fix
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
              
              div = (dble(2*3.14159265)**(dble(nmes(i,yk)/2)))*sqrt(exp(det))
              
              vrais_Y = vrais_Y * exp(-Y4/2.d0)/div

           end if
        end if
        
        sumMesYk = sumMesYk+nmes(i,yk)
        sumntr = sumntr + ntr(yk)
     end do ! fin boucle yk

!if(i.lt.4) print*,"avant survie, vrais_Y",vrais_Y
     ! partie survie

     if (nbevt.ne.0) then

        ! calcul de brisq en chaque composante et risq, surv et surv0 pour chaque evt
        risq=0.d0
        surv=0.d0
        surv0=0.d0
        survint=0.d0

        sumnrisq=0
        do ke=1,nbevt

           brisq=0.d0
           if (logspecif.eq.1) then
              do k=1,nprisq(ke)
                 brisq(k)=exp(b1(sumnrisq+k))
              end do

           else
              do k=1,nprisq(ke)
                 brisq(k)=b1(sumnrisq+k)*b1(sumnrisq+k)
              end do
           end if

           call fct_risq_irtsre(i,ke,brisq,risq,surv,surv0,survint)

           sumnrisq = sumnrisq + nprisq(ke)
        end do


        ! variables explicatives de la survie
        Xevt=0.d0
        bevt=0.d0
        bevtint=0.d0
        if (nxevt.ne.0) then

           m=0
           do ke=1,nbevt
              nevtxcurr=0
              do k=1,nv 

                 if (idtdv(k).ne.1) then

                    if (idsurv(k).eq.1) then  
                       m=m+1
                       bevt(m)=b1(nrisqtot+nevtxcurr+1)
                       Xevt(m)=X(nmescur+1,k)
                    else
                       if (idsurv(k).eq.2) then   
                          m=m+1
                          bevt(m)=b1(nrisqtot+nevtxcurr+ke)
                          Xevt(m)=X(nmescur+1,k)
                       end if
                    end if

                 else ! i.e timedepvar

                    if (idsurv(k).eq.1) then  
                       bevtint(ke)=b1(nrisqtot+nevtxcurr+1)
                    else 
                       if (idsurv(k).eq.2) then
                          bevtint(ke)=b1(nrisqtot+nevtxcurr+ke)
                       end if
                    end if
                 end if
                 nevtxcurr=nevtxcurr+nevtparx(k)
              end do
           end do

        end if

        ! vrais survie
        surv_glob=0.d0
        surv0_glob=0.d0
        varexpsurv=0.d0
        nxevtcurr=0
        fevt=0.d0
        m=0
        do ke=1,nbevt

           ! calculer Xevt * bevt
           varexpsurv=0.d0
           if (nxevt.ne.0) then
              varexpsurv=DOT_PRODUCT(Xevt((nxevtcurr+1):(nxevtcurr+nxevt))&
                   ,bevt((nxevtcurr+1):(nxevtcurr+nxevt)))
           end if

           ! effets aleatoires partages
           easurv=0.d0
           if(nea.gt.0) then
              easurv=DOT_PRODUCT(ui,b1((nrisqtot+nvarxevt+m+1)&
                   :(nrisqtot+nvarxevt+m+nea)))
              m = m+nea
           end if

           ! avoir evt au temps Ti si Devt=1
           if (Devt(i).eq.ke) then
              fevt=risq(ke)*exp(varexpsurv+easurv)
              if (ind_survint(i).eq.1) then
                 fevt=fevt*exp(bevtint(ke))
              end if
           end if

           ! risque cumule jusque Ti
           Surv_glob=surv_glob + survint(ke)*exp(varexpsurv+easurv) + &
                exp(bevtint(ke)+varexpsurv+easurv)*(surv(ke)-survint(ke))

           ! troncature : risque cumule au temps T0
           surv0_glob=surv0_glob+surv0(ke)*exp(varexpsurv+easurv)

           nxevtcurr=nxevtcurr+nxevt
        end do

        ! vraisemblance de la partie survie
        vrais_surv = exp(-Surv_glob)

        if(Devt(i).eq.1) vrais_surv = vrais_surv * fevt

        if (idtrunc.eq.1) then
           ! vrais_surv = vrais_surv / exp(-surv0_glob)
           trunc = trunc + exp(-surv0_glob)
        end if

        ! vrais totale 
        som = som + vrais_Y * vrais_surv
        !print*, "vrais_Y=",vrais_Y, " vrais_surv=", vrais_surv

     else !  pas de survie

        som = som + vrais_Y

     end if

     !if(l.lt.4) print*,"l=", l, " som=",som
  end do ! fin boucle MC

  vrais_irtsre_i = log(som) - log(dble(nMC)) + jacobien
  
  if (idtrunc.eq.1) then
     vrais_irtsre_i = vrais_irtsre_i - log(trunc) + log(dble(nMC))
  end if
  !print*,"i=",i," vrais_Y=",vrais_Y,"jac=",jacobien," vrais_surv=",vrais_surv," vrais_irtsre_i=",vrais_irtsre_i
!  if(expectancy.eq.1) print*,"vrais_surv= ", vrais_surv, "vrais_irtsre_i=",vrais_irtsre_i
654 continue

  return

end function vrais_irtsre_i





double precision function vrais_irtsre(b,m,id,thi,jd,thj)


  use modirtsre,only:ns,nmes,nmescur
  
  implicit none

  integer::m,i,id,jd
  double precision::thi,thj,vrais_irtsre_i,temp
  double precision,dimension(m)::b

  nmescur=0
  vrais_irtsre=0.d0
  do i=1,ns
     temp=vrais_irtsre_i(b,m,id,thi,jd,thj,i)          
     vrais_irtsre = vrais_irtsre + temp
     if (temp.eq.-1.d9 .or. temp/temp.ne.1) then 
        !if (temp/temp.ne.1) write(*,*)"i=",i,"vrais= ",temp
        !if (temp.eq.-1.d9) then 
        vrais_irtsre = -1.d9
        !print*,"dans vrais i=",i," vrais_irtsre=",vrais_irtsre," m=",m," b=",b
        ! if(verbose==1) write(*,*)"i=",i,"vrais= ",temp
        goto 541
     end if
     nmescur = nmescur + sum(nmes(i,:))
  end do
541 continue
  return

end function vrais_irtsre




! =============================================
! subroutine de creation de design matrix
! =============================================




subroutine design_splines_irtsre (ier)

  use modirtsre

  implicit none

  integer ::jj,l,k,ier,yk,q,sumnval,qq,sumqqval
  double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n,hn,hht

  ier=0
  jj=0
  l=0
  q=0
  qq=0
  sumqqval=0
  sumnval=0
  do yk=1,ny
     if (idlink(yk).eq.2) then 
        q=q+1
        do jj=1,nvalSPL(q)      !     ou se trouve la valeur de zi

           do k = 2,ntr(yk)-2
              if ((uniqueY(sumqqval+sumnval+jj).ge.zitr(k-1,q)).and.(uniqueY(sumqqval+sumnval+jj).lt.zitr(k,q))) then
                 l=k-1
              end if
           End do


           if (uniqueY(sumqqval+sumnval+jj).eq.zitr(ntr(yk)-2,q)) then
              l=ntr(yk)-3
           end if

           ht2 = zitr(l+1,q)-uniqueY(sumqqval+sumnval+jj)
           htm= uniqueY(sumqqval+sumnval+jj)-zitr(l-1,q)
           ht = uniqueY(sumqqval+sumnval+jj)-zitr(l,q)
           ht3 = zitr(l+2,q)-uniqueY(sumqqval+sumnval+jj)
           hht = uniqueY(sumqqval+sumnval+jj)-zitr(l-2,q)
           h = zitr(l+1,q)-zitr(l,q)
           hh= zitr(l+1,q)-zitr(l-1,q)
           hn= zitr(l+1,q)-zitr(l-2,q)
           h2n=zitr(l+2,q)-zitr(l-1,q)
           h2= zitr(l+2,q)-zitr(l,q)
           h3= zitr(l+3,q)-zitr(l,q)

           if (uniqueY(sumqqval+sumnval+jj).ne.zitr(ntr(yk)-2,q)) then
              mm2(sumnval+jj) = (3.d0*ht2*ht2)/(hh*h*hn)
              mm1(sumnval+jj) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
              mm(sumnval+jj)  = (3.d0*ht*ht)/(h3*h2*h)

           end if
           if (uniqueY(sumqqval+sumnval+jj).eq.zitr(ntr(yk)-2,q)) then
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

     if(idlink(yk).eq.3) then
        qq = qq + 1
        sumqqval = sumqqval + nvalORD(qq)
     end if

  end do


765 continue

end subroutine design_splines_irtsre
!fin design_splines





subroutine splines_irtsre(k)
  use modirtsre
  implicit none

  integer::k
  integer::i,kk,n,l
  double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n, &
       hn,hh2,hh3



  l=0
  Tmm=0.d0
  Tmm1=0.d0
  Tmm2=0.d0
  Tmm3=0.d0
  Tim=0.d0
  Tim1=0.d0
  Tim2=0.d0
  Tim3=0.d0
  Tmm0=0.d0
  Tmm01=0.d0
  Tmm02=0.d0
  Tmm03=0.d0
  Tim0=0.d0
  Tim01=0.d0
  Tim02=0.d0
  Tim03=0.d0
  Tmmt=0.d0
  Tmmt1=0.d0
  Tmmt2=0.d0
  Tmmt3=0.d0
  Timt=0.d0
  Timt1=0.d0
  Timt2=0.d0
  Timt3=0.d0


  zi(-2,k)=zi(1,k)
  zi(-1,k)=zi(1,k)
  zi(0,k)=zi(1,k)
  zi(nz(k)+1,k)=zi(nz(k),k)
  zi(nz(k)+2,k)=zi(nz(k),k)
  zi(nz(k)+3,k)=zi(nz(k),k)


  n=nz(k)+2
  !------------------- Tsurv ---------------------------
  Do i=1,ns

     do kk=2,n-2
        if ((Tsurv(i).ge.zi(kk-1,k)).and.  &
             Tsurv(i).lt.zi(kk,k)) then
           l=kk-1
        end if
     end do

     if (Tsurv(i).eq.zi(n-2,k)) then
        l=n-3
     end if

     ht = Tsurv(i)-zi(l,k)
     htm = Tsurv(i)-zi(l-1,k)
     h2t = Tsurv(i)-zi(l+2,k)
     ht2 = zi(l+1,k)-Tsurv(i)
     ht3 = zi(l+3,k)-Tsurv(i)
     hht = Tsurv(i)-zi(l-2,k)
     h = zi(l+1,k)-zi(l,k)
     hh = zi(l+1,k)-zi(l-1,k)
     h2 = zi(l+2,k)-zi(l,k)
     h3 = zi(l+3,k)-zi(l,k)
     h4 = zi(l+4,k)-zi(l,k)
     h3m = zi(l+3,k)-zi(l-1,k)
     h2n = zi(l+2,k)-zi(l-1,k)
     hn = zi(l+1,k)-zi(l-2,k)
     hh3 = zi(l+1,k)-zi(l-3,k)
     hh2 = zi(l+2,k)-zi(l-2,k)

     if (Tsurv(i).ne.zi(n-2,k)) then

        Tmm3(ns*(k-1)+i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        Tmm2(ns*(k-1)+i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
             +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))  &
             +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        Tmm1(ns*(k-1)+i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
             +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))   &
             +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        Tmm(ns*(k-1)+i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

     end if

     if (Tsurv(i).eq.zi(n-2,k)) then

        Tmm3(ns*(k-1)+i) = 0.d0
        Tmm2(ns*(k-1)+i) = 0.d0
        Tmm1(ns*(k-1)+i) = 0.d0
        Tmm(ns*(k-1)+i) = 4.d0/h

     end if

     Tim3(ns*(k-1)+i) = (0.25d0*(Tsurv(i)-zi(l-3,k))*Tmm3(ns*(k-1)+i)) &
          +(0.25d0*hh2*Tmm2(ns*(k-1)+i))        &
          +(0.25d0*h3m*Tmm1(ns*(k-1)+i))+(0.25d0*h4*Tmm(ns*(k-1)+i))
     Tim2(ns*(k-1)+i) = (0.25d0*hht*Tmm2(ns*(k-1)+i))  &
          +(h3m*Tmm1(ns*(k-1)+i)*0.25d0)+(h4*Tmm(ns*(k-1)+i)*0.25d0)
     Tim1(ns*(k-1)+i) = (htm*Tmm1(ns*(k-1)+i)*0.25d0)+(h4*Tmm(ns*(k-1)+i)*0.25d0)
     Tim(ns*(k-1)+i) = ht*Tmm(ns*(k-1)+i)*0.25d0

     !------------------- Tsurv0 --------------------------

     if (idtrunc.eq.1) then

        do kk=2,n-2
           if ((Tsurv0(i).ge.zi(kk-1,k)).and.   &
                Tsurv0(i).lt.zi(kk,k)) then
              l=kk-1
           end if
        end do

        if (Tsurv0(i).eq.zi(n-2,k)) then
           l=n-3
        end if

        ht = Tsurv0(i)-zi(l,k)
        htm = Tsurv0(i)-zi(l-1,k)
        h2t = Tsurv0(i)-zi(l+2,k)
        ht2 = zi(l+1,k)-Tsurv0(i)
        ht3 = zi(l+3,k)-Tsurv0(i)
        hht = Tsurv0(i)-zi(l-2,k)
        h = zi(l+1,k)-zi(l,k)
        hh = zi(l+1,k)-zi(l-1,k)
        h2 = zi(l+2,k)-zi(l,k)
        h3 = zi(l+3,k)-zi(l,k)
        h4 = zi(l+4,k)-zi(l,k)
        h3m = zi(l+3,k)-zi(l-1,k)
        h2n = zi(l+2,k)-zi(l-1,k)
        hn = zi(l+1,k)-zi(l-2,k)
        hh3 = zi(l+1,k)-zi(l-3,k)
        hh2 = zi(l+2,k)-zi(l-2,k)

        if (Tsurv0(i).ne.zi(nz(k)-2,k)) then

           Tmm03(ns*(k-1)+i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))

           Tmm02(ns*(k-1)+i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))   &
                +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmm01(ns*(k-1)+i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))    &
                +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmm0(ns*(k-1)+i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

        end if

        if (Tsurv0(i).eq.zi(n-2,k)) then

           Tmm03(ns*(k-1)+i) = 0.d0
           Tmm02(ns*(k-1)+i) = 0.d0
           Tmm01(ns*(k-1)+i) = 0.d0
           Tmm0(ns*(k-1)+i) = 4.d0/h

        end if

        Tim03(ns*(k-1)+i) = (0.25d0*(Tsurv0(i)-zi(l-3,k))*Tmm03(ns*(k-1)+i))  &
             +(0.25d0*hh2*Tmm02(ns*(k-1)+i))           &
             +(0.25d0*h3m*Tmm01(ns*(k-1)+i))+(0.25d0*h4*Tmm0(ns*(k-1)+i))
        Tim02(ns*(k-1)+i) = (0.25d0*hht*Tmm02(ns*(k-1)+i))                  &
             +(h3m*Tmm01(ns*(k-1)+i)*0.25d0)+(h4*Tmm0(ns*(k-1)+i)*0.25d0)
        Tim01(ns*(k-1)+i) = (htm*Tmm01(ns*(k-1)+i)*0.25d0)+(h4*Tmm0(ns*(k-1)+i)*0.25d0)
        Tim0(ns*(k-1)+i) = ht*Tmm0(ns*(k-1)+i)*0.25d0

     end if


     !------------------- Tsurvint --------------------------
     if (ind_survint(i).eq.1) then
        do kk=2,n-2
           if ((Tsurvint(i).ge.zi(kk-1,k)).and. &
                Tsurvint(i).lt.zi(kk,k)) then
              l=kk-1
           end if
        end do

        if (Tsurvint(i).eq.zi(nz(k)-2,k)) then
           l=n-3
        end if

        ht = Tsurvint(i)-zi(l,k)
        htm = Tsurvint(i)-zi(l-1,k)
        h2t = Tsurvint(i)-zi(l+2,k)
        ht2 = zi(l+1,k)-Tsurvint(i)
        ht3 = zi(l+3,k)-Tsurvint(i)
        hht = Tsurvint(i)-zi(l-2,k)
        h = zi(l+1,k)-zi(l,k)
        hh = zi(l+1,k)-zi(l-1,k)
        h2 = zi(l+2,k)-zi(l,k)
        h3 = zi(l+3,k)-zi(l,k)
        h4 = zi(l+4,k)-zi(l,k)
        h3m = zi(l+3,k)-zi(l-1,k)
        h2n = zi(l+2,k)-zi(l-1,k)
        hn = zi(l+1,k)-zi(l-2,k)
        hh3 = zi(l+1,k)-zi(l-3,k)
        hh2 = zi(l+2,k)-zi(l-2,k)

        if (Tsurvint(i).ne.zi(nz(k)-2,k)) then

           Tmmt3(ns*(k-1)+i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
           Tmmt2(ns*(k-1)+i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn)) &
                +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))     &
                +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmmt1(ns*(k-1)+i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))       &
                +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmmt(ns*(k-1)+i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

        end if

        if (Tsurvint(i).eq.zi(nz(k)-2,k)) then

           Tmmt3(ns*(k-1)+i) = 0.d0
           Tmmt2(ns*(k-1)+i) = 0.d0
           Tmmt1(ns*(k-1)+i) = 0.d0
           Tmmt(ns*(k-1)+i) = 4.d0/h

        end if

        Timt3(ns*(k-1)+i) = (0.25d0*(Tsurvint(i)-zi(l-3,k))*Tmmt3(ns*(k-1)+i)) &
             +(0.25d0*hh2*Tmmt2(ns*(k-1)+i))               &
             +(0.25d0*h3m*Tmmt1(ns*(k-1)+i))+(0.25d0*h4*Tmmt(ns*(k-1)+i))
        Timt2(ns*(k-1)+i) = (0.25d0*hht*Tmmt2(ns*(k-1)+i))                   &
             +(h3m*Tmmt1(ns*(k-1)+i)*0.25d0)+(h4*Tmmt(ns*(k-1)+i)*0.25d0)
        Timt1(ns*(k-1)+i) = (htm*Tmmt1(ns*(k-1)+i)*0.25d0)+(h4*Tmmt(ns*(k-1)+i)*0.25d0)
        Timt(ns*(k-1)+i) = ht*Tmmt(ns*(k-1)+i)*0.25d0

     else
        Timt3(ns*(k-1)+i) =Tim3(ns*(k-1)+i)
        Timt2(ns*(k-1)+i) =Tim2(ns*(k-1)+i)
        Timt1(ns*(k-1)+i) =Tim1(ns*(k-1)+i)
        Timt(ns*(k-1)+i) =Tim(ns*(k-1)+i)
     end if

  End Do

end subroutine splines_irtsre





subroutine fct_risq_irtsre(i,k,brisq,risq,surv,surv0,survint)

  use modirtsre
        
  implicit none
  
  integer::i,k,g
  double precision,dimension(nprisq(k))::brisq
  double precision,dimension(nbevt)::risq,surv,surv0,survint
  
  integer::j,l,ll,kk,ii
  double precision::som
  
  if (typrisq(k).eq.2.and.logspecif.eq.1) then
     
     surv(k)=brisq(1)*(tsurv(i)-zi(1,k))**brisq(2)
     
     risq(k)=brisq(1)*brisq(2)*(tsurv(i)-zi(1,k))**(brisq(2)-1)
     if (idtrunc.eq.1) then
        surv0(k)=brisq(1)*(tsurv0(i)-zi(1,k))**brisq(2)
     end if
     if (ind_survint(i).eq.1) then
        survint(k)=brisq(1)*(tsurvint(i)-zi(1,k))**brisq(2)
     else
        survint(k)=surv(k)
     end if
     
  end if
  if (typrisq(k).eq.2.and.logspecif.eq.0) then
     
     surv(k)=(brisq(1)*(tsurv(i)-zi(1,k)))**brisq(2)
     
     risq(k)=brisq(1)*brisq(2)*(brisq(1)*(tsurv(i)-zi(1,k)))**(brisq(2)-1)
     if (idtrunc.eq.1) then
        surv0(k)=(brisq(1)*(tsurv0(i)-zi(1,k)))**brisq(2)
     end if
     if (ind_survint(i).eq.1) then
        survint(k)=(brisq(1)*(tsurvint(i)-zi(1,k)))**brisq(2)
     else
        survint(k)=surv(k)
     end if
     
  end if
 
  if (typrisq(k).eq.1) then
     do j=1,nz(k)-1
        som=0.d0
        do l=1,j-1
           som=som+brisq(l)*(zi(l+1,k)-zi(l,k))
        end do
        if (idtrunc.eq.1) then
           if (Tsurv0(i).ge.zi(j,k).and.Tsurv0(i).le.zi(j+1,k)) then
              surv0(k)=som+brisq(j)*(Tsurv0(i)-zi(j,k))
           end if
        end if
        if (Tsurv(i).ge.zi(j,k).and.Tsurv(i).le.zi(j+1,k)) then
           surv(k)=som+brisq(j)*(Tsurv(i)-zi(j,k))
           risq(k)=brisq(j)
        end if
        if (ind_survint(i).eq.1) then
           if (Tsurvint(i).ge.zi(j,k).and.Tsurvint(i).le.zi(j+1,k)) &
                then
              survint(k)=som+brisq(j)*(Tsurvint(i)-zi(j,k))
           end if
        end if
     end do
     if (ind_survint(i).eq.0) then
        survint(k)=surv(k)
     end if
  end if
  
  

  if (typrisq(k).eq.3) then
     !------------ survie et risq pour Tsurv ----------------
     ll=0
     if (Tsurv(i).eq.zi(nz(k),k)) then
        ll=nz(k)-1
     end if
     som=0.d0
     do kk=2,nz(k)
        if ((Tsurv(i).ge.zi(kk-1,k)).and.(Tsurv(i).lt.zi(kk,k))) &
             then
           ll=kk-1
        end if
     end do
     if (ll.gt.1) then
        do ii=1,ll-1
           som=som+brisq(ii)
        end do
     end if
     
     surv(k)=som+brisq(ll)*Tim3(ns*(k-1)+i)+brisq(ll+1)*Tim2(ns*(k-1)+i) &
          +brisq(ll+2)*Tim1(ns*(k-1)+i)+brisq(ll+3)*Tim(ns*(k-1)+i)
     risq(k)=brisq(ll)*Tmm3(ns*(k-1)+i)+brisq(ll+1)*Tmm2(ns*(k-1)+i)     &
          +brisq(ll+2)*Tmm1(ns*(k-1)+i)+brisq(ll+3)*Tmm(ns*(k-1)+i)
     
     !------------ survie et risq pour Tsurv0 ----------------
     
     if (idtrunc.eq.1) then
        ll=0
        if (Tsurv0(i).eq.zi(nz(k),k)) then
           ll=nz(k)-1
        end if
        som=0.d0
        do kk=2,nz(k)
           if ((Tsurv0(i).ge.zi(kk-1,k)).and.(Tsurv0(i).lt.zi(kk,k))) &
                then
              ll=kk-1
           end if
        end do
        !               if (ll.lt.1.or.ll.gt.nz-1) then
        !                  write(*,*) 'probleme dans fct_risq splines'
        !                  write(*,*) 'll=',ll,'T=',Tsurv0(i)
        !                  stop
        !               end if
        if (ll.gt.1) then
           do ii=1,ll-1
              som=som+brisq(ii)
           end do
        end if
        
        surv0(k)=som+brisq(ll)*Tim03(ns*(k-1)+i)+brisq(ll+1)*Tim02(ns*(k-1)+i) &
             +brisq(ll+2)*Tim01(ns*(k-1)+i)+brisq(ll+3)*Tim0(ns*(k-1)+i)
        
     end if
     
     !------------ survie et risq pour Tsurvint ----------------


     if (ind_survint(i).eq.1) then

        !               write(*,*)'i',i,tsurvint(i),ind_survint(i),tsurv(i),tsurv0(i)
!               write(*,*)timt3(i),Timt2(i),timt1(i)

        ll=0
        if (Tsurvint(i).eq.zi(nz(k),k)) then
           ll=nz(k)-1
        end if
        som=0.d0
        do kk=2,nz(k)
           if((Tsurvint(i).ge.zi(kk-1,k)).and.(Tsurvint(i).lt.zi(kk,k))) &
                then
              ll=kk-1
           end if
        end do
        !               if (ll.lt.1.or.ll.gt.nz-1) then
        !                  write(*,*) 'probleme dans fct_risq splines'
        !                  write(*,*) 'll=',ll,'T=',Tsurvint(i)
        !                  stop
        !               end if
        if (ll.gt.1) then
           do ii=1,ll-1
              som=som+brisq(ii)
           end do
        end if
               
        survint(k)=som+brisq(ll)*Timt3(ns*(k-1)+i)+brisq(ll+1)*Timt2(ns*(k-1)+i) &
             +brisq(ll+2)*Timt1(ns*(k-1)+i)+brisq(ll+3)*Timt(ns*(k-1)+i)
        
     else
        survint(k)=surv(k)
     end if
     
  end if
  
  
end subroutine fct_risq_irtsre



!!!! postfit !!!

!! fct_risq_estime
subroutine fct_risq_estime_irtsre(k,brisq,time,nsim,risq,surv)

  use modirtsre
        
  implicit none
  
  integer::i,k,nsim
  double precision,dimension(nprisq(k))::brisq
  double precision,dimension(nsim)::time
  double precision,dimension(nsim,nbevt)::risq,surv
  
  integer::j,l,ll,kk,ii
  double precision::som
  double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n,hn,hh3,hh2
  double precision::Tmm3_est,Tmm2_est,Tmm1_est,Tmm_est,Tim3_est,Tim2_est,Tim1_est,Tim_est

  if(typrisq(k).eq.3) then
     zi(-2,k)=zi(1,k)
     zi(-1,k)=zi(1,k)
     zi(0,k)=zi(1,k)
     zi(nz(k)+1,k)=zi(nz(k),k)
     zi(nz(k)+2,k)=zi(nz(k),k)
     zi(nz(k)+3,k)=zi(nz(k),k)
  end if
        
  l=0
  do i=1,nsim
     
     if (typrisq(k).eq.2.and.logspecif.eq.1) then
        
        surv(i,k)=brisq(1)*(time(i)-zi(1,k))**brisq(2)     
        risq(i,k)=brisq(1)*brisq(2)*(time(i)-zi(1,k))**(brisq(2)-1)
        
     end if
     if (typrisq(k).eq.2.and.logspecif.eq.0) then
     
        surv(i,k)=(brisq(1)*(time(i)-zi(1,k)))**brisq(2)
        risq(i,k)=brisq(1)*brisq(2)*(brisq(1)*(time(i)-zi(1,k)))**(brisq(2)-1)
     
     end if
 

     if (typrisq(k).eq.1) then
        do j=1,nz(k)-1
           som=0.d0
           do l=1,j-1
              som=som+brisq(l)*(zi(l+1,k)-zi(l,k))
           end do
           if (time(i).ge.zi(j,k).and.time(i).le.zi(j+1,k)) then
              surv(i,k)=som+brisq(j)*(time(i)-zi(j,k))
              risq(i,k)=brisq(j)
           end if
        end do
     end if
  

     if (typrisq(k).eq.3) then
        ll=0
        if (time(i).eq.zi(nz(k),k)) then
           ll=nz(k)-1
        end if
        som=0.d0
        do kk=2,nz(k)
           if ((time(i).ge.zi(kk-1,k)).and.(time(i).lt.zi(kk,k))) &
                then
              ll=kk-1
           end if
        end do
        if (ll.gt.1) then
           do ii=1,ll-1
              som=som+brisq(ii)
           end do
        end if

        ht = time(i)-zi(ll,k)
        htm = time(i)-zi(ll-1,k)
        h2t = time(i)-zi(ll+2,k)
        ht2 = zi(ll+1,k)-time(i)
        ht3 = zi(ll+3,k)-time(i)
        hht = time(i)-zi(ll-2,k)
        h = zi(ll+1,k)-zi(ll,k)
        hh = zi(ll+1,k)-zi(ll-1,k)
        h2 = zi(ll+2,k)-zi(ll,k)
        h3 = zi(ll+3,k)-zi(ll,k)
        h4 = zi(ll+4,k)-zi(ll,k)
        h3m = zi(ll+3,k)-zi(ll-1,k)
        h2n = zi(ll+2,k)-zi(ll-1,k)
        hn = zi(ll+1,k)-zi(ll-2,k)
        hh3 = zi(ll+1,k)-zi(ll-3,k)
        hh2 = zi(ll+2,k)-zi(ll-2,k)
        
     
        if (time(i).ne.zi(nz(k),k)) then
           
           Tmm3_est = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
           Tmm2_est = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))  &
                +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmm1_est = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))   &
                +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmm_est = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
           
        end if
        
        if (time(i).eq.zi(nz(k),k)) then
           
           Tmm3_est = 0.d0
           Tmm2_est = 0.d0
           Tmm1_est = 0.d0
           Tmm_est = 4.d0/h
           
        end if
        
        Tim3_est = (0.25d0*(time(i)-zi(ll-3,k))*Tmm3_est) &
             +(0.25d0*hh2*Tmm2_est)        &
             +(0.25d0*h3m*Tmm1_est)+(0.25d0*h4*Tmm_est)
        Tim2_est = (0.25d0*hht*Tmm2_est)  &
             +(h3m*Tmm1_est*0.25d0)+(h4*Tmm_est*0.25d0)
        Tim1_est = (htm*Tmm1_est*0.25d0)+(h4*Tmm_est*0.25d0)
        Tim_est = ht*Tmm_est*0.25d0
        
        surv(i,k)=som+brisq(ll)*Tim3_est+brisq(ll+1)*Tim2_est &
             +brisq(ll+2)*Tim1_est+brisq(ll+3)*Tim_est
        risq(i,k)=brisq(ll)*Tmm3_est+brisq(ll+1)*Tmm2_est     &
             +brisq(ll+2)*Tmm1_est+brisq(ll+3)*Tmm_est
        
     end if

  end do
  
end subroutine fct_risq_estime_irtsre





subroutine transfos_estimees_irtsre(b,npm,nsim,marker,transfY)
  
  use modirtsre
  
  implicit none
  
  integer::kk,nsim,npm,j,k,yk,sumntr,numSPL,l,avantntr
  double precision,dimension(nsim*ny)::marker,transfY
  double precision,dimension(maxval(ntr))::splaa
  double precision,dimension(maxval(ntr))::Xspl
  double precision,dimension(nsim)::mmm,mmm1,mmm2,iim,iim1,iim2
  double precision::aa1,bb1,dd1,aa,bb,betai,eps,pas,ytemp,cc1
  double precision, dimension(npm)::b,b1
  double precision ::ht,htm,ht2,ht3,hht,h,hh,h2,h3,h2n,hn
  
  
  b1=0.d0
  eps=1.D-20
  do k=1,npm
     b1(k)=b(k)
  end do

  marker=0.d0
  transfY=0.d0

  avantntr = nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor
  
  sumntr = 0
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

        splaa(1)=b1(avantntr+sumntr+1)
        do kk=2,ntr(yk)
           splaa(kk)=b1(avantntr+sumntr+kk)*b1(avantntr+sumntr+kk)
        end do
        
        !calcul de H(y) 
        do j=1,nsim
           ! ou se trouve la valeur
           l=0
           
           do k = 2,ntr(yk)-2
              if ((marker((yk-1)*nsim+j).ge.zitr(k-1,numSPL)).and.(marker((yk-1)*nsim+j).lt.zitr(k,numSPL))) then
                 l=k-1
              end if
           end do
           
           if (marker((yk-1)*nsim+j).eq.zitr(ntr(yk)-2,numSPL)) then
              l=ntr(yk)-3
           end if
           
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
           
           if (marker((yk-1)*nsim+j).ne.zitr(ntr(yk)-2,numSPL)) then
              mmm2(j) = (3.d0*ht2*ht2)/(hh*h*hn)
              mmm1(j) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
              mmm(j)  = (3.d0*ht*ht)/(h3*h2*h)
           end if
           if (marker((yk-1)*nsim+j).eq.zitr(ntr(yk)-2,numSPL)) then
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

     else if (idlink(yk).eq.1) then
        
        aa1=exp(b1(avantntr+sumntr+1))/ &
             (1+exp(b1(avantntr+sumntr+1)))
        bb1=exp(b1(avantntr+sumntr+2))/ &
             (1+exp(b1(avantntr+sumntr+2)))
        bb1=aa1*(1.d0-aa1)*bb1
        cc1=b1(avantntr+sumntr+3)
        dd1=abs(b1(avantntr+sumntr+4))
        
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
           transfY((yk-1)*nsim+j)=(marker((yk-1)*nsim+j)-b1(avantntr+sumntr+1)) &
                /abs(b1(avantntr+sumntr+2))
        end do

     end if

     !! cas idlink = 3 fait dans R
     
     sumntr = sumntr + ntr(yk)
     
  end do

end subroutine transfos_estimees_irtsre












!!! pour calculer les esperances de vie !!!
!! -> on se sert de la fonction de vraisemblance, eventuellement seulement sur une partie des items, et avec expectancy=1


subroutine proba_irtsre(Y0,X0,Tentr0,Tevt0,Devt0,ind_survint0 &
     ,idea0,idg0,idcor0,idcontr0,idsurv0,idtdv0 &
     ,typrisq0,nz0,zi0,nbevt0,idtrunc0,logspecif0 &
     ,ny0,nv0,nobs0,nea0,nmes0,idiag0,ncor0,nalea0&
     ,epsY0,idlink0,nbzitr0,zitr0,uniqueY0,indiceY0 &
     ,nvalSPLORD0,fix0,methInteg0,nMC0,dimMC0,seqMC0,npmtot0,btot,proba)

  use modirtsre
  use optim

  IMPLICIT NONE

  !Declaration des variables en entree
  integer,intent(in)::nv0,ny0,nMC0,methInteg0,dimMC0
  integer, intent(in)::nobs0,idiag0,npmtot0,nea0,ncor0,nalea0
  integer,intent(in)::idtrunc0,logspecif0,nbevt0
  double precision,dimension(1),intent(in)::Tentr0,Tevt0
  integer,dimension(1),intent(in)::ind_survint0,Devt0
  integer, dimension(nv0),intent(in)::idtdv0,idsurv0
  integer,dimension(nbevt0),intent(in)::typrisq0,nz0
  double precision,dimension(maxval(nz0),nbevt0),intent(in)::zi0    
  double precision,dimension(ny0),intent(in)::epsY0
  integer, dimension(ny0),intent(in)::idlink0,nbzitr0,nvalSPLORD0
  double precision,dimension(maxval(nbzitr0),ny0),intent(in)::zitr0
  integer,dimension(nobs0),intent(in)::indiceY0
  double precision,dimension(sum(nvalSPLORD0(:))),intent(in)::uniqueY0
  integer, dimension(nv0),intent(in)::idea0,idg0,idcor0,idcontr0
  integer,dimension(1,ny0)::nmes0   
  double precision,dimension(nobs0),intent(in)::Y0
  double precision,dimension(nobs0*nv0),intent(in)::X0
  integer,dimension(npmtot0),intent(in)::fix0
  double precision,dimension(dimMC0*nMC0),intent(in)::seqMC0
  !!double precision,dimension(ny0),intent(in)::item0
  double precision, dimension(npmtot0), intent(inout) :: btot

  !Declaration des variables en sortie
  double precision,intent(out)::proba

  !Variables locales
  integer::jtemp,i,j,npm,ier,k,ktemp,yk,k1,k2,mi,nbfix
  integer::ke,sumnrisq,it,ns0,id,jd
  double precision::eps,thi,thj
  double precision,dimension(npmtot0)::mvc,b
  double precision,external::vrais_irtsre_i

  !! on a un seul sujet donc maxmes=ni
  maxmes=sum(nmes0(1,:))
  nmescur=0

    allocate(rangeY(ny0),minY(ny0),maxY(ny0),idlink(ny0),ntr(ny0),epsY(ny0))

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

  methInteg = methInteg0
  nMC = nMC0

!print*,"ok2"," nvalsplord=",nvalsplord0, " nvalord=",nvalord
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
     if (idlink(k).eq.0) ntr(k)=2     
     if (idlink(k).eq.1) ntr(k)=4
     if (idlink(k).eq.2) then
        k1=k1+1
        ntr(k)=nbzitr0(k)+2

        zitr(1:nbzitr0(k),k1)=zitr0(1:nbzitr0(k),k)
        zitr(-1,k1)=zitr(1,k1)
        zitr(0,k1)=zitr(1,k1)
        zitr(ntr(k)-1,k1)=zitr(ntr(k)-2,k1)
        zitr(ntr(k),k1)=zitr(ntr(k)-1,k1)
     end if
     if (idlink(k).eq.3) then
        k2 = k2+1
        ntr(k) = nvalORD(k2)-1
     end if
  end do

  ntrtot = sum(ntr)
  
  !if (verbose==1)       write(*,*)'zitr',zitr

  ns0 = 1
  allocate(Y(nobs0),X(nobs0,nv0),uniqueY(ntotvalSPL+ntotvalORD) &
       ,idea(nv0),idg(nv0),idcor(nv0),idcontr(nv0),nmes(ns0,ny0) &
       ,indiceY(nobs0))

  allocate(Tsurv0(ns0),Tsurv(ns0),Tsurvint(ns0),ind_survint(ns0),Devt(ns0))
  allocate(typrisq(nbevt0),nz(nbevt0),nprisq(nbevt0),nevtparx(nv0),nxcurr(nv0))
  allocate(idsurv(nv0),idtdv(nv0))

  ! zi : contient noeuds pour hazard (ou min et max si Weibull)
  if(any(typrisq0.eq.3)) then
     allocate(zi(-2:maxval(nz0)+3,nbevt0))
  else
     allocate(zi(maxval(nz0),nbevt0))
  end if


  eps=1.d-20

  ! enregistrement pour les modules
  nbevt=nbevt0    
  typrisq=typrisq0
  idtrunc=idtrunc0
  Tsurv0=Tentr0   
  Tsurv=Tevt0    
  devt=devt0    
  ind_survint=ind_survint0
  logspecif=logspecif0 
  Tsurvint=Tsurv
  ny=ny0
  ns=1
  nv=nv0
  nobs=nobs0
  ncor=ncor0
  nalea=nalea0
  idiag=idiag0

  !     if (verbose==1) write(*,*)'ntotvalSPL',ntotvalSPL

  if (ntotvalSPL+ntotvalORD.gt.0) uniqueY(1:ntotvalSPL+ntotvalORD)=uniqueY0(1:ntotvalSPL+ntotvalORD)

  nmes=0
  Y=0.d0
  X=0.d0
  idea=0
  idg=0
  idcor=0
  idcontr=0
  idsurv=0
  idtdv=0
  ktemp=0

  do k=1,nv
     idsurv(k)=idsurv0(k)
     idtdv(k)=idtdv0(k)
     idea(k)=idea0(k)
     idg(k)=idg0(k)
     idcor(k)=idcor0(k)
     idcontr(k)=idcontr0(k)

     jtemp=0
     DO i=1,ns
        do yk=1,ny            
           if (k.eq.1) then
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
          ! write(*,*)'X k:',X(1:50,2)

  ! definition Tsurvint 
  nvdepsurv=0
  if(sum(ind_survint).gt.0) then
     nvdepsurv=1
     do k=1,nv
        if (idtdv(k).eq.1) then
           it=0
           do i=1,ns
              Tsurvint(i)=X(it+1,k)
              it=it+maxval(nmes(i,:))
           end do
        end if
     end do
  end if


  ! prm fixes
  allocate(fix(npmtot0))
  fix=0
  fix(1:npmtot0)=fix0(1:npmtot0)
  nbfix=sum(fix)
  if(nbfix.eq.0) then
     allocate(bfix(1))
  else
     allocate(bfix(nbfix))
  end if
  bfix=0.d0

  ! creation des parametres

  nprisq=0
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

     nrisqtot = nrisqtot+nprisq(ke)  ! nb total de prm pour hazards
     zi(1:nz(ke),ke)=zi0(1:nz(ke),ke)
  end do

  ! nvarxevt = nombre total de coef pour survie (sans prm hazard)
  nxevt=0
  nevtparx=0
  do j=1,nv

     if(idtdv(j).ne.1) then

        if(idsurv(j).eq.1) then
           nevtparx(j) = 1
           nxevt = nxevt + 1
        end if
        if(idsurv(j).eq.2) then 
           nevtparx(j) = nbevt
           nxevt = nxevt + 1
        end if
     end if

  end do

  nvarxevt = sum(nevtparx) + nvdepsurv

  nea=0
  nef=0
  ncontr=0
  do k=1,nv
     if (idg(k).eq.1) then
        nef = nef + 1
     end if
     nea=nea+idea(k)
     ncontr=ncontr+idcontr(k)*(ny-1) 
  end do
  nef = nef - 1 !intercept pas estime

  nasso = nbevt*nea


  if (idiag.eq.1) then
     nvc=nea-1
  else if(idiag.eq.0) then
     nvc=(nea+1)*nea/2-1
  end if

  npmtot = nrisqtot+nvarxevt+nasso+nef+ncontr+nvc+ncor+ntrtot+nalea+ny
  !print*,nrisqtot,nvarxevt,nef,ncontr,nvc,ncor,ntrtot,nalea,ny,nasso
  !print*,"npmtot=",npmtot
  !print*,"btot=",btot
  ! if (verbose==1 ) write(*,*)"npm0=",npm0,"npm=",npm

  !  write(*,*)'idlink',idlink
  !  write(*,*)'idea',idea
  !  write(*,*)'NVC',nvc


  ! points qmc
  if(methInteg.ne.3) then 
     allocate(seqMC(1))
  else
     allocate(seqMC(dimMC0*nMC))
     seqMC = seqMC0(1:dimMC0*nMC) 
  end if

  ! creation du vecteur b avec slt les prm a estimer
  b=0.d0
  npm=0
  k=0
  do j=1,npmtot
     if(fix0(j).eq.0) then
        npm=npm+1
        b(npm)=btot(j)
     end if
     if(fix0(j).eq.1) then
        k=k+1
        bfix(k)=btot(j)
     end if
  end do


  ! creer base de splines si au moins un hazard splines
  if(any(typrisq.eq.3)) then
     allocate(Tmm(ns*nbevt),Tmm1(ns*nbevt),Tmm2(ns*nbevt),Tmm3(ns*nbevt),Tim(ns*nbevt)         &
          ,Tim1(ns*nbevt),Tim2(ns*nbevt),Tim3(ns*nbevt),Tmm0(ns*nbevt),Tmm01(ns*nbevt),Tmm02(ns*nbevt)  &
          ,Tmm03(ns*nbevt),Tim0(ns*nbevt),Tim01(ns*nbevt),Tim02(ns*nbevt),Tim03(ns*nbevt),          &
          Tmmt(ns*nbevt),Tmmt1(ns*nbevt),Tmmt2(ns*nbevt),Tmmt3(ns*nbevt),Timt(ns*nbevt),Timt1(ns*nbevt) &
          ,Timt2(ns*nbevt),Timt3(ns*nbevt))

     do ke=1,nbevt
        if(typrisq(ke).eq.3) then
           call splines_irtsre(ke)
        end if
     end do

  end if



  ! base de splines transfos
  if (any(idlink.eq.2)) then 
     call design_splines_irtsre(ier)
     if (ier.eq.-1) then
        proba=-1.d9
        go to 1589
     end if
  end if

  ! indicateur vraisemblance ou esperance de vie
  expectancy = 1 ! on va calculer une esperance de vie
  !! donc une valeur j pour un outcome ordinal veut dire y <= j

  ! calculer la proba
  id=0
  jd=0
  thi=0.d0
  thj=0.d0
  i=1

  !print*,"dans proba Y=",Y, " nobs=", nobs, " nmes=", nmes
  !print*," indiceY=", indiceY, " minY=",minY, " maxY=", maxY, " nvalORD=", nvalORD
  
  proba = vrais_irtsre_i(b,npm,id,thi,jd,thj,i)
  
1589 continue

  if (any(typrisq.eq.3)) then
     deallocate(Tmm,Tmm1,Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,    &
          Tmm01,Tmm02,Tmm03,Tim0,Tim01,Tim02,Tim03,Tmmt,Tmmt1,     &
          Tmmt2,Tmmt3,Timt,Timt1,Timt2,Timt3)
  endif

  deallocate(Tsurv0,Tsurv,Tsurvint &
       ,ind_survint,zi,devt,typrisq,nz,nprisq,idsurv,idtdv &
       ,nevtparx,nxcurr)

  deallocate(Y,X,idea,idg,idcor,idcontr,nmes,uniqueY,indiceY,ntr)


  deallocate(zitr,mm,mm1,mm2,im,im1,im2,minY,maxY,rangeY,idlink,nvalSPL,nvalORD,epsY)

  deallocate(fix,bfix,seqMC)


  !write(*,*)'fin'
  return
end subroutine proba_irtsre


