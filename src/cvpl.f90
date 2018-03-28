


!--------------------------------------------------------------------
!                  commun_cvpl MODULE
!--------------------------------------------------------------------



module commun_cvpl
  integer,save::ns_s,nmes_curr_s
  double precision,dimension(:),allocatable,save::Y_s,time_cvpl
  double precision,dimension(:,:),allocatable,save::X_s
  integer,dimension(:),allocatable,save::nmes_s
  double precision,dimension(:),allocatable,save::Tsurv0_s,Tsurv_s,Tsurvint_s
  integer,dimension(:),allocatable,save::Devt_s
  integer,dimension(:),allocatable,save::ind_survint_s
  double precision,dimension(:),allocatable,save::Tmm_valt,Tmm1_valt, &
       Tmm2_valt,Tmm3_valt,Tim_valt,Tim1_valt,Tim2_valt,Tim3_valt
  integer,dimension(:),allocatable::indT
end module commun_cvpl





!--------------------------------------------------------------------
!                  Function CVPL
!--------------------------------------------------------------------


subroutine cvpl(Y0,X0,Prior0,idprob0,idea0,idg0,idcor0 &
     ,idcom0,idspecif0,idtdv0 &
     ,ns0,ng0,ncor0,nv0,nobs0,nmes0,idiag0,nwg0,npmtot0,time0,typrisq0 &
     ,idtrunc0,risqcom0,nz0,zi0,Tentr0,Tevt0 &
     ,devt0,ind_survint0,vopt &
     ,nT,valT,btot,epoir,rl_cond,ns_vect,nevt_vect &
     ,contribt,logspecif0 &
     ,fix0,npm0)


  use commun_cvpl
  use commun_comp
  use optim


  implicit none


  !        Declaration des variables en entree
  integer,intent(in)::nv0,logspecif0,npm0
  integer,intent(in)::ns0,ng0,nobs0,idiag0,nwg0,npmtot0,ncor0
  integer,dimension(nv0),intent(in)::idea0,idg0,idprob0,idcor0
  integer,dimension(nv0),intent(in)::idcom0,idspecif0,idtdv0
  integer,dimension(ns0),intent(in)::nmes0,prior0
  double precision,dimension(ns0)::Tentr0,Tevt0
  double precision,dimension(nobs0),intent(in)::Y0,time0
  double precision,dimension(nobs0*nv0),intent(in)::X0
  integer,intent(in)::risqcom0,typrisq0,nz0,idtrunc0
  integer,dimension(ns0)::Devt0,ind_survint0
  integer,dimension(npmtot0),intent(in)::fix0
  double precision,dimension(nz0),intent(in)::zi0


  !        Declaration des variable en entree et sortie
  double precision, dimension(npmtot0), intent(inout) :: btot


  !-------------------> argument propres a epoir
  integer,intent(in)::nt
  double precision,dimension(nt),intent(out)::rl_cond,epoir
  double precision,dimension(nt),intent(in)::valT
  integer,dimension(nt)::ns_vect,nevt_vect
  double precision,dimension(ns0*nt)::contribt
  double precision, dimension(npmtot0*(npmtot0+1)/2),intent(in)::vopt
  !-------------------> variables locales joint
  integer :: i,j,k,ktemp,nmestot
  double precision::eps
  integer::npm
  !-------------------> variables locales joint
  double precision::trace3,rl_condt
  double precision,dimension(ns0)::rlindiv
  integer::t,nbevt_s,nmescur,nmescur_s,jtemp
  double precision,dimension(npm0,npm0)::J_cond,mat3
  double precision,dimension(npm0,npm0)::H_1
  double precision, dimension(npmtot0) :: mvc,b
  integer::ier,ke,it,jj,kk,nbfix
  !-------------> initialisation des sorties

  ! verification des donnees en entree

  eps=1.d-20
  maxmes=0
  do i=1,ns0
     if (nmes0(i).gt.maxmes) then
        maxmes=nmes0(i)
     end if
  end do

  !-------------> allocation generale

  !----------------- allocation generale ---------------------------
  allocate(Y(nobs0),idprob(nv0),X(nobs0,nv0)   &
       ,idea(nv0),idg(nv0),idcor(nv0),nmes(ns0),Tsurv0(ns0),Tsurv(ns0)    &
       ,Tsurvint(ns0),ind_survint(ns0),idcom(nv0),idspecif(nv0)   &
       ,idtdv(nv0),devt(ns0),prior(ns0),time_cvpl(nobs0))


  nbevt=1 !pas de competitif!
  idlink=-1 !link=NULL!

  allocate(risqcom(nbevt),typrisq(nbevt),nz(nbevt),nprisq(nbevt) &
       ,nrisq(nbevt),nxevtspec(nbevt),nevtparx(nv0),nxcurr(nv0))


  Tsurv0=Tentr0
  Tsurv=Tevt0
  devt=devt0
  ind_survint=ind_survint0
  logspecif=logspecif0


  typrisq(1)=typrisq0
  risqcom(1)=risqcom0
  idtrunc=idtrunc0


  Y=0.d0
  X=0.d0
  idprob=0
  idea=0
  idg=0
  idcor=0
  idcom=0
  idspecif=0
  idtdv=0
  prior=0
  nmes=0
  time_cvpl=0.d0


  nz(1)=nz0
  select case (typrisq(1))
  case (1)
     allocate(zi(nz0,1))
     zi(1:nz0,1)=zi0(1:nz0)
  case (2)
     allocate(zi(nz0,1))
     zi(1:nz0,1)=zi0(1:nz0)
  case (3)
     allocate(zi(-2:nz0+3,1))
     do i=1,nz(1)
        zi(i,1)=zi0(i)
     end do
  end select



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

  ns=ns0
  ng=ng0
  nv=nv0
  ncor=ncor0
  nobs=nobs0

  if (nwg0.eq.0) then
     nwg=0
  else
     nwg=ng-1
  end if
  idiag=idiag0
  !------------------------
  if (typrisq(1).eq.1) then
     nprisq(1)=nz(1)-1
  end if
  if (typrisq(1).eq.2) then
     nprisq(1)=2
  end if
  if (typrisq(1).eq.3) then
     nprisq(1)=nz(1)+2
  end if
  !--------------------
  nmes=nmes0
  nmestot=0
  ktemp=0
  do k=1,nv
     idprob(k)=idprob0(k)
     idea(k)=idea0(k)
     idg(k)=idg0(k)
     idcor(k)=idcor0(k)
     idcom(k)=idcom0(k)
     idtdv(k)=idtdv0(k)
     idspecif(k)=idspecif0(k)
     jtemp=0
     DO i=1,ns
        if (k.eq.1) then
           nmes(i)=nmes0(i)
           prior(i)=prior0(i)
           do j=1,nmes(i)
              nmestot=nmestot+1
              Y(nmestot)=Y0(nmestot)
              time_cvpl(nmestot)=time0(nmestot)
           end do
        end if
        do j=1,nmes(i)
           ktemp=ktemp+1
           jtemp=jtemp+1
           X(jtemp,k)=X0(ktemp)
        end do
     end do
  end do

  nrisq=0
  if (risqcom(1).eq.1) then
     nrisq(1)=nprisq(1)
  end if
  if (risqcom(1).eq.2) then
     nrisq(1)=nprisq(1)+ng-1
  end if
  if (risqcom(1).eq.0) then
     nrisq(1)=nprisq(1)*ng
  end if
  nrisqtot = sum(nrisq)

  ! definition Tsurvint 
  tsurvint=tsurv
  nvdepsurv=0
  if(sum(ind_survint).gt.0) then
     nvdepsurv=1
     do k=1,nv
        if (idtdv(k).eq.1) then
           it=0
           do i=1,ns
              Tsurvint(i)=X(it+1,k)
              it=it+nmes(i)
           end do
        end if
     end do
  end if


  ! parms fixes pour les vexp dans la modelisation de l'evenement
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




  !-----------------------------------------------------------
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
        ncg=ncg+1      ! nb var. dans melange
     end if
     nea=nea+idea(k)
     nprob=nprob+(idprob(k))*(ng-1)
     nvarprob=nvarprob+idprob(k)
  end do


  !       if((ng.eq.1.and.ncg.gt.0).or.(ng.eq.1.and.nprob.gt.0)) then
  !    istop=12
  !    go to 1588
  !       end if


  !  nb effets fixes = nb effets fixes sans melange
  !                  + ng fois le nb de var dans melange

  if (idiag.eq.1) then
     nvc=nea
  else if(idiag.eq.0) then
     nvc=nea*(nea+1)/2
  end if

  nef=ncssg+ncg*ng !
  npmtot=nprob+nvarxevt+nrisqtot+nef+nvc+nwg+ncor+1

  ntrtot=1

  if (idiag.eq.1.and.nvc.gt.0) then
     DO j=1,nvc
        btot(nprob+nvarxevt+nrisqtot+nef+j)=dsqrt(abs(btot(nprob+nvarxevt+nrisqtot+nef+j)))
     END DO
  end if
  ! si idiag=0, on met dans le vecteur des parms, les parms
  ! de la transformee de Cholesky

  if (idiag.eq.0.and.nvc.gt.0) then

     DO j=1,nvc
        mvc(j)=btot(nprob+nvarxevt+nrisqtot+nef+j)
     END DO
     CALL DMFSD(mvc,nea,EPS,IER)
     DO j=1,nvc
        btot(nprob+nvarxevt+nrisqtot+nef+j)=mvc(j)
     END DO
  end if

  if (nwg.gt.0) then
     do i=1,nwg
        btot(nprob+nvarxevt+nrisqtot+nef+nvc+i)=abs(btot(nprob+nvarxevt+nrisqtot+nef+nvc+i))
     end do
  end if

  if (typrisq(1).eq.3) then

     allocate(Tmm(ns0),Tmm1(ns0),Tmm2(ns0),Tmm3(ns0),Tim(ns0)         &
          ,Tim1(ns0),Tim2(ns0),Tim3(ns0),Tmm0(ns0),Tmm01(ns0),Tmm02(ns0)  &
          ,Tmm03(ns0),Tim0(ns0),Tim01(ns0),Tim02(ns0),Tim03(ns0),          &
          Tmmt(ns0),Tmmt1(ns0),Tmmt2(ns0),Tmmt3(ns0),Timt(ns0),Timt1(ns0) &
          ,Timt2(ns0),Timt3(ns0))

     allocate(Tmm_valt(nt),Tmm1_valt(nt),Tmm2_valt(nt),Tmm3_valt(nt), &
          Tim_valt(nt),Tim1_valt(nt),Tim2_valt(nt),Tim3_valt(nt))
     call splines_cvpl(nt,valt)
  end if

  !######################################################################################
  !------------------------------------- Debut cvpl -------------------------------------
  !######################################################################################
  !-------> amis en matrice de H_1


  !        write(*,*)'debut CVPL nt:',nt
  !        write(*,*)'debut CVPL nt*ns:',nt*ns

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

  H_1 = 0.d0
  jj=0
  do j = 1,npmtot
     if(fix(j).eq.0) then
        jj=jj+1
        kk=jj-1
        do k =j,npmtot
           if(fix(k).eq.0) then
              kk=kk+1
              H_1(jj,kk) = vopt(j+k*(k-1)/2)
              H_1(kk,jj) = vopt(j+k*(k-1)/2)
           end if
        end do
     end if
  end do

  eps=1.d-20

  contribt=0.d0
  ns_vect=0
  nevt_vect=0

  !        write(*,*)'debut CVPL 2'
  !-------> allocation tableau cvpl
  allocate(ind_survint_s(ns),nmes_s(ns),Y_s(ns0*maxmes),X_s(ns0*maxmes,nv0))
  allocate(Tsurv0_s(ns),Tsurv_s(ns),Tsurvint_s(ns),Devt_s(ns),indT(ns))
  !-------> fin allocation


  do t=1,nT

     !--------> initialisation


    !            write(*,*)' boucle', t


     ns_s = 0
     ind_survint_s = 0
     nmes_s = 0
     Y_s = 0.d0
     X_s = 0.d0
     Tsurv0_s = 0.d0
     Tsurv_s = 0.d0
     Tsurvint_s = 0.d0
     Devt_s = 0
     indT = 0
     nbevt_s=0
     !---------------> boucle sur les individu
     nmescur=0
     nmescur_s=0
     ktemp=0


     do i=1,ns

        !          write(*,*)'i',i,Tsurvint(i),valT(t),time_cvpl(nmescur+1)

        if (Tsurvint(i).ge.valT(t).and.(time_cvpl(nmescur+1).le.valT(t))) then

           indT(i)=1
           ns_s=ns_s+1

           do j=1,nmes(i)
              if (time_cvpl(nmescur+j).le.valT(t)) then
                 nmes_s(ns_s)=nmes_s(ns_s)+1
                 nmescur_s=nmescur_s+1
                 Y_s(nmescur_s)=Y(nmescur+j)
                 do k=1,nv
                    X_s(nmescur_s,k) = X(nmescur+j,k)
                 end do
              end if
           end do

           Tsurv_s(ns_s)=Tsurv(i)
           Tsurv0_s(ns_s)=Tsurv0(i)

           do k=1,nvdepsurv
              Tsurvint_s(ns_s)=Tsurvint(i)
              ind_survint_s(ns_s)=ind_survint(i)
           end do
           Devt_s(ns_s)=Devt(i)
           nbevt_s=nbevt_s+Devt_s(ns_s)

           !              else
           !                 write(*,*)'i exclu',i,Tsurv(i),time_cvpl(nmescur+1),X(nmescur+1,2),valt(t)
           !
        end if
        nmescur=nmescur+nmes(i)
     end do

    ! write(*,*)'apres boucle sujets',ns_s," ",ns

     ns_vect(t)=ns_s
     nevt_vect(t)=nbevt_s

     J_cond=0.d0
     rlindiv=0.d0

     !          write(*,*)'avant derivc'
     !         write(*,*)'parms=',(b(j),j=1,npm)
     !         write(*,*) 'bfix=',bfix

     call derivc_condT(b,npm,J_cond,rlindiv,t,valT(t))

     !           write(*,*)'apres derivc'

     rl_condt=0.d0
     do i=1,ns
        contribt(ns*(t-1)+i)=rlindiv(i)
        if (rlindiv(i).eq.-1.d9) then
           rl_cond(t)=-1.d9
           epoir(t)=1.d9
           goto 5289
        end if
        rl_condt=rl_condt+rlindiv(i)
     end do

     MAT3=MATMUL(H_1,J_cond)
     trace3=0.d0
     do k=1,npm
        trace3=trace3+MAT3(k,k)
     end do

     epoir(t)=-rl_condt/dble(ns_s)+(trace3*dble(ns)/(dble(ns_s)*dble(ns-1)))
     rl_cond(t)=rl_condt/dble(ns_s)

     if (epoir(t).ne.epoir(t)) then
        epoir(t)=1.d9
     end if
     if (rl_cond(t).ne.rl_cond(t)) then
        rl_cond(t)=-1.d9
     end if

5289 continue


  end do

!  write(*,*)'avant deallocate 1'

  deallocate(Y,idprob,X,idea,idg,idcor,nmes,Tsurv0,Tsurv,Tsurvint, &
       ind_survint,idcom,idspecif,idtdv,devt,prior,zi,time_cvpl)



!  write(*,*)'avant deallocate 2'

  deallocate(ind_survint_s,nmes_s,Y_s,X_s,Tsurv0_s,Tsurv_s,Tsurvint_s,Devt_s,indT)



!  write(*,*)'avant deallocate 3'

  if (typrisq(1).eq.3) then
     deallocate(Tmm,Tmm1,Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,    &
          Tmm01,Tmm02,Tmm03,Tim0,Tim01,Tim02,Tim03,Tmmt,Tmmt1,     &
          Tmmt2,Tmmt3,Timt,Timt1,Timt2,Timt3)

     deallocate(Tmm_valt,Tmm1_valt,Tmm2_valt,Tmm3_valt,Tim_valt, &
          Tim1_valt,Tim2_valt,Tim3_valt)
  endif


!  write(*,*)'avant deallocate 4'
  deallocate(risqcom,typrisq,nz,nprisq,nrisq,nxevtspec,nevtparx,nxcurr)

!  write(*,*)'avant deallocate 5'


  deallocate(fix,bfix)

  !       write(*,*)"dans CVPL : valeur de la vraisemblance"
  !       do t=1,nT
  !          write(*,*)valt(t),epoir(t),rl_cond(t),ns_vect(t),nevt_vect(t)
  !       end do
  !        write(*,*)"fin du programme"

!  write(*,*)'avant return'



  return



end subroutine cvpl
























!-----------------------------------------------------------
!                        derivc_condt
!------------------------------------------------------------
subroutine derivc_condt(b,m,v,rlindiv,t1,valt)

  ! Calcul du gradient et de la matrice remplacant la hessienne
  ! par Fisher scoring empirique
  use commun_comp
  use commun_cvpl

  IMPLICIT NONE

  integer::m,i,k,ii,id,t1
  double precision::funcpi_condt,thn,th,z,valt,vrais_comp_i,temp1,temp2
  double precision,dimension(m,1)::Uscore, Uscore2
  double precision,dimension(m)::b
  double precision,dimension(m,m)::v
  double precision,dimension(ns)::rlindiv

  V=0.d0
  z=0.d0
  rlindiv=0.d0
  id=0
  ! Calcul des gradients par sujets et totaux
  ii=0
  nmes_curr=0
  nmes_curr_s=0

  DO i=1,ns
     Uscore=0.d0
     Uscore2=0.d0
      !          rlindiv(i)=funcpij(b,m,id,z,id,z,i)
     if (indT(i).eq.1) then
        ii=ii+1
        rlindiv(i)=funcpi_condt(b,m,id,z,id,z,ii,t1,valt)
        if (rlindiv(i).eq.-1.d9) then 
           V=0.d0
           rlindiv=-1.d9
           goto 777
        end if
        do k=1,m
           th=DMAX1(1.d-6, 1.d-4 * DABS(b(k)))
           thn=-1.D0*th
           temp1=funcpi_condt(b,m,k,th,id,z,ii,t1,valt)
           temp2=funcpi_condt(b,m,k,thn,id,z,ii,t1,valt)
           if (temp1.eq.-1.d9.or.temp2.eq.-1.d9) then 
              V=0.d0
              rlindiv=-1.d9
              goto 777
           end if
           Uscore(k,1)=-(temp1-temp2)/(2.d0*th)
           temp1=vrais_comp_i(b,m,k,th,id,z,i)
           temp2=vrais_comp_i(b,m,k,thn,id,z,i)
           if (temp1.eq.-1.d9.or.temp2.eq.-1.d9) then 
              V=0.d0
              rlindiv=-1.d9
              goto 777
           end if

           Uscore2(k,1)=-(temp1-temp2)/(2.d0*th)
        END DO
        nmes_curr_s = nmes_curr_s + nmes_s(ii)
     end if
      V=V+MATMUL(Uscore,transpose(Uscore2))

     nmes_curr = nmes_curr + nmes(i)
  end do

777 continue

  return

end subroutine derivc_condt



!-----------------------------------------------------------
!                        FUNCPA_CONDT
!------------------------------------------------------------

double precision function funcpa_condt(b,m,id,thi,jd,thj,t1,valt)
  use commun_comp
  use commun_cvpl

  implicit none

  integer::m,i,id,jd,t1
  double precision::rl,thi,thj,funcpi_condt
  double precision,dimension(m)::b
  double precision::valt

  rl=0.d0
  nmes_curr_s=0
  do i=1,ns_s
     rl = rl + funcpi_condt(b,m,id,thi,jd,thj,i,t1,valt)
     nmes_curr_s = nmes_curr_s + nmes_s(i)
  end do
  funcpa_condt =rl
  return

end function funcpa_condt

!-----------------------------------------------------------
!                        FUNCPI_CONDT
!------------------------------------------------------------


double precision function funcpi_condt(b,npm,id,thi,jd,thj,i,t1,valt)

  use commun_comp
  use commun_cvpl
  use optim
  IMPLICIT NONE

  integer ::i,j,k,l,m,g,l2,m2,id,jd,jj,npm,t1
  integer ::ier,nmoins,j1,j2
  double precision,dimension(maxmes,nv) ::Z,P,X00,X2
  double precision,dimension(nvarprob) ::Xprob,bprob
  double precision,dimension(nv,nv) ::Ut,Ut1
  double precision,dimension(maxmes,maxmes) ::VC,Se
  double precision,dimension(npm) :: b
  double precision,dimension(npmtot) :: b1
  double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
  double precision,dimension(nv) :: b0,b2
  double precision :: vrais,eps,det
  double precision ::thi,thj,temp
  double precision ::Y4,expo
  double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
  double precision,dimension(ng) :: pi
  double precision,dimension(nrisq(1))::brisq
  double precision::retard,entretard,vrais_survie
  double precision::vrais_long,expo_long
  double precision,dimension(nxevt)::bevt,Xevt
  double precision,dimension(nbevt)::bevtint
  double precision,dimension(ng)::risq,surv,surv0,survint,survT
  double precision::valt
  logical::is_nan
  integer::nevtxcurr,sumnrisq,ke

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

  Ut=0.d0
  If (idiag.eq.1) then
     do j=1,nea
        do k=1,nea
           if (j.eq.k) then
              Ut(j,k)=b1(nprob+nvarxevt+nrisqtot+nef+j)
           else
              Ut(j,k)=0.d0
           end if
        end do
     end do
  end if


  If (idiag.eq.0) then
     do j=1,nea
        do k=1,j
           Ut(j,k)=b1(nprob+nvarxevt+nrisqtot+nef+k+j*(j-1)/2)
        end do
     end do
  end if

  ! calcul de brisq e chaque composante et risq, surv et surv0 pour tous les i et tous les g

  risq=0.d0
  surv=0.d0
  survT=0.d0
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

        call fct_risq_it(i,t1,brisq,g,risq,surv,surv0,survint,survt,valt)


        if (risqcom(ke).eq.2.and.ng.gt.1.and.g.lt.ng) then
           risq(g)=risq(g)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
           surv(g)=surv(g)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
           survint(g)=survint(g)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
           surv0(g)=surv0(g)*exp(b1(nprob+sumnrisq+nprisq(ke)+g)) 
           survt(g)=survt(g)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
        end if

     end do
     sumnrisq = sumnrisq + nrisq(ke)
  end do




  ! do g=1,ng

  !    brisq=0.d0
  !    if (logspecif.eq.1) then
  !       if (risqcom(1).eq.0) then
  !          do k=1,nprisq(1)
  !             brisq(k)=exp(b1(nprob+nprisq(1)*(g-1)+k))
  !          end do
  !       elseif (risqcom.eq.1) then
  !          do k=1,nprisq
  !             brisq(k)=exp(b1(nprob+k))
  !          end do
  !       elseif (risqcom.eq.2) then
  !          do k=1,nprisq
  !             brisq(k)=exp(b1(nprob+k))
  !          end do
  !       end if
  !    else
  !       if (risqcom.eq.0) then
  !          do k=1,nprisq
  !             brisq(k)=b1(nprob+nprisq*(g-1)+k)*b1(nprob+nprisq*(g-1)+k)
  !          end do
  !       elseif (risqcom.eq.1) then
  !          do k=1,nprisq
  !             brisq(k)=b1(nprob+k)*b1(nprob+k)
  !          end do
  !       elseif (risqcom.eq.2) then
  !          do k=1,nprisq
  !             brisq(k)=b1(nprob+k)*b1(nprob+k)
  !          end do
  !       end if
  !    end if
  !    call fct_risq_it(i,t1,brisq,g,risq,surv,surv0,survint,survT,valt)


  !    if (risqcom.eq.2.and.ng.gt.1.and.g.lt.ng) then
  !       risq(g)=risq(g)*exp(b1(nprob+nprisq+g))
  !       surv(g)=surv(g)*exp(b1(nprob+nprisq+g))
  !       survint(g)=survint(g)*exp(b1(nprob+nprisq+g))
  !       surv0(g)=surv0(g)*exp(b1(nprob+nprisq+g))
  !       survT(g)=survT(g)*exp(b1(nprob+nprisq+g))
  !    end if
  ! end do


  ! ----------- boucle sur les individus -------------
  entretard=0.d0
  vrais_survie=0.d0
  vrais=0.d0
  vrais_long=0.d0


  !------------------  debut boucle sujet ---------------------



  ! -------- creation de Vi = ZiGZi'+se*seIni+Corr ----------
  ! creation de Zi

  Z=0.d0
  l=0
  do k=1,nv
     if (idea(k).eq.1) then
        l=l+1
        do j=1,nmes_s(i)
           Z(j,l)=dble(X_s(nmes_curr_s+j,k))
        end do
     end if

  end do
  ! creation de Se=s2*I+Corr et Y1

  Se=0.d0
  tcor=0.d0
  if (ncor.gt.0) then
     do k=1,nv
        if (idcor(k).eq.1) then
           do j=1,nmes_s(i)
              tcor(j) = X_s(nmes_curr_s+j,k)
           end do
        end if
     end do
  end if
  do j1=1,nmes_s(i)
     do j2=1,nmes_s(i)
        if (j1.eq.j2) Se(j1,j2) = b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)**2
        if (ncor.eq.1) then
           Se(j1,j2) = Se(j1,j2)+b1(nprob+nrisqtot+nvarxevt &
                       +nef+nvc+nwg+ncor)**2*min(tcor(j1),tcor(j2))
        else if (ncor.eq.2) then
           Se(j1,j2) = Se(j1,j2)+b1(nprob+nrisqtot+nvarxevt &
                       +nef+nvc+nwg+ncor)**2*exp(-b1(nprob+nrisqtot+nvarxevt &
                       +nef+nvc+nwg+1)*abs(tcor(j1)-tcor(j2)))
        end if
     end do
  end do


  Y1=0.d0
  do j=1,nmes_s(i)
     Y1(j)=dble(Y_s(nmes_curr_s + j))
  end do

  ! ! creer Xevt:
  !       Xevt=0.d0

  !          l=0
  !          do k=1,nv
  !             if (idxevt(k).eq.1.or.idxevt(k).eq.2) then
  !                l=l+1
  !                Xevt(l)=X_s(nmes_curr_s+1,k)
  !             end if
  !          end do

  ! creation de P=Zi*Ut et V=P*P' que si non spec aux classes

  IF (nwg.eq.0.OR.NG.EQ.1) then


     P=0.d0
     P=MATMUL(Z,Ut)
     VC=0.d0
     VC=MATMUL(P,transpose(P))+Se

     ! Vi en vecteur
     jj=0
     Vi=0.d0
     do j=1,nmes_s(i)
        do k=j,nmes_s(i)
           jj=j+k*(k-1)/2
           Vi(jj)=VC(j,k)
        end do
     end do
     CALL DSINV(Vi,nmes_s(i),eps,ier,det)
     if (ier.eq.-1) then
        funcpi_condt=-1.d9
        goto 654
     end if

     !     retransformation du vecteur Vi en matrice :

     VC=0.d0
     do j=1,nmes_s(i)
        do k=1,nmes_s(i)
           if (k.ge.j) then
              VC(j,k)=Vi(j+k*(k-1)/2)
           else
              VC(j,k)=Vi(k+j*(j-1)/2)
           end if
        end do
     end do
  END IF

  !     debut du calcul de la vraisemblance
  vrais=vrais-nmes_s(i)*dlog(dble(2*3.14159265))
  vrais_long=vrais_long-nmes_s(i)*dlog(dble(2*3.14159265))

  ! Cas 1 : ng=1
  IF (ng.eq.1) then


     Xevt=0.d0
     bevt=0.d0
     bevtint=0.d0
     if (nxevt.ne.0) then
        l=0
        do ke=1,nbevt
           nevtxcurr=0
           do k=1,nv 

              if (idtdv(k).ne.1) then

                 if (idcom(k).eq.1) then  
                    l=l+1
                    bevt(l)=b1(nprob+nrisqtot+nevtxcurr+1)
                    Xevt(l)=X(nmes_curr+1,k)
                 else
                    if (idspecif((ke-1)*nv+k).eq.1) then   
                       l=l+1
                       bevt(l)=b1(nprob+nrisqtot+nevtxcurr+ke)
                       Xevt(l)=X(nmes_curr+1,k)
                    end if
                 end if

              else ! i.e timedepvar

                 if (idcom(k).eq.1) then  
                    bevtint(ke)=b1(nprob+nrisqtot+nevtxcurr+1)
                 else 
                    if (idspecif((ke-1)*nv+k).eq.1) then
                       bevtint(ke)=b1(nprob+nrisqtot+nevtxcurr+ke)
                    end if
                 end if
              end if
              nevtxcurr=nevtxcurr+nevtparx(k)
           end do
        end do
     end if


     ! bevt=0.d0

     ! bevtint=0.d0
     ! l=1
     ! do k=1,nv!warning nvar
     !    if (idxevt(k).eq.1) then
     !       bevt(l)=b1(nprob+nrisq+l)
     !       l=l+1
     !    end if
     ! end do

     ! if (l-1.ne.nvarxevt-nvdepsurv) then
     !    !              write(*,*)'probleme nvarxevt'
     !    funcpi_condt=-1.d9
     !    goto 654
     ! end if

     ! if (nvdepsurv.ne.0) then
     !    bevtint(1)=b1(nprob+nrisq+nvarxevt)
     ! end if


     !--------------------------------------------------------------

     vrais=vrais-det
     vrais_long=vrais_long-det

     b0=0.d0
     l=0
     X00=0.d0
     do k=1,nv
        if (idg(k).ne.0) then
           l=l+1
           do j=1,nmes_s(i)
              X00(j,l)=dble(X_s(nmes_curr_s+j,k))
           end do
           b0(l)=b1(nprob+nrisqtot+nvarxevt+l)
        end if
     end do

     mu=0.d0
     y2=0.d0
     y3=0.d0
     y4=0.d0

     mu=matmul(X00,b0)

     Y2=Y1-mu
     Y3=matmul(VC,Y2)
     Y4=DOT_PRODUCT(Y2,Y3)
     vrais=vrais-Y4
     vrais_long=vrais_long-Y4


     !------------------------------------------------------------------
     if (ind_survint_s(i).eq.1.and.tsurvint_s(i).lt.valt) then
        vrais_long=vrais_long-2*exp(DOT_PRODUCT(Xevt,bevt))*(survint(1)+exp(bevtint(1))       &
             *(survT(1)-survint(1)))
     else
        vrais_long=vrais_long-2*exp(DOT_PRODUCT(Xevt,bevt))*(survt(1))
     end if


     if (Devt_s(i).eq.1) then
        if (ind_survint_s(i).eq.1) then
           risq(1)=risq(1)*exp(bevtint(1))
        end if
        if (risq(1).le.0.or.is_nan(risq(1))) then
           funcpi_condt=-1.d9
           goto 654
        end if

        vrais=vrais+2*(log(risq(1))+DOT_PRODUCT(Xevt,bevt) &
             -exp(DOT_PRODUCT(Xevt,bevt))*(survint(1)+     &
             exp(bevtint(1))*(surv(1)-survint(1))))


        vrais_survie=vrais_survie+2*(log(risq(1))        &
             +DOT_PRODUCT(Xevt,bevt)                        &
             -exp(DOT_PRODUCT(Xevt,bevt))*(survint(1)+  &
             exp(bevtint(1))*(surv(1)-survint(1))))


     end if
     if (Devt_s(i).eq.0) then
        vrais=vrais-2*exp(DOT_PRODUCT(Xevt,bevt)) &
             *(survint(1)+exp(bevtint(1))       &
             *(surv(1)-survint(1)))

        vrais_survie=vrais_survie                 &
             -2*exp(DOT_PRODUCT(Xevt,bevt))      &
             *(survint(1)+exp(bevtint(1))      &
             *(surv(1)-survint(1)))

     end if
     entretard=entretard-surv0(1)*exp(DOT_PRODUCT(Xevt,bevt))


  ELSE



     ! cas 2 :  ng>1  composantes
     ! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))

     if (prior(i).ne.0) then
        pi=0.d0
        pi(prior(i))=1.d0
     else
        Xprob=0.d0
        Xprob(1)=1
        l=0
        do k=1,nv
           if (idprob(k).eq.1) then
              l=l+1
              Xprob(1+l)=X_s(nmes_curr_s+1,k)
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
     retard=0.d0
     ! creation des vecteurs de variables explicatives
     l=0
     m=0
     X00=0.d0
     X2=0.d0
     do k=1,nv
        if (idg(k).eq.2) then
           l=l+1
           do j=1,nmes_s(i)
              X2(j,l)=dble(X_s(nmes_curr_s+j,k))
           end do
        else if (idg(k).eq.1) then
           m=m+1
           do j=1,nmes_s(i)
              X00(j,m)=dble(X_s(nmes_curr_s+j,k))
           end do
        end if
     end do

     b2=0.d0
     b0=0.d0
     expo=0.d0
     expo_long=0.d0

     Do g=1,ng
        !-------------------------------------------------------


        if (nxevt.ne.0) then
           bevt=0.d0
           bevtint=0.d0
           Xevt=0.d0
           l=0
           ! for each covariate, number of previous effects before ke
           nxcurr=0
           do ke=1,nbevt     
              ! number of previous effects in b1 before k    
              nevtxcurr=0
              do k=1,nv                  
                 if (idtdv(k).ne.1) then
                    if (idcom(k).eq.1) then  
                       if (idspecif(k).eq.2) then 
                          l=l+1
                          bevt(l)=b1(nprob+nrisqtot+nevtxcurr+g)
                          Xevt(l)=X(nmes_curr+1,k)
                       else 
                          l=l+1
                          bevt(l)=b1(nprob+nrisqtot+nevtxcurr+1)
                          Xevt(l)=X(nmes_curr+1,k)
                       end if

                    else 

                       if (idspecif((ke-1)*nv+k).eq.1) then   
                          l=l+1
                          bevt(l)=b1(nprob+nrisqtot+nevtxcurr+nxcurr(k)+1)
                          Xevt(l)=X(nmes_curr+1,k)
                          nxcurr(k)=nxcurr(k)+1
                       end if
                       if (idspecif((ke-1)*nv+k).eq.2) then  
                          l=l+1
                          bevt(l)=b1(nprob+nrisqtot+nevtxcurr+nxcurr(k)+g)
                          Xevt(l)=X(nmes_curr+1,k)
                          nxcurr(k)=nxcurr(k)+ng
                       end if
                    end if

                 else ! i.e timedepvar

                    if (idcom(k).eq.1) then  
                       if(idspecif(k).eq.1) then
                          bevtint(ke)=b1(nprob+nrisqtot+nevtxcurr+1)
                       end if
                       if(idspecif(k).eq.2) then
                          bevtint(ke)=b1(nprob+nrisqtot+nevtxcurr+g)
                       end if

                    else 
                       if (idspecif((ke-1)*nv+k).eq.1) then 
                          bevtint(ke)=b1(nprob+nrisqtot+nevtxcurr+nxcurr(k)+1)
                          nxcurr(k)=nxcurr(k)+1
                       end if
                       if (idspecif((ke-1)*nv+k).eq.2) then
                          bevtint(ke)=b1(nprob+nrisqtot+nevtxcurr+nxcurr(k)+g)
                          nxcurr(k)=nxcurr(k)+ng
                       end if
                    end if
                 end if

                 nevtxcurr=nevtxcurr+nevtparx(k)
              end do
           end Do

        end if



        !     write(*,*)'g=',g,'bevt=',(bevt(k),k=1,l)
        !-------------------------------------------------------

        nmoins=0
        l2=0
        m2=0
        do k=1,nv
           if (idg(k).eq.1) then
              m2=m2+1
              b0(m2)=b1(nprob+nrisqtot+nvarxevt+nmoins+1)
              nmoins=nmoins+1
           else if (idg(k).eq.2) then
              l2=l2+1
              b2(l2)=b1(nprob+nrisqtot+nvarxevt+nmoins+g)
              nmoins=nmoins+ng
           end if
        end do

        ! variance covariance si spec aux classes :
        Ut1=0.d0
        Ut1=Ut

        if (nwg.ne.0) then

           Ut1=0.d0
           if (g.eq.ng) then
              Ut1=Ut
           else
              Ut1=Ut*abs(b1(nprob+nvarxevt+nrisqtot+nef+nvc+g))
           end if


           P=0.d0
           P=MATMUL(Z,Ut1)
           VC=0.d0
           VC=MATMUL(P,transpose(P))+Se

           ! Vi en vecteur
           Vi=0.d0
           jj=0
           do j=1,nmes_s(i)
              do k=j,nmes_s(i)
                 jj=j+k*(k-1)/2
                 Vi(jj)=VC(j,k)
              end do
           end do

           CALL DSINV(Vi,nmes_s(i),eps,ier,det)
           if (ier.eq.-1) then
              !               write(*,*)'probleme dsinv'
              funcpi_condt=-1.d9
              goto 654
           end if

           !     retransformation du vecteur Vi en matrice :

           VC=0.d0
           do j=1,nmes_s(i)
              do k=1,nmes_s(i)
                 if (k.ge.j) then
                    VC(j,k)=Vi(j+k*(k-1)/2)
                 else
                    VC(j,k)=Vi(k+j*(j-1)/2)
                 end if
              end do
           end do

        end if


        mu=0.d0
        y2=0.d0
        y3=0.d0
        y4=0.d0
        mu=matmul(X00,b0)+matmul(X2,b2)
        Y2=Y1-mu
        Y3=Matmul(VC,Y2)
        Y4=DOT_PRODUCT(Y2,Y3)


        if (ind_survint_s(i).eq.1) then
           risq(g)=risq(g)*exp(bevtint(1))
        end if

        if (Devt_s(i).eq.1) then
           expo=expo+pi(g)*risq(g)*                        &
                exp(DOT_PRODUCT(Xevt,bevt)                    &
                +(-det-Y4)/2.d0-exp(DOT_PRODUCT(Xevt,bevt))*  &
                (survint(g)+exp(bevtint(1))*(surv(g)      &
                -survint(g))))
        end if
        if (Devt_s(i).eq.0) then
           expo=expo+pi(g)*exp((-det-Y4)/2.d0             &
                -exp(DOT_PRODUCT(Xevt,bevt))*                 &
                (survint(g)+exp(bevtint(1))*(surv(g)      &
                -survint(g))))
        end if

        ! dans expolong, j'ajoute la survie en T


        if (ind_survint_s(i).eq.1.and.tsurvint_s(i).lt.valt) then
           expo_long=expo_long+pi(g)*exp((-det-Y4)/2.d0-exp(DOT_PRODUCT(Xevt,bevt))* &
                (survint(g)+exp(bevtint(g))*(survT(g)-survint(g))))
        else
           expo_long=expo_long+pi(g)*exp((-det-Y4)/2.d0-exp(DOT_PRODUCT(Xevt,bevt))* &
                (survT(g)))
        end if

        retard=retard+pi(g)*exp(-surv0(g) &
             *exp(DOT_PRODUCT(Xevt,bevt)))

     end do



     if (expo.le.0.or.is_nan(expo)) then
        funcpi_condt=-1.d9
        goto 654
     end if

     entretard=entretard+log(retard)

     vrais=vrais+2*log(expo)
     vrais_long=vrais_long+2*log(expo_long)
  END IF

  if (idtrunc.eq.0) then
     entretard=0.d0
  end if


  if (ng.eq.1) then
     vrais_survie=vrais_survie/2.d0-entretard
  end if

  !      if (vrais.lt.-1.d9) then
  !         write(*,*)'vrais inf -10^9',vrais
  !         funcpi_condt=-1.d9
  !         goto 654
  !       end if


  ! FIN BOUCLE SUJET


  funcpi_condt=vrais/2.D0-entretard - vrais_long/2.D0



654 continue


  !      if (is_nan(funcpi_condt).or.abs(funcpi_condt).gt.1.d30) then
  if (is_nan(funcpi_condt)) then
     funcpi_condt=-1.d9
  end if


  return

end function funcpi_condt

!---------------------------------------------------------------
!
!SPLINES
!
!---------------------------------------------------------------


subroutine splines_cvpl(nt,valt)

  use commun_comp
  use commun_cvpl

  implicit none

  integer::i,k,n,l,nt
  double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n, &
       hn,hh2,hh3
  double precision,dimension(nt)::valt


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

  Tmm_valt=0.d0
  Tmm1_valt=0.d0
  Tmm2_valt=0.d0
  Tmm3_valt=0.d0

  Tim_valt=0.d0
  Tim1_valt=0.d0
  Tim2_valt=0.d0
  Tim3_valt=0.d0



  zi(-2,1)=zi(1,1)
  zi(-1,1)=zi(1,1)
  zi(0,1)=zi(1,1)
  zi(nz(1)+1,1)=zi(nz(1),1)
  zi(nz(1)+2,1)=zi(nz(1),1)
  zi(nz(1)+3,1)=zi(nz(1),1)
  n=nz(1)+2

  !------------------- Tsurv ---------------------------
  Do i=1,ns
     l=0
     do k=2,n-2
        if ((Tsurv(i).ge.zi(k-1,1)).and.  &
             Tsurv(i).lt.zi(k,1)) then
           l=k-1
        end if
     end do

     if (Tsurv(i).eq.zi(n-2,1)) then
        l=n-3
     end if

     ht = Tsurv(i)-zi(l,1)
     htm = Tsurv(i)-zi(l-1,1)
     h2t = Tsurv(i)-zi(l+2,1)
     ht2 = zi(l+1,1)-Tsurv(i)
     ht3 = zi(l+3,1)-Tsurv(i)
     hht = Tsurv(i)-zi(l-2,1)
     h = zi(l+1,1)-zi(l,1)
     hh = zi(l+1,1)-zi(l-1,1)
     h2 = zi(l+2,1)-zi(l,1)
     h3 = zi(l+3,1)-zi(l,1)
     h4 = zi(l+4,1)-zi(l,1)
     h3m = zi(l+3,1)-zi(l-1,1)
     h2n = zi(l+2,1)-zi(l-1,1)
     hn = zi(l+1,1)-zi(l-2,1)
     hh3 = zi(l+1,1)-zi(l-3,1)
     hh2 = zi(l+2,1)-zi(l-2,1)

     if (Tsurv(i).ne.zi(n-2,1)) then
        Tmm3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        Tmm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
             +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))  &
             +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        Tmm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
             +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))   &
             +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        Tmm(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
     end if

     if (Tsurv(i).eq.zi(n-2,1)) then
        Tmm3(i) = 0.d0
        Tmm2(i) = 0.d0
        Tmm1(i) = 0.d0
        Tmm(i) = 4.d0/h
     end if

     Tim3(i) = (0.25d0*(Tsurv(i)-zi(l-3,1))*Tmm3(i)) &
          +(0.25d0*hh2*Tmm2(i))        &
          +(0.25d0*h3m*Tmm1(i))+(0.25d0*h4*Tmm(i))
     Tim2(i) = (0.25d0*hht*Tmm2(i))  &
          +(h3m*Tmm1(i)*0.25d0)+(h4*Tmm(i)*0.25d0)
     Tim1(i) = (htm*Tmm1(i)*0.25d0)+(h4*Tmm(i)*0.25d0)
     Tim(i) = ht*Tmm(i)*0.25d0

     !------------------- Tsurv0 --------------------------

     if (idtrunc.eq.1) then
        l=0
        do k=2,n-2      
           if ((Tsurv0(i).ge.zi(k-1,1)).and.   &
                Tsurv0(i).lt.zi(k,1)) then
              l=k-1
           end if
        end do

        if (Tsurv0(i).eq.zi(n-2,1)) then
           l=n-3
        end if

        ht = Tsurv0(i)-zi(l,1)
        htm = Tsurv0(i)-zi(l-1,1)
        h2t = Tsurv0(i)-zi(l+2,1)
        ht2 = zi(l+1,1)-Tsurv0(i)
        ht3 = zi(l+3,1)-Tsurv0(i)
        hht = Tsurv0(i)-zi(l-2,1)
        h = zi(l+1,1)-zi(l,1)
        hh = zi(l+1,1)-zi(l-1,1)
        h2 = zi(l+2,1)-zi(l,1)
        h3 = zi(l+3,1)-zi(l,1)
        h4 = zi(l+4,1)-zi(l,1)
        h3m = zi(l+3,1)-zi(l-1,1)
        h2n = zi(l+2,1)-zi(l-1,1)
        hn = zi(l+1,1)-zi(l-2,1)
        hh3 = zi(l+1,1)-zi(l-3,1)
        hh2 = zi(l+2,1)-zi(l-2,1)

        if (Tsurv0(i).ne.zi(nz(1)-2,1)) then

           Tmm03(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))

           Tmm02(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))   &
                +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmm01(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))    &
                +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmm0(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

        end if

        if (Tsurv0(i).eq.zi(n-2,1)) then

           Tmm03(i) = 0.d0
           Tmm02(i) = 0.d0
           Tmm01(i) = 0.d0
           Tmm0(i) = 4.d0/h

        end if

        Tim03(i) = (0.25d0*(Tsurv0(i)-zi(l-3,1))*Tmm03(i))  &
             +(0.25d0*hh2*Tmm02(i))           &
             +(0.25d0*h3m*Tmm01(i))+(0.25d0*h4*Tmm0(i))
        Tim02(i) = (0.25d0*hht*Tmm02(i))                  &
             +(h3m*Tmm01(i)*0.25d0)+(h4*Tmm0(i)*0.25d0)
        Tim01(i) = (htm*Tmm01(i)*0.25d0)+(h4*Tmm0(i)*0.25d0)
        Tim0(i) = ht*Tmm0(i)*0.25d0

     end if


     !------------------- Tsurvint --------------------------
     if (ind_survint(i).eq.1) then
        l=0
        do k=2,n-2
           if ((Tsurvint(i).ge.zi(k-1,1)).and. &
                Tsurvint(i).lt.zi(k,1)) then
              l=k-1
           end if
        end do

        if (Tsurvint(i).eq.zi(nz(1)-2,1)) then
           l=n-3
        end if

        ht = Tsurvint(i)-zi(l,1)
        htm = Tsurvint(i)-zi(l-1,1)
        h2t = Tsurvint(i)-zi(l+2,1)
        ht2 = zi(l+1,1)-Tsurvint(i)
        ht3 = zi(l+3,1)-Tsurvint(i)
        hht = Tsurvint(i)-zi(l-2,1)
        h = zi(l+1,1)-zi(l,1)
        hh = zi(l+1,1)-zi(l-1,1)
        h2 = zi(l+2,1)-zi(l,1)
        h3 = zi(l+3,1)-zi(l,1)
        h4 = zi(l+4,1)-zi(l,1)
        h3m = zi(l+3,1)-zi(l-1,1)
        h2n = zi(l+2,1)-zi(l-1,1)
        hn = zi(l+1,1)-zi(l-2,1)
        hh3 = zi(l+1,1)-zi(l-3,1)
        hh2 = zi(l+2,1)-zi(l-2,1)

        if (Tsurvint(i).ne.zi(nz(1)-2,1)) then
           Tmmt3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
           Tmmt2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn)) &
                +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))     &
                +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmmt1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))       &
                +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmmt(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
        end if

        if (Tsurvint(i).eq.zi(nz(1)-2,1)) then
           Tmmt3(i) = 0.d0
           Tmmt2(i) = 0.d0
           Tmmt1(i) = 0.d0
           Tmmt(i) = 4.d0/h
        end if

        Timt3(i) = (0.25d0*(Tsurvint(i)-zi(l-3,1))*Tmmt3(i)) &
             +(0.25d0*hh2*Tmmt2(i))               &
             +(0.25d0*h3m*Tmmt1(i))+(0.25d0*h4*Tmmt(i))
        Timt2(i) = (0.25d0*hht*Tmmt2(i))                   &
             +(h3m*Tmmt1(i)*0.25d0)+(h4*Tmmt(i)*0.25d0)
        Timt1(i) = (htm*Tmmt1(i)*0.25d0)+(h4*Tmmt(i)*0.25d0)
        Timt(i) = ht*Tmmt(i)*0.25d0
     else
        Timt3(i) =Tim3(i)
        Timt2(i) =Tim2(i)
        Timt1(i) =Tim1(i)
        Timt(i) =Tim(i)
     end if
  End Do


  Do i=1,nt
     l=0
     do k=2,n-2 
        if ((ValT(i).ge.zi(k-1,1)).and.  &
             ValT(i).lt.zi(k,1)) then
           l=k-1
        end if
     end do

     if (ValT(i).eq.zi(n-2,1)) then
        l=n-3
     end if

     ht = ValT(i)-zi(l,1)
     htm = ValT(i)-zi(l-1,1)
     h2t = ValT(i)-zi(l+2,1)
     ht2 = zi(l+1,1)-ValT(i)
     ht3 = zi(l+3,1)-ValT(i)
     hht = ValT(i)-zi(l-2,1)
     h = zi(l+1,1)-zi(l,1)
     hh = zi(l+1,1)-zi(l-1,1)
     h2 = zi(l+2,1)-zi(l,1)
     h3 = zi(l+3,1)-zi(l,1)
     h4 = zi(l+4,1)-zi(l,1)
     h3m = zi(l+3,1)-zi(l-1,1)
     h2n = zi(l+2,1)-zi(l-1,1)
     hn = zi(l+1,1)-zi(l-2,1)
     hh3 = zi(l+1,1)-zi(l-3,1)
     hh2 = zi(l+2,1)-zi(l-2,1)

     if (ValT(i).ne.zi(n-2,1)) then
        Tmm3_valt(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        Tmm2_valt(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
             +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))  &
             +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        Tmm1_valt(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
             +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))   &
             +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        Tmm_valt(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)
     end if

     if (ValT(i).eq.zi(n-2,1)) then
        Tmm3_valt(i) = 0.d0
        Tmm2_valt(i) = 0.d0
        Tmm1_valt(i) = 0.d0
        Tmm_valt(i) = 4.d0/h
     end if

     Tim3_valt(i)=(0.25d0*(ValT(i)-zi(l-3,1))*Tmm3_valt(i)) &
          +(0.25d0*hh2*Tmm2_valt(i))+(0.25d0*h3m*Tmm1_valt(i))+(0.25d0*h4*Tmm_valt(i))
     Tim2_valt(i)=(0.25d0*hht*Tmm2_valt(i))+(h3m*Tmm1_valt(i)*0.25d0)+(h4*Tmm_valt(i)*0.25d0)
     Tim1_valt(i)=(htm*Tmm1_valt(i)*0.25d0)+(h4*Tmm_valt(i)*0.25d0)
     Tim_valt(i)=ht*Tmm_valt(i)*0.25d0
  end do

end subroutine splines_cvpl

!---------------------------------------------------------------
!
!FCT_RISQ_IT
!
!---------------------------------------------------------------


subroutine fct_risq_it(i,t1,brisq,g,risq,surv,surv0,survint,survT,valT)

  use commun_cvpl
  use commun_comp

  implicit none
  ! risq contient le risque instantane d'evenement
  ! surv contient le risq cumule d'evenement et non la fonction de survie directe


  double precision, dimension(nprisq(1))::brisq
  integer ::i,g,ll,kk,ii,l,j,t1
  double precision,dimension(ng)::risq,surv,surv0,survint,survT
  double precision::som,valT



  !   do i=1,ns
!change CPL 20022015

  if (typrisq(1).eq.2.and.logspecif.eq.1) then

     surv(g)=brisq(1)*((tsurv_s(i)))**brisq(2)

     risq(g)=brisq(1)*brisq(2) &
          *(Tsurv_s(i))**(brisq(2)-1)
     if (idtrunc.eq.1) then
        surv0(g)=brisq(1)*((tsurv0_s(i)))**brisq(2)
     end if
     if (ind_survint_s(i).eq.1) then
        survint(g)=brisq(1)*((tsurvint_s(i))) &
             **brisq(2)
     else
        survint(g)=surv(g)
     end if

     survT(g)=brisq(1)*(valt)**brisq(2)

  end if
  if (typrisq(1).eq.2.and.logspecif.eq.0) then

     surv(g)=(brisq(1)*tsurv_s(i))**brisq(2)

     risq(g)=brisq(1)*brisq(2) &
          *(brisq(1)*Tsurv_s(i))**(brisq(2)-1)
     if (idtrunc.eq.1) then
        surv0(g)=(brisq(1)*tsurv0_s(i))**brisq(2)
     end if
     if (ind_survint_s(i).eq.1) then
        survint(g)=(brisq(1)*tsurvint_s(i)) &
             **brisq(2)
     else
        survint(g)=surv(g)
     end if

     survT(g)=(brisq(1)*valt)**brisq(2)

  end if
! end change CPL 20022015 
  if (typrisq(1).eq.1) then
     do j=1,nz(1)-1
        som=0.d0
        do l=1,j-1
           som=som+brisq(l)*(zi(l+1,1)-zi(l,1))
        end do

        if (idtrunc.eq.1) then
           if (Tsurv0_s(i).ge.zi(j,1).and.Tsurv0_s(i).le.zi(j+1,1)) then
              surv0(g)=som+brisq(j)*(Tsurv0_s(i)-zi(j,1))
           end if
        end if

        if (Tsurv_s(i).ge.zi(j,1).and.Tsurv_s(i).le.zi(j+1,1)) then
           surv(g)=som+brisq(j)*(Tsurv_s(i)-zi(j,1))
           risq(g)=brisq(j)
        end if

        if (ind_survint_s(i).eq.1) then
           if (Tsurvint_s(i).ge.zi(j,1).and.Tsurvint_s(i).le.zi(j+1,1)) then
              survint(g)=som+brisq(j)*(Tsurvint_s(i)-zi(j,1))
           end if
        end if

        if (valT.ge.zi(j,1).and.valT.le.zi(j+1,1)) then
           survT(g)=som+brisq(j)*(valT-zi(j,1))
        end if

     end do
     if (ind_survint_s(i).eq.0) then
        survint(g)=surv(g)
     end if
  end if

  if (typrisq(1).eq.3) then
     !------------ survie et risq pour Tsurv ----------------
     ll=0
     if (Tsurv_s(i).eq.zi(nz(1),1)) then
        ll=nz(1)-1
     end if
     som=0.d0

     do kk=2,nz(1)
        if ((Tsurv_s(i).ge.zi(kk-1,1)).and.(Tsurv_s(i).lt.zi(kk,1))) then
           ll=kk-1
        end if
     end do

     if (ll.gt.1) then
        do ii=1,ll-1
           som=som+brisq(ii)
        end do
     end if

     surv(g)=som+brisq(ll)*Tim3(i)+brisq(ll+1)*Tim2(i) &
          +brisq(ll+2)*Tim1(i)+brisq(ll+3)*Tim(i)
     risq(g)=brisq(ll)*Tmm3(i)+brisq(ll+1)*Tmm2(i)     &
          +brisq(ll+2)*Tmm1(i)+brisq(ll+3)*Tmm(i)

     !----------------------- survT ---------------------------
     ll=0
     if (valT.eq.zi(nz(1),1)) then
        ll=nz(1)-1
     end if
     som=0.d0

     do kk=2,nz(1)
        if ((valT.ge.zi(kk-1,1)).and.(valT.lt.zi(kk,1))) then
           ll=kk-1
        end if
     end do

     if (ll.gt.1) then
        do ii=1,ll-1
           som=som+brisq(ii)
        end do
     end if

     survT(g)=som+brisq(ll)*Tim3_valT(t1)+brisq(ll+1)*Tim2_valT(t1) &
          +brisq(ll+2)*Tim1_valT(t1)+brisq(ll+3)*Tim_valT(t1)


     !------------ survie et risq pour Tsurv0 ----------------

     if (idtrunc.eq.1) then
        ll=0
        if (Tsurv0_s(i).eq.zi(nz(1),1)) then
           ll=nz(1)-1
        end if
        som=0.d0
        do kk=2,nz(1)
           if ((Tsurv0_s(i).ge.zi(kk-1,1)).and.(Tsurv0_s(i).lt.zi(kk,1))) &
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

        surv0(g)=som+brisq(ll)*Tim03(i)+brisq(ll+1)*Tim02(i) &
             +brisq(ll+2)*Tim01(i)+brisq(ll+3)*Tim0(i)

     end if

     !------------ survie et risq pour Tsurvint ----------------


     if (ind_survint_s(i).eq.1) then
        ll=0
        if (Tsurvint_s(i).eq.zi(nz(1),1)) then
           ll=nz(1)-1
        end if
        som=0.d0
        do kk=2,nz(1)
           if((Tsurvint_s(i).ge.zi(kk-1,1)).and.(Tsurvint_s(i).lt.zi(kk,1))) &
                then
              ll=kk-1
           end if
        end do
        !               if (ll.lt.1.or.ll.gt.nz-1) then
        !                  stop
        !               end if
        if (ll.gt.1) then
           do ii=1,ll-1
              som=som+brisq(ii)
           end do
        end if

        survint(g)=som+brisq(ll)*Timt3(i)+brisq(ll+1)*Timt2(i) &
             +brisq(ll+2)*Timt1(i)+brisq(ll+3)*Timt(i)

     else
        survint(g)=surv(g)
     end if

  end if

end subroutine fct_risq_it




!-------------------------------------------------------------------
! IS NAN
!-------------------------------------------------------------------

logical function is_nan(x)

  implicit none

  double precision,intent(in)::x

  if (x .ne. x) then
     is_nan=.true.
  else
     is_nan=.false.
  end if

  return

end function is_nan
