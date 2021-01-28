

    !!! signification des variables de commun_mpjord !!!

!!$  nbK : nombre de processus latents
!!$  ny : nombre de marqueurs par processus
!!$  nbevt : nombre d evenements concurrents
!!$  ns : nombre de sujet dans l'echantillon
!!$  ng : nombre de classes latentes
!!$  nv : nombre total de covariables pour les MM
!!$  nv2 : nombre total de covariables pour probas classes et survie
!!$  idiag : 1 si matrice varcov des effets aleatoires est diagonale; 0 sinon
!!$  ncssg : nb de prm pas en mixture dans modele mixte (MM)sca
!!$  nvc : nb de prm a estimer dans varcov effets aleatoires
!!$  nea : nb de variable en effet aleatoire dans MM
!!$  ncg : nb de prm en mixture dans MM
!!$  nwg : nb de prm pr proportionnalite de varcov entre classes
!!$  ncor : nb de prm pr correlation des erreurs
!!$  nprob : nb de prm dans modele pr appartenance aux classes
!!$  nvarprob : nb de variable dans modele appartenance aux classes
!!$  maxmes : nb max de mesures par sujet
!!$  nobs : nb d observations
!!$  ntr : nb de prm pr transformation
!!$  nerr : nombre d erreurs de mesures (=0 si lcmm, =ny sinon)
!!$  nalea : nombre d intercept aleatoires specifiques aux marqueurs
!!$  nvarxevt : nb total de prm pr survie (sans prm hazard)
!!$  nef : nb d effets fixes a estimer dans MM
!!$  ncontr : nb de contrastes a estimer
!!$  idlink : type de transformation (0=linear,1=beta,2=splines,3=ordinal)
!!$  typrisq : 1=piecewise 2=Weibull 3=splines                   
!!$  risqcom : 0=specifique a classe 1=commun 2=proportionnel    
!!$  idtrunc : 1 si donnees tronquees a gauche; 0 sinon          
!!$  Tsurv0 : temps d'entree si troncature; 0 sinon             
!!$  Tsurv : temps d'evenement ou temps de censure           
!!$  Tsurvint : temps traitement ou idem Tsurv                  
!!$  devt : 0 si censure; 1 si evt 1; 2 si evt 2            
!!$  ind_survint : 1 si traitement; 0 sinon (pour chaque suje
!!$  logspecif : 1 si parametrisation exp pour brisq; 0 si 
!!$  nvdepsurv : 1 si effet traitement dans modele; 0 sinon
!!$  contrainte : type de parametrisation (0=hlme, 1=lcmm, 2=multlcmm)

!!! contenu de mpjord.f0 :
!! subroutine mpjord l.
!! function vrais_mpjord_i l.
!! function vrais_mpjord l.
!! subroutine design_splines_mpjord l.
!! subroutine transfo_estimee_mpjord l.

!! a faire dans multlcmm :
!! $nbmod -> de longueur ny, mais slt rempli pour ordinal (comme nbnodes)



!! ajout ordinal :
! idlink=3
! mettre les modalites observees dans uniqueY0
! faire indiceY0
! mettre slt minY et maxY dans zitr0
! methInteg : 1=MCO 2=MCA 3=MCQ
! seqMC : sequence de pts pour QMC (halton, sobol)
! besoin de alnorm defini dans hetmixOrd
! besoin de bgos defini dans predictCont

     module commun_mpjord

      implicit none
      integer,save::nbK,ns,ng,nv2,ncssg,ncg,  &
      nvarprob,maxmes,nobs,nvarxevt,nySPL,ntotvalSPL,nyORD,ntotvalORD, &
      nbevt,logspecif,idtrunc,nvdepsurv,nrisqtot,nxevt,npmtot
      integer,save::contraint,ntrtot,nprob,nmescur,nMC,methInteg
      integer,dimension(:),allocatable,save::nef,ncontr,nvc,nw,ncor,nerr, &
           nalea,ntr,ntrK,ny,nea,idiag,nv,needMC,anyord
      integer,dimension(:),allocatable,save::risqcom,typrisq,nz,nprisq,nrisq,nxevtspec,nevtparx,nxcurr
      double precision,dimension(:),allocatable,save::Y,uniqueY,minY,maxY,epsY,seqMC
      integer,dimension(:),allocatable,save::indiceY
      double precision,dimension(:,:),allocatable,save ::X,Xns,zi,zitr
      double precision,dimension(:),allocatable,save::Tsurv0,Tsurv,Tsurvint
      integer,dimension(:),allocatable,save::Devt,ind_survint
      integer,dimension(:),allocatable,save ::idea,idg,idcor,idcontr
      integer,dimension(:),allocatable,save::idcom,idspecif,idtdv,idprob,idlink
      integer,dimension(:,:),allocatable,save::nmes,nmesparK
      integer,dimension(:),allocatable,save::prior,nvalSPL,nvalORD
      double precision,dimension(:),allocatable,save::mm,mm1,mm2,im,im1,im2
      double precision,dimension(:),allocatable,save::Tmm,Tmm1,&
           Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,Tmm01,Tmm02,Tmm03,Tim0, &
           Tim01,Tim02,Tim03,Tmmt,Tmmt1,Tmmt2,Tmmt3,Timt,Timt1,&
           Timt2,Timt3
      double precision,dimension(:),allocatable,save::Tmm_est,Tmm1_est &
          ,Tmm2_est,Tmm3_est,Tim_est,Tim1_est,Tim2_est,Tim3_est
      double precision,dimension(:),allocatable,save::brisq_est
      !double precision,save::vrais_surv
      integer,dimension(:),allocatable,save::fix
      double precision,dimension(:),allocatable,save::bfix
      !double precision::vrais_S,vrais_Y1,vrais_Y2,jac1,jac2,vrais_Y1_G1,vrais_Y2_G1,vrais_Y1_G2,vrais_Y2_G2
    end module commun_mpjord

  


    subroutine mpjord(K0,ny0,nbevt0,ng0,ns0,Y0,nobs0,X0,nv0,Xns0,nv20, &
         Prior0,Tentr0,Tevt0,Devt0,ind_survint0, &
         idnv0,idnv20,idspecif0,idlink0, &
         epsY0,nbzitr0,zitr0,uniqueY0,nvalSPLORD0,indiceY0, &
         haz0,nz0,zi0, &
         nmes0,nea0,nw0,ncor0,nalea0,idiag0,int4, &
         npmtot0,btot,Vopt,vrais,ni,istop,gconv,ppi0,ppitest0,resid_m, &
         resid_ss,pred_m_g,pred_ss_g,pred_RE,pred_RE_Y,convBLG, &
         time,risq_est,risqcum_est,marker,transfY,nsim,Yobs,statscoretest, &
         pbH0,fix0,MC,maxdim0,seqMC0)
      

      use parameters
      use commun_mpjord
      use optim
      !use sobol

      IMPLICIT NONE

      !Declaration des variables en entree
      integer,intent(in)::K0,nbevt0,ns0,ng0,nv20,nobs0,npmtot0,nsim,maxdim0
      integer,dimension(K0)::ny0,nv0,nea0,nw0,ncor0,nalea0,idiag0
      integer,dimension(ns0,sum(ny0(:)))::nmes0   
      integer,dimension(4),intent(in)::int4 ! int4=c(idtrunc0,logspecif0,maxiter0,contrainte0)
      double precision,dimension(nobs0),intent(in)::Y0
      integer,dimension(nobs0),intent(in)::indiceY0
      double precision,dimension(nobs0*sum(nv0(:))),intent(in)::X0
      double precision,dimension(ns0*nv20),intent(in)::Xns0
      double precision, dimension(ns0),intent(in)::Tentr0,Tevt0
      integer, dimension(ns0),intent(in)::prior0,Devt0,ind_survint0
      integer,dimension(4*sum(nv0(:))),intent(in)::idnv0
      integer,dimension(3*nv20),intent(in)::idnv20
      integer, dimension(nv20*nbevt0),intent(in)::idspecif0
      integer, dimension(sum(ny0(:))),intent(in)::idlink0,nbzitr0,nvalSPLORD0
      double precision, dimension(sum(ny0(:))),intent(in)::epsY0
      double precision,dimension(maxval(nbzitr0),sum(ny0(:))),intent(in)::zitr0
      double precision,dimension(sum(nvalSPLORD0(:))),intent(in)::uniqueY0 
      integer,dimension(2*nbevt0),intent(in)::haz0 !haz0=c(typrisq0,risqcom0)
      integer,dimension(nbevt0),intent(in)::nz0
      double precision,dimension(maxval(nz0),nbevt0),intent(in)::zi0
      double precision,dimension(3),intent(in)::convBLG
      integer,dimension(npmtot0),intent(in)::pbH0,fix0
      integer,dimension(2),intent(in)::MC
      double precision,dimension(maxdim0*MC(2)),intent(in)::seqMC0
      
      !Declaration des variable en entree et sortie
      double precision,dimension(npmtot0), intent(inout) :: btot

      !Declaration des variables en sortie
      double precision,intent(out)::vrais
      double precision,dimension(3),intent(out)::gconv
      double precision,dimension(ns0*ng0),intent(out)::ppi0,ppitest0
      double precision,dimension(nobs0),intent(out)::resid_m,resid_ss,Yobs
      double precision,dimension(nobs0*ng0),intent(out)::pred_m_g
      double precision,dimension(nobs0*ng0),intent(out)::pred_ss_g
      double precision,dimension(ns0*sum(nea0(:))),intent(out)::pred_RE
      double precision,dimension(ns0*sum(nalea0(:))),intent(out)::pred_RE_Y
      double precision,dimension(nsim*sum(ny0)),intent(out)::marker,transfY
      double precision,dimension(nsim),intent(out)::time
      double precision,dimension(nsim*ng0,nbevt0),intent(out)::risq_est,risqcum_est
      double precision,dimension(npmtot0*(npmtot0+1)/2),intent(out)::Vopt
      integer,intent(out)::ni,istop
      double precision,dimension(1+nbevt0),intent(out)::statscoretest

      !Variables locales
      integer::jtemp,i,g,j,npm,ier,k,ktemp,ig,ke,sumnrisq,dejaspl,nbfix,nytot
      integer::mi,k1,k2,l,neftot,nvctot,ncontrtot,nwtot,ncortot,nerrtot,naleatot
      integer::tmp,sumny,m,jj,ll,yk
      double precision::eps,ca,cb,dd
      double precision,dimension(ns0,ng0)::PPI,ppiy
      double precision,dimension(npmtot0)::mvc,b
      double precision,dimension(npmtot0*(npmtot0+3)/2)::V
      double precision,external::vrais_mpjord


      ! sorties initialisees
      ppi0=0.d0
      Vopt=0.d0
      gconv=0.d0
      pred_ss_g=0.d0
      pred_m_g=0.d0
      pred_RE=0.d0
      pred_RE_Y=0.d0
      marker=0.d0
      transfY=0.d0
      resid_m=0.d0
      resid_ss=0.d0
      vrais=0.d0
      ni=0
      Yobs=0.d0
      statscoretest=0.d0

      ! nmbre max de mesures pour 1 sujet
      allocate(nmesparK(ns0,K0))
      nmesparK=0
      maxmes=0
      do i=1,ns0
         sumny=0
         do k=1,K0
            mi=0
            do m=1,ny0(k)
               mi=mi+nmes0(i,sumny+m)
            end do
            sumny = sumny + ny0(k)

            nmesparK(i,k) = mi
            
            if (mi.gt.maxmes) then
               maxmes=mi
            end if
         end do
      end do

      ! prm pour algo estimation
      epsa=convBLG(1)
      epsb=convBLG(2)
      epsd=convBLG(3)
      maxiter=int4(3)

      ! alloc pour partie modele mixte
      contraint=int4(4)
      nbK=K0
      nytot = sum(ny0(:))
      allocate(minY(nytot),maxY(nytot),idlink(nytot),ntr(nytot),epsY(nytot))
      allocate(nef(nbK),ncontr(nbK),nea(nbK),nvc(nbK),nw(nbK),ncor(nbK),nerr(nbK),nalea(nbK),ntrK(nbK),ny(nbK),idiag(nbK),nv(nbK))
      allocate(needMC(nbK),anyord(nbK))

      ny=ny0
      epsY=epsY0
      nySPL=0 
      nyORD=0 
      do k=1,nytot
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
      
      ntotvalSPL=0
      if(nySPL.gt.0) then 
         allocate(nvalSPL(nySPL))
      else
         allocate(nvalSPL(1))
         nvalSPL(1) = 0
      end if

      ntotvalORD=0
      if(nyORD.gt.0) then 
         allocate(nvalORD(nyORD))         
      else
         allocate(nvalORD(1))
         nvalORD(1) = 0
      end if

      k1=0
      k2=0
      do k=1,nytot
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

      methInteg = MC(1)
      nMC = MC(2)
      
      if(all(idlink/=2)) then
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
         allocate(mm(ntotvalSPL+ntotvalORD),mm1(ntotvalSPL+ntotvalORD),mm2(ntotvalSPL+ntotvalORD), &
          im(ntotvalSPL+ntotvalORD),im1(ntotvalSPL+ntotvalORD),im2(ntotvalSPL+ntotvalORD))
      end if

     
      zitr=0.d0  
      k1=0
      k2=0
      ntr=0
      do k=1,nytot
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
            ntr(k)=nvalORD(k2)-1
         end if
      end do

      ntrK=0
      sumny=0
      do k=1,nbK
         do yk=1,ny(k)
            ntrK(k) = ntrK(k) + ntr(sumny+yk)
         end do
         sumny = sumny + ny(k)
      end do
    


      allocate(Y(nobs0),X(nobs0,sum(nv0(:))),Xns(ns0,nv20),prior(ns0),  &
           idea(sum(nv0(:))),idg(sum(nv0(:))),idcontr(sum(nv0(:))),&
           idcor(sum(nv0(:))), nmes(ns0,nytot),&
           idprob(nv20),idcom(nv20),idspecif(nv20*nbevt0),idtdv(nv20), &
           uniqueY(ntotvalSPL+ntotvalORD),indiceY(nobs0))
      
      ns=ns0
      ng=ng0
      nv=nv0
      nv2=nv20
      nobs=nobs0
      nbevt=nbevt0

      nmes=nmes0
      Y=0.d0
      X=0.d0
      Xns=0.d0
      idprob=0  !classmb
      idea=0    !random
      idg=0     !fixed
      idcontr=0 !contrast
      idcor=0   !cor
      idcom=0     !survival
      idtdv=0     !timedepvar
      idspecif=0  !survival

      do j=1,sum(nv(:))
         idg(j)=idnv0(j)
         idcontr(j)=idnv0(sum(nv(:))+j)
         idea(j)=idnv0(2*sum(nv(:))+j)
         idcor(j)=idnv0(3*sum(nv(:))+j)
      end do
      
      jj=0
      do j=1,nv2
         idprob(j)=idnv20(j)
         idcom(j)=idnv20(nv2+j)
         idtdv(j)=idnv20(2*nv2+j)
         if(nbevt.gt.0) then
            do ke=1,nbevt
               jj=jj+1
               idspecif(jj)=idspecif0(jj)
            end do
         end if
      end do
 
      if (ntotvalSPL+ntotvalORD.gt.0) uniqueY(1:(ntotvalSPL+ntotvalORD))=uniqueY0(1:(ntotvalSPL+ntotvalORD))
      !uniqueY=uniqueY0
      indiceY=indiceY0

      prior=prior0
      Y=Y0

      ! X en matrice
      ktemp=0
      do l=1,sum(nv(:))
         jtemp=0
         do i=1,ns
            do k=1,nytot
               do j=1,nmes(i,k)
                  ktemp=ktemp+1
                  jtemp=jtemp+1
                  X(jtemp,l)=X0(ktemp)
               end do
            end do
         end do
      end do
      
      ! Xns en matrice
      ktemp=0
      do l=1,nv2
         jtemp=0
         do i=1,ns
            ktemp=ktemp+1
            jtemp=jtemp+1
            Xns(jtemp,l)=Xns0(ktemp)
         end do
      end do
      

      !  alloc pour partie survie

      if(nbevt.gt.0) then
         allocate(Tsurv0(ns0),Tsurv(ns0),Tsurvint(ns0),ind_survint(ns0),Devt(ns0))
         allocate(risqcom(nbevt0),typrisq(nbevt0),nz(nbevt0),nprisq(nbevt0), &
              nrisq(nbevt0),nxevtspec(nbevt0),nevtparx(nv20),nxcurr(nv20))     
         
         typrisq=0
         risqcom=0
         nz=0
         do ke=1,nbevt
            typrisq(ke) = haz0(ke)
            risqcom(ke) = haz0(nbevt+ke)
            nz(ke)=nz0(ke) ! nb de noeuds pour hazard (ou 2 pr Weibull)
         end do

         idtrunc=int4(1)
         Tsurv0=Tentr0   
         Tsurv=Tevt0    
         devt=devt0    
         ind_survint=ind_survint0
         logspecif=int4(2)
         Tsurvint=Tsurv ! Tsurvint defini plus loin
         
         ! zi : contient noeuds pour hazard (ou min et max si Weibull)
         if(any(typrisq.eq.3)) then
            allocate(zi(-2:maxval(nz0)+3,nbevt))
         else
            allocate(zi(maxval(nz0),nbevt))
         end if
         
         
         ! definition Tsurvint 
         nvdepsurv=0
         if(sum(ind_survint).gt.0) then
            nvdepsurv=1
            do k=1,nv2
               if (idtdv(k).eq.1) then
                  do i=1,ns
                     Tsurvint(i)=Xns(i,k)
                  end do
               end if
            end do
         end if
         
         
         nprisq=0
         nrisq=0
         nrisqtot=0
         
         do ke=1,nbevt
            
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


         ! nvarxevt = nombre total de coef pour survie (sans prm hazard)
         nxevt=0
         nxevtspec=0
         nevtparx=0
         do j=1,nv2
            if (idtdv(j).eq.1) then 
               if(idcom(j).eq.1.and.idspecif(j).eq.1) then
                  nevtparx(j)=1
               end if
               if(idcom(j).eq.1.and.idspecif(j).eq.2) then 
                  nevtparx(j)=ng
               end if
               if(idcom(j).eq.0) then
                  do ke=1,nbevt
                     if(idspecif((ke-1)*nv2+j).eq.1) then
                        nevtparx(j)=nevtparx(j)+1
                     end if
                     if(idspecif((ke-1)*nv2+j).eq.2) then
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
                     if(idspecif((ke-1)*nv2+j).eq.1) then
                        nevtparx(j)=nevtparx(j)+1
                        nxevtspec(ke)=nxevtspec(ke)+1
                     end if
                     if(idspecif((ke-1)*nv2+j).eq.2) then
                        nevtparx(j)=nevtparx(j)+ng
                        nxevtspec(ke)=nxevtspec(ke)+1
                     end if
                  end do
               end if
            end if
            nvarxevt = sum(nevtparx)
            nxevt=sum(nxevtspec)

         end do

      else ! pas de survie
         
         nrisqtot=0
         nvarxevt=0
         allocate(nrisq(1)) ! besoin pour brisq dans vrais
         nrisq=0
         
      end if


      nvarprob=sum(idprob)
      nprob=(ng-1)*nvarprob
      
      if(ng.eq.1.and.nprob.gt.0) then
         istop=12
         go to 1236
      end if

      idiag=idiag0
      nw=nw0
      ncor=ncor0
      nalea=nalea0
      
      nef=0
      ncontr=0
      nea=0
      nvc=0
      nerr=0
      ll=0
      do k=1,nbK
         ncssg=0
         ncg=0
         do l=1,nv(k)
            if (idg(ll+l).eq.1) then
               ncssg=ncssg+1      ! nb var. sans melange
            else if (idg(ll+l).eq.2) then
               ncg=ncg+1      ! nb var. dans melange
            end if
            
            if(ng.eq.1.and.ncg.gt.0) then
               istop=12
               go to 1236
            end if

            nef(k)=ncssg+ng*ncg
            ncontr(k)=ncontr(k)+idcontr(ll+l)*(ny(k)-1)
            if(contraint.ne.0) nef(k)=nef(k)-1
            nea(k)=nea(k)+idea(ll+l)
            if(idiag(k).eq.1) then
               nvc(k)=nea(k)
            else
               nvc(k)=nea(k)*(nea(k)+1)/2
            end if
            if(contraint.eq.2) nvc(k)=nvc(k)-1
            if(contraint.ne.1) nerr(k)=ny(k)
         end do
         ll = ll+nv(k)
      end do

      
      neftot=sum(nef(:))
      ncontrtot=sum(ncontr(:))
      nvctot=sum(nvc(:))
      nwtot=sum(nw(:))
      ncortot=sum(ncor(:))
      nerrtot=sum(nerr(:))
      naleatot=sum(nalea(:))
      ntrtot=sum(ntr(:))
      
      npmtot=nprob+nrisqtot+nvarxevt+ &
           neftot+ncontrtot+nvctot+nwtot+ncortot+ &
           nerrtot+naleatot+ntrtot
      
      
      if (npmtot.ne.npmtot0) then
         print*,"npmtot",npmtot,"npmtot0",npmtot0
         print*,nprob,nrisqtot,nvarxevt,nef,ncontr,nvc,nw,ncor,nerr,nalea,ntrtot
         istop=4
         goto 1236
      end if
      


      ! creer base de splines si transfo splines
      if (any(idlink.eq.2)) then 
         call design_splines_mpjord(ier)
         if (ier.eq.-1) then          
            istop=9
            go to 1236
         end if
      end if
      
      ! creer base de splines si au moins un hazard splines
      if(nbevt.gt.0) then
         dejaspl=0
         if(any(typrisq.eq.3)) then
            allocate(Tmm(ns),Tmm1(ns),Tmm2(ns),Tmm3(ns),Tim(ns)         &
                 ,Tim1(ns),Tim2(ns),Tim3(ns),Tmm0(ns),Tmm01(ns),Tmm02(ns)  &
                 ,Tmm03(ns),Tim0(ns),Tim01(ns),Tim02(ns),Tim03(ns),          &
                 Tmmt(ns),Tmmt1(ns),Tmmt2(ns),Tmmt3(ns),Timt(ns),Timt1(ns) &
                 ,Timt2(ns),Timt3(ns))
            
            do ke=1,nbevt
               if(typrisq(ke).eq.3 .and. dejaspl.eq.0) then
                  call splines_mpj(ke)
                  dejaspl=1
               end if
            end do
            
         end if
      end if
      

      ! indicateurs calcul en ordinal et integration MC
      needMC=0
      anyord=0
      do k=1,nbK
         do m=1,ny(k)
            if(idlink(m).eq.3) then
               anyord(k)=1
               if(nea(k).gt.0 .or. ncor(k).gt.0 .or. nalea(k).gt.0) needMC(k)=1
            end if
         end do
      end do

      if(all(needMC.ne.1).or.methInteg.ne.3) then
         allocate(seqMC(1))
      else
         allocate(seqMC(maxdim0*nMC))
         seqMC = seqMC0(1:maxdim0*nMC) !sobol1(1:nMC) 
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


      ! changer prm de B par prm a estimer
      tmp=0
      sumny=0
      do k=1,nbK

         if(nvc(k).gt.0) then
            if (idiag(k).eq.1) then
               do j=1,nvc(k)
                  btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j)= &
                       dsqrt(abs(btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j)))
               end do
            else

               if(contraint.eq.2) then
                  mvc(1)=1.d0
                  do j=1,nvc(k)
                     mvc(1+j)=btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j)
                  end do
               else
                  do j=1,nvc(k)
                     mvc(j)=btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+j)
                  end do
               end if
               
               eps=1.d-20
               call DMFSD(mvc,nea(k),eps,ier)
               
               if(contraint.eq.2) then
                  do j=1,nvc(k)
                     btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j)=mvc(1+j)
                  end do
               else
                  do j=1,nvc(k)
                     btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+j)=mvc(j)
                  end do
               end if
            end if
         end if
         
         if (nw(k).gt.0) then
            do j=1,nw(k)
               btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+j)= &
                    abs(btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+j))
            end do
         end if
            
         tmp = tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+ &
              sum(ntr(sumny+1:ny(k)))
         sumny = sumny+ny(k)
      end do


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


      ! optimisation
      ca=0.d0
      cb=0.d0
      dd=0.d0
      v=0.d0
      !write(*,*)'B avant marq98',npm,(b(i),i=1,npm)

      !vrais_surv=0.d0
      !vrais_S=0.d0
      !vrais_Y1=0.d0
      !vrais_Y2=0.d0
      !vrais_Y1_G1=0.d0
      !vrais_Y2_G1=0.d0
      !vrais_Y1_G2=0.d0
      !vrais_Y2_G2=0.d0
      !jac1=0.d0
      !jac2=0.d0
      call marq98(b,npm,ni,v,vrais,ier,istop,ca,cb,dd,vrais_mpjord)

      !print*,"fin estimation"
      !print*,"vrais_S=",vrais_S
      !print*,"vrais_Y1=",vrais_Y1/2.d0
      !print*,"vrais_Y2=",vrais_Y2/2.d0
      !print*,"vrais_Y1_G1=",vrais_Y1_G1
      !print*,"vrais_Y2_G1=",vrais_Y2_G1
      !print*,"vrais_Y1_G2=",vrais_Y1_G2
      !print*,"vrais_Y2_G2=",vrais_Y2_G2
      !print*,"jac=",jac1,jac2
      ! valeur des criteres de convergence
      gconv=0.d0
      gconv(1)=ca
      gconv(2)=cb
      gconv(3)=dd
      
      ! variance des parametres
      vopt=0.d0
      vopt(1:npm*(npm+1)/2)=v(1:npm*(npm+1)/2)

      !  injecter le b estime dans btot
      k=0
      do j=1,npmtot
         if(fix0(j).eq.0) then
            k=k+1
            btot(j)=b(k)
         end if
      end do


      if(istop.eq.4 .or. istop.eq.12) then
         goto 1589
      end if

      ! calculs post-estimation

      if(istop.eq.1 .or. istop.eq.2 .or. istop.eq.3) then 

      ! proba d'appartenance aux classes
         ppi=0.d0
         ppiy=0.d0
         if(ng>1) then
           ! call postprob_mpj(btot,npmtot,ppi,ppiy)
         end if

      ! residus et predictions des effets aleatoires
         !call residuals_mpj(btot,npmtot,ppi,resid_m,pred_m_g,resid_ss,pred_ss_g,pred_RE,pred_RE_Y,Yobs)

!print*,"fin residuals"
      ! recopier pour sortie
         ig=0
         do i=1,ns
            do g=1,ng0
               ig=ig+1
               ppi0(ig)=PPI(i,g)
               ppitest0(ig)=ppiy(i,g)
            end do
         end do     
   

      ! test independance conditionnelle
         if(sum(nea(:)).gt.0) then
           ! call scoretest_mpj(btot,npmtot,statscoretest)
            statscoretest=9999.d0
         else
            statscoretest=9999.d0
         end if


      ! grille de valeur pour les risques de base
         if(nbevt.gt.0) then
            allocate(brisq_est(maxval(nprisq)))
            
            dejaspl=0
            if(any(typrisq.eq.3)) then
               allocate(Tmm_est(nsim),Tmm1_est(nsim),Tmm2_est(nsim),Tmm3_est(nsim), &
                    Tim_est(nsim),Tim1_est(nsim),Tim2_est(nsim),Tim3_est(nsim))
               
               do ke=1,nbevt
                  if(typrisq(ke).eq.3 .and. dejaspl.eq.0) then
                     call splines_estime_mpj(time,nsim,ke)
                     dejaspl=1
                  endif
               end do
            end if
            
            sumnrisq=0
            do ke=1,nbevt
               
               brisq_est=0.d0
               
               do g=1,ng
                  
                  if (logspecif.eq.1) then
                     if (risqcom(ke).eq.0) then
                        do k=1,nprisq(ke)
                           brisq_est(k)=exp(btot(nprob+sumnrisq+nprisq(ke)*(g-1)+k))
                           if(g.eq.ng.and.k.eq.nprisq(ke))  sumnrisq = sumnrisq+nrisq(ke)
                        end do
                     elseif(risqcom(ke).eq.1) then
                        do k=1,nprisq(ke)
                           brisq_est(k)=exp(btot(nprob+sumnrisq+k))
                           if(g.eq.ng.and.k.eq.nprisq(ke))  sumnrisq = sumnrisq+nrisq(ke)
                        end do
                     elseif (risqcom(ke).eq.2) then
                        do k=1,nprisq(ke)
                           brisq_est(k)=exp(btot(nprob+sumnrisq+k))
                           if(g.eq.ng.and.k.eq.nprisq(ke))  sumnrisq = sumnrisq+nrisq(ke)
                        end do
                     end if
                     
                  else
                     if (risqcom(ke).eq.0) then
                        do k=1,nprisq(ke)
                           brisq_est(k)=btot(nprob+sumnrisq+nprisq(ke)*(g-1)+k) &
                                *btot(nprob+sumnrisq+nprisq(ke)*(g-1)+k)
                           if(g.eq.ng.and.k.eq.nprisq(ke))  sumnrisq = sumnrisq+nrisq(ke)
                        end do
                     elseif(risqcom(ke).eq.1) then
                        do k=1,nprisq(ke)
                           brisq_est(k)=btot(nprob+sumnrisq+k)*btot(nprob+sumnrisq+k)
                           if(g.eq.ng.and.k.eq.nprisq(ke))  sumnrisq = sumnrisq+nrisq(ke)
                        end do
                     elseif (risqcom(ke).eq.2) then
                        do k=1,nprisq(ke)
                           brisq_est(k)=btot(nprob+sumnrisq+k)*btot(nprob+sumnrisq+k)
                           if(g.eq.ng.and.k.eq.nprisq(ke))  sumnrisq = sumnrisq+nrisq(ke)
                        end do
                     end if
                  end if
                  
                  call fct_risq_estime_mpj(ke,brisq_est,time,nsim,g,risq_est,risqcum_est)
                  
                  if (risqcom(ke).eq.2.and.ng.gt.1.and.g.lt.ng) then
                     do i=1,nsim
                        risq_est((g-1)*nsim+i,ke)=risq_est((g-1)*nsim+i,ke)* &
                             exp(btot(nprob+sumnrisq+nprisq(ke)+g))
                        risqcum_est((g-1)*nsim+i,ke)=risqcum_est((g-1)*nsim+i,ke)* & 
                             exp(btot(nprob+sumnrisq+nprisq(ke)+g))
                     end do
                  end if
                  
               end do
               
            end do
            
            deallocate(brisq_est)
            if(any(typrisq.eq.3)) then 
               deallocate(Tmm_est,Tmm1_est,Tmm2_est,Tmm3_est,Tim_est &
                    ,Tim1_est,Tim2_est,Tim3_est)
            end if
         end if

!print*,"fin risq"

         ! grille de valeurs pour la transfo
         if(any(idlink.ne.-1)) then
            call transfo_estimee_mpjord(btot,npmtot,nsim,marker,transfY)
         end if
         
!print*,"fin transfo"



      end if

 1589 continue



!--------------- liberation de la memoire ------------------------


      if(nbevt.gt.0) then
         if (any(typrisq.eq.3)) then
            deallocate(Tmm,Tmm1,Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,    &
                 Tmm01,Tmm02,Tmm03,Tim0,Tim01,Tim02,Tim03,Tmmt,Tmmt1,     &
                 Tmmt2,Tmmt3,Timt,Timt1,Timt2,Timt3)
         endif
      end if
      
      deallocate(pbH,bfix,fix,seqMC)
      
1236  continue
      
      deallocate(nmesparK,minY,maxY,idlink,ntr,ntrK,epsY,nef,ncontr,nea,nvc,nw,ncor,nerr,&
           nalea,ny,idiag,nvalSPL,nvalORD,nv)
      
      deallocate(Y,X,Xns,idprob,idea,idg,idcontr,idcor,nmes,prior,zitr,mm,mm1,mm2,im,im1,&
           im2,indiceY,uniqueY,idcom,idspecif,idtdv,needMC,anyord)
      
      
      if(nbevt.gt.0) then
         deallocate(Tsurv0,Tsurv,Tsurvint,ind_survint,zi,devt,risqcom,typrisq, &
              nz,nprisq,nrisq,nxevtspec,nevtparx,nxcurr)
      else
         deallocate(nrisq)
      end if
      

      return


end subroutine mpjord



! ===================================================================
! ===================  vraisemblance sujet i  =======================
! ===================================================================

double precision function vrais_mpjord_i(b,npm,id,thi,jd,thj,i)
  
  use parameters
  use commun_mpjord
  use optim

  implicit none

  integer::i,id,jd,npm
  double precision::thi,thj
  double precision,dimension(npm)::b

  integer::k,l,m,sumnrisq,ke,g,j,yk,nevtxcurr,nxevtcurr,sumMesYk,sumntr,numSPL,l1
  integer::sumny,tmp,kk,ll,ii,j1,j2,jj,ier,nmoins,q,sumparK,sumnv,sumnea,verb,npts
  double precision::eps,temp,fevt,entretard,varexpsurv,surv_glob,surv0_glob,jacobien
  double precision,dimension(npmtot)::b1
  double precision,dimension(ng)::pi
  double precision,dimension(nvarprob)::Xprob,bprob 
  double precision,dimension(maxval(nrisq))::brisq
  double precision,dimension(ng,nbevt)::risq,surv,surv0,survint
  double precision,dimension(nxevt)::Xevt,bevt
  double precision,dimension(nbevt)::bevtint
  double precision,dimension(-1:maxval(ntr)-3)::splaa
  double precision::aa1,bb1,dd1,aa,bb,betai,cc1,ytemp,som,eta0,vrais_Y,det,Y4
  double precision::vrais_surv,expo,retard,beta_densite
  double precision,dimension(sum(nea),sum(nea))::Utot
  double precision,dimension(maxval(nea),maxval(nea))::Ut
  double precision,dimension(maxmes,maxval(nea))::Z,P
  double precision,dimension(maxmes,maxmes)::VC
  double precision,dimension(nbK*maxmes,nbK*maxmes)::Corr,chCorr !??
  double precision,dimension(nbK*maxmes)::tcor
  double precision,dimension(maxmes)::mu,Y2,Y3
  double precision,dimension(nbK*maxmes)::Y1
  double precision,dimension(nbK*maxmes*(nbK*maxmes+1)/2)::Vi
  double precision,dimension(maxmes,(maxval(ncontr)+sum(idcontr)))::X01
  double precision,dimension(maxval(ncontr)+sum(idcontr))::b01
  double precision,dimension(maxmes,maxval(nv))::X00,X2
  double precision,dimension(maxval(nv))::b0,b2
  double precision,dimension(maxval(nea))::ui,usim
 !  double precision,dimension(sum(nea))::usim
  double precision,dimension(maxmes)::wi,wsim
 !  double precision,dimension(nbK*maxmes)::wsim
 !  double precision,dimension(sum(ny))::asim
  double precision::ai,binf,bsup,asim
  double precision::SX,x22,vrais_K,vrais_l,div
  double precision,external::alnorm
  
  vrais_mpjord_i=0.d0

  verb=0
  !if(i.eq.493.or.i.eq.495) verb=1
  !if(i.lt.3) verb=1

  if(verb.eq.1) print*,"vrais pour i=",i

  
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

  
  if (nbevt.gt.0) then
     ! calcul de brisq et risq, surv et surv0 pour chaque evt et tous les g
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

           call fct_risq_mpj_i(i,ke,brisq,g,risq,surv,surv0,survint)


           if (risqcom(ke).eq.2.and.ng.gt.1.and.g.lt.ng) then
              risq(g,ke)=risq(g,ke)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
              surv(g,ke)=surv(g,ke)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
              survint(g,ke)=survint(g,ke)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
              surv0(g,ke)=surv0(g,ke)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
           end if


        end do
        sumnrisq = sumnrisq + nrisq(ke)
     end do

  end if
 ! print*,"risq ok"
  !! transformation des Y
  Y1=0.d0
  jacobien=0.d0
  
  splaa=0.d0
  numSPL=0
  tmp=0
  sumny=0
  sumparK=0
  do k=1,nbK
     sumntr=0
     sumMesYk=0
     do yk=1,ny(k)
        
        if (idlink(sumny+yk).eq.0) then  ! Linear link

           do j=1,nmes(i,sumny+yk)
              Y1(sumparK+sumMesYk+j)=(dble(Y(nmescur+sumMesYk+sumparK+j))-b1(nprob+nrisqtot+nvarxevt&
                   +tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+ &
                   sumntr+1)) &
                   /abs(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+ &
                   nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr+2))
              jacobien = jacobien - log(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ &
                   ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr+2))
           end do

        else if (idlink(sumny+yk).eq.1) then  ! Beta link
           aa1=exp(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                ncor(k)+nerr(k)+nalea(k)+sumntr+1))/ &
                (1+exp(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+&
                ncor(k)+nerr(k)+nalea(k)+sumntr+1)))
           
           bb1=exp(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                ncor(k)+nerr(k)+nalea(k)+sumntr+2))/ &
                (1+exp(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+&
                ncor(k)+nerr(k)+nalea(k)+sumntr+2)))

           cc1=abs(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                ncor(k)+nerr(k)+nalea(k)+sumntr+3))

           dd1=abs(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                ncor(k)+nerr(k)+nalea(k)+sumntr+4))

           bb1=aa1*(1.d0-aa1)*bb1
           aa=aa1*aa1*(1-aa1)/bb1-aa1
           bb=aa*(1-aa1)/aa1

           do j=1,nmes(i,sumny+yk)

              ytemp=(dble(Y(nmescur+sumMesYk+sumparK+j))-minY(sumny+yk)+epsY(sumny+yk)) &
                   /(maxY(sumny+yk)-minY(sumny+yk)+2*epsY(sumny+yk))
              Y1(sumMesYk+sumparK+j)=(betai(aa,bb,ytemp)-cc1)/dd1

              if (abs(Y1(sumMesYk+sumparK+j)-999.d0).lt.1.d-6) then
                 vrais_mpjord_i=-1.d9
                 print*,"-1.d9 Y1=999"
                 goto 654
              end if

              jacobien = jacobien + log(abs(beta_densite(ytemp,aa,bb))/dd1)
              jacobien=jacobien-log(abs(maxY(sumny+yk)-minY(sumny+yk)+2*epsY(sumny+yk)))
           end do

        else if (idlink(sumny+yk).eq.2) then ! Splines link
           numSPL=numSPL+1

           splaa=0.d0
           eta0=0.d0
           eta0=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                ncor(k)+nerr(k)+nalea(k)+sumntr+1)

           do kk=2,ntr(sumny+yk)
              splaa(kk-3)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+&
                   nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr+kk)*&
                   b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                   ncor(k)+nerr(k)+nalea(k)+sumntr+kk)
           end do
           !if(i<=5 .and. id==0 .and. jd==0) print*,"eta0=",eta0,"splaa=",sqrt(splaa)
           do j=1,nmes(i,sumny+yk)
              ll=0
              !if(i<=5 .and. id==0 .and. jd==0) print*,"Y=",Y(nmescur+sumMesYk+sumparK+j)
              if (abs(Y(nmescur+sumMesYk+sumparK+j)-zitr(ntr(sumny+yk)-2,numSPL)).lt.1.d-6) then
                 ll=ntr(sumny+yk)-3
              end if

              som=0.d0
              do kk = 2,ntr(sumny+yk)-2
                 if ((Y(nmescur+sumMesYk+sumparK+j).ge.zitr(kk-1,numSPL)).and. &
                      (Y(nmescur+sumMesYk+sumparK+j).lt.zitr(kk,numSPL))) then
                    ll=kk-1
                 end if
              end do
           
              if (ll.lt.1.or.ll.gt.ntr(sumny+yk)-3) then          
                 vrais_mpjord_i=-1.d9
                 print*,"-1.d9 ll<1 ou ll>ntrtot-3"
                 goto 654
              end if
              if (ll.gt.1) then
                 do ii=2,ll
                    som=som+splaa(ii-3)
                 end do
              end if

              Y1(sumMesYk+sumparK+j)=eta0+som +splaa(ll-2)*im2(indiceY(nmescur+sumMesYk+sumparK+j)) &
                   +splaa(ll-1)*im1(indiceY(nmescur+sumMesYk+sumparK+j))&
                   +splaa(ll)*im(indiceY(nmescur+sumMesYk+sumparK+j))

              jacobien = jacobien + log(splaa(ll-2)*mm2(indiceY(nmescur+sumMesYk+sumparK+j))&
                   +splaa(ll-1)*mm1(indiceY(nmescur+sumMesYk+sumparK+j))&
                   +splaa(ll)*mm(indiceY(nmescur+sumMesYk+sumparK+j)))
             ! if(i<=5 .and. id==0 .and. jd==0) print*,"Ytrans=",Y1(sumMesYk+sumparK+j)
           end do
        else if(idlink(sumny+yk).eq.-1.or.idlink(sumny+yk).eq.3) then
           do j=1,nmes(i,sumny+yk)
              Y1(sumMesYk+sumparK+j)=Y(nmescur+sumMesYk+sumparK+j)
           end do
        end if
        sumntr = sumntr + ntr(sumny+yk)
        sumMesYk = sumMesYk + nmes(i,sumny+yk)
     end do
     
     sumny = sumny + ny(k)
     tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr
     sumparK = sumparK + nmesparK(i,k)
  end do
  
!  print*,"transfo ok"

  
  !! matrices Ut et Corr
  Utot=0.d0
  Corr=0.d0
  sumnea=0
  sumnv=0
  tmp=0
  sumparK=0
  do k=1,nbK

     !! effets aleatoires
     if (nea(k).gt.0) then

        if(contraint.eq.2) then
           Utot(sumnea+1,sumnea+1)=1.d0
           if (nea(k).gt.1) then 

              If (idiag(k).eq.1) then
                 do j=2,nea(k)
                    do kk=2,nea(k)
                       if (j.eq.kk) then
                          Utot(sumnea+j,sumnea+kk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j-1)
                       else
                          Utot(sumnea+j,sumnea+kk)=0.d0
                       end if
                    end do
                 end do
              end if

              if (idiag(k).eq.0) then
                 do j=2,nea(k)
                    do kk=1,j
                       Utot(sumnea+j,sumnea+kk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+kk-1+j*(j-1)/2)
                    end do
                 end do
              end if

           end if

        else ! pas de variance contrainte a 1

           if (idiag(k).eq.1) then
              do j=1,nea(k)
                 do kk=1,nea(k)
                    if (j.eq.kk) then
                       Utot(sumnea+j,sumnea+kk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j)
                    else
                       Utot(sumnea+j,sumnea+kk)=0.d0
                    end if
                 end do
              end do
           end if

           if (idiag(k).eq.0) then
              do j=1,nea(k)
                 do kk=1,j
                    Utot(sumnea+j,sumnea+kk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+kk+j*(j-1)/2)
                 end do
              end do
           end if
        end if
     end if

     !! BM ou AR
     Corr=0.d0
     tcor=0.d0
     if(ncor(k).gt.0) then

        do kk=1,nv(k)
           if (idcor(sumnv+kk).eq.1) then
              do j=1,nmesparK(i,k)
                 tcor(sumparK+j) = X(nmescur+sumparK+j,sumnv+kk)
              end do
           end if
        end do

        do j1=1,nmesparK(i,k)
           do j2=1,nmesparK(i,k)
              if (ncor(k).eq.1) then 
                 Corr(sumparK+j1,sumparK+j2) = Corr(sumparK+j1,sumparK+j2)+&
                      b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k))* &
                      b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+ &
                      nw(k)+ncor(k)) * min(tcor(j1),tcor(j2))
              else if (ncor(k).eq.2) then
                 Corr(sumparK+j1,sumparK+j2) = Corr(sumparK+j1,sumparK+j2)+ &
                      b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k))* &
                      b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+ &
                      nw(k)+ncor(k))* &
                      exp(-b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                      nvc(k)+nw(k)+1) * abs(tcor(j1)-tcor(j2)))
              end if
           end do
        end do

        !! passer en cholesky si MC
        if(needMC(k).eq.1) then
           jj=0
           Vi=0.d0
           do j=1,nmesparK(i,k)
              do kk=j,nmesparK(i,k)
                 jj=j+kk*(kk-1)/2
                 Vi(jj)=Corr(sumparK+j,sumparK+kk)
              end do
           end do

           CALL DMFSD(Vi,nmesparK(i,k),EPS,IER)

           Corr=0.d0
           do j=1,nmesparK(i,k)
              do kk=1,j
                 Corr(sumparK+j,sumparK+kk)=Vi(kk+j*(j-1)/2)
              end do
           end do
        end if

     end if

     sumnea = sumnea + nea(k)
     tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+ntrK(k)
     sumparK = sumparK + nmesparK(i,k) 
  end do

!if(verb.eq.1) print*,"U=",Utot

  
  if(ng.eq.1) then

     ! !! vrais survie !!
     vrais_surv=0.d0
     entretard=0.d0
     if(nbevt.gt.0) then

        Xevt=0.d0
        bevt=0.d0
        bevtint=0.d0
        if (nxevt.ne.0) then
           l=0
           do ke=1,nbevt
              nevtxcurr=0
              do k=1,nv2 
                 
                 if (idtdv(k).ne.1) then
                    
                    if (idcom(k).eq.1) then  
                       l=l+1
                       bevt(l)=b1(nprob+nrisqtot+nevtxcurr+1)
                       Xevt(l)=Xns(i,k)
                    else
                       if (idspecif((ke-1)*nv2+k).eq.1) then   
                          l=l+1
                          bevt(l)=b1(nprob+nrisqtot+nevtxcurr+ke)
                          Xevt(l)=Xns(i,k)
                       end if
                    end if
                    
                 else ! i.e timedepvar
                    
                    if (idcom(k).eq.1) then  
                       bevtint(ke)=b1(nprob+nrisqtot+nevtxcurr+1)
                    else 
                       if (idspecif((ke-1)*nv2+k).eq.1) then
                          bevtint(ke)=b1(nprob+nrisqtot+nevtxcurr+ke)
                       end if
                    end if
                 end if
                 nevtxcurr=nevtxcurr+nevtparx(k)
              end do
           end do
           
        end if
        
        surv_glob=0.d0
        surv0_glob=0.d0
        varexpsurv=0.d0
        nxevtcurr=0
        fevt=0.d0       
        do ke=1,nbevt
           
           varexpsurv=0.d0
           if (nxevt.ne.0) then 
              varexpsurv=DOT_PRODUCT(Xevt((nxevtcurr+1):(nxevtcurr+nxevtspec(ke)))&
                   ,bevt((nxevtcurr+1):(nxevtcurr+nxevtspec(ke))))
           end if
           
           if (Devt(i).eq.ke) then
              fevt=log(risq(1,ke))+varexpsurv
              if (ind_survint(i).eq.1) then
                 fevt=fevt+bevtint(ke)
              end if
           end if
           
           Surv_glob=surv_glob+ survint(1,ke)*exp(varexpsurv)+ &
                exp(bevtint(ke)+varexpsurv)*(surv(1,ke)-survint(1,ke))
           surv0_glob=surv0_glob+surv0(1,ke)*exp(varexpsurv)
           nxevtcurr=nxevtcurr+nxevtspec(ke)
        end do
        
        vrais_surv=vrais_surv+(fevt-surv_glob)  
        entretard=entretard-surv0_glob
        
     end if ! fin survie

    !vrais_S = vrais_S+vrais_surv


   
     ! !! vrais longitudinal !!
     vrais_Y=0.d0
     tmp=0
     sumny=0
     sumnv=0


     ui=0.d0
     wi=0.d0
     sumnea=0
     sumparK=0
     sumnv=0
     sumny=0
     do k=1,nbK
        vrais_K = 0.d0

        !! variables en effet aleatoire
        Z=0.d0
        ll=0
        do kk=1,nv(k)
           if (idea(sumnv+kk).eq.1) then
             ll=ll+1
             do j=1,nmesparK(i,k)
                Z(j,ll)=dble(X(nmescur+sumparK+j,sumnv+kk))
             end do
           end if
        end do
        
        ! variables en effet fixe
        b0=0.d0
        b01=0.d0
        l1=0
        ll=0
        m=0
        X00=0.d0
        X01=0.d0
        do kk=1,nv(k)
           if (idg(sumnv+kk).ne.0) then
             l1=l1+1
             do j=1,nmesparK(i,k)
                X00(j,l1)=dble(X(nmescur+sumparK+j,sumnv+kk))  
             end do
             if (kk.eq.1.and.contraint.ne.0) then
                b0(l1)=0.d0 ! intercept fixe a 0
             else
                ll=ll+1
                b0(l1)=b1(nprob+nrisqtot+nvarxevt+tmp+ll)
             end if
           end if

           !contrast : 
           if (idcontr(sumnv+kk).ne.0) then
              m=m+1
              sumMesYk=0
              do yk=1,ny(k)
                ! creation matrice design des contrastes: X01
                do j=1,nmes(i,sumny+yk)
                   X01(sumMesYk+j,(m-1)*ny(k)+yk) = dble(X(nmescur+sumparK+sumMesYk+j,sumnv+kk))
                end do
                sumMesYk=sumMesYk+nmes(i,sumny+yk)
                ! creation vecteur parms des contrastes: b01
                if (yk<ny(k)) then
                   b01((m-1)*ny(k)+yk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+(m-1)*(ny(k)-1)+yk)
                else
                   b01((m-1)*ny(k)+ny(k)) =-sum(b1(nprob+nrisqtot+nvarxevt+tmp+ &
                        nef(k)+(m-1)*(ny(k)-1)+1 &
                        :nprob+nrisqtot+nvarxevt+tmp+nef(k)+(m-1)*(ny(k)-1)+ny(k)-1))
                end if
              end do
           end if
        end do
        
        !if(verb.eq.1) print*,"X ok", " nmescur=",nmescur," sumparK=",sumparK," j=",j," sumnv=",sumnv," kk=",kk, " b0=",b0
        if(anyord(k).eq.0) then 
           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! !!!!!!!!!!!!!!!!!!!   ng=1 et que du continu   !!!!!!!!!!!!!!!!!!!
           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           VC=0.d0
           sumMesYk = 0
           do yk=1,ny(k)
              do j1=1,nmes(i,sumny+yk)
                 if(contraint.eq.1) then
                    VC(sumMesYk+j1,sumMesYk+j1) =  Corr(sumparK+sumMesYk+j1,sumparK+sumMesYk+j1) + 1
                 else
                    VC(sumMesYk+j1,sumMesYk+j1) =  Corr(sumparK+sumMesYk+j1,sumparK+sumMesYk+j1) + &
                         b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                         nvc(k)+nw(k)+ncor(k)+yk)**2 !variance de l'erreur yk
                 end if
                 if (nalea(k).eq.ny(k)) then ! intercept aleatoire de yk
                    do j2=1,nmes(i,sumny+yk)
                       VC(sumMesYk+j1,sumMesYk+j2) = VC(sumMesYk+j1,sumMesYk+j2) + &
                            b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                            nvc(k)+nw(k)+ncor(k)+nerr(k)+yk)**2
                    end do
                 end if
              end do
              sumMesYk = sumMesYk + nmes(i,sumny+yk)
           end do

           ! Ut = sous matrice de Utot
           Ut=0.d0
           do j1=1,nea(k)
              do j2=1,j1
                 Ut(j1,j2) = Utot(sumnea+j1,sumnea+j2)
              end do
           end do

           ! matrice de variance de Yik
           P=0.d0
           P=MATMUL(Z,Ut)
           VC=MATMUL(P,transpose(P))+VC

           ! Vi en vecteur
           jj=0
           Vi=0.d0
           do j1=1,nmesparK(i,k)
              do j2=j1,nmesparK(i,k)
                 jj=j1+j2*(j2-1)/2
                 Vi(jj)=VC(j1,j2)
              end do
           end do

           ! inversion
           CALL dsinv(Vi,nmesparK(i,k),eps,ier,det) ! det=log(determinant(Vi))
           if (ier.eq.-1) then
              vrais_mpjord_i=-1.d9
              print*,"-1.d9 dsinv"
              print*,"b=",b
              !print*,"bfix=",bfix
              !print*,"id=",id,"thi=",thi
              !print*,"jd=",jd,"thj=",thj
              !print*,"fix=",fix
              goto 654
           end if

           ! retransformation du vecteur Vi en matrice :
           VC=0.d0
           do j1=1,nmesparK(i,k)
              do j2=1,nmesparK(i,k)
                if (j2.ge.j1) then
                   VC(j1,j2)=Vi(j1+j2*(j2-1)/2)
                else
                   VC(j1,j2)=Vi(j2+j1*(j1-1)/2)
                end if
              end do
           end do


           ! calcul de la vrais_Y
           mu=0.d0
           Y2=0.d0
           Y3=0.d0
           Y4=0.d0
           mu=matmul(X00,b0)+matmul(X01,b01)
           do j=1,nmesparK(i,k)
               Y2(j) = Y1(sumparK+j)-mu(j)
           end do
           Y3=matmul(VC,Y2)
           Y4=DOT_PRODUCT(Y2,Y3) 

           vrais_K = vrais_K-(nmesparK(i,k)*dlog(dble(2*3.14159265))+det+Y4)/2.d0
           !if(k.eq.1) vrais_Y1 = vrais_Y1 -nmesparK(i,k)*dlog(dble(2*3.14159265))-det-Y4
           !if(k.eq.2) vrais_Y2 = vrais_Y2 -nmesparK(i,k)*dlog(dble(2*3.14159265))-det-Y4
        
        else
           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! !!!!!!!!!!!!!!!!!!! ng=1 et au moins 1 ordinal !!!!!!!!!!!!!!!!!!!
           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           npts=needMC(k)*(nMC-1)+1
           ! boucle pour les nMC simulations du Monte Carlo
           !if(verb.eq.1)print*,"nMC=",nMC
           !print*,"npts=",npts, " needMC=",needMC(k)
           som=0.d0
           do l=1,npts
             !print*," MC l=",l
             vrais_l=1.d0
           
              if(methInteg.eq.1) then 
                 ! !!!!!!!!!!!!! MCO !!!!!!!!!!!!!
              
                 ! simuler les effets aleatoires
                 if(nea(k).gt.0) then
                    usim=0.d0
                    x22=0.d0
                    SX=1.d0
                    do j=1,nea(k)
                       call bgos(SX,0,usim(j),x22,0.d0)
                       !if(ui(j).lt.-3.or.ui(j).gt.3) print*,"usim=",ui(j)
                    end do
                  end if
              
                 ! simuler le BM ou AR
                 if(ncor(k).gt.0) then
                    wsim=0.d0
                    x22=0.d0
                    SX=1.d0
                    do j=1,nmesparK(i,k)
                       call bgos(SX,0,wsim(j),x22,0.d0)
                    end do
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
                   if(nea(k).gt.0) then
                      usim=0.d0
                      x22=0.d0
                      SX=1.d0
                      do j=1,nea(k)
                         call bgos(SX,0,usim(j),x22,0.d0)
                      end do
                   end if
                 
                   ! simuler le BM ou AR
                   if(ncor(k).gt.0) then
                      wsim=0.d0
                      x22=0.d0
                      SX=1.d0
                      do j=1,nmesparK(i,k)
                        call bgos(SX,0,wsim(j),x22,0.d0)
                      end do
                   end if
                 
                  end if
              
               else 
                  ! !!!!!!!!!!!!! QMC !!!!!!!!!!!!!
              
                   ! simuler les effets aleatoires
                   if(nea(k).gt.0) then
                      usim=0.d0
                      do j=1,nea(k)
                         usim(j)=seqMC(nMC*(j-1)+l)
                      end do
                   end if
                 
                   ! simuler le BM ou AR
                   if(ncor(k).gt.0) then
                      wsim=0.d0
                      do j=1,nmesparK(i,k)
                        wsim(j)=seqMC(nMC*(nea(k)+j-1)+l)
                      end do
                   end if
              
               end if
           

             if(verb.eq.1) print*,"MC ok"

             if(verb.eq.1) then
              !print*, "k=",k, "  anyord(k)=",anyord(k), "  needMC(k)=",needMC(k)
              !print*,"nmescur=",nmescur," sumparK=",sumparK," nmesparK=",nmesparK(i,k), " sumnv=",sumnv
             end if
          
          
             !print*,"i=",i
                            
             ! esperance conditionnelle
             mu=0.d0
             mu=matmul(X00,b0)+matmul(X01,b01)
             !if(verb.eq.1)print*,"Xb=",mu            
              
             if(nea(k).gt.0) then

                ! Ut = sous matrice de Utot
                Ut=0.d0
                do j1=1,nea(k)
                !if(verb.eq.1) print*,"usim=",usim(sumnea+j1)
                   do j2=j1,nea(k)
                      Ut(j2,j1) = Utot(sumnea+j2,sumnea+j1)
                   end do
                end do
                ui = matmul(Ut,usim)
                !if(i.lt.2.and.l.lt.1) print*,"Ut=",Ut, "  ui=",ui
                if(verb.eq.1) print*,"ui=",ui

                mu=mu+matmul(Z,ui)
               
              end if

              if(ncor(k).gt.0) mu=mu+wi ! wi a faire !!

              !if(verb.eq.1) print*,"mu=",mu

              sumMesYk=0
              sumntr=0
              q=0
              ai=0.d0
              do yk =1,ny(k)

                 if(idlink(yk).eq.3) then
                    !! yk est ordinal
                    q=q+1
                    !print*,"q-eme outcome ordinal, q=",q, "  niq=",nmes(i,sumny+yk)
                    !! on simule l'EA specifique au test
                    asim=0.d0
                    if(nalea(k).gt.0) then
                       if(methInteg.eq.1) then
                          !! MCO
                          call bgos(SX,0,asim,x22,0.d0)
                          ai = b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                               nvc(k)+nw(k)+ncor(k)+nerr(k)+yk)*asim
                       else if(methInteg.eq.2) then
                          !! MCA
                          if(mod(l,2).eq.0) then
                             ! si l est pair on prend l'oppose du precedent
                             ai = -ai
                          else
                             call bgos(SX,0,asim,x22,0.d0)
                             ai = b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                                  nvc(k)+nw(k)+ncor(k)+nerr(k)+yk)*asim
                          end if
                       else
                          !! QMC
                          asim = seqMC(nMC*(nea(k)+nmesparK(i,k))+l)
                          ai = b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                                  nvc(k)+nw(k)+ncor(k)+nerr(k)+yk)*asim
                       end if
                       do j=1,nmes(i,sumny+yk)
                         !! on ajoute ai a mu
                         mu(sumparK+sumMesYk+j) = mu(sumparK+sumMesYk+j)+ai
                       end do
                    end if
                    
                    do j=1,nmes(i,sumny+yk)
                       !! trouver binf et bsup tq binf < lambda + epsilon < bsup :
                       
                       !! initialiser au premier seuil
                       binf = b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                            nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr+1)
                       bsup = binf
                      if(verb.eq.1) print*,"Y=",Y(nmescur+sumparK+sumMesYk+j)
                       !! si Y>minY ajouter b1(..)^2
                       if(indiceY(nmescur+sumparK+sumMesYk+j).gt.1) then
                          do ll=2,min(indiceY(nmescur+sumparK+sumMesYk+j),ntr(sumny+yk))
                             bsup = bsup + b1(nprob+nrisqtot+nvarxevt+tmp+ &
                                  nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+ &
                                  nalea(k)+sumntr+ll)**2
                             if(ll.lt.indiceY(nmescur+sumparK+sumMesYk+j)) then
                                binf = binf + b1(nprob+nrisqtot+nvarxevt+tmp+ &
                                     nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+ &
                                     nalea(k)+sumntr+ll)**2
                             end if
                          end do
                       end if
                       
                       if(verb.eq.1) print*,"seuils=",binf, " ",bsup, " sumntr=",sumntr
                       !if(verb.eq.1) print*,"mu=",mu, " sumMesYk=",sumMesYk
                       if(verb.eq.1) print*,"std=",b1(nprob+nrisqtot+nvarxevt+tmp+&
                               nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+yk)
                       !! centrer et standardiser
                       binf = binf - mu(sumparK+sumMesYk+j)
                       bsup = bsup - mu(sumparK+sumMesYk+j)
                       if(nerr(k).gt.0) then
                          binf = binf/b1(nprob+nrisqtot+nvarxevt+tmp+&
                               nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+yk)
                          bsup = bsup/b1(nprob+nrisqtot+nvarxevt+tmp+&
                               nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+yk)
                       end if
                       if(verb.eq.1) print*,"binf=",binf, " bsup=",bsup

                       if(indiceY(nmescur+sumparK+sumMesYk+j).eq.1) then
                          !! si Y=minY
                          vrais_l = vrais_l*alnorm(binf,.false.)
                          if(verb.eq.1) print*,"P(Y=binf)=", alnorm(binf,.false.) 
                       else if(indiceY(nmescur+sumparK+sumMesYk+j).eq.nvalORD(q)) then
                          !! si Y=maxY
                          vrais_l = vrais_l*(1.d0-alnorm(bsup,.false.))
                          if(verb.eq.1) print*,"P(Y=bsup)=", 1-alnorm(bsup,.false.)
                          
                       else
                          !! minY < Y < maxY
                          vrais_l = vrais_l*(alnorm(bsup,.false.)-alnorm(binf,.false.))
                          if(verb.eq.1) print*,"P(binf<Y<bsup)=", alnorm(bsup,.false.),"-",alnorm(binf,.false.)
                       end if
                       if(verb.eq.1) print*,'vrais_l=',vrais_l
                    end do                    
                 else
                    !! yk est continu
                    
                    !! variance de Y|ui,wi
                    VC=0.d0
                    do j1=1,nmes(i,sumny+yk)
                       if(contraint.eq.1) then
                          VC(j1,j1) =  1
                       else
                          VC(j1,j1) = b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                               nvc(k)+nw(k)+ncor(k)+yk)**2 !variance de l'erreur yk
                       end if
                       if (nalea(k).eq.ny(k)) then ! intercept aleatoire de yk
                          do j2=1,nmes(i,sumny+yk)
                             VC(j1,j2) = VC(j1,j2) + &
                                  b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                                  nvc(k)+nw(k)+ncor(k)+nerr(k)+yk)**2
                          end do
                       end if
                    end do

                    ! Vi en vecteur
                    jj=0
                    Vi=0.d0
                    do j1=1,nmes(i,sumny+yk)
                       do j2=j1,nmes(i,sumny+yk)
                          jj=j1+j2*(j2-1)/2
                          Vi(jj)=VC(j1,j2)
                       end do
                    end do

                    ! inversion
                    CALL dsinv(Vi,nmes(i,sumny+yk),eps,ier,det)
                    if (ier.eq.-1) then
                       vrais_mpjord_i=-1.d9
                       print*,"-1.d9 dsinv"
                       print*,"b=",b
                       !print*,"bfix=",bfix
                       !print*,"id=",id,"thi=",thi
                       !print*,"jd=",jd,"thj=",thj
                       !print*,"fix=",fix
                       goto 654
                    end if

                    ! retransformation du vecteur Vi en matrice :
                    VC=0.d0
                    do j1=1,nmes(i,sumny+yk)
                       do j2=1,nmes(i,sumny+yk)
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
                    do j=1,nmes(i,sumny+yk)
                       Y2(j) = Y1(sumparK+sumMesYk+j)-mu(sumMesYk+j)
                    end do
                    Y3=matmul(VC,Y2)
                    Y4=DOT_PRODUCT(Y2,Y3) 

                    div = (dble(2*3.14159265)**(nmes(i,sumny+yk)/2))*sqrt(exp(det))
                    vrais_l = vrais_l*exp(-Y4/2.d0)/div   

                    if(verb.eq.1) print*,"vrais continue = ",exp(-Y4/2.d0)/div                 
                 end if
                 
                 sumMesYk = sumMesYk+nmes(i,sumny+yk)
                 sumntr = sumntr + ntr(sumny+yk)
              end do ! fin boucle yk  

              som = som + vrais_l/dble(npts)   

           end do ! fin boucle l MC

           vrais_K = log(som)

        end if ! fin anyord

        vrais_Y = vrais_Y + vrais_K

        sumny = sumny + ny(k)
        tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+ntrK(k)
        sumparK = sumparK + nmesparK(i,k)
        sumnv = sumnv + nv(k)
        sumnea = sumnea + nea(k)
     end do ! fin boucle nbK

      if(verb.eq.1) print*,"vrais_Y=",vrais_Y

     ! !! log-vraisemblance totale
     !if(verb.eq.1) then
     !print*,"vrais_surv=",vrais_surv
     !print*,"vrais_Y=",vrais_Y
     !print*,"jacobien=",jacobien
     !print*,"entretard=",entretard
     !end if
     vrais_mpjord_i = vrais_surv + vrais_Y + jacobien
     
     ! entree retardee
     if(idtrunc.eq.1) then 
        vrais_mpjord_i = vrais_mpjord_i - entretard
     end if
     
   else ! ng>1
     
     ! !! proba appartenance aux classes !!
     if (prior(i).ne.0) then
        
        pi=0.d0
        pi(prior(i))=1.d0
        
     else
        
        ! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
        Xprob=0.d0
        l=0
        do k=1,nv2
           if (idprob(k).eq.1) then
              l=l+1
              Xprob(l)=Xns(i,k)
           end if
        end do
        pi=0.d0
        temp=0.d0
        do g=1,ng-1
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
     
     ! boucle sur les classes
     expo=0.d0
     vrais_mpjord_i=0.d0
     !vrais_Y=1.d0
     retard=0.d0
     do g=1,ng
        
        ! !! partie survie de la classe g !!
        vrais_surv=1.d0
        if(nbevt.gt.0) then
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
                 do k=1,nv2
                    if (idtdv(k).ne.1) then
                       if (idcom(k).eq.1) then  
                          if (idspecif(k).eq.2) then 
                             l=l+1
                             bevt(l)=b1(nprob+nrisqtot+nevtxcurr+g)
                             Xevt(l)=Xns(i,k)
                          else 
                             l=l+1
                             bevt(l)=b1(nprob+nrisqtot+nevtxcurr+1)
                             Xevt(l)=Xns(i,k)
                          end if

                       else 

                          if (idspecif((ke-1)*nv2+k).eq.1) then   
                             l=l+1
                             bevt(l)=b1(nprob+nrisqtot+nevtxcurr+nxcurr(k)+1)
                             Xevt(l)=X(i,k)
                             nxcurr(k)=nxcurr(k)+1
                          end if
                          if (idspecif((ke-1)*nv2+k).eq.2) then  
                             l=l+1
                             bevt(l)=b1(nprob+nrisqtot+nevtxcurr+nxcurr(k)+g)
                             Xevt(l)=Xns(i,k)
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
                          if (idspecif((ke-1)*nv2+k).eq.1) then 
                             bevtint(ke)=b1(nprob+nrisqtot+nevtxcurr+nxcurr(k)+1)
                             nxcurr(k)=nxcurr(k)+1
                          end if
                          if (idspecif((ke-1)*nv2+k).eq.2) then
                             bevtint(ke)=b1(nprob+nrisqtot+nevtxcurr+nxcurr(k)+g)
                             nxcurr(k)=nxcurr(k)+ng
                          end if
                       end if
                    end if

                    nevtxcurr=nevtxcurr+nevtparx(k)
                 end do
              end do

           end if

          
           surv_glob=0.d0
           surv0_glob=0.d0
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
              surv0_glob=surv0_glob+surv0(g,ke)*exp(varexpsurv)

              if (Devt(i).eq.ke) then
                 fevt=risq(g,ke)*exp(varexpsurv)                     
                 if (ind_survint(i).eq.1) then
                    fevt=fevt*exp(bevtint(ke))
                 end if
              end if
              
              nxevtcurr=nxevtcurr+nxevtspec(ke)
           end do

           ! vrais survie de la classe g
           if (Devt(i).ne.0) then
              vrais_surv = fevt*exp(-surv_glob)
           else
              vrais_surv = exp(-surv_glob)
           end if

           retard = retard + pi(g)*exp(-surv0_glob)
           !print*,"i=",i," g=",g;
           !print*,"Y4=",Y4," surv=",surv_glob;
          ! print*,"Tsurv=",Tsurv(i);
           
        end if ! fin survie

        
        ! !! partie longitudinale de la classe g !!
        sumny=0
        tmp=0
        sumparK=0
        sumnv=0
        vrais_Y=1.d0
        do k=1,nbK
      !     print*,"k=",k
           sumntr=0

           !! variance des effets aleatoires
           Ut=0.d0
           if (nea(k)>0) then

              if(contraint.eq.2) then
                 Ut=0.d0
                 Ut(1,1)=1.d0
                 if (nea(k)>1) then 

                    if (idiag(k).eq.1) then
                       do j=2,nea(k)
                          do kk=2,nea(k)
                             if (j.eq.kk) then
                                Ut(j,kk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j-1)
                             else
                                Ut(j,kk)=0.d0
                             end if
                          end do
                       end do
                    end if

                    if (idiag(k).eq.0) then
                       do j=2,nea(k)
                          do kk=1,j
                             Ut(j,kk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+kk-1+j*(j-1)/2)
                          end do
                       end do
                    end if

                 end if

              else ! pas de variance contrainte a 1

                 if (idiag(k).eq.1) then
                    do j=1,nea(k)
                       do kk=1,nea(k)
                          if (j.eq.kk) then
                             Ut(j,kk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j)
                          else
                             Ut(j,kk)=0.d0
                          end if
                       end do
                    end do
                 end if

                 if (idiag(k).eq.0) then
                    do j=1,nea(k)
                       do kk=1,j
                          Ut(j,kk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+kk+j*(j-1)/2)
                       end do
                    end do
                 end if
              end if

              ! si variance specifique a la classe
              if(nw(k).ne.0.and.g.lt.ng) then
                 Ut = abs(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+g))*Ut
              end if
           end if
         !  print*,"Ut ok"
           !! variables en effet aleatoire
           Z=0.d0
           l=0
           do kk=1,nv(k)
              if (idea(sumnv+kk).eq.1) then
                 l=l+1
                 do j=1,nmesparK(i,k)
                    Z(j,l)=dble(X(nmescur+sumparK+j,sumnv+kk))
                 end do
              end if
           end do

           !print*,"Z ok"
           
           ! variables en effet fixe
           l=0
           m=0
           q=0
           b2=0.d0
           b0=0.d0
           nmoins=0
           X00=0.d0
           X2=0.d0
           X01=0.d0
           b01=0.d0
           do kk=1,nv(k)
              
              ! sans mixture
              if (idg(sumnv+kk).eq.1) then
                 m=m+1
                 do j=1,nmesparK(i,k)
                    X00(j,m)=dble(X(nmescur+sumparK+j,sumnv+kk))
                 end do
                 if (kk.eq.1.and.contraint.ne.0) then
                    b0(m)=0.d0
                 else
                    b0(m)=b1(nprob+nrisqtot+nvarxevt+tmp+nmoins+1)
                    nmoins=nmoins+1
                 end if

                 ! avec mixture
              else if (idg(sumnv+kk).eq.2) then
                 l=l+1
                 do j=1,nmesparK(i,k)
                    X2(j,l)=dble(X(nmescur+sumparK+j,sumnv+kk))
                 end do
                 if (kk.eq.1.and.contraint.ne.0) then
                    if (g.eq.1) then
                       b2(l)=0.d0
                       nmoins=nmoins+ng-1
                    else
                       b2(l)=b1(nprob+nrisqtot+nvarxevt+tmp+nmoins+g-1)
                       nmoins=nmoins+ng-1
                    end if
                 else
                    b2(l)=b1(nprob+nrisqtot+nvarxevt+tmp+nmoins+g)
                    nmoins=nmoins+ng
                 end if
              end if

              !contrast : 
              if (idcontr(sumnv+kk).ne.0) then
                 q=q+1
                 sumMesYk=0
                 do yk=1,ny(k)
                    ! creation matrice design des contrastes: X01
                    do j=1,nmes(i,sumny+yk)
                       X01(sumMesYk+j,(q-1)*ny(k)+yk) = dble(X(nmescur+sumparK+sumMesYk+j,sumnv+kk))
                    end do
                    sumMesYk=sumMesYk+nmes(i,sumny+yk)
                    ! creation vecteur parms des contrastes: b01
                    if (yk<ny(k)) then
                       b01((q-1)*ny(k)+yk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ &
                            (q-1)*(ny(k)-1)+yk)
                    else
                       b01((q-1)*ny+ny) =-sum(b1((nprob+nrisqtot+nvarxevt+tmp+ &
                            nef(k)+(q-1)*(ny(k)-1)+1) &
                            :(nprob+nrisqtot+nvarxevt+tmp+ &
                            nef(k)+(q-1)*(ny(k)-1)+ny(k)-1)))
                    end if
                 end do
              end if
           end do

           
           !! matrice Corr (cor+err+alea)
           Corr=0.d0
           tcor=0.d0
           if(ncor(k).gt.0) then

              do kk=1,nv(k)
                 if (idcor(sumnv+kk).eq.1) then
                    do j=1,nmesparK(i,k)
                       tcor(j) = X(nmescur+sumparK+j,sumnv+kk)
                    end do
                 end if
              end do
              do j1=1,nmesparK(i,k)
                 do j2=1,nmesparK(i,k)
                    if (ncor(k).eq.1) then 
                       Corr(j1,j2) = Corr(j1,j2)+b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+&
                            ncontr(k)+nvc(k)+nw(k)+ncor(k))* &
                            b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+ &
                            nw(k)+ncor(k)) * min(tcor(j1),tcor(j2))
                    else if (ncor(k).eq.2) then
                       Corr(j1,j2) = Corr(j1,j2)+b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+&
                            ncontr(k)+nvc(k)+nw(k)+ncor(k))* &
                            b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+ &
                            nw(k)+ncor(k))* &
                            exp(-b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                            nvc(k)+nw(k)+1) * abs(tcor(j1)-tcor(j2)))
                    end if
                 end do
              end do
           end if

           
        if(anyord(k).eq.0) then 
           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! !!!!!!!!!!!!!!!!!!!   ng>1 et que du continu   !!!!!!!!!!!!!!!!!!!
           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
           sumMesYk=0
           do yk=1,ny(k)
              do j1=1,nmes(i,sumny+yk)
                 if(contraint.eq.1) then
                    Corr(sumMesYk+j1,sumMesYk+j1) =  Corr(sumMesYk+j1,sumMesYk+j1) + 1
                 else
                    Corr(sumMesYk+j1,sumMesYk+j1) =  Corr(sumMesYk+j1,sumMesYk+j1) + &
                         b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                         nvc(k)+nw(k)+ncor(k)+yk)**2 !variance de l'erreur yk
                 end if
                 if (nalea(k).eq.ny(k)) then ! intercept aleatoire de yk
                    do j2=1,nmes(i,sumny+yk)
                       Corr(sumMesYk+j1,sumMesYk+j2) = Corr(sumMesYk+j1,sumMesYk+j2) + &
                            b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                            nvc(k)+nw(k)+ncor(k)+nerr(k)+yk)**2
                    end do
                 end if
              end do
        !   print*,"Corr",yk," ok"
              sumMesYk = sumMesYk + nmes(i,sumny+yk)
           end do

           ! matrice de variance de Yik
           P=0.d0
           P=MATMUL(Z,Ut)
           VC=0.d0
           VC=MATMUL(P,transpose(P))+Corr

           ! Vi en vecteur
           jj=0
           Vi=0.d0
           do j1=1,nmesparK(i,k)
              do j2=j1,nmesparK(i,k)
                 jj=j1+j2*(j2-1)/2
                 Vi(jj)=VC(j1,j2)
              end do
           end do
           
           CALL dsinv(Vi,nmesparK(i,k),eps,ier,det)
           if (ier.eq.-1) then
              vrais_mpjord_i=-1.d9
              print*,"-1.d9 dsinv"
             ! print*,"i=",i, "nmespark=",nmesparK(i,k)
             ! print*,"Vi=",Vi
             ! print*,"b=",b
             ! print*,"bfix=",bfix
             ! print*,"id=",id,"thi=",thi
             ! print*,"jd=",jd,"thj=",thj
             ! print*,"fix=",fix
              goto 654
           end if

           ! retransformation du vecteur Vi en matrice :
           VC=0.d0
           do j1=1,nmesparK(i,k)
              do j2=1,nmesparK(i,k)
                 if (j2.ge.j1) then
                    VC(j1,j2)=Vi(j1+j2*(j2-1)/2)
                 else
                    VC(j1,j2)=Vi(j2+j1*(j1-1)/2)
                 end if
              end do
           end do

           
           ! calcul de la vrais pour Yik
           mu=0.d0
           y2=0.d0
           y3=0.d0
           y4=0.d0
           mu=matmul(X00,b0)+matmul(X2,b2)+matmul(X01,b01)
           
           do j=1,nmesparK(i,k)
              Y2(j) = Y1(sumparK+j)-mu(j)
           end do
           !Y2=Y1-mu
           Y3=Matmul(VC,Y2)
           Y4=0.d0
           Y4=DOT_PRODUCT(Y2,Y3)
           
           vrais_Y = vrais_Y*exp((-det-Y4)/2.d0)

          ! if(i.le.3) then
          !    print*,"i=",i," g=",g, "pi(g)=",pi(g)
          !    print*,"det=",det, " Y4=", Y4
          !   ! print*,"b0=",b0, " b2=",b2
          !   ! print*,'mu=',mu
          !    print*,"vrais_Y = ",vrais_Y
          ! end if
           
!           if(g.eq.1) then
!              if(k.eq.1) vrais_Y1_G1 = vrais_Y1_G1 +((-det-Y4)/2)
!              if(k.eq.2) vrais_Y2_G1 = vrais_Y2_G1 + ((-det-Y4)/2)
!           end if
!           if(g.eq.2) then
!              if(k.eq.1) vrais_Y1_G2 = vrais_Y1_G2 +((-det-Y4)/2)
!              if(k.eq.2) vrais_Y2_G2 = vrais_Y2_G2 + ((-det-Y4)/2)
!           end if

           
              
  !          print*,"b0",b0
!            print*,"b2",b2
!            print*,"mu=",mu
!            print*," det=",det
!            print*,"Y4",Y4
!            print*,"Y1",Y1
!            print*,"Y2",Y2
           ! print*,"vrais_Y ok", exp((-det-Y4)/2)


           
        else
           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! !!!!!!!!!!!!!!!!!!! ng>1 et au moins 1 ordinal !!!!!!!!!!!!!!!!!!!
           ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           
           !cholesky de Corr (variance de AR ou BM)
           if(needMC(k).eq.1 .and. ncor(k).gt.0) then
              jj=0
              Vi=0.d0
              do j=1,maxmes
                 do kk=j,maxmes
                    jj=j+kk*(kk-1)/2
                    Vi(jj)=Corr(j,kk)
                 end do
              end do
              CALL DMFSD(Vi,maxmes,EPS,IER)
              chCorr=0.d0
              do j=1,maxmes
                 do kk=1,j
                    chCorr(j,kk)=Vi(kk+j*(j-1)/2)
                 end do
              end do
           end if

           ! boucle pour les nMC simulations du Monte Carlo
           !nMC=needMC(k)*(nptsMC-1)+1
           !print*,"nMC=",nMC
           som=0.d0
           do l=1,nMC
              expo=1.d0
              ui=0.d0
              wi=0.d0
              ai=0.d0
              !print*," MC l=",l
              if(needMC(k).eq.1) then

                 if(g.eq.1.or.nbK.gt.1) then !** a voir pour avoir les memes points MC d'une classe a l'autre 
                 if(methInteg.eq.1) then 
                    ! !!!!!!!!!!!!! MCO !!!!!!!!!!!!!
                    
                    ! simuler les effets aleatoires
                    if(nea(k).gt.0) then
                       x22=0.d0
                       usim=0.d0
                       SX=1.d0
                       do j=1,nea(k)
                          call bgos(SX,0,usim(j),x22,0.d0)
                          !print*,"usim=",usim(j)
                       end do
                       ui=0.d0
                       ui=matmul(Ut,usim)
                       !print*," ui=",ui, " Ut=",Ut, "  nea=",nea(k)
                    end if

                    ! simuler le BM ou AR
                    if(ncor(k).gt.0) then
                       x22=0.d0
                       usim=0.d0
                       SX=1.d0
                       do j=1,maxmes
                          call bgos(SX,0,wsim(j),x22,0.d0)
                       end do
                       wi=0.d0
                       wi=matmul(chCorr,wsim)
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
                       if(nea(k).gt.0) then
                          x22=0.d0
                          usim=0.d0
                          SX=1.d0
                          do j=1,nea(k)
                             call bgos(SX,0,usim(j),x22,0.d0)
                          end do
                          ui=0.d0
                          ui=matmul(Ut,usim)
                       end if
                       
                       ! simuler le BM ou AR
                       if(ncor(k).gt.0) then
                          x22=0.d0
                          usim=0.d0
                          SX=1.d0
                          do j=1,maxmes
                             call bgos(SX,0,wsim(j),x22,0.d0)
                          end do
                          wi=0.d0
                          wi=matmul(chCorr,wsim)
                       end if
                       
                    end if

                 else 
                    ! !!!!!!!!!!!!! QMC !!!!!!!!!!!!!

                    ! a voir si on peut utiliser LowDiscrepancy.f
                    ! ou si on cree les points dans R avec randtoolbox
                    

                 end if
              end if !**

              else
                 !! pas de MC a faire, donc pas d'EA
                 
                 ui=0.d0
                 wi=0.d0

              end if

              ! esperance conditionnelle
              mu=matmul(X00,b0)+matmul(X01,b01)+matmul(X2,b2)
              mu=mu+matmul(Z,ui)+wi

              
              sumMesYk=0
              sumntr=0
              q=0
              do yk =1,ny(k)

                 if(idlink(yk).eq.3) then
                    !! yk est ordinal
                    q=q+1
                    ai=0.d0
                    if(nalea(k).gt.0) then
                       if(g.eq.1.or.nbK.gt.1) then !**
                       if(methInteg.eq.1) then
                          !! MCO
                          call bgos(SX,0,ai,x22,0.d0)
                          ai = b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                               nvc(k)+nw(k)+ncor(k)+nerr(k)+yk)*ai
                       else if(methInteg.eq.2) then
                          !! MCA
                          if(mod(l,2).eq.0) then
                             ! si l est pair on prend l'oppose du precedent
                             ai = -ai
                          else
                             call bgos(SX,0,ai,x22,0.d0)
                             ai = b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                                  nvc(k)+nw(k)+ncor(k)+nerr(k)+yk)*ai
                          end if
                       else
                          !! QMC
                       end if
                    end if !**
                    end if
                    
                    do j=1,nmes(i,sumny+yk)
                       
                       !! on ajoute ai a mu
                       mu(sumparK+sumMesYk+j) = mu(sumparK+sumMesYk+j)+ai
                       
                       !! trouver binf et bsup tq binf < lambda + epsilon < bsup :
                       
                       !! initialiser au premier seuil
                       binf = b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                            nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr+1)
                       bsup = binf
                      
                       !! si Y>minY ajouter b1(..)^2
                       if(indiceY(nmescur+sumparK+sumMesYk+j).gt.1) then
                          do ll=2,indiceY(nmescur+sumparK+sumMesYk+j)
                             bsup = bsup + b1(nprob+nrisqtot+nvarxevt+tmp+ &
                                  nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+ &
                                  nalea(k)+sumntr+ll)**2
                             if(ll.lt.indiceY(nmescur+sumparK+sumMesYk+j)) then
                                binf = binf + b1(nprob+nrisqtot+nvarxevt+tmp+ &
                                     nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+ &
                                     nalea(k)+sumntr+ll)**2
                             end if
                          end do
                       end if
                       
                       !! centrer et standardiser
                       binf = binf - mu(sumparK+sumMesYk+j)
                       bsup = bsup - mu(sumparK+sumMesYk+j)
                       if(nerr(k).gt.0) then
                          binf = binf/b1(nprob+nrisqtot+nvarxevt+tmp+&
                               nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+yk)
                          bsup = bsup/b1(nprob+nrisqtot+nvarxevt+tmp+&
                               nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+yk)
                       end if

                       if(indiceY(nmescur+sumparK+sumMesYk+j).eq.1) then
                          !! si Y=minY
                          expo = expo*alnorm(binf,.false.)
                          
                       else if(indiceY(nmescur+sumparK+sumMesYk+j).eq.nvalORD(q)) then
                          !! si Y=maxY
                          expo = expo*(1.d0-alnorm(bsup,.false.))
                          
                       else
                          !! minY < Y < maxY
                          expo = expo*(alnorm(bsup,.false.)-alnorm(binf,.false.))
                          
                       end if
                      ! print*,"expo=",expo
                    end do
                    
                 else
                    !! yk est continu
                    
                    !! variance de Y|ui,wi a calculer une seule fois
                    if(l.eq.1) then
                       VC=0.d0
                       do j1=1,nmes(i,sumny+yk)
                          if(contraint.eq.1) then
                             VC(j1,j1) =  1
                          else
                             VC(j1,j1) = b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                                  nvc(k)+nw(k)+ncor(k)+yk)**2 !variance de l'erreur yk
                          end if
                          if (nalea(k).eq.ny(k)) then ! intercept aleatoire de yk
                             do j2=1,nmes(i,sumny+yk)
                                VC(j1,j2) = VC(j1,j2) + &
                                     b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+ &
                                     nvc(k)+nw(k)+ncor(k)+nerr(k)+yk)**2
                             end do
                          end if
                       end do

                       ! Vi en vecteur
                       jj=0
                       Vi=0.d0
                       do j1=1,nmes(i,sumny+yk)
                          do j2=j1,nmes(i,sumny+yk)
                             jj=j1+j2*(j2-1)/2
                             Vi(jj)=VC(j1,j2)
                          end do
                       end do

                       ! inversion
                       CALL dsinv(Vi,nmes(i,sumny+yk),eps,ier,det)
                       if (ier.eq.-1) then
                          vrais_mpjord_i=-1.d9
                          print*,"-1.d9 dsinv"
                          print*,"b=",b
                          !print*,"bfix=",bfix
                          !print*,"id=",id,"thi=",thi
                          !print*,"jd=",jd,"thj=",thj
                          !print*,"fix=",fix
                          goto 654
                       end if

                       ! retransformation du vecteur Vi en matrice :
                       VC=0.d0
                       do j1=1,nmes(i,sumny+yk)
                          do j2=1,nmes(i,sumny+yk)
                             if (j2.ge.j1) then
                                VC(j1,j2)=Vi(j1+j2*(j2-1)/2)
                             else
                                VC(j1,j2)=Vi(j2+j1*(j1-1)/2)
                             end if
                          end do
                       end do
                    end if

                    ! calcul de la vrais
                    Y2=0.d0
                    Y3=0.d0
                    Y4=0.d0
                    do j=1,nmes(i,sumny+yk)
                       Y2(j) = Y1(sumparK+sumMesYk+j)-mu(sumMesYk+j)
                    end do
                    Y3=matmul(VC,Y2)
                    Y4=DOT_PRODUCT(Y2,Y3) 

                    !! partie a sommer sur l
                    expo = expo * exp((-det-Y4)/2.d0)
                      ! print*,"expo=",expo

                                       
                 end if
                 
                 sumMesYk = sumMesYk+nmes(i,sumny+yk)
                 sumntr = sumntr + ntr(sumny+yk)
              end do ! fin boucle yk

              !! somme MC
              som = som + expo
             ! print*,"som=",som
           end do ! fin boucle l MC

           !! passer en log
          ! print*,"vrais_Y avant ajout som ",vrais_Y
           vrais_Y = vrais_Y*som/nMC

           
        end if
           
           sumny = sumny + ny(k)
           tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+ntrK(k)
           sumparK=sumparK+nmesparK(i,k)
           sumnv = sumnv +nv(k)
        end do ! fin boucle K

        ! sommer sur les classes
        expo = expo + pi(g)*vrais_surv*vrais_Y
        
     end do ! fin boucle classes
     
     ! log-vraisemblance totale
     vrais_mpjord_i = -sum(nmesparK(i,:))/2.d0*dlog(dble(2*3.14159265)) + log(expo) + jacobien
     ! if(i.le.3) then
     !    print*,"log(expo)=",log(expo)
     !    print*,"jac=",jacobien
     !    print*,"nmes=",sum(nmesparK(i,:))
     !    print*,"vrais_mpjord_i=",vrais_mpjord_i
     ! endif
     
     ! entree retardee
     if(idtrunc.eq.1) then 
        vrais_mpjord_i = vrais_mpjord_i - log(retard)
     end if
     
    ! print*,"apres idtrunc"
  end if

  
654 continue
 ! print*,"fin vrais i"
  return
  
end function vrais_mpjord_i






double precision function vrais_mpjord(b,m,id,thi,jd,thj)


  use commun_mpjord,only:ns,nmes,nmescur

  implicit none

  integer::m,i,id,jd
  double precision::thi,thj,vrais_mpjord_i,temp
  double precision,dimension(m)::b

  
  nmescur=0
  vrais_mpjord=0.d0
  do i=1,ns
     !call intpr("dans vrais tot i=",-1,i,1)
     !print*,vrais_mpjord_i(b,m,id,thi,jd,thj,i) 
     temp=vrais_mpjord_i(b,m,id,thi,jd,thj,i)  
     !write(*,*)"i=",i,"vrais= ",temp," id=",id," jd=",jd
     !call dblepr("vrais=",-1,temp,1) 
     !write(*,*)"vrais= ",temp        
     vrais_mpjord = vrais_mpjord + temp
     !print*,"vrais_mpjord=",vrais_mpjord
     !if (temp.ne.-1.d9 .or. temp/temp.ne.1) then 
     if (abs(temp+1.d9).lt.1.d-6 .or. temp.ne.temp) then 
        write(*,*)"i=",i,"vrais= ",temp," id=",id," jd=",jd
        vrais_mpjord = -1.d9
        ! if(verbose==1) write(*,*)"i=",i,"vrais= ",temp
        goto 541
     end if
     nmescur = nmescur + sum(nmes(i,:))
    ! print*,"nmescur=",nmescur
    ! print*,"i=",i," vrais_mpjord=",vrais_mpjord
     !if(i.eq.2) goto 541
  end do
541 continue
     !print*,"vrais_mpjord=",vrais_mpjord
  return

end function vrais_mpjord

          







! =============================================
! subroutine de creation de design matrix
! =============================================

      subroutine design_splines_mpjord (ier)

      use commun_mpjord

      implicit none

      integer ::jj,l,k,ier,yk,q,sumnval,nytot
      double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n,hn,hht

      nytot=sum(ny(:))
      ier=0
      jj=0
      l=0
      q=0
      sumnval=0
      do yk=1,nytot
         if (idlink(yk).eq.2) then 
            q=q+1
            do jj=1,nvalSPL(q)      !     ou se trouve la valeur de zi

               do k = 2,ntr(yk)-2
                  if ((uniqueY(sumnval+jj).ge.zitr(k-1,q)).and.(uniqueY(sumnval+jj).lt.zitr(k,q))) then
                     l=k-1
                  end if
                End do  

            if (abs(uniqueY(sumnval+jj)-zitr(ntr(yk)-2,q)).lt.1.d-6) then
               l=ntr(yk)-3
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
            
            if (abs(uniqueY(sumnval+jj)-zitr(ntr(yk)-2,q)).gt.1.d-6) then
               mm2(sumnval+jj) = (3.d0*ht2*ht2)/(hh*h*hn)
               mm1(sumnval+jj) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
               mm(sumnval+jj)  = (3.d0*ht*ht)/(h3*h2*h)

            end if
            if (abs(uniqueY(sumnval+jj)-zitr(ntr(yk)-2,q)).lt.1.d-6) then
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
        
     else if (idlink(yk).eq.3) then        
        ! incrementer sumnval aussi si ordinal car uniqueY contient les modalites
        sumnval = sumnval + nvalORD(q)
        
        
     end if

end do


765     continue

end subroutine design_splines_mpjord







subroutine transfo_estimee_mpjord(b,npm,nsim,marker,transfY)

  use commun_mpjord

  implicit none

  integer::kk,nsim,npm,j,k,yk,sumntr,numSPL,l,sumny,tmp
  double precision,dimension(nsim*sum(ny(:)))::marker,transfY
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

  !     print*,"dans trasnfos"
  !      write(*,*)'infos',minY,maxY,nsim,npm
  !      write(*,*)'b',(b1(j),j=1,npm)


  marker=0.d0
  transfY=0.d0

  sumntr = 0
  sumny = 0
  tmp = 0
  numSPL =0
  do k=1,nbK

     sumntr = 0

     do yk=1,ny(k)

        pas=(maxY(sumny+yk)-minY(sumny+yk))/dble(nsim-1)
        j=1
        marker((sumny+yk-1)*nsim+1)=minY(sumny+yk)
        do while(j.lt.nsim)
           j=j+1
           marker((sumny+yk-1)*nsim+j)=marker((sumny+yk-1)*nsim+j-1)+pas
        end do
        marker((sumny+yk)*nsim)=maxY(sumny+yk)


        if (idlink(sumny+yk).eq.2) then
           numSPL = numSPL+1
           splaa=0.d0

           splaa(1)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                ncor(k)+nerr(k)+nalea(k)+sumntr+1)
           do kk=2,ntr(sumny+yk)
              splaa(kk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                   ncor(k)+nerr(k)+nalea(k)+sumntr+kk)* &
                   b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                   ncor(k)+nerr(k)+nalea(k)+sumntr+kk)
           end do

           !calcul de H(y)
           do j=1,nsim
              ! ou se trouve la valeur
              l=0

              do kk = 2,ntr(sumny+yk)-2
                 if ((marker((sumny+yk-1)*nsim+j).ge.zitr(kk-1,numSPL)).and. &
                      (marker((sumny+yk-1)*nsim+j).lt.zitr(kk,numSPL))) then
                    l=kk-1
                 end if
              end do

              if (marker((sumny+yk-1)*nsim+j).eq.zitr(ntr(sumny+yk)-2,numSPL)) then
                 l=ntr(sumny+yk)-3
              end if

              !         if (l.lt.1.or.l.gt.ntrtot-1) then
              !            write(*,*)'probleme estim splines',l
              !            write(*,*)'j=',j,'test(j)',test(j)
              !            stop
              !         end if


              ht2 = zitr(l+1,numSPL)-marker((sumny+yk-1)*nsim+j)
              htm= marker((sumny+yk-1)*nsim+j)-zitr(l-1,numSPL)
              ht = marker((sumny+yk-1)*nsim+j)-zitr(l,numSPL)
              ht3 = zitr(l+2,numSPL)-marker((sumny+yk-1)*nsim+j)
              hht = marker((sumny+yk-1)*nsim+j)-zitr(l-2,numSPL)
              h = zitr(l+1,numSPL)-zitr(l,numSPL)
              hh= zitr(l+1,numSPL)-zitr(l-1,numSPL)
              hn= zitr(l+1,numSPL)-zitr(l-2,numSPL)
              h2n=zitr(l+2,numSPL)-zitr(l-1,numSPL)
              h2= zitr(l+2,numSPL)-zitr(l,numSPL)
              h3= zitr(l+3,numSPL)-zitr(l,numSPL)

              if (marker((sumny+yk-1)*nsim+j).ne.zitr(ntr(sumny+yk)-2,numSPL)) then
                 mmm2(j) = (3.d0*ht2*ht2)/(hh*h*hn)
                 mmm1(j) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
                 mmm(j)  = (3.d0*ht*ht)/(h3*h2*h)
              end if
              if (marker((sumny+yk-1)*nsim+j).eq.zitr(ntr(sumny+yk)-2,numSPL)) then
                 mmm2(j) = 0.d0
                 mmm1(j) = 0.d0
                 mmm(j)  = 3.d0/h
              end if

              iim2(j) = hht*mmm2(j)/(3.d0) + h2n*mmm1(j)/(3.d0) + h3*mmm(j)/(3.d0)
              iim1(j) = htm*mmm1(j)/(3.d0) + h3*mmm(j)/(3.d0)

              iim(j) = ht*mmm(j)/(3.d0)

              !-------- transformation et IC de la transformation :

              Xspl=0.d0
              Xspl(1)=1
              do kk=2,l
                 Xspl(kk)=1
              end do
              Xspl(l+1)=iim2(j)
              Xspl(l+2)=iim1(j)
              Xspl(l+3)=iim(j)
              transfY((sumny+yk-1)*nsim+j)= dot_product(Xspl,splaa)
           end do
           !fin H(y)


        else if (idlink(sumny+yk).eq.1) then

           aa1=exp(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                ncor(k)+nerr(k)+nalea(k)+sumntr+1))/ &
                (1+exp(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+&
                ncor(k)+nerr(k)+nalea(k)+sumntr+1)))

           bb1=exp(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                ncor(k)+nerr(k)+nalea(k)+sumntr+2))/ &
                (1+exp(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+&
                ncor(k)+nerr(k)+nalea(k)+sumntr+2)))

           cc1=abs(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                ncor(k)+nerr(k)+nalea(k)+sumntr+3))

           dd1=abs(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ &
                ncor(k)+nerr(k)+nalea(k)+sumntr+4))

           bb1=aa1*(1.d0-aa1)*bb1
           aa=aa1*aa1*(1-aa1)/bb1-aa1
           bb=aa*(1-aa1)/aa1

           do j=1,nsim
              ytemp=(marker((sumny+yk-1)*nsim+j)-minY(sumny+yk)+epsY(sumny+yk))/(maxY(sumny+yk)-minY(sumny+yk)+2*epsY(sumny+yk))
              transfY((sumny+yk-1)*nsim+j)=(betai(aa,bb,ytemp)-cc1)/dd1
              if (transfY((sumny+yk-1)*nsim+j).eq.999.d0) then
                 !                    write(*,*)'problem'
              end if

           end do


        else if (idlink(sumny+yk).eq.0) then

           do j=1,nsim
              transfY((sumny+yk-1)*nsim+j) = (marker((sumny+yk-1)*nsim+j) - &
                   b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+ &
                   nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr+1)) &
                   /abs(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+&
                   nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr+2))
           end do
        end if

        sumntr = sumntr + ntr(sumny+yk)

     end do ! fin boucle yk

     sumny = sumny + ny(k)
     tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr

  end do ! fin boucle k

  !     write(*,*)(marker(j),j=1,ny*nsim)
  !     write(*,*)(transfY(j),j=1,ny*nsim)          
end subroutine transfo_estimee_mpjord
