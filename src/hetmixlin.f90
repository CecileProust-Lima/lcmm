!==========================================================
!
!      Program for estimating linear mixed models
!         with a mixture of distribution for
!                  the random effects
!
!        Cecile Proust, Helene Jacqmin-Gadda
!
!
!       Corresponding author :
!       C cile Proust, INSERM U897, ISPED,
!       146  rue L\'eo Saignat,
!       33076 Bordeaux cedex, France.
!       Tel: (33) 5 57 57 45 79; Fax: (33) 5 56 24 00 81;
!       e-mail : cecile.proust@isped.u-bordeaux2.fr
!
!                                        10/12/08
!===========================================================
!
!
!--------- MAIN VARIABLES  ---------------------------------
!
! nv : number of covariates in the model
! ng : number of components in the model
! ns : number of units
! Numero(i) : identification number of unit i
! Y(i,j) : dependent variable for unit i at time j
! X(i,j,k) : covariate k for unit i at time j
! idea(k) : indicator of random effect on covariate k :
!                =0 if no random effect
!                =1 if random effect
! idg(k) : indicator of mixture on covariate k :
!                =0 if no effect
!                =1 if overall effect
!                =2 if class-specific effect
! idprob(k) : indicator that covariate k is in the class-membership model :
!                =0 if no effect
!                =1 if overall effect
!                =2 if class-specific effect
! idiag : indicator of the random-effect covariance matrix
!         structure
!                =0 if unstructured matrix
!                =1 if diagonal matrix
!
! maxiter : number of iteration
!
!
!c --------- REMARKS ---------------------------------------
!
! Transformation for the covariance matrix V : Cholesky matrix U
!     => V=U'U
!
! multinomial logistic model for class probabilities
!
!Algorithm of optimisation : modified Marquardt
!(if the loglikelihood at iteration k is not improved
!then the step for the update is modified to ensure
!that the loglikelihood is improved)
!
!
!----------- PROGRAM COMPOSITION --------------------------
!
!    - SUBOURTINE HETMIXLIN
!    - function FUNCPA
!    - subroutine POSTPROB
!
!----------------------------------------------------------





      module commun

      implicit none
      integer,save ::ns,ng,nv,idiag,ncssg,nvc,nea,ncg,nwg,ncor
      integer,save ::nprob,nvarprob,maxmes,nobs,npmtot
      double precision,dimension(:),allocatable,save::Y
      double precision,dimension(:,:),allocatable,save ::X
      integer,dimension(:),allocatable,save ::idea,idg,idprob,idcor
      integer,dimension(:),allocatable,save :: nmes,prior
      integer,dimension(:),allocatable,save::fix
      double precision,dimension(:),allocatable,save::bfix
      end module commun




!*************************************************************

!================ SUBROUTINES ================================

      subroutine hetmixlin(Y0,X0,Prior0, idprob0,idea0,idg0,idcor0  &
          ,ns0,ng0,nv0,nobs0,nea0,nmes0,idiag0,nwg0,ncor0   &
          ,npmtot0,btot,vopt,vrais,ni,istop,gconv,ppi0,resid_m0,resid_ss0 &
          ,pred_m_g0,pred_ss_g0,pred_RE,convB,convL,convG,maxiter0,fix0,pbH0)

      use commun
      use parameters
      use optim

      IMPLICIT NONE


        !D claration des variables en entree
      integer,intent(in):: nv0,maxiter0,nea0
      integer, intent(in) :: ns0, ng0, nobs0, idiag0, nwg0, npmtot0,ncor0
      integer, dimension(nv0), intent(in) :: idea0,idg0,idprob0,idcor0
      integer, dimension(ns0), intent(in) :: nmes0,Prior0
      double precision, dimension(nobs0), intent(in) :: Y0
      double precision, dimension(nobs0*nv0), intent(in) :: X0
      double precision, intent(in) :: convB, convL, convG
      integer,dimension(npmtot0),intent(in)::fix0,pbH0

        !D claration des variable en entr e et sortie
      double precision, dimension(npmtot0), intent(inout) :: btot

        !D claration des variables en sortie
      double precision, intent(out) :: vrais
      double precision, dimension(3), intent(out) :: gconv
      double precision, dimension(ns0*ng0), intent(out) :: ppi0
      double precision, dimension(nobs0), intent(out) :: resid_m0
      double precision, dimension(nobs0), intent(out) :: resid_ss0
      double precision, dimension(nobs0*ng0), intent(out) :: pred_m_g0
      double precision, dimension(nobs0*ng0), intent(out) :: pred_ss_g0
      double precision, dimension(npmtot0*(npmtot0+1)/2),intent(out) :: vopt
      double precision, dimension(npmtot0*(npmtot0+3)/2) :: v
      integer, intent(out) :: ni, istop

        !Variables locales
      integer :: jtemp,nef,i,g,j,ij,npm,ier,k,ktemp,ig,nmestot,it,nbfix
      double precision :: eps, ca, cb, dd
      double precision, dimension(ns0,ng0) :: PPI
      double precision, dimension(npmtot0) :: mvc,b
      double precision, dimension(ns0*nea0), intent(out)::pred_RE
      double precision, dimension(nobs0) :: resid_m, resid_ss
      double precision, dimension(nobs0,ng0):: pred_m_g, pred_ss_g
      double precision,external::funcpa


      !write(*,*)'B',npm0,b
      !write(*,*)'ncor0',ncor0


      maxmes=0
      do i=1,ns0
         if (nmes0(i).gt.maxmes) then
            maxmes=nmes0(i)
         end if
      end do

      epsa=convB
      epsb=convL
      epsd=convG
      maxiter=maxiter0

      ! pas de H restreint pr hlme
         
      allocate(Y(ns0*maxmes),idprob(nv0),X(ns0*maxmes,nv0)    &
     ,idea(nv0),idg(nv0),nmes(ns0),prior(ns0),idcor(nv0))



      ppi0=1.d0
      resid_m0=0.d0
      resid_ss0=0.d0
      pred_m_g0=0.d0
      pred_ss_g0=0.d0
      pred_re=0.d0

      eps=1.d-20
!enrigstrement pour les modules
      ns=ns0
      ng=ng0
      nv=nv0
      nobs=nobs0
      if (nwg0.eq.0) then
         nwg=0
      else
         nwg=ng-1
      end if
          ncor=ncor0

      idiag=idiag0
      prior=0
      nmes=0
      Y=0.d0
      X=0.d0
      idprob=0
      idea=0
      idg=0
          idcor=0
      nmestot=0
      ktemp=0
      do k=1,nv
         idprob(k)=idprob0(k)
         idea(k)=idea0(k)
         idg(k)=idg0(k)
             idcor(k) = idcor0(k)
         jtemp=0
         it=0
         DO i=1,ns
            if (k.eq.1) then
               nmes(i)=nmes0(i)
                   prior(i)=prior0(i)
               do j=1,nmes(i)
                  nmestot=nmestot+1
                  jtemp=jtemp+1
                          Y(jtemp)=Y0(jtemp)
               end do
            end if

            do j=1,nmes(i)
               ktemp=ktemp+1
                       it=it+1
               X(it,k)=X0(ktemp)
            end do
         end do
      end do

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
      

!creation des parametres

      nea=nea0
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
         nprob=nprob+(idprob(k))*(ng-1)
         nvarprob=nvarprob+idprob(k)
      end do

      if((ng.eq.1.and.ncg.gt.0).or.(ng.eq.1.and.nprob.gt.0)) then
         istop=12
         go to 1589
      end if


! nb effets fixes = nb effets fixes sans melange
!                 + ng fois le nb de var dans melange


      if (idiag.eq.1) then
         nvc=nea
      else if(idiag.eq.0) then
         nvc=nea*(nea+1)/2
      end if

      nef=nprob+ncssg+ncg*ng
      npmtot=nef+nvc+nwg+ncor+1

      if (idiag.eq.1) then
         DO j=1,nvc
            btot(nef+j)=dsqrt(abs(btot(nef+j)))
         END DO
      end if

!si idiag=0, on met dans le vecteur des parms, les parms
!de la transformee de Cholesky

      if (idiag.eq.0) then


         DO j=1,nvc
            mvc(j)=btot(nef+j)
         END DO

         CALL DMFSD(mvc,nea,EPS,IER)
         DO j=1,nvc
            btot(nef+j)=mvc(j)
         END DO
      end if
      if (nwg.gt.0) then
         do i=1,nwg
            btot(nef+nvc+i)=abs(btot(nef+nvc+i))
         end do
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
      pbH(1:npm)=0 
      k=0
      do j=1,npmtot
         if(fix0(j).eq.0) then
            k=k+1
            pbH(k)=pbH0(j)
         end if
      end do

!lancement de l'optimisation

        IF (npm.eq.1) then
          istop=12

        else


         ca=0.d0
         cb=0.d0
         dd=0.d0

         call marq98(b,npm,ni,V,vrais,ier,istop,ca,cb,dd,funcpa)

!        write(*,*)
!        write(*,*)'    FIN OPTIMISATION  ..... '
!        write(*,*)'istop',istop,'vrais',vrais

         gconv=0.d0
         gconv(1)=ca
         gconv(2)=cb
         gconv(3)=dd


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
                 

!probas posteriori

         if (istop.eq.1.or.istop.eq.2.or.istop.eq.3) then
            if (ng.gt.1) then
               call postprob(btot,npmtot,PPI)
            end if
                        
            call residuals(btot,npmtot,ppi,resid_m,pred_m_g,resid_ss  &
                ,pred_ss_g,pred_RE)

            ig=0
            ij=0
            do i=1,ns
                if (ng.gt.1) then
                    do g=1,ng0
                        ig=ig+1
                        ppi0(ig)=PPI(i,g)
                    end do
                end if
               do j=1,nmes(i)
                  ij=ij+1
                  resid_ss0(ij)=resid_ss(ij)
                  resid_m0(ij)=resid_m(ij)
                  do g=1,ng0
                     pred_ss_g0(ij+nmestot*(g-1))=pred_ss_g(ij,g)
                     pred_m_g0(ij+nmestot*(g-1))=pred_m_g(ij,g)
                  end do
               end do
            end do

         else
            ig=0
            ij=0
            do i=1,ns
               do g=1,ng0
                  ig=ig+1
                  ppi0(ig)=0.d0
               end do
               do j=1,nmes(i)
                  ij=ij+1
                  resid_ss0(ij)=0.d0
                  resid_m0(ij)=0.d0
                  do g=1,ng0
                     pred_ss_g0(ij+nmestot*(g-1))=0.d0
                     pred_m_g0(ij+nmestot*(g-1))=0.d0
                  end do
               end do
            end do

         end if

      end if

      deallocate(pbH)

 1589 continue

      deallocate(Y,X,idprob,idea,idg,nmes,prior,idcor)
      deallocate(fix,bfix)

      return
      end subroutine hetmixlin



!-----------------------------------------------------------
!                       FUNCPA
!------------------------------------------------------------


      double precision function funcpa(b,npm,id,thi,jd,thj)

      use commun
      use optim

      IMPLICIT NONE

      integer ::i,j,k,l,m,g,l2,m2,id,jd,jj,npm,nef,it,j1,j2
      integer ::ier,nmoins,kk
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(maxmes,nv) ::X00,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nea,nea) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Corr
      double precision,dimension(npm) :: b
      double precision,dimension(npmtot) :: b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: vrais,eps,det
      double precision ::thi,thj,temp
      double precision ::Y4,expo
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
      double precision,dimension(ng) :: pi

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


      
     ! write(*,*)'npm',npm,ncor
      !write(*,*)'b1',b1

!c----------- rappel des parametres utilises ---------

      nef=nprob+ncssg+ncg*ng

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

!----------- boucle sur les individus -------------
       kk=0
       vrais=0.d0
       it=0
       do i=1,ns

!          write(*,*)'i',i,vrais

!-------- creation de Vi = ZiGZi'+Ri+se*seIni ----------

!creation de Zi

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
        
        !matrice Ci=Ri+s2*I
        
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
               if (j1.eq.j2) Corr(j1,j2) = b1(npmtot)*b1(npmtot)
               if (ncor.eq.1) then 
                  Corr(j1,j2) = Corr(j1,j2)+b1(nef+nvc+nwg+ncor)*b1(nef+nvc+nwg+ncor)*min(tcor(j1),tcor(j2))
               else if (ncor.eq.2) then
                  Corr(j1,j2) = Corr(j1,j2)+b1(nef+nvc+nwg+2)*b1(nef+nvc+nwg+2)*exp(-b1(nef+nvc+nwg+1)*abs(tcor(j1)-tcor(j2)))
               end if
            end do
         end do        

        !creation de  et Y1
        
        Y1=0.d0
        do j=1,nmes(i)
           kk=kk+1
           Y1(j)=dble(Y(kk))
        end do
        
!creation de P=Zi*Ut et V=P*P' que si non spec aux classes

         if (nwg.eq.0.OR.NG.EQ.1) then
            P=0.d0
            P=MATMUL(Z,Ut)
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Corr

! Vi en vecteur

            jj=0
            Vi=0.d0
            do j=1,nmes(i)
               do k=j,nmes(i)
                  jj=j+k*(k-1)/2
                  Vi(jj)=VC(j,k)
               end do
            end do

            CALL DSINV(Vi,nmes(i),eps,ier,det)
            if (ier.eq.-1) then
               funcpa=-1.d9
               goto 654
            end if

!     retransformation du vecteur Vi en matrice :

            VC=0.d0
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
!     debut du calcul de la vraisemblance

       vrais=vrais-nmes(i)*dlog(dble(2*3.14159265))


! contribution individuelle a la vraisemblance

! cas 1 : ng=1

       if (ng.eq.1) then

            vrais=vrais-det
            b0=0.d0
            l=0
            X00=0.d0
            do k=1,nv
               if (idg(k).ne.0) then
                  l=l+1
                  do j=1,nmes(i)
                        X00(j,l)=dble(X(it+j,k))
                  end do
                  b0(l)=b1(nprob+l)
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

! cas 2 :  ng>1  composantes
         else
            if (prior(i).ne.0) then
                pi=0.d0
                pi(prior(i))=1.d0
            else
! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
               Xprob=0.d0
               Xprob(1)=1
               l=0
               do k=1,nv
                 if (idprob(k).eq.1) then
                   l=l+1
                   Xprob(1+l)=X(it+1,k)
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
            X00=0.d0
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
                     X00(j,m)=dble(X(it+j,k))
                  end do
               end if
            end do

            b2=0.d0
            b0=0.d0
            expo=0.d0
            do g=1,ng
               nmoins=0
               l2=0
               m2=0
               do k=1,nv
                  if (idg(k).eq.1) then
                     m2=m2+1
                     b0(m2)=b1(nprob+nmoins+1)
                     nmoins=nmoins+1
                  else if (idg(k).eq.2) then
                     l2=l2+1
                     b2(l2)=b1(nprob+nmoins+g)
                     nmoins=nmoins+ng
                  end if
               end do

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
                  VC=MATMUL(P,transpose(P))+Corr

! Vi en vecteur
                  Vi=0.d0
                  jj=0
                  do j=1,nmes(i)
                     do k=j,nmes(i)
                        jj=j+k*(k-1)/2
                        Vi(jj)=VC(j,k)
                     end do
                  end do

                  CALL DSINV(Vi,nmes(i),eps,ier,det)
                  if (ier.eq.-1) then
                     funcpa=-1.d9
                     goto 654
                  end if

!     retransformation du vecteur Vi en matrice :

                  VC=0.d0
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
               y2=0.d0
               y3=0.d0
               y4=0.d0
               mu=matmul(X00,b0)+matmul(X2,b2)
               Y2=Y1-mu
               Y3=Matmul(VC,Y2)
               Y4=0.d0
               Y4=DOT_PRODUCT(Y2,Y3)

               expo = expo+pi(g)*exp((-det-Y4)/2.d0)
            end do

            if (expo.le.0.d0) then
                funcpa=-1.d9
                goto 654
            end if
            vrais=vrais+2*log(expo)
         end if
         it=it+nmes(i)
      end do

! FIN BOUCLE SUJET

      funcpa=vrais/2.D0

 654  continue

      return

      end function funcpa


!------------------------------------------------------------
!                      POSTPROB
!------------------------------------------------------------

!-------------------------------------------------------------
!
!          Subroutine pour calculer les
!      probabilites a posteriori de suivre chacune
!      des composantes g pour chacun des sujets i
!
!                                 25/02/03
!-------------------------------------------------------------

      subroutine postprob(b,npm,PPI)


      use commun
       use optim


      implicit none

      integer ::i,j,k,l,m,g,l2,m2,jj,it,npm,nef,ier,nmoins,kk,j1,j2
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(maxmes,nv) ::X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nea,nea) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Corr
      double precision,dimension(npm) :: b,b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: eps,det
      double precision ::temp
      double precision ::Y4,f
      double precision,dimension(ng) ::fi,pi
      double precision,dimension(ns,ng) ::PPI
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor



      eps=1.D-20

      PPI=0.D0

      do k=1,npm
         b1(k)=b(k)
      end do


!----------- rappel des parametres utilises ---------

      nef=nprob+ncssg+ncg*ng

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
               Ut(j,k)=b1(nef+k+j*(j-1)/2)
            end do
         end do
      end if


! ----------- boucle sur les individus -------------
      kk=0
        it=0
      do i=1,ns

! -------- creation de Vi = ZiGZi'+R+se*seIni ----------

! creation de Zi

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



! creation de Y1

                 Y1=0.d0
         do j=1,nmes(i)
             kk=kk+1
             Y1(j)=dble(Y(kk))
         end do
                 

        !matrice Ci=Ri+s2*I
        
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
               if (j1.eq.j2) Corr(j1,j2) = b1(npm)*b1(npm)
               if (ncor.eq.1) then 
                  Corr(j1,j2) = Corr(j1,j2)+b1(npm-1)*b1(npm-1)*min(tcor(j1),tcor(j2))
               else if (ncor.eq.2) then
                  Corr(j1,j2) = Corr(j1,j2)+b1(npm-1)*b1(npm-1)*exp(-b1(npm-2)*abs(tcor(j1)-tcor(j2)))
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
            PPI=-1.d0
            go to 147
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

       if (prior(i).ne.0) then
          pi=0.d0
          pi(prior(i))=1.d0
       else
! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
       Xprob=0.d0
       Xprob(1)=1
       l=0
       do k=1,nv
          if (idprob(k).eq.1) then
             l=l+1
             Xprob(1+l)=X(it+1,k)
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

       pi(ng)=1/(1+temp)

       do g=1,ng-1
          pi(g)=pi(g)*pi(ng)
       end do
      endif
!     write(*,*)'pi',(pi(g),g=1,ng)

!     creation des vecteurs de variables explicatives
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


!     calcul de la vraisemblance par composante
       f=0.d0
       fi=0.d0
            b2=0.d0
            b0=0.d0
       do g=1,ng
          nmoins=0
          l2=0
          m2=0
          do k=1,nv
             if (idg(k).eq.1) then
                m2=m2+1
                b0(m2)=b1(nprob+nmoins+1)
                nmoins=nmoins+1
             else if (idg(k).eq.2) then
                l2=l2+1
                b2(l2)=b1(nprob+nmoins+g)
                nmoins=nmoins+ng
             end if
          end do

!     variance covariance si spec aux classes :

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
                PPI=-1.d0
                goto 147
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
          fi(g)=fi(g)- nmes(i)*log(dble(2*3.14159265))
          fi(g)=fi(g) -det
          fi(g)=fi(g) - Y4
          fi(g)=fi(g)/(2.d0)
          fi(g)=exp(fi(g))

       end do
       f=DOT_PRODUCT(pi,fi)
       do g=1,ng
          PPI(i,g)=pi(g)*fi(g)/f
       end do

       it=it+nmes(i)

      end do

 147  continue
      return

      end subroutine postprob

!------------------------------------------------------------
!                      RESIDUALS
!------------------------------------------------------------


      subroutine residuals(b1,npm,ppi,resid_m,pred_m_g,resid_ss, &
        pred_ss_g,pred_RE)

      use commun
      use optim

      implicit none
      integer ::i,j,k,l,m,g,l2,m2,jj,npm,nef,j1,j2
      integer ::ier,nmoins,nmes_cur,n2,nmoins2,kk
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(maxmes,nv) ::X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nea) :: err2
      double precision,dimension(nea,nea) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Corr,VC1,SigmaE,CovDev
      double precision,dimension(npm) ::b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob,b3
      double precision :: eps,det
      double precision ::temp
      double precision,dimension(nea,maxmes)::Valea
      double precision,dimension(maxmes) :: mu,Y1,Y2,pred1,err1,tcor
      double precision,dimension(ng) :: pi
      double precision,dimension(nobs)::resid_m &
       ,pred_m,resid_ss,pred_ss
      double precision,dimension(nobs,ng)::     &
    pred_m_g,pred_ss_g
      double precision,dimension(ns,ng) ::PPI
      double precision,dimension(ns*nea)::pred_RE






      eps=1.D-20

!c----------- rappel des parametres utilises ---------

      nef=nprob+ncssg+ncg*ng



!      write(*,*)'nvc',nvc,'nea',nea,'nwg',nwg,'nef',nef
! creation de Ut, decomposition de cholesky pour G
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
               Ut(j,k)=b1(nef+k+j*(j-1)/2)
            end do
         end do
      end if

!         write(*,*)'Ut',Ut

! ----------- boucle sur les individus -------------

      resid_m=0.d0
      pred_m=0.d0
      pred_m_g=0.d0
      pred_ss=0.d0
      resid_ss=0.d0
      pred_ss_g=0.d0
      pred_re=0.d0
      nmes_cur=0
      kk=0
      do i=1,ns

! -------- creation de Vi = ZiGZi'+se*seIni ----------

! creation de Zi

         Z=0.d0
         l=0
         do k=1,nv
            if (idea(k).eq.1) then
               l=l+1
               do j=1,nmes(i)
                  Z(j,l)=dble(X(nmes_cur+j,k))
               end do
            end if

         end do

        !matrice Corr et sigmaE
        
        Corr=0.d0
        tcor=0.d0
        SigmaE=0.d0
        if (ncor.gt.0) then
           do k=1,nv
              if (idcor(k).eq.1) then
                 do j=1,nmes(i)
                    tcor(j) = X(nmes_cur+j,k)
                 end do
              end if
           end do
         end if
         do j1=1,nmes(i)
            do j2=1,nmes(i)
               if (j1.eq.j2) sigmaE(j1,j2) = b1(npm)*b1(npm)
               if (ncor.eq.1) then 
                  Corr(j1,j2) = Corr(j1,j2)+b1(npm-1)*b1(npm-1)*min(tcor(j1),tcor(j2))
               else if (ncor.eq.2) then
                  Corr(j1,j2) = Corr(j1,j2)+b1(npm-1)*b1(npm-1)*exp(-b1(npm-2)*abs(tcor(j1)-tcor(j2)))
               end if
            end do
         end do  

!     creation de Y1

         Y1=0.d0
         do j=1,nmes(i)
            kk=kk+1
            Y1(j)=dble(Y(kk))
         end do

!     creation de P=Zi*Ut et V=P*P' que si non spec aux classes


!     debut du calcul de la vraisemblance

!     cas 1 : ng=1

         if (ng.eq.1) then


            Valea=0.d0
            VC=0.d0
            P=0.d0
            Covdev=0.d0


            P=MATMUL(Z,Ut)
            Valea=MATMUL(Ut,transpose(P))
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Corr+SigmaE
            CovDev = MATMUL(P,transpose(P))+Corr

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
               do j=1,nmes(i)
                  resid_m(nmes_cur+j)=9999.d0
                  pred_m(nmes_cur+j)=9999.d0
                  pred_m_g(nmes_cur+j,1)=9999.d0

                  resid_ss(nmes_cur+j)=9999.d0
                  pred_ss(nmes_cur+j)=9999.d0
                  pred_ss_g(nmes_cur+j,1)=9999.d0
               end do
               do k=1,nea
                  pred_RE((i-1)*nea+k)=9999.d0
               end do
               goto 654
            end if

!     retransformation du vecteur Vi en matrice :
            VC1=0.d0
            do j=1,nmes(i)
               do k=1,nmes(i)
                  if (k.ge.j) then
                     VC1(j,k)=Vi(j+k*(k-1)/2)
                  else
                     VC1(j,k)=Vi(k+j*(j-1)/2)
                  end if
               end do
            end do


            b0=0.d0
            l=0
            X0=0.d0
            do k=1,nv
               if (idg(k).ne.0) then
                  l=l+1
                  do j=1,nmes(i)
                     X0(j,l)=dble(X(nmes_cur+j,k))
                  end do
                  b0(l)=b1(nprob+l)
               end if
            end do

            mu=matmul(X0,b0)
            Y2=Y1-mu

            err1=0.d0
            err1=MATMUL(VC1,Y2)
            err2=0.d0
            err2=MATMUL(Valea,err1)
            pred1=0.d0
            pred1=mu+MATMUL(CovDev,err1)
            do j=1,nmes(i)
               resid_m(nmes_cur+j)=Y2(j)
               pred_m(nmes_cur+j)=mu(j)
               pred_m_g(nmes_cur+j,1)=mu(j)

               resid_ss(nmes_cur+j)=Y1(j)-pred1(j)
               pred_ss(nmes_cur+j)=pred1(j)
               pred_ss_g(nmes_cur+j,1)=pred_ss(nmes_cur+j)

            end do


            do k=1,nea
               pred_RE((i-1)*nea+k)=err2(k)
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

!     write(*,*)'g=',g,'nvarprob',nvarprob,'bprob='
!     &,(bprob(k),k=1,nvarprob)

               temp=temp+exp(DOT_PRODUCT(bprob,Xprob))
               pi(g)=exp(DOT_PRODUCT(bprob,Xprob))

            end do
            pi(ng)=1.d0/(1.d0+temp)
            do g=1,ng-1
               pi(g)=pi(g)*pi(ng)
            end do
          end if
!     write(*,*)'pi',(pi(g),g=1,ng)

!     creation des vecteurs de variables explicatives
            l=0
            m=0
            X0=0.d0
            X2=0.d0
            do k=1,nv
               if (idg(k).eq.2) then
                  l=l+1
                  do j=1,nmes(i)
                     X2(j,l)=dble(X(nmes_cur+j,k))
                  end do
               else if (idg(k).eq.1) then
                  m=m+1
                  do j=1,nmes(i)
                     X0(j,m)=dble(X(nmes_cur+j,k))
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
                     m2=m2+1
                     b0(m2)=b1(nprob+nmoins+1)
                     nmoins=nmoins+1
                  else if (idg(k).eq.2) then
                     l2=l2+1
                     b2(l2)=b1(nprob+nmoins+g)
                     nmoins=nmoins+ng
                  end if
                  IF (IDEA(k).EQ.1.and.idg(k).eq.1) THEN
                     n2=n2+1
                     b3(n2)=b1(nprob+nmoins2+1)
                     nmoins2=nmoins2+1
                  else if(IDEA(k).EQ.1.and.idg(k).eq.2) THEN
                     n2=n2+1
                     b3(n2)=b1(nprob+nmoins2+g)
                     nmoins2=nmoins2+1
                  end if
               end do

               mu=0.d0
               mu=matmul(X0,b0)+matmul(X2,b2)
                VC=0.d0
               P=0.d0
               Ut1=Ut
               if (nwg.ne.0) then
                  if (g.eq.ng) then
                     Ut1=Ut
                  else
                     Ut1=Ut*b1(nef+nvc+g)
                  end if
               end if
               P=0.d0
               Valea=0.d0
               VC=0.d0
               P=MATMUL(Z,Ut1)
               Valea=MATMUL(Ut1,transpose(P))
               VC=MATMUL(P,transpose(P))+Corr+SigmaE
               CovDev = MATMUL(P,transpose(P))+Corr


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
                  do j=1,nmes(i)
                     resid_m(nmes_cur+j)=9999.d0
                     pred_m(nmes_cur+j)=9999.d0

                     resid_ss(nmes_cur+j)=9999.d0
                     pred_ss(nmes_cur+j)=9999.d0
                     do l=1,ng
                        pred_ss_g(nmes_cur+j,l)=9999.d0
                        pred_m_g(nmes_cur+j,l)=9999.d0
                     end do
                  end do
                  do k=1,nea
                    pred_RE((i-1)*nea+k)=9999.d0
                  end do
                  goto 654
               end if

!     retransformation du vecteur Vi en matrice :
               VC1=0.d0
               do j=1,nmes(i)
                  do k=1,nmes(i)
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
               err2=MATMUL(Valea,err1)
               pred1=0.d0
               pred1=mu+MATMUL(CovDev,err1)

               do j=1,nmes(i)
                  pred_m_g(nmes_cur+j,g)=mu(j)
                  pred_ss_g(nmes_cur+j,g)=pred1(j)

                  resid_ss(nmes_cur+j)=resid_ss(nmes_cur+j)  &
                 +ppi(i,g)*(Y1(j)-pred1(j))
                  pred_ss(nmes_cur+j)=pred_ss(nmes_cur+j)+ppi(i,g)   &
               *pred1(j)
                  pred_m(nmes_cur+j)=pred_m(nmes_cur+j)+pi(g)*mu(j)
                  resid_m(nmes_cur+j)=resid_m(nmes_cur+j)+pi(g)*(Y2(j))

               end do

               do k=1,nea
                  pred_RE((i-1)*nea+k)=pred_RE((i-1)*nea+k)+ppi(i,g)*  &
                  err2(k)
               end do

            end do

         end if

 654  continue

          nmes_cur=nmes_cur+nmes(i)
       enddo

! FIN BOUCLE SUJET


      end subroutine residuals








