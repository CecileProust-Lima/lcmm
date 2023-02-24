


 
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------

      subroutine loglikhlmevarhetero(Y0,X0,prior0,idprob0,idea0,idg0,idcor0,iderr0,ns0,ng0,nv0,&
           nobs0,nea0,nmes0,idiag0,nwg0,ncor0,nerr0,npm0,b0,nfix0,bfix0,fix0,loglik)

      use commun
      use optim

      IMPLICIT NONE


        !D claration des variables en entree
      integer,intent(in):: nv0,nea0,nfix0
!AB add
      integer, intent(in) :: ns0, ng0, nobs0, idiag0, nwg0, npm0,ncor0,nerr0
      integer, dimension(nv0), intent(in) :: idea0,idg0,idprob0
      integer, dimension(ns0), intent(in) :: nmes0,Prior0,idcor0,iderr0
      double precision, dimension(nobs0), intent(in) :: Y0
      double precision, dimension(nobs0*nv0), intent(in) :: X0
      integer,dimension(npm0+nfix0),intent(in)::fix0
      double precision,dimension(nfix0),intent(in)::bfix0

        !D claration des variable en entr e et sortie
      double precision, dimension(npm0), intent(inout) :: b0

        !D claration des variables en sortie
      double precision,intent(out)::loglik

        !Variables locales
      integer :: jtemp,nef,i,j,k,ktemp,nmestot,it,nbfix,npmtot0
      double precision,external::funcpahlme

      !write(*,*)'B',npm0,b
      !write(*,*)'ncor0',ncor0

      npmtot0 = npm0+nfix0

      maxmes=0
      do i=1,ns0
         if (nmes0(i).gt.maxmes) then
            maxmes=nmes0(i)
         end if
      end do
      
      ! pas de H restreint pr hlme
         ! AB add 
      allocate(Y(ns0*maxmes),idprob(nv0),X(ns0*maxmes,nv0)    &
     ,idea(nv0),idg(nv0),nmes(ns0),prior(ns0),idcor(nv0),iderr(nv0))

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
!AB add
      nerr=nerr0  
      prior=0
      nmes=0
      Y=0.d0
      X=0.d0
      idprob=0
      idea=0
      idg=0
      idcor=0
!AB add
      iderr=0
      nmestot=0
      ktemp=0
      do k=1,nv
         idprob(k)=idprob0(k)
         idea(k)=idea0(k)
         idg(k)=idg0(k)
         idcor(k) = idcor0(k)
!AB add
         iderr(k)=iderr0(k)
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
      do k=1,npmtot0
         fix(k)=fix0(k)
      end do
      nbfix=sum(fix)
      if(nbfix.eq.0) then
         allocate(bfix(1))
      else
         allocate(bfix(nbfix))
      end if
      bfix(1:nbfix)=bfix0(1:nbfix)
      

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
!AB replace
      npmtot=nef+nvc+nwg+ncor+nerr


  ! calcul de la vrais !  
  loglik=funcpahlme(b0,npm0)
 

 1589 continue
!AB add
      deallocate(Y,X,idprob,idea,idg,nmes,prior,idcor,iderr)
      deallocate(fix,bfix)

      return
      end subroutine loglikhlmevarhetero



!-----------------------------------------------------------
!                       FUNCPAHLME
!------------------------------------------------------------


      double precision function funcpahlme(b,npm)

      use commun
      use optim

      IMPLICIT NONE

      integer ::i,j,k,l,m,g,l2,m2,jj,npm,nef,it,j1,j2
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
      double precision ::temp
      double precision ::Y4,expo
      !AB add error that will contain for a subject i the variable var_error of i at all times, need to be an int as ncor 
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
      integer,dimension(maxmes) :: error
      double precision,dimension(ng) :: pi

      b1=0.d0
      eps=1.D-20
      l=0
      m=0
      do k=1,npmtot
         if(fix(k).eq.0) then
            l=l+1
            b1(k)=b(l)
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
               Ut(j,k)=b1(nef+k+j*(j-1)/2)
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
        !AB added for subject i the value of error

        error=1
        if(iderr(nv)==1) then
           do j=1,nmes(i)
              error(j) = int(X(it+j,nv))
           end do
        end if

                 ! AB replaced : 
        ! do j1=1,nmes(i)
           ! do j2=1,nmes(i)
              ! if (j1.eq.j2) Corr(j1,j2) = b1(npmtot)*b1(npmtot)
               !if (ncor.eq.1) then 
                !  Corr(j1,j2) = Corr(j1,j2)+b1(npmtot-1)*b1(npmtot-1)*min(tcor(j1),tcor(j2))
               !else if (ncor.eq.2) then
                !  Corr(j1,j2) = Corr(j1,j2)+b1(npmtot-1)*b1(npmtot-1)*exp(-b1(npmtot-2)*abs(tcor(j1)-tcor(j2)))
               !end if
           ! end do
        !end do
        ! by :
        
         do j1=1,nmes(i)
            do j2=1,nmes(i)
               if (j1.eq.j2) Corr(j1,j2) = b1(nef+nvc+nwg+ncor+error(j1))*b1(nef+nvc+nwg+ncor+error(j1))
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

            CALL dsinv(Vi,nmes(i),eps,ier,det)
            if (ier.eq.-1) then
               funcpahlme=-1.d9
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

                  CALL dsinv(Vi,nmes(i),eps,ier,det)
                  if (ier.eq.-1) then
                     funcpahlme=-1.d9
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
                funcpahlme=-1.d9
                goto 654
            end if
            vrais=vrais+2*log(expo)
         end if
         it=it+nmes(i)
      end do

! FIN BOUCLE SUJET

      funcpahlme=vrais/2.D0

 654  continue

      return

      end function funcpahlme
