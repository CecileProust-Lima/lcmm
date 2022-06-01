!==========================================================
!
!      Latent class mixed model for curvilinear outcome
!           using  Beta CDF or quadratic I-Splines
!
!        Cecile Proust, Helene Jacqmin-Gadda
!
!
!       Corresponding author :
!       Cecile Proust, INSERM U897, ISPED,
!       146  rue L\'eo Saignat,
!       33076 Bordeaux cedex, France.
!       Tel: (33) 5 57 57 45 79; Fax: (33) 5 56 24 00 81;
!       e-mail : cecile.proust@isped.u-bordeaux2.fr
!
!                                       21/12/2010
!===========================================================
! - Version fortran 90
!






!----------------------------------------------------------
!
!- Module COMMUN avec les donnees dynamiques
!
!----------------------------------------------------------




      module communc

      implicit none
      integer,save ::ns,ng,nv,idiag,ncssg,nvc,nea,ncg,nwg,ncor  &
      ,nprob,nvarprob,maxmes,nobs,ntrtot,nrisq,nvarxevt,nef &
      ,idlink,npmtot
      double precision,dimension(:),allocatable,save::Y
      double precision,dimension(:,:),allocatable,save ::X
      integer,dimension(:),allocatable,save ::idea,idg,idprob,idcor
      integer,dimension(:),allocatable,save :: nmes,prior
      double precision,save :: minY,maxY,epsY
      double precision,dimension(:),allocatable,save::zitr
      double precision,dimension(:),allocatable,save::mm,mm1,mm2,im,im1,im2
      integer,save::rangeY
      integer,dimension(:),allocatable,save::fix
      double precision,dimension(:),allocatable,save::bfix
      end module communc





      module donnees_indivc

      implicit none
      double precision,dimension(:,:),allocatable::Ut1
      double precision,dimension(:),allocatable::mu
      double precision,dimension(:,:),allocatable::Z
      integer::numpat,nmescur
      integer,parameter ::nf=1
      double precision,dimension(:),allocatable,save::seuils
      end module donnees_indivc









!===========================================================
!
      !SUBROUTINES
!
!===========================================================



!-----------------------------------------------------------
!                        FUNCPA
!------------------------------------------------------------


      double precision function funcpac(b,npm,id,thi,jd,thj)

      use communc
!      use optim
      IMPLICIT NONE
      integer ::i,j,k,l,m,g,l2,m2,id,jd,jj,npm,it,nmestot,ll,ii,j1,j2
      integer ::ier,nmoins,kk
      double precision,dimension(maxmes,nv) ::Z,P,X00,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nv,nv) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Corr
      double precision,dimension(npm) :: b
      double precision,dimension(npmtot) :: b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: vrais,eps,det,som,thi,thj,temp
      double precision ::Y4,expo,jacobien,beta_densite,ytemp
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
      double precision,dimension(ng) :: pi
      double precision,dimension(-1:(ntrtot-3))::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,cc1

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



!            print*,"dans funcpa"
!----------- rappel des parametres utilises ---------

  !if(id.eq.0.and.jd.eq.0) then
  ! write(*,*) "avant de remplir Ut :"
  ! write(*,*) "b1=",b1
  ! write(*,*) "nef=",nef
  ! write(*,*) "npm=",npm
  ! write(*,*) "nea=",nea
  ! write(*,*) "nprob+ncssg+ncg*ng+nvarxevt+nrisq-1=",nprob,"+",ncssg,"+",ncg,"*",ng,"+",nvarxevt,"+",nrisq,"-",1
  !end if
   
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
      nmestot=0
      vrais=0.d0
      jacobien=0.d0


      it=0
      do i=1,ns

! -------- creation de Vi = ZiGZi'+Ri+se*seIni = ZiGZi'+ Corri ----------
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
                 
                 
!matrice Corri=Ri+s2*I
        
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
               if (j1.eq.j2) Corr(j1,j2) = 1
               if (ncor.eq.1) then 
                  Corr(j1,j2) = Corr(j1,j2)+b1(npmtot)*b1(npmtot)*min(tcor(j1),tcor(j2))
               else if (ncor.eq.2) then
                  Corr(j1,j2) = Corr(j1,j2)+b1(npmtot)*b1(npmtot)*exp(-b1(npmtot-1)*abs(tcor(j1)-tcor(j2)))
               end if
            end do
         end do    
          ! print*,"Corr ok",Corr
                 
! creation de Y1

         Y1=0.d0

         if (idlink.eq.0) then  ! Linear link

            do j=1,nmes(i)
               nmestot=nmestot+1
               Y1(j)=(dble(Y(nmestot))-b1(nef+nvc+nwg+1))/abs(b1(nef+nvc+nwg+2))
               jacobien = jacobien - log(b1(nef+nvc+nwg+2))
            end do

         else if (idlink.eq.1) then  ! Beta link


            aa1=exp(b1(nef+nvc+nwg+1))/ &
             (1+exp(b1(nef+nvc+nwg+1)))
            bb1=exp(b1(nef+nvc+nwg+2))/ &
             (1+exp(b1(nef+nvc+nwg+2)))
            bb1=aa1*(1.d0-aa1)*bb1

            cc1=abs(b1(nef+nvc+nwg+3))

            dd1=abs(b1(nef+nvc+nwg+4))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1

            do j=1,nmes(i)

               nmestot=nmestot+1

               ytemp=(dble(Y(nmestot))-minY+epsY)/(maxY-minY+2*epsY)
               Y1(j)=(betai(aa,bb,ytemp)-cc1)/dd1


               if (Y1(j).eq.999.d0) then
                  funcpac=-1.d9
                  goto 654
               end if

               jacobien = jacobien + log(abs(beta_densite(ytemp,aa,bb))/dd1)
               jacobien=jacobien-log(abs(maxY-minY+2*epsY))
            end do

         else if (idlink.eq.2) then ! Splines link



            bb=b1(nef+nvc+nwg+1)

            do kk=2,ntrtot
            splaa(kk-3)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
            end do

            do j=1,nmes(i)

               nmestot=nmestot+1

               ll=0
               if (Y(nmestot).eq.zitr(ntrtot-2)) then
                  ll=ntrtot-3
               end if
               som=0.d0
               do kk = 2,ntrtot-2
                  if ((Y(nmestot).ge.zitr(kk-1)).and. &
                      (Y(nmestot).lt.zitr(kk))) then
                     ll=kk-1
                  end if
               end do
               if (ll.lt.1.or.ll.gt.ntrtot-3) then
                funcpac=-1.d9
                goto 654
               end if
               if (ll.gt.1) then
                  do ii=2,ll
                     som=som+splaa(ii-3)
                  end do
               end if

               Y1(j)=bb+ som +splaa(ll-2)*im2(nmestot)  &
                +splaa(ll-1)*im1(nmestot)+splaa(ll)*im(nmestot)

               jacobien = jacobien + log(splaa(ll-2)*mm2(nmestot) &
                 +splaa(ll-1)*mm1(nmestot)+splaa(ll)*mm(nmestot))
            end do
         end if


!         if (i.lt.3)then
!            write(*,*)'nmes',nmes(i),b1((nef+nvc+nwg+1):npm),nef
!            write(*,*)'Y1',Y1
!            write(*,*)'Y',Y(nmestot-nmes(i)+1:nmestot)
!         end if


! creation de P=Zi*Ut et V=P*P' que si non spec aux classes

         if (nwg.eq.0.OR.NG.EQ.1) then
            P=0.d0
            P=MATMUL(Z,Ut)
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Corr
               !print*,"variance ok"
! Vi en vecteur

            jj=0
            Vi=0.d0
            do j=1,nmes(i)
               do k=j,nmes(i)
                  jj=j+k*(k-1)/2
                  Vi(jj)=VC(j,k)
               end do
            end do
                  !print*,"avant dsinv"
            CALL dsinv(Vi,nmes(i),eps,ier,det)
            !print*,"apres dsinv"
            if (ier.eq.-1) then
               funcpac=-1.d9
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

        ! print*,"v inverse",VC
!     debut du calcul de la vraisemblance
       vrais=vrais-nmes(i)*dlog(dble(2*3.14159265))
! contribution individuelle a la vraisemblance
          !print*,"debut du calcul de la vraisemblance"
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
! idg ne 0 pour l'intercept forcement donc on met le parm a 0
                  if (k.eq.1) then
                     b0(l)=0.d0
                  else
                     b0(l)=b1(nprob+l-1)
                  end if
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
! if(i<9) then
! print*, "i=",i
! print*, "ni=",nmes(i)," log(det)=",det," jac=",jacobien
! print*, "mu=",mu
! print*, "Zi=",Z
! print*, "B=",MATMUL(Ut,transpose(Ut))
! endif

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


!               if (thi.eq.0.and.thj.eq.0) then
!                if (i.eq.1) then
!                    write(*,*)'g',g,b2
!                    write(*,*)'g',g,b0
!                    stop
!                    end if
!                    end if


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
                     funcpac=-1.d9
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
            vrais=vrais+2*log(expo)

         end if

         it=it+nmes(i)
       end do

! FIN BOUCLE SUJET

      funcpac=vrais/2.D0+jacobien
          !print*,"fin funcpa"
 654  continue

      return

      end function funcpac


!------------------------------------------------------------
!                      POSTPROB
!------------------------------------------------------------

!-------------------------------------------------------------
!
!          Subroutine pour calculer les
!      probabilites a posteriori de suivre chacune
!      des composantes g pour chacun des sujets i
!
!-------------------------------------------------------------

      subroutine postprobc(b,npm,PPI)
      use communc
!      use optim
      implicit none

      integer ::i,j,k,l,m,g,l2,m2,jj,it,npm,ier,nmoins,kk,nmestot,ii,ll,j1,j2
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(maxmes,nv) ::X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nea,nea) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Corr
      double precision,dimension(npm) :: b,b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: eps,det,som,temp,Y4,f,ytemp
      double precision,dimension(ng) ::fi,pi
      double precision,dimension(ns,ng) ::PPI
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
      double precision,dimension(-1:(ntrtot-3))::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,cc1


      eps=1.D-20

      PPI=0.D0

      do k=1,npm
         b1(k)=b(k)
      end do


!----------- rappel des parametres utilises ---------


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
      nmestot=0
      it=0
      do i=1,ns

! -------- creation de Vi = ZiGZi'+se*seIni+Ri ----------

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

!matrice Corri=Ri+s2*I
        
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
               if (j1.eq.j2) Corr(j1,j2) = 1
               if (ncor.eq.1) then 
                  Corr(j1,j2) = Corr(j1,j2)+b1(npm)*b1(npm)*min(tcor(j1),tcor(j2))
               else if (ncor.eq.2) then
                  Corr(j1,j2) = Corr(j1,j2)+b1(npm)*b1(npm)*exp(-b1(npm-1)*abs(tcor(j1)-tcor(j2)))
               end if
            end do
         end do     

! creation de Y1
         Y1=0.d0
         if (idlink.eq.0) then  ! Linear link

            do j=1,nmes(i)
               nmestot=nmestot+1
               Y1(j)=dble(Y(nmestot)-b1(nef+nvc+nwg+1))/abs(b1(nef+nvc+nwg+2))
            end do

         else if (idlink.eq.1) then  ! Beta link


            aa1=exp(b1(nef+nvc+nwg+1))/ &
              (1+exp(b1(nef+nvc+nwg+1)))
            bb1=exp(b1(nef+nvc+nwg+2))/ &
              (1+exp(b1(nef+nvc+nwg+2)))
            bb1=aa1*(1.d0-aa1)*bb1
            cc1=b1(nef+nvc+nwg+3)
            dd1=abs(b1(npm-ncor))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1

            do j=1,nmes(i)

               nmestot=nmestot+1

               ytemp=(dble(Y(nmestot))-minY+epsY)/(maxY-minY+2*epsY)
               Y1(j)=(betai(aa,bb,ytemp)-cc1)/dd1


               if (Y1(j).eq.999.d0) then
                PPI=-1.d0
                go to 147
               end if

            end do

         else if (idlink.eq.2) then ! Splines link

            bb=b1(nef+nvc+nwg+1)


            do kk=2,ntrtot
               splaa(kk-3)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
            end do

            do j=1,nmes(i)

               nmestot=nmestot+1

               ll=0
               if (Y(nmestot).eq.zitr(ntrtot-2)) then
                  ll=ntrtot-3
               end if
               som=0.d0
               do kk = 2,ntrtot-2
                  if ((Y(nmestot).ge.zitr(kk-1)).and.(Y(nmestot).lt.zitr(kk))) then
                     ll=kk-1
                  end if
               end do
               if (ll.lt.1.or.ll.gt.ntrtot-3) then
                PPI=-1.d0
                go to 147
               end if
               if (ll.gt.1) then
                  do ii=2,ll
                     som=som+splaa(ii-3)
                  end do
               end if

               Y1(j)=bb+som +splaa(ll-2)*im2(nmestot) &
               +splaa(ll-1)*im1(nmestot)+splaa(ll)*im(nmestot)
            end do
        end if

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

         CALL dsinv(Vi,nmes(i),eps,ier,det)
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

             CALL dsinv(Vi,nmes(i),eps,ier,det)
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

      end subroutine postprobc

!------------------------------------------------------------
!                      RESIDUALS
!------------------------------------------------------------


      subroutine residualsc(b1,npm,ppi,resid_m,pred_m_g,resid_ss,pred_ss_g,pred_RE,Yobs)

      use communc

!      use optim

      implicit none
      integer ::i,j,k,l,m,g,l2,m2,jj,npm,nmestot,ll,ii,j1,j2
      integer ::ier,nmoins,nmes_cur,n2,nmoins2,kk
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(maxmes,nv) ::X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nea) ::err2
      double precision,dimension(nea,nea) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Corr,VC1,SigmaE,CovDev
      double precision,dimension(npm) ::b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob,b3
      double precision :: eps,det,temp
      double precision,dimension(nea,maxmes)::Valea
      double precision,dimension(maxmes) :: mu,Y1,Y2,pred1,err1,tcor
      double precision,dimension(ng) :: pi
      double precision,dimension(nobs)::resid_m,resid_ss,Yobs
      double precision,dimension(nobs*ng)::pred_m_g,pred_ss_g
      double precision,dimension(ns,ng) ::PPI
      double precision,dimension(ns*nea)::pred_RE
      double precision,dimension(-1:(ntrtot-3))::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,ytemp,som,cc1

      eps=1.D-20

!----------- rappel des parametres utilises ---------

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

! ----------- boucle sur les individus -------------

      resid_m=0.d0
      pred_m_g=0.d0
      resid_ss=0.d0
      pred_ss_g=0.d0
      Yobs=0.d0
      nmes_cur=0
      kk=0
      nmestot=0
      pred_RE=0.d0
      do i=1,ns

!     -------- creation de Vi = ZiGZi'+Ri+se*seIni ----------
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
               if (j1.eq.j2) SigmaE(j1,j2) = 1
               if (ncor.eq.1) then 
                  Corr(j1,j2) = Corr(j1,j2)+b1(npm)*b1(npm)*min(tcor(j1),tcor(j2))
               else if (ncor.eq.2) then
                  Corr(j1,j2) = Corr(j1,j2)+b1(npm)*b1(npm)*exp(-b1(npm-1)*abs(tcor(j1)-tcor(j2)))
               end if
            end do
         end do     
                 
! creation de Y1
         Y1=0.d0
         if (idlink.eq.0) then  ! Linear link

            do j=1,nmes(i)
               nmestot=nmestot+1
               Y1(j)=dble(Y(nmestot)-b1(nef+nvc+nwg+1))/abs(b1(nef+nvc+nwg+2))
               Yobs(nmes_cur+j)=Y1(j)
            end do

         elseif (idlink.eq.1) then  ! Beta link


            aa1=exp(b1(nef+nvc+nwg+1))/ &
             (1+exp(b1(nef+nvc+nwg+1)))
            bb1=exp(b1(nef+nvc+nwg+2))/ &
             (1+exp(b1(nef+nvc+nwg+2)))
            bb1=aa1*(1.d0-aa1)*bb1
            cc1=b1(nef+nvc+nwg+3)
            dd1=abs(b1(npm-ncor))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1


            do j=1,nmes(i)

               nmestot=nmestot+1

               ytemp=(dble(Y(nmestot))-minY+epsY)/(maxY-minY+2*epsY)
               Y1(j)=(betai(aa,bb,ytemp)-cc1)/dd1

               if (Y1(j).eq.999.d0) then
                  do k=1,nmes(i)
                     resid_m(nmes_cur+k)=9999.d0
                     pred_m_g(nmes_cur+k)=9999.d0
                     Yobs(nmes_cur+k)=9999.d0
                     resid_ss(nmes_cur+k)=9999.d0
                     pred_ss_g(nmes_cur+k)=9999.d0
                  end do
                  do k=1,nea
                     pred_RE((i-1)*nea+k)=9999.d0
                  end do
                  goto 654
               else
                  Yobs(nmes_cur+j)=Y1(j)
               end if

            end do

         elseif (idlink.eq.2) then ! Splines link


            bb=b1(nef+nvc+nwg+1)

            do kk=2,ntrtot
               splaa(kk-3)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
            end do

            do j=1,nmes(i)

               nmestot=nmestot+1

               ll=0
               if (Y(nmestot).eq.zitr(ntrtot-2)) then
                  ll=ntrtot-3
               end if
               som=0.d0
               do kk = 2,ntrtot-2
                  if ((Y(nmestot).ge.zitr(kk-1)).and.(Y(nmestot).lt.zitr(kk))) then
                     ll=kk-1
                  end if
               end do
               if (ll.lt.1.or.ll.gt.ntrtot-3) then
               do k=1,nmes(i)
                  resid_m(nmes_cur+k)=9999.d0
                  pred_m_g(nmes_cur+k)=9999.d0
                  Yobs(nmes_cur+k)=9999.d0
                  resid_ss(nmes_cur+k)=9999.d0
                  pred_ss_g(nmes_cur+k)=9999.d0
               end do
               do k=1,nea
                    pred_RE((i-1)*nea+k)=9999.d0
                end do
               goto 654
               end if
               if (ll.gt.1) then
                  do ii=2,ll
                     som=som+splaa(ii-3)
                  end do
               end if

               Y1(j)=bb+som +splaa(ll-2)*im2(nmestot) &
                 +splaa(ll-1)*im1(nmestot)+splaa(ll)*im(nmestot)

               Yobs(nmes_cur+j)=Y1(j)


            end do

         end if
!     debut du calcul de la vraisemblance

!     cas 1 : ng=1

         if (ng.eq.1) then


            Valea=0.d0
            VC=0.d0
            P=0.d0
            CovDev=0.d0

            P=MATMUL(Z,Ut)
            Valea=MATMUL(Ut,transpose(P))
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Corr+SigmaE

!covDev covariance matrix for individual deviation
            CovDev=MATMUL(P,transpose(P))+Corr

!     Vi en vecteur

            jj=0
            do j=1,nmes(i)
               do k=j,nmes(i)
                  jj=j+k*(k-1)/2
                  Vi(jj)=VC(j,k)
               end do
            end do

            CALL dsinv(Vi,nmes(i),eps,ier,det)
            if (ier.eq.-1) then
               do j=1,nmes(i)
                  resid_m(nmes_cur+j)=9999.d0
                  pred_m_g(nmes_cur+j)=9999.d0
                  resid_ss(nmes_cur+j)=9999.d0
                  pred_ss_g(nmes_cur+j)=9999.d0
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
                  if (k.gt.1) then
                     b0(l)=b1(nprob+l-1)
                  end if
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
               pred_m_g(nmes_cur+j)=mu(j)

               resid_ss(nmes_cur+j)=Y1(j)-pred1(j)
               pred_ss_g(nmes_cur+j)=pred1(j)
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
               temp=temp+exp(DOT_PRODUCT(bprob,Xprob))
               pi(g)=exp(DOT_PRODUCT(bprob,Xprob))

            end do
            pi(ng)=1.d0/(1.d0+temp)
            do g=1,ng-1
               pi(g)=pi(g)*pi(ng)
            end do


           end if
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

! covDev covariance matrix for individual deviation

            CovDev=MATMUL(P,transpose(P))+Corr


!     Vi en vecteur
               jj=0
               do j=1,nmes(i)
                  do k=j,nmes(i)
                     jj=j+k*(k-1)/2
                     Vi(jj)=VC(j,k)
                  end do
               end do
               CALL dsinv(Vi,nmes(i),eps,ier,det)
               if (ier.eq.-1) then
                  do j=1,nmes(i)
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
                  pred_m_g((g-1)*nobs+nmes_cur+j)=mu(j)
                  pred_ss_g((g-1)*nobs+nmes_cur+j)=pred1(j)

                  resid_ss(nmes_cur+j)=resid_ss(nmes_cur+j) &
                      +ppi(i,g)*(Y1(j)-pred1(j))
                  resid_m(nmes_cur+j)=resid_m(nmes_cur+j)+pi(g)*(Y2(j))
               end do
               do k=1,nea
                  pred_RE((i-1)*nea+k)=pred_RE((i-1)*nea+k)+ppi(i,g)*err2(k)
               end do

            end do

         end if

 654  continue

         nmes_cur=nmes_cur+nmes(i)
      end do

! FIN BOUCLE SUJET


      end subroutine residualsc





! =============================================
! subroutine de creation de design matrix
! =============================================




      subroutine design_splines (ier)

      use communc

      implicit none

      integer ::j,i,jj,l,k,ier
      double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n,hn,hht

      ier=0
      jj=0
      l=0
      do i=1,ns
         do j=1,nmes(i)
            jj=jj+1
!     ou se trouve la valeur de zi

            do k = 2,ntrtot-2
               if ((Y(jj).ge.zitr(k-1)).and.(Y(jj).lt.zitr(k))) then
                  l=k-1
               endif

            end do


            if (Y(jj).eq.zitr(ntrtot-2)) then
               l=ntrtot-3
            end if

            ht2 = zitr(l+1)-Y(jj)
            htm= Y(jj)-zitr(l-1)
            ht = Y(jj)-zitr(l)
            ht3 = zitr(l+2)-Y(jj)
            hht = Y(jj)-zitr(l-2)
            h = zitr(l+1)-zitr(l)
            hh= zitr(l+1)-zitr(l-1)
            hn= zitr(l+1)-zitr(l-2)
            h2n=zitr(l+2)-zitr(l-1)
            h2= zitr(l+2)-zitr(l)
            h3= zitr(l+3)-zitr(l)

            if (Y(jj).ne.zitr(ntrtot-2)) then
               mm2(jj) = (3.d0*ht2*ht2)/(hh*h*hn)
               mm1(jj) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
               mm(jj)  = (3.d0*ht*ht)/(h3*h2*h)

            end if
            if (Y(jj).eq.zitr(ntrtot-2)) then
               mm2(jj) = 0.d0
               mm1(jj) = 0.d0
               mm(jj)  = 3.d0/h
            end if

            if (mm2(jj).lt.0.or.mm1(jj).lt.0.or.mm(jj).lt.0) then
                ier=-1
                goto 765
            end if

            im2(jj)=hht*mm2(jj)/(3.d0)+ h2n*mm1(jj)/(3.d0) &
             +h3*mm(jj)/(3.d0)
            im1(jj)=htm*mm1(jj)/(3.d0)+h3*mm(jj)/(3.d0)
            im(jj)=ht*mm(jj)/(3.d0)

         end do
      end do


765     continue

      end subroutine design_splines









!---------------------------------------------------------------

      !transfo estimee

!---------------------------------------------------------------

      subroutine transfo_estimee(b,npm,nsim,marker,transfY)

      use communc

      implicit none

      integer::kk,nsim,npm,j,k
      double precision,dimension(nsim)::marker,transfY
      double precision,dimension(ntrtot)::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,eps,pas,ytemp,cc1
      double precision, dimension(npm)::b,b1



      b1=0.d0
      eps=1.D-20
      do k=1,npm
         b1(k)=b(k)
      end do


!       write(*,*)'infos',minY,maxY,nsim,npm
!       write(*,*)'b',(b1(j),j=1,npm)



       marker=0.d0
       transfY=0.d0

       pas=(maxY-minY)/dble(nsim-1)
       j=1
       marker(1)=minY
       do while(j.lt.nsim)
           j=j+1
           marker(j)=marker(j-1)+pas
       end do
       marker(nsim)=maxY

!       write(*,*)(marker(j),j=1,nsim)

       if (idlink.eq.2) then
        splaa=0.d0


        splaa(1)=b1(nef+nvc+nwg+1)
        do kk=2,ntrtot
           splaa(kk)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
             end do

!           write(*,*)'ntr',(b1(nprob+nrisq
!     &       +nvarxevt+ncssg+ncg*ng+kk),kk=1,ntrtot)


            call   estim_splines_ssstd(nsim,splaa,marker,transfY)

        else if (idlink.eq.1) then

            aa1=exp(b1(nef+nvc+nwg+1))/ &
             (1+exp(b1(nef+nvc+nwg+1)))
            bb1=exp(b1(nef+nvc+nwg+2))/ &
             (1+exp(b1(nef+nvc+nwg+2)))
            bb1=aa1*(1.d0-aa1)*bb1
            cc1=b1(nef+nvc+nwg+3)
            dd1=abs(b1(npm))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1

            do j=1,nsim
                  ytemp=(marker(j)-minY+epsY)/(maxY-minY+2*epsY)
                  transfY(j)=(betai(aa,bb,ytemp)-cc1)/dd1
                  if (transfY(j).eq.999.d0) then
!                    write(*,*)'problem'
                  end if

               end do


        else if (idlink.eq.0) then

                 do j=1,nsim
                    transfY(j)=(marker(j)-b1(nef+nvc+nwg+1))/abs(b1(nef+nvc+nwg+2))
                 end do
        end if
        end subroutine transfo_estimee





!---------------------------------------------------------------

!      TRANSFOS SPLINES SANS STD ERR

!---------------------------------------------------------------





        subroutine estim_splines_ssstd(nsim,aa,test,transf)

       use communc

       implicit none

       integer::nsim,j,k,l
       double precision,dimension(nsim)::mmm,mmm1,mmm2 &
      ,iim,iim1,iim2,transf,test
       double precision,dimension(ntrtot)::aa,Xspl
       double precision ::ht,htm,ht2,ht3,hht,h,hh,h2,h3,h2n,hn


! matrice de transition pour delta-metho (carre des parms 2,..,ntr)

       do j=1,nsim
! ou se trouve la valeur
         l=0

         do k = 2,ntrtot-2
               if ((test(j).ge.zitr(k-1)).and.(test(j).lt.zitr(k)))then
                  l=k-1
               endif
            end do

           if (test(j).eq.zitr(ntrtot-2)) then
               l=ntrtot-3
            end if

!         if (l.lt.1.or.l.gt.ntrtot-1) then
!            write(*,*)'probleme estim splines',l
!            write(*,*)'j=',j,'test(j)',test(j)
!            stop
!         end if


               ht2 = zitr(l+1)-test(j)
               htm= test(j)-zitr(l-1)
               ht = test(j)-zitr(l)
               ht3 = zitr(l+2)-test(j)
               hht = test(j)-zitr(l-2)
               h = zitr(l+1)-zitr(l)
               hh= zitr(l+1)-zitr(l-1)
               hn= zitr(l+1)-zitr(l-2)
               h2n=zitr(l+2)-zitr(l-1)
               h2= zitr(l+2)-zitr(l)
               h3= zitr(l+3)-zitr(l)

               if (test(j).ne.zitr(ntrtot-2)) then
                  mmm2(j) = (3.d0*ht2*ht2)/(hh*h*hn)
                  mmm1(j) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
                  mmm(j)  = (3.d0*ht*ht)/(h3*h2*h)
               end if
               if (test(j).eq.zitr(ntrtot-2)) then
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
            transf(j)= dot_product(Xspl,aa)
      end do



      end subroutine estim_splines_ssstd























!=======================================================================
!============= AJOUT 1/02/2012 :
!               VRAIS ORDINALE POUR OUTCOME DISCRET
!=======================================================================




!-----------------------------------------------------------
!                   VRAIS_CONT
!------------------------------------------------------------




      subroutine vrais_cont(b,m,id,thi,jd,thj,vraiscont)


        use communc,only:ns,nmes
        use donnees_indivc,only:nmescur

        implicit none

        integer::m,i,id,jd
        double precision::thi,thj,vrais_cont_i,vraiscont
        double precision,dimension(m)::b

        vraiscont=0.d0
        nmescur=0
        do i=1,ns
           vraiscont = vraiscont + vrais_cont_i(b,m,id,thi,jd,thj,i)
           nmescur = nmescur + nmes(i)
        end do

        return

      end subroutine vrais_cont




!-----------------------------------------------------------
!                   VRAIS_DISCRET
!------------------------------------------------------------


      subroutine vrais_discret(b,m,id,thi,jd,thj,vraisdiscret)


        use communc,only:ns,nmes
        use donnees_indivc,only:nmescur

        implicit none

        integer::m,i,id,jd
        double precision::thi,thj,vrais_discret_i,vraisdiscret
        double precision,dimension(m)::b

        vraisdiscret=0.d0
        nmescur=0
        do i=1,ns
           vraisdiscret = vraisdiscret + vrais_discret_i(b,m,id,thi,jd,thj,i)
           nmescur = nmescur + nmes(i)
        end do

        return

      end subroutine vrais_discret



      double precision function vrais_discret_i(b,npm,id,thi,jd,thj,i)

        use communc
!        use optim
        use donnees_indivc

        IMPLICIT NONE
        integer ::i,j,k,l,m,g,l2,m2,npm,id,jd
        integer ::nmoins,nspl,ier,ii,kk,ll
        double precision,dimension(maxmes,nv) ::X00,X2
        double precision,dimension(nv) ::Xprob
        double precision,dimension(nea,nea) ::Ut
        double precision,dimension(nea) ::Xea
        double precision,dimension(npm) :: b
        double precision,dimension(npmtot) :: b1
        double precision,dimension(nv) :: b0,b2,bprob
        double precision :: expo,temp,vraisobsc,vraisdiscret,ytemp,ytemp1, &
             aa1,bb1,cc1,dd1,aa,bb,thi,thj,som,betai
        double precision,dimension(ng) :: pi
!        double precision,dimension(rangeY+1) :: valY
!        double precision,dimension(rangeY+1,3) :: imseuil
        double precision,dimension(rangeY) :: valY
        double precision,dimension(rangeY,3) :: imseuil
        double precision,dimension(-1:ntrtot-3)::splaa


        allocate(Ut1(nea,nea),mu(maxmes),Z(maxmes,nea),seuils(rangeY))

        b1=0.d0
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

!        write(*,*)'b1',b1

        ! ============ creation des seuils pour le calcul de vrais


        seuils=0.d0
        if (idlink.eq.0) then  ! Linear link
           do k=1,rangeY
              ytemp=(minY+k-b1(nef+nvc+nwg+1))/abs(b1(nef+nvc+nwg+2))
              ytemp1=(minY+k-1-b1(nef+nvc+nwg+1))/abs(b1(nef+nvc+nwg+2))
              seuils(k)=(ytemp+ytemp1)/2
           end do

        else if (idlink.eq.1) then  ! Beta link

           aa1=exp(b1(nef+nvc+nwg+1))/ &
                (1+exp(b1(nef+nvc+nwg+1)))
           bb1=exp(b1(nef+nvc+nwg+2))/ &
                (1+exp(b1(nef+nvc+nwg+2)))
           bb1=aa1*(1.d0-aa1)*bb1

           cc1=abs(b1(nef+nvc+nwg+3))

           dd1=abs(b1(nef+nvc+nwg+4))

           aa=aa1*aa1*(1-aa1)/bb1-aa1
           bb=aa*(1-aa1)/aa1

           do k=1,rangeY
!      creation des seuils
! soit milieu echelle Lambda + ATTENTION CHANGEMENT TAILLE VALY
!              ytemp=(minY+k-minY+epsY)/(maxY-minY+2*epsY)
!              ytemp1=(minY+k-1-minY+epsY)/(maxY-minY+2*epsY)
!              seuils(k)=(betai(aa,bb,ytemp1)+betai(aa,bb,ytemp)-2*cc1)/dble(2*dd1)

! soit milieu echelle Y
              ytemp=dble(minY+k-5.d-1-minY+epsY)/dble(maxY-minY+2*epsY)
              seuils(k)=(betai(aa,bb,ytemp)-cc1)/dble(dd1)
           end do

        else if (idlink.eq.2) then ! Splines link

          bb=b1(nef+nvc+nwg+1)
          splaa=0.d0
          do kk=2,ntrtot
             splaa(kk-3)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
          end do


          valY=minY
          imseuil=0.d0

! creation des seuils : soit ValY sur le discret:
!          nspl=rangeY+1
!          do k=1,rangeY+1
!             valY(k)=k-1
!          end do

         ! creation des seuils : soit ValY au milieu +0.5
          nspl=rangeY
          do k=1,rangeY
             valY(k)=minY+k-5.d-1
          end do

          call splines_seuils(nspl,valY,imseuil,ier)
          if(ier.eq.-1) then
!             write(*,*)'ier eq -1 vrais_i'
             vraisdiscret=-1.d9
             goto 2587
          endif

          Do k=1,rangeY

!      creation des seuils
! Si milieu de Lambda : changer def de Val(Y) au dessus + de commentariser la deuxième partie
! soit milieu echelle Y, alors mettre  "valY(k)=k-5.d-1" et garder que la première partie ci-dessous
             ll=0
             if (ValY(k).eq.zitr(ntrtot-2)) then
                ll=ntrtot-3
             end if
             som=0.d0
             do kk = 2,ntrtot-2
                if ((ValY(k).ge.zitr(kk-1)).and. &
                     (ValY(k).lt.zitr(kk))) then
                   ll=kk-1
                end if
             end do
!             if (ll.lt.1.or.ll.gt.ntrtot-3) then
!                write(*,*) 'probleme dans seuilcont splines'
!                write(*,*) 'll=',ll,'Y+1=',valY(k)
!                stop
!             end if
             if (ll.gt.1) then
                do ii=2,ll
                   som=som+splaa(ii-3)
                end do
             end if
             seuils(k)=bb+ som +splaa(ll-2)*imseuil(k,3)  &
                  +splaa(ll-1)*imseuil(k,2)+splaa(ll)*imseuil(k,1)

! Que si milieu intervalle echelle Lambda
!             ! puis h(k+1)
!
!             ll=0
!             if (ValY(k+1).eq.zitr(ntrtot-2)) then
!                ll=ntrtot-3
!             end if
!             som=0.d0
!             do kk = 2,ntrtot-2
!                if ((ValY(k+1).ge.zitr(kk-1)).and. &
!                     (ValY(k+1).lt.zitr(kk)))then
!                   ll=kk-1
!                end if
!             end do
!             if (ll.lt.1.or.ll.gt.ntrtot-3) then
!             end if
!             if (ll.gt.1) then
!                do ii=2,ll
!                   som=som+splaa(ii-3)
!                end do
!             end if
!
!             seuils(k)=seuils(k)+ bb+ som +splaa(ll-2)*imseuil(k+1,3)  &
!                  +splaa(ll-1)*imseuil(k+1,2)+splaa(ll)*imseuil(k+1,1)
!
!             seuils(k)=seuils(k)/2.d0

 
          end do


!      if (i.eq.1) then
!         do k=1,nspl
!            write(*,*)'imseuil',(imseuil(k,l),l=1,3)
!         end do
!         write(*,*)'ValY',(valY(k),k=1,nspl)
!         write(*,*)"seuils",(seuils(k),k=1,rangeY)
!       end if


       end if



!----------- rappel des parametres utilises ---------


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


! ----------- boucle sur les individus -------------
      vraisdiscret=0.d0
      numpat=i

      ! -------- creation de Vi = ZiGZi'+se*seIni ----------
! creation de Zi

      Z=0.d0
      l=0
      do k=1,nv
         if (idea(k).eq.1) then
            l=l+1
            do j=1,nmes(i)
               Z(j,l)=dble(X(nmescur+j,k))
            end do
         end if
      end do

      ! cas 1 : ng=1

      if (ng.eq.1) then


         b0=0.d0
         l=0
         X00=0.d0
         do k=1,nv
            if (idg(k).ne.0) then
               l=l+1
               do j=1,nmes(i)
                  X00(j,l)=dble(X(nmescur+j,k))
               end do
               ! idg ne 0 pour l'intercept forcement donc on met le parm a 0
               if (k.eq.1) then
                  b0(l)=0.d0
               else
                  b0(l)=b1(nprob+l-1)
               end if
            end if
         end do


         mu=0.d0
         mu=matmul(X00,b0)

!         write(*,*)'mu',mu

         if (nea.gt.0) then
            temp=vraisobsc()
         else
            Xea=0.d0
            call vraistotc(nea,Xea,nf,temp)
         end if


         if(temp.lt.1.d-300) then
!            write(*,*)'temp a 1.d-300',i
            temp=1.d-300
         end if

         vraisdiscret=vraisdiscret+log(temp)


         !            if(thi.eq.0.and.thj.eq.0) then
         !               write(*,*)'i',i,temp,log(temp),vrais
         !            end if

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
                  Xprob(1+l)=X(nmescur+1,k)
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
                  X2(j,l)=dble(X(nmescur+j,k))
               end do
            else if (idg(k).eq.1) then
               m=m+1
               do j=1,nmes(i)
                  X00(j,m)=dble(X(nmescur+j,k))
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


            !               if (thi.eq.0.and.thj.eq.0) then
            !                if (i.eq.1) then
            !                    write(*,*)'g',g,b2
            !                    write(*,*)'g',g,b0
            !                    stop
            !                    end if
            !                    end if


            ! variance covariance si spec aux classes :

            Ut1=Ut
            if (nwg.ne.0.and.g.ne.ng) then
               Ut1=Ut*abs(b1(nef+nvc+g))
            end if

            mu=0.d0
            mu=matmul(X00,b0)+matmul(X2,b2)

            if (nea.gt.0) then
               temp=vraisobsc()
            else
               Xea=0.d0
               call vraistotc(nea,Xea,nf,temp)
            end if

            if (temp.lt.1.d-300) then

!               write(*,*)'temp a 1.d-300',i
               temp=1.d-300
            end if
            expo = expo+pi(g)*temp


         end do
         vraisdiscret=vraisdiscret+log(expo)

      end if
      ! FIN BOUCLE SUJET

 2587 continue
      deallocate(Ut1,mu,Z,seuils)
      !            if(thi.eq.0.and.thj.eq.0) then
      !               write(*,*)'funcpao',funcpao
      !            end if

      vrais_discret_i=vraisdiscret
      return

    end function vrais_discret_i


!=============== FIN VRAIS_DISCRET_i ======================
!==========================================================


!==========================================================
!======= pour les splines : definition des bases splines


      subroutine splines_seuils(nspl,valY,imseuil,ier)


      use communc

      implicit none
      integer ::k,l,nspl,ier,i
      double precision,dimension(nspl)::valY
      double precision,dimension(nspl,3)::imseuil
      double precision,dimension(3)::mmseuil
      double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n,hn,hht




      imseuil=0.d0
      ier=0
      l=0
      do i=1,nspl
         mmseuil=0.d0
         do k = 2,ntrtot-2
            if ((ValY(i).ge.zitr(k-1)).and.(ValY(i).lt.zitr(k))) then
               l=k-1
            endif
         end do

         if (ValY(i).eq.zitr(ntrtot-2)) then
            l=ntrtot-3
         end if

         ht2 = zitr(l+1)-ValY(i)
         htm= ValY(i)-zitr(l-1)
         ht = ValY(i)-zitr(l)
         ht3 = zitr(l+2)-ValY(i)
         hht = ValY(i)-zitr(l-2)
         h = zitr(l+1)-zitr(l)
         hh= zitr(l+1)-zitr(l-1)
         hn= zitr(l+1)-zitr(l-2)
         h2n=zitr(l+2)-zitr(l-1)
         h2= zitr(l+2)-zitr(l)
         h3= zitr(l+3)-zitr(l)

         if (ValY(i).ne.zitr(ntrtot-2)) then
            mmseuil(3) = (3.d0*ht2*ht2)/(hh*h*hn)
            mmseuil(2) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
            mmseuil(1)  = (3.d0*ht*ht)/(h3*h2*h)
         end if
         if (ValY(i).eq.zitr(ntrtot-2)) then
            mmseuil(3) = 0.d0
            mmseuil(2) = 0.d0
            mmseuil(1)  = 3.d0/h
         end if

         if (mmseuil(3).lt.0.or.mmseuil(2).lt.0.or.mmseuil(1).lt.0) then
            ier=-1
!            write(*,*)'mmseuil',mmseuil(3),mmseuil(2),mmseuil(1)
            goto 765
         end if

         imseuil(i,3)=hht*mmseuil(3)/(3.d0)+ h2n*mmseuil(2)/(3.d0) &
              +h3*mmseuil(1)/(3.d0)
         imseuil(i,2)=htm*mmseuil(2)/(3.d0)+h3*mmseuil(1)/(3.d0)
         imseuil(i,1)=ht*mmseuil(1)/(3.d0)

      end do
765     continue

      end subroutine splines_seuils



! ================================================================
!
!
!     Integration numerique : calcul de funcpa pour 1 sujet
!     dans vraisobs()
!
! =================================================================

      double precision function vraisobsc()

!        use lois_normalesc
        use donnees_indivc
        use communc


        implicit none
        integer::j
        external ::vraistotc
        INTEGER ::NDIM2, NF2, MINPTS, MAXPTS, RESTAR
        INTEGER ::NEVAL, IFAIL
        DOUBLE PRECISION ::EPSABS, EPSREL
        DOUBLE PRECISION,dimension(2) ::RESULT
        DOUBLE PRECISION,dimension(2) :: ABSERR2
        DOUBLE PRECISION,dimension(1000) ::WORK
        double precision,dimension(1)::Xea
        integer :: npg
        !     fin de declaration des variables de mvndst
!        integer :: npt
!        double precision::deter,rl,uniran,phid
!        integer ::line,nloop,iquad,iprint,ifault
        double precision :: funvls

        double precision, dimension(2,51)::Gauss

        ndim2=nea
        nf2=nf
        result=0.d0

        if (nea.gt.1) then

           MINPTS=30
           MAXPTS=500               !  100 !500
           EPSABS=1.d-100
           EPSREL=1.d-100
           RESTAR=0
           call  hrmsym( NDIM2, NF2, MINPTS, MAXPTS, vraistotc, EPSABS, &
                EPSREL, RESTAR, RESULT, ABSERR2, NEVAL, IFAIL, WORK)

           if (result(1).le. 1.d-300) then
              result(1)=1.d-300
           end if
           vraisobsc=result(1)

           return

        else
              npg=30
              ! on definit les points
              call gausshermite(gauss,npg)

              ! boucle pour faire l'integration multiple
              do j=1,npg
                 Xea(1)=gauss(1,j)
                 call vraistotc(nea,Xea,nf2,funvls)
                 result(1)=result(1)+funvls*gauss(2,j)
              end do

              if (result(1).le. 1.d-300) then
                 result(1)=1.d-300
              end if

              vraisobsc=result(1)
              return

        end if




      end function vraisobsc






! ===============================================================
!
!     Vraisemblance sachant les effets aleatoires
! ===============================================================




      subroutine vraistotc(NDIM2,Xea,nf2,FUNVLS)

      use donnees_indivc
      use communc
!      use lois_normales

      implicit none
      double precision ::vraisind,gamma1,gamma0,alnorm
      double precision :: funvls,sup,inf
      integer::ndim2,nf2,i,k,j,ind
      logical ::upper
      double precision,dimension(NDIM2)::Xea
      double precision,dimension(NDIM2)::ui
      double precision, dimension(maxmes)::mu1
      i=numpat
      nf2=nf2




      mu1=0.d0
      ui=0.d0
      if (NDIM2.ge.1) then
         ui=MATMUL(Ut1,Xea)
         mu1=mu+MATMUL(Z,ui)
      else
         mu1=mu
      end if

!      if (i.lt.3) then
!    write(*,*)'seuils',seuils
!         write(*,*)'mu1',mu1
!         write(*,*)'mu',mu
!         write(*,*)'Xea',Xea
!      end if



      vraisind=1.d0

      upper=.false.

      do j=1,nmes(i)
         ind=0
         if (Y(nmescur+j).eq.minY) then
            gamma0=(seuils(1)-mu1(j))
 !           write(*,*)'gamma0',gamma0
            vraisind=vraisind*(alnorm(gamma0,upper))
            ind=1
         else
            sup=seuils(1)
            inf=sup
            do k=1,rangeY-1
               sup=seuils(k+1)
               if (Y(nmescur+j).eq.dble(minY+k)) then
                  gamma1=(sup-mu1(j))
                  gamma0=(inf-mu1(j))
 !           write(*,*)'gamma0',gamma0,gamma1
                  vraisind=vraisind*(alnorm(gamma1,upper) &
                 -alnorm(gamma0,upper))
                  ind=1
               end if
               inf=sup
            end do
            if (Y(nmescur+j).eq.maxY) then
               gamma0=(sup-mu1(j))
   !         write(*,*)'gamma0',gamma0
               vraisind=vraisind*(1.d0-alnorm(gamma0,upper))
               ind=1
            end if

         end if

      end do

!      if (i.lt.3) then
!         write(*,*)'vraistotc',vraisind
!      end if

      FUNVLS=vraisind

      end subroutine vraistotc




!===========================================================
!      DERIV pour UACV
!===========================================================


      subroutine computUACV(b,m,rlindiv,vopt,UACV)

! Calcul du gradient et de la matrice remplacant la hessienne
! par Fisher scoring empirique
        use communc
        use donnees_indivc,only:nmescur

        IMPLICIT NONE

        integer::m,i,k,id,j
        double precision::vrais_discret_i,vrais_cont_i, &
             rldiscret,UACV,trace,th0,thn,th,rlcont
        double precision,dimension(m,1)::Uscore, Uscore2
        double precision,dimension(m)::b
        double precision,dimension(m*(m+1)/2)::vopt
        double precision,dimension(m,m)::J_cond,H_1,MAT
        double precision,dimension(ns)::rlindiv

        J_cond=0.d0
        th0=0.d0
        rlindiv=0.d0
        id=0



! Calcul des gradients par sujets et totaux
        nmescur=0
        rldiscret=0.d0
        rlcont=0.d0
        DO i=1,ns
           Uscore=0.d0
           Uscore2=0.d0
           rlindiv(i)=vrais_discret_i(b,m,id,th0,id,th0,i)
           rldiscret=rldiscret+rlindiv(i)
           rlcont=rlcont+vrais_cont_i(b,m,id,th0,id,th0,i)
           do k=1,m
              th=DMAX1(1.d-6, 1.d-4 * DABS(b(k)))
              thn=-1.D0*th
              Uscore(k,1)=-(vrais_discret_i(b,m,k,th,id,th0,i) &
                   - vrais_discret_i(b,m,k,thn,id,th0,i))/(2.d0*th)
              Uscore2(k,1)=-(vrais_cont_i(b,m,k,th,id,th0,i) &
                   - vrais_cont_i(b,m,k,thn,id,th0,i))/(2.d0*th)
           END DO

!           write(*,*)'Uscore',Uscore
!           write(*,*)'Uscore2',Uscore2

           J_cond=J_cond+MATMUL(Uscore,transpose(Uscore2))
           nmescur = nmescur + nmes(i)
        end do

        H_1 = 0.d0
        do j = 1,m
           do k =j,m
              H_1(j,k) = vopt(j+k*(k-1)/2)
              H_1(k,j) = vopt(j+k*(k-1)/2)
           end do
        end do

        MAT=MATMUL(H_1,J_cond)
        trace=0.d0
        do k=1,m
           trace=trace+MAT(k,k)
        end do

        rldiscret=rldiscret/dble(ns)
        UACV=-rldiscret+trace/dble(ns-1)

!        write(*,*)'rldiscret',rldiscret,trace
!        write(*,*)'rlcont',ns,rlcont

!        write(*,*)'MAT',MAT
 !       write(*,*)'J_cond',J_cond
  !      write(*,*)'H_1',H_1

        return

        end subroutine computUACV







!=================== VRAIS_CONT_i : utilise pour l'instant par computUACV mais pourrait etre dans FUNCPA





!-----------------------------------------------------------
!                       VRAIS_CONT_i
!------------------------------------------------------------


      double precision function vrais_cont_i(b,npm,id,thi,jd,thj,i) 

      use communc
!      use optim
      use donnees_indivc,only:nmescur
      IMPLICIT NONE
      integer ::i,j,k,l,m,g,l2,m2,id,jd,jj,npm,nmestot,ll,ii
      integer ::ier,nmoins,kk
      double precision,dimension(maxmes,nv) ::Z,P,X00,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nv,nv) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Se
      double precision,dimension(npm) :: b
      double precision,dimension(npmtot) :: b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: vrais,eps,det,som,thi,thj,temp
      double precision ::Y4,expo,jacobien,beta_densite,ytemp
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3
      double precision,dimension(ng) :: pi
      double precision,dimension(-1:(ntrtot-3))::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,cc1


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

! ----------- boucle sur les individus -------------
      nmestot=nmescur
      vrais=0.d0
      jacobien=0.d0
! -------- creation de Vi = ZiGZi'+se*seIni ----------
! creation de Zi

         Z=0.d0
         l=0
         do k=1,nv
            if (idea(k).eq.1) then
               l=l+1
               do j=1,nmes(i)
                  Z(j,l)=dble(X(nmescur+j,k))
               end do
            end if
         end do
! creation de s2*I et Y1

         Se=0.d0
         Y1=0.d0

         if (idlink.eq.0) then  ! Linear link

            do j=1,nmes(i)
               nmestot=nmestot+1
               Se(j,j)=1 ! erreur de mesure = parm de transfo
               Y1(j)=(dble(Y(nmestot))-b1(nef+nvc+nwg+1))/abs(b1(nef+nvc+nwg+2))

               jacobien = jacobien - log(b1(nef+nvc+nwg+2))
            end do

         else if (idlink.eq.1) then  ! Beta link


            aa1=exp(b1(nef+nvc+nwg+1))/ &
             (1+exp(b1(nef+nvc+nwg+1)))
            bb1=exp(b1(nef+nvc+nwg+2))/ &
             (1+exp(b1(nef+nvc+nwg+2)))
            bb1=aa1*(1.d0-aa1)*bb1

            cc1=abs(b1(nef+nvc+nwg+3))

            dd1=abs(b1(nef+nvc+nwg+4))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1

            do j=1,nmes(i)

               nmestot=nmestot+1
               Se(j,j)=1.d0

               ytemp=(dble(Y(nmestot))-minY+epsY)/(maxY-minY+2*epsY)
               Y1(j)=(betai(aa,bb,ytemp)-cc1)/dd1


               if (Y1(j).eq.999.d0) then
                  vrais_cont_i=-1.d9
                  goto 654
               end if

               jacobien = jacobien + log(abs(beta_densite(ytemp,aa,bb))/dd1)
               jacobien=jacobien-log(abs(maxY-minY+2*epsY))
            end do

         else if (idlink.eq.2) then ! Splines link



            bb=b1(nef+nvc+nwg+1)

            do kk=2,ntrtot
            splaa(kk-3)=b1(nef+nvc+nwg+kk)*b1(nef+nvc+nwg+kk)
            end do

            do j=1,nmes(i)

               nmestot=nmestot+1
               Se(j,j)=1.d0

               ll=0
               if (Y(nmestot).eq.zitr(ntrtot-2)) then
                  ll=ntrtot-3
               end if
               som=0.d0
               do kk = 2,ntrtot-2
                  if ((Y(nmestot).ge.zitr(kk-1)).and. &
                      (Y(nmestot).lt.zitr(kk))) then
                     ll=kk-1
                  end if
               end do
               if (ll.lt.1.or.ll.gt.ntrtot-3) then
                vrais_cont_i=-1.d9
                goto 654
               end if
               if (ll.gt.1) then
                  do ii=2,ll
                     som=som+splaa(ii-3)
                  end do
               end if

               Y1(j)=bb+ som +splaa(ll-2)*im2(nmestot)  &
                +splaa(ll-1)*im1(nmestot)+splaa(ll)*im(nmestot)

               jacobien = jacobien + log(splaa(ll-2)*mm2(nmestot) &
                 +splaa(ll-1)*mm1(nmestot)+splaa(ll)*mm(nmestot))
            end do
         end if


!         if (i.lt.3)then
!            write(*,*)'nmes',nmes(i),b1((nef+nvc+nwg+1):npm),nef
!            write(*,*)'Y1',Y1
!            write(*,*)'Y',Y(nmestot-nmes(i)+1:nmestot)
!         end if


! creation de P=Zi*Ut et V=P*P' que si non spec aux classes

         if (nwg.eq.0.OR.NG.EQ.1) then
            P=0.d0
            P=MATMUL(Z,Ut)
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Se

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
               vrais_cont_i=-1.d9
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
                     X00(j,l)=dble(X(nmescur+j,k))
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
                    Xprob(1+l)=X(nmescur+1,k)
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
                        X2(j,l)=dble(X(nmescur+j,k))
                  end do
               else if (idg(k).eq.1) then
                  m=m+1
                  do j=1,nmes(i)
                     X00(j,m)=dble(X(nmescur+j,k))
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


!               if (thi.eq.0.and.thj.eq.0) then
!                if (i.eq.1) then
!                    write(*,*)'g',g,b2
!                    write(*,*)'g',g,b0
!                    stop
!                    end if
!                    end if


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
                  VC=MATMUL(P,transpose(P))+Se

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
                     vrais_cont_i=-1.d9
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
            vrais=vrais+2*log(expo)

         end if

! FIN BOUCLE SUJET

      vrais_cont_i=vrais/2.D0+jacobien

 654  continue

      return

      end function vrais_cont_i




      subroutine logliklcmmcont(Y0,X0,Prior0,idprob0,idea0,idg0,idcor0,ns0,ng0,nv0, &
           nobs0,nea0,nmes0,idiag0,nwg0,ncor0,npm0,b0, &
           ppi0,resid_m ,resid_ss,pred_m_g,pred_ss_g,pred_RE, &
           epsY0,idlink0,nbzitr0,zitr0,marker,transfY,nsim0,Yobs,Ydiscret,vraisdiscret,&
           UACV,rlindiv,V,fix0,nfix0,bfix0,estim0,loglik)
        
!      use parameters
      use communc
!      use optim

      IMPLICIT NONE

        !Declaration des variables en entree
      integer,intent(in)::nv0,Ydiscret,estim0,nfix0
      double precision,intent(in)::epsY0
      integer,intent(in)::idlink0,nbzitr0
      double precision,dimension(nbzitr0),intent(in)::zitr0

      integer, intent(in)::ns0,ng0,nobs0,idiag0,nwg0,npm0,nea0,nsim0,ncor0
      integer, dimension(nv0),intent(in)::idea0,idg0,idprob0,idcor0
      integer, dimension(ns0),intent(in)::nmes0,prior0
      double precision,dimension(nobs0),intent(in)::Y0
      double precision,dimension(nobs0*nv0),intent(in)::X0
      integer,dimension(npm0+nfix0),intent(in)::fix0
      double precision,dimension(nfix0),intent(in)::bfix0
          
        !Declaration des variable en entree et sortie
      double precision, dimension(npm0+nfix0), intent(in) :: b0
      double precision,dimension(npm0*(npm0+3)/2)::V ! pour computeUACV
        !Declaration des variables en sortie
      double precision,intent(out)::loglik,vraisdiscret,UACV
      double precision,dimension(ns0*ng0),intent(out)::ppi0
      double precision,dimension(ns0),intent(out)::rlindiv
      double precision,dimension(nobs0),intent(out)::resid_m,resid_ss,Yobs
      double precision,dimension(nobs0*ng0),intent(out)::pred_m_g
      double precision,dimension(nobs0*ng0),intent(out)::pred_ss_g
      double precision,dimension(ns0*nea0),intent(out)::pred_RE
      double precision,dimension(nsim0),intent(out)::marker,transfY
        !Variables locales
      integer::jtemp,i,g,j,ij,ier,k,ktemp,ig,nmestot,it,nbfix
      integer::k2,id,jd,npmtot0
      double precision::thi,thj
      double precision,dimension(ns0,ng0)::PPI
      double precision,dimension(npm0+nfix0)::btot
      double precision,external::funcpac

         !print*,"debut fortran"
! sorties initialisees

       ppi0=0.d0
       pred_ss_g=0.d0
       pred_m_g=0.d0
       pred_RE=0.d0
       marker=0.d0
       transfY=0.d0
       resid_m=0.d0
       resid_ss=0.d0
       loglik=0.d0
       Yobs=0.d0
       vraisdiscret=0.d0
       UACV=0.d0
       rlindiv=0.d0



! en prevision de l'extension au conjoint
      nrisq=0
      nvarxevt=0
! fin en prevision

      maxmes=0
      do i=1,ns0
         if (nmes0(i).gt.maxmes) then
            maxmes=nmes0(i)
         end if
      end do


      minY=zitr0(1)
      maxY=zitr0(nbzitr0)

      rangeY=0
      if (Ydiscret.eq.1) rangeY=INT(maxY)-INT(minY)


      epsY=epsY0
      idlink=idlink0
      if (idlink.eq.0) ntrtot=2
      if (idlink.eq.1) ntrtot=4
      if (idlink.eq.3) ntrtot=INT(zitr0(nbzitr0))-INT(zitr0(1))
      if (idlink.eq.2) then
         ntrtot=nbzitr0+2
         allocate(zitr(-1:(ntrtot)))

         allocate(mm(nobs0),mm1(nobs0),mm2(nobs0),im(nobs0),im1(nobs0),im2(nobs0))

         zitr(1:nbzitr0)=zitr0(1:nbzitr0)
         zitr(-1)=zitr(1)
         zitr(0)=zitr(1)
         zitr(ntrtot-1)=zitr(ntrtot-2)
         zitr(ntrtot)=zitr(ntrtot-1)
        else
         allocate(zitr(1))
         allocate(mm(1),mm1(1),mm2(1),im(1),im1(1),im2(1))
      end if


      allocate(Y(nobs0),idprob(nv0),X(nobs0,nv0) &
      ,idea(nv0),idg(nv0),idcor(nv0),nmes(ns0),prior(ns0))

             !print*,"apres alloc"
! enrigstrement pour les modules
      ns=ns0
      ng=ng0
      nv=nv0
      nobs=nobs0
          ncor=ncor0
      if (nwg0.eq.0) then
         nwg=0
      else
         nwg=ng-1
      end if

      idiag=idiag0

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
         idcor(k)=idcor0(k)

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
         
!         write(*,*)'X k:',X(1:50,k)
      end do


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

      if((ng.eq.1.and.ncg.gt.0).or.(ng.eq.1.and.nprob.gt.0)) then
         go to 1589
      end if


!  nb effets fixes = nb effets fixes sans melange
!                  + ng fois le nb de var dans melange


      if (idiag.eq.1) then
         nvc=nea
      else if(idiag.eq.0) then
         nvc=nea*(nea+1)/2
      end if

      nef=nprob+ncssg+ncg*ng+nvarxevt+nrisq-1
      npmtot=nef+nvc+nwg+ntrtot+ncor

             
    !  if (nwg.gt.0) then
    !     do i=1,nwg
    !        btot(nef+nvc+i)=abs(btot(nef+nvc+i))
    !     end do
    !  end if


      if (idlink.eq.2) then
         call design_splines(ier)
         if (ier.eq.-1) then
            go to 1589
         end if
      end if



! lancement de l'optimisation

      IF (estim0.eq.1) then

         id=0
         jd=0
         thi=0.d0
         thj=0.d0
         
         loglik=funcpac(b0,npm0,id,thi,jd,thj)

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
           call transfo_estimee(btot,npmtot,nsim0,marker,transfY)
!         end if 

! probas posteriori

!     write(*,*)'avant postprobc'


         !if (istop.eq.1) then
            if (ng.gt.1) then
               call postprobc(btot,npmtot,PPI)
            end if



!            write(*,*)'avant residuals'

            call residualsc(btot,npmtot,ppi,resid_m,pred_m_g,resid_ss &
          ,pred_ss_g,pred_RE,Yobs)
            ig=0
            ij=0
            do i=1,ns
               do g=1,ng0
                  ig=ig+1
                  ppi0(ig)=PPI(i,g)
               end do
            end do


            if (Ydiscret.eq.1.and.ncor.eq.0) then
               id=0
               thi=0.d0
               call vrais_discret(b0,npm0,id,thi,id,thi,vraisdiscret)

!               write(*,*)'avant UACV'
               call computUACV(b0,npm0,rlindiv,v,UACV)
!               write(*,*)'arpes UACV'
            end if

         end if




!      write(*,*)'avant deallocate'

 1589 continue

      deallocate(Y,X,idprob,idea,idg,idcor,nmes,prior)

      deallocate(zitr,mm,mm1,mm2,im,im1,im2)

      deallocate(fix,bfix)


!      write(*,*)'fin'
      return
    end subroutine logliklcmmCont


