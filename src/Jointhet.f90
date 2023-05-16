

    !!! signification des variables de commun_comp !!!

!!$  ns : nombre de sujet dans l'echantillon
!!$  ng : nombre de classes latentes
!!$  nv : nombre total de covariables
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
!!$  ntrtot : nb de prm pr transformation
!!$  nvarxevt : nb total de prm pr survie (sans prm hazard)
!!$  nef : nb d effets fixes a estimer dans MM
!!$  idlink : type de transformation (0=linear,1=beta,2=splines)
!!$  nbevt : nombre d'evenements concurrents                   
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



     module commun_comp

      implicit none
      integer,save ::ns,ng,nv,idiag,ncssg,nvc,nea,ncg,nwg,ncor  &
      ,nprob,nvarprob,maxmes,nobs,ntrtot,nvarxevt,nef &
      ,idlink,nbevt,logspecif,nmes_curr,idtrunc,nvdepsurv,nrisqtot,nxevt,npmtot
      integer,dimension(:),allocatable,save::risqcom,typrisq,nz,nprisq,nrisq,nxevtspec,nevtparx,nxcurr
      double precision,dimension(:),allocatable,save::Y,uniqueY
      integer,dimension(:),allocatable,save::indiceY
      double precision,dimension(:,:),allocatable,save ::X,zi
      double precision,dimension(:),allocatable,save::Tsurv0,Tsurv,Tsurvint
      integer,dimension(:),allocatable,save::Devt,ind_survint
      integer,dimension(:),allocatable,save ::idea,idg,idprob,idcor
      integer,dimension(:),allocatable,save::idcom,idspecif,idtdv
      integer,dimension(:),allocatable,save :: nmes,prior
      double precision,dimension(:,:),allocatable,save ::pprior
      double precision,save :: minY,maxY,epsY
      double precision,dimension(:),allocatable,save::zitr
      double precision,dimension(:),allocatable,save::mm,mm1,mm2,im,im1,im2
      integer,save::rangeY,nvalSPL
      double precision,dimension(:),allocatable,save::Tmm,Tmm1,&
           Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,Tmm01,Tmm02,Tmm03,Tim0, &
           Tim01,Tim02,Tim03,Tmmt,Tmmt1,Tmmt2,Tmmt3,Timt,Timt1,&
           Timt2,Timt3
      double precision,dimension(:),allocatable,save::Tmm_est,Tmm1_est &
          ,Tmm2_est,Tmm3_est,Tim_est,Tim1_est,Tim2_est,Tim3_est
      double precision,dimension(:),allocatable,save::brisq_est
      double precision,save::vrais_surv
      integer,dimension(:),allocatable,save::fix
      double precision,dimension(:),allocatable,save::bfix
    end module commun_comp





!-----------------------------------------------------------
!                       VRAIS_COMP_i
!------------------------------------------------------------


      double precision function vrais_comp_i(b,npm,id,thi,jd,thj,i)

!      use parameters
      use commun_comp
!      use optim
     ! use donnees_indivm,only:nmescur

      IMPLICIT NONE
      integer ::i,j,k,l,m,g,l2,m2,id,jd,jj,npm,ll,ii
      integer ::ier,nmoins,kk,j1,j2,q
      integer::ke,sumnrisq,nevtxcurr,nxevtcurr
      double precision,dimension(maxmes,nv) ::X00,X2
      double precision,dimension(maxmes,nea) ::Z,P
      double precision,dimension(nvarprob) ::Xprob,bprob 
      double precision,dimension(nea,nea) ::Ut,Ut1
      double precision,dimension(maxmes,maxmes) ::VC,Corr
      double precision,dimension(npm) :: b
      double precision,dimension(npmtot) :: b1
      double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
      double precision,dimension(nv) :: b0,b2
      double precision,dimension(nxevt)::Xevt,bevt
      double precision,dimension(nbevt)::bevtint
      double precision,dimension(maxval(nrisq))::brisq

      double precision :: vrais,eps,det,som,thi,thj,temp,eta0,fevt,varexpsurv
      double precision ::Y4,expo,jacobien,beta_densite,ytemp &
           ,surv_glob,surv0_glob
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
      double precision,dimension(ng) :: pi
      double precision,dimension(-1:ntrtot-3)::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,cc1,entretard,retard
      double precision,dimension(ng,nbevt)::risq,surv,surv0,survint



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
               surv0(g,ke)=surv0(g,ke)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
            end if


        end do
        sumnrisq = sumnrisq + nrisq(ke)
     end do


   
  
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


      vrais=0.d0
      jacobien=0.d0

      entretard=0.d0
      retard=0.d0

! -------- creation de Vi = ZiGZi'+se*seIni ----------
! creation de Zi

         Z=0.d0
         l=0
         do k=1,nv
            if (idea(k).eq.1) then
               l=l+1
               do j=1,nmes(i)
                  Z(j,l)=dble(X(nmes_curr+j,k))
               end do
            end if
         end do

!matrice Corr=Ri+s2*I

        Corr=0.d0
        tcor=0.d0
        if (ncor.gt.0) then
           do k=1,nv
              if (idcor(k).eq.1) then
                 do j=1,nmes(i)
                    tcor(j) = X(nmes_curr+j,k)
                 end do
              end if
           end do
         do j1=1,nmes(i)
            do j2=1,nmes(i)
               if (ncor.eq.1) then
                  Corr(j1,j2) = Corr(j1,j2)+b1(nprob+nrisqtot+nvarxevt &
                       +nef+nvc+nwg+ncor)**2*min(tcor(j1),tcor(j2))
               else if (ncor.eq.2) then
                  Corr(j1,j2) = Corr(j1,j2)+b1(nprob+nrisqtot+nvarxevt &
                       +nef+nvc+nwg+ncor)**2*exp(-b1(nprob+nrisqtot+nvarxevt &
                       +nef+nvc+nwg+1)*abs(tcor(j1)-tcor(j2)))
               end if
               if(j1.eq.j2) then
                  if(idlink.eq.-1) then 
                     Corr(j1,j2) = Corr(j1,j2)+b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)**2
                  else
                     Corr(j1,j2) = Corr(j1,j2) + 1 !variance de l'erreur
                  end if
               end if
            end do
         end do
         end if

         if(ncor.eq.0) then
            do j1=1,nmes(i)
               if(idlink.ne.-1) Corr(j1,j1) = 1 ! variance de l'erreur
               if(idlink.eq.-1) Corr(j1,j1) = b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)**2
            end do
         end if



         ! creation de Y1
         Y1=0.d0
         splaa=0.d0

         if (idlink.eq.0) then  ! Linear link

            do j=1,nmes(i)
               Y1(j)=(dble(Y(nmes_curr+j))-b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg &
                    +ncor+1))/abs(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2))
               jacobien = jacobien-log(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2))
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

               ytemp=(dble(Y(nmes_curr+j))-minY+epsY)/(maxY-minY+2*epsY)
               Y1(j)=(betai(aa,bb,ytemp)-cc1)/dd1


               if (Y1(j).eq.999.d0) then
                  !print*,"pb betai"
                  vrais_comp_i=-1.d9
                  goto 654
               end if

               jacobien = jacobien + log(abs(beta_densite(ytemp,aa,bb))/dd1)
               jacobien=jacobien-log(abs(maxY-minY+2*epsY))

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
               
               if (Y(nmes_curr+j).eq.zitr(ntrtot-2)) then
                  ll=ntrtot-3
               end if

               som=0.d0
               do kk = 2,ntrtot-2
                  if ((Y(nmes_curr+j).ge.zitr(kk-1)).and. &
                      (Y(nmes_curr+j).lt.zitr(kk))) then
                     ll=kk-1
                  end if
               end do

               if (ll.lt.1.or.ll.gt.ntrtot-3) then
                  !print*,"pb ll"
                vrais_comp_i=-1.d9
                goto 654
               end if
               if (ll.gt.1) then
                  do ii=2,ll
                     som=som+splaa(ii-3)
                  end do
               end if



               Y1(j)=eta0+som +splaa(ll-2)*im2(indiceY(nmes_curr+j)) &
                +splaa(ll-1)*im1(indiceY(nmes_curr+j))&
                + splaa(ll)*im(indiceY(nmes_curr+j))

               jacobien = jacobien + log(splaa(ll-2)*mm2(indiceY(nmes_curr+j)) &
                 +splaa(ll-1)*mm1(indiceY(nmes_curr+j))&
                 +splaa(ll)*mm(indiceY(nmes_curr+j)))

!                write(*,*)'Y',Y1(j),j,jacobien
            end do
         else if (idlink.eq.-1) then
            do j=1,nmes(i)
               Y1(j)=dble(Y(nmes_curr+j))
            end do            
         end if ! fin transfos




! creation de P=Zi*Ut et V=P*P' que si non spec aux classes

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
               vrais_comp_i=-1.d9
               !print*,"pb inv Vi"
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



          vrais=vrais-det
          b0=0.d0
          l=0
          m=0
          X00=0.d0
          do k=1,nv
             if (idg(k).ne.0) then
                l=l+1
                do j=1,nmes(i)
                   X00(j,l)=dble(X(nmes_curr+j,k))
                end do
                ! idg ne 0 pour l'intercept forcement donc on met le parm a 0
                if (k.eq.1.and.idlink.ne.-1) then
                   b0(l)=0.d0
                else
                   m = m+1
                   b0(l)=b1(nprob+nrisqtot+nvarxevt+m)
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

          vrais=vrais+2*(fevt-surv_glob)  
          entretard=entretard-surv0_glob


! cas 2 :  ng>1  composantes
       else
            
            
            if (prior(i).ne.0) then
               pi=0.d0
               pi(prior(i))=1.d0
            else
                pi=1.d0   
                if(nprob.gt.0) then
                   ! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
                   Xprob=0.d0
                   !Xprob(1)=1
                   l=0
                   do k=1,nv
                      if (idprob(k).eq.1) then
                         l=l+1
                         Xprob(l)=X(nmes_curr+1,k)
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
                
                do g=1,ng
                   pi(g) = pi(g)*pprior(i,g)
                end do
                
             end if


! creation des vecteurs de variables explicatives
            l=0
            m=0
            q=0
            X00=0.d0
            X2=0.d0
            do k=1,nv
               if (idg(k).eq.2) then
                  l=l+1
                  do j=1,nmes(i)
                        X2(j,l)=dble(X(nmes_curr+j,k))
                  end do
               else if (idg(k).eq.1) then
                  m=m+1
                  do j=1,nmes(i)
                     X00(j,m)=dble(X(nmes_curr+j,k))
                  end do
               end if
            end do

            b2=0.d0
            b0=0.d0


            expo=0.d0
            do g=1,ng
            
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



! variance covariance si spec aux classes :

               Ut1=Ut
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
                     !print*,"pb dsinv ier=-1"
                     vrais_comp_i=-1.d9
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
               Y4=DOT_PRODUCT(Y2,Y3)

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
                  end IF
                  nxevtcurr=nxevtcurr+nxevtspec(ke)
               end do

               if (Devt(i).ne.0) then
                  expo=expo+pi(g)*fevt*exp(&
                       +(-det-Y4)/2.d0-surv_glob)
               ELSE
                  expo=expo+pi(g)*exp((-det-Y4)/2.d0-surv_glob)
               end if

               retard = retard+pi(g)*exp(-surv0_glob)

            
            end do ! fin boucle classe
         
            entretard=entretard+log(retard)

 
            vrais=vrais+2*log(expo)


         end if

      vrais_comp_i=vrais/2.D0+jacobien


      if(idtrunc.eq.1) then 
         vrais_comp_i = vrais_comp_i - entretard
      end if


 654  continue


      return

      end function vrais_comp_i





      double precision function vrais_comp(b,m,id,thi,jd,thj)


        use commun_comp,only:ns,nmes,nmes_curr,vrais_surv

        implicit none

        integer::m,i,id,jd
        double precision::thi,thj,vrais_comp_i,temp
        double precision,dimension(m)::b

        vrais_surv=0.d0

        nmes_curr=0
        vrais_comp=0.d0
        do i=1,ns

           temp=vrais_comp_i(b,m,id,thi,jd,thj,i)

           vrais_comp = vrais_comp + temp
           if (temp.eq.-1.d9 .or. temp/temp.ne.1) then
              vrais_comp = -1.d9
              !print*,"temp",temp
              goto 541
           end iF
           nmes_curr = nmes_curr + nmes(i)
        end do

 !print*,'id=',id,'jd=',jd,'vrais',vrais_comp
        541 continue

        return

      end function vrais_comp




      subroutine fct_risq_comp_i(i,k,brisq,g,risq,surv,surv0,survint)

        use commun_comp
        
        implicit none

        integer::i,k,g
        double precision,dimension(nprisq(k))::brisq
        double precision,dimension(ng,nbevt)::risq,surv,surv0,survint

        integer::j,l,ll,kk,ii
        double precision::som

! CPL change 20022015
        if (typrisq(k).eq.2.and.logspecif.eq.1) then

            surv(g,k)=brisq(1)*(tsurv(i))**brisq(2)

            risq(g,k)=brisq(1)*brisq(2)*(tsurv(i))**(brisq(2)-1)
            if (idtrunc.eq.1) then
               surv0(g,k)=brisq(1)*(tsurv0(i))**brisq(2)
            end if
            if (ind_survint(i).eq.1) then
               survint(g,k)=brisq(1)*(tsurvint(i))**brisq(2)
            else
               survint(g,k)=surv(g,k)
            end if

         end if
        if (typrisq(k).eq.2.and.logspecif.eq.0) then

            surv(g,k)=(brisq(1)*tsurv(i))**brisq(2)

            risq(g,k)=brisq(1)*brisq(2)*(brisq(1)*tsurv(i))**(brisq(2)-1)
            if (idtrunc.eq.1) then
               surv0(g,k)=(brisq(1)*tsurv0(i))**brisq(2)
            end if
            if (ind_survint(i).eq.1) then
               survint(g,k)=(brisq(1)*tsurvint(i))**brisq(2)
            else
               survint(g,k)=surv(g,k)
            end if

         end if
! end CPL change 20022015

         if (typrisq(k).eq.1) then
            do j=1,nz(k)-1
               som=0.d0
               do l=1,j-1
                  som=som+brisq(l)*(zi(l+1,k)-zi(l,k))
               end do
               if (idtrunc.eq.1) then
                  if (Tsurv0(i).ge.zi(j,k).and.Tsurv0(i).le.zi(j+1,k)) then
                     surv0(g,k)=som+brisq(j)*(Tsurv0(i)-zi(j,k))
                  end if
               end if
               if (Tsurv(i).ge.zi(j,k).and.Tsurv(i).le.zi(j+1,k)) then
                  surv(g,k)=som+brisq(j)*(Tsurv(i)-zi(j,k))
                  risq(g,k)=brisq(j)
               end if
               if (ind_survint(i).eq.1) then
                  if (Tsurvint(i).ge.zi(j,k).and.Tsurvint(i).le.zi(j+1,k)) &
                       then
                     survint(g,k)=som+brisq(j)*(Tsurvint(i)-zi(j,k))
                  end if
               end if
            end do
            if (ind_survint(i).eq.0) then
               survint(g,k)=surv(g,k)
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

            surv(g,k)=som+brisq(ll)*Tim3(i)+brisq(ll+1)*Tim2(i) &
                 +brisq(ll+2)*Tim1(i)+brisq(ll+3)*Tim(i)
            risq(g,k)=brisq(ll)*Tmm3(i)+brisq(ll+1)*Tmm2(i)     &
                 +brisq(ll+2)*Tmm1(i)+brisq(ll+3)*Tmm(i)

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

               surv0(g,k)=som+brisq(ll)*Tim03(i)+brisq(ll+1)*Tim02(i) &
                    +brisq(ll+2)*Tim01(i)+brisq(ll+3)*Tim0(i)

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

               survint(g,k)=som+brisq(ll)*Timt3(i)+brisq(ll+1)*Timt2(i) &
                    +brisq(ll+2)*Timt1(i)+brisq(ll+3)*Timt(i)

            else
               survint(g,k)=surv(g,k)
            end if

         end if


      end subroutine fct_risq_comp_i




      subroutine splines(k)
        use commun_comp
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

           Tmm3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
           Tmm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                  +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))  &
                  +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                  +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))   &
                  +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmm(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

        end if

        if (Tsurv(i).eq.zi(n-2,k)) then

           Tmm3(i) = 0.d0
           Tmm2(i) = 0.d0
           Tmm1(i) = 0.d0
           Tmm(i) = 4.d0/h

        end if

        Tim3(i) = (0.25d0*(Tsurv(i)-zi(l-3,k))*Tmm3(i)) &
                         +(0.25d0*hh2*Tmm2(i))        &
            +(0.25d0*h3m*Tmm1(i))+(0.25d0*h4*Tmm(i))
        Tim2(i) = (0.25d0*hht*Tmm2(i))  &
            +(h3m*Tmm1(i)*0.25d0)+(h4*Tmm(i)*0.25d0)
        Tim1(i) = (htm*Tmm1(i)*0.25d0)+(h4*Tmm(i)*0.25d0)
        Tim(i) = ht*Tmm(i)*0.25d0

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

              Tmm03(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))

              Tmm02(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                   +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))   &
                   +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
              Tmm01(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                   +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))    &
                   +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
              Tmm0(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

           end if

           if (Tsurv0(i).eq.zi(n-2,k)) then

              Tmm03(i) = 0.d0
              Tmm02(i) = 0.d0
              Tmm01(i) = 0.d0
              Tmm0(i) = 4.d0/h

           end if

           Tim03(i) = (0.25d0*(Tsurv0(i)-zi(l-3,k))*Tmm03(i))  &
                +(0.25d0*hh2*Tmm02(i))           &
                +(0.25d0*h3m*Tmm01(i))+(0.25d0*h4*Tmm0(i))
           Tim02(i) = (0.25d0*hht*Tmm02(i))                  &
                +(h3m*Tmm01(i)*0.25d0)+(h4*Tmm0(i)*0.25d0)
           Tim01(i) = (htm*Tmm01(i)*0.25d0)+(h4*Tmm0(i)*0.25d0)
           Tim0(i) = ht*Tmm0(i)*0.25d0

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

              Tmmt3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
              Tmmt2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn)) &
                  +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))     &
                  +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
              Tmmt1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                  +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))       &
                  +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
              Tmmt(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

           end if

           if (Tsurvint(i).eq.zi(nz(k)-2,k)) then

              Tmmt3(i) = 0.d0
              Tmmt2(i) = 0.d0
              Tmmt1(i) = 0.d0
              Tmmt(i) = 4.d0/h

           end if

           Timt3(i) = (0.25d0*(Tsurvint(i)-zi(l-3,k))*Tmmt3(i)) &
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

      end subroutine splines




! =============================================
! subroutine de creation de design matrix
! =============================================




      subroutine design_splines_comp (ier)

      use commun_comp

      implicit none

      integer ::jj,l,k,ier
      double precision ::ht,htm,ht2,ht3,h,hh,h2,h3,h2n,hn,hht

      ier=0
      jj=0
      l=0
    
      do jj=1,nvalSPL      !     ou se trouve la valeur de zi
         do k = 2,ntrtot-2
            if ((uniqueY(jj).ge.zitr(k-1)).and.(uniqueY(jj).lt.zitr(k))) then
               l=k-1
            end if
         End do


         if (uniqueY(jj).eq.zitr(ntrtot-2)) then
            l=ntrtot-3
         end if

         ht2 = zitr(l+1)-uniqueY(jj)
         htm= uniqueY(jj)-zitr(l-1)
         ht = uniqueY(jj)-zitr(l)
         ht3 = zitr(l+2)-uniqueY(jj)
         hht = uniqueY(jj)-zitr(l-2)
         h = zitr(l+1)-zitr(l)
         hh= zitr(l+1)-zitr(l-1)
         hn= zitr(l+1)-zitr(l-2)
         h2n=zitr(l+2)-zitr(l-1)
         h2= zitr(l+2)-zitr(l)
         h3= zitr(l+3)-zitr(l)
            
         if (uniqueY(jj).ne.zitr(ntrtot-2)) then
            mm2(jj) = (3.d0*ht2*ht2)/(hh*h*hn)
            mm1(jj) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
            mm(jj)  = (3.d0*ht*ht)/(h3*h2*h)
         end if

         if (uniqueY(jj).eq.zitr(ntrtot-2)) then
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
      
      
765   continue

    end subroutine design_splines_comp





    subroutine postprob_comp(b,npm,ppi,ppiy)


      use commun_comp
!      use optim

      implicit none

      integer ::i,j,k,l,m,g,l2,m2,jj,it,npm,ier,nmoins,kk,j1,j2,ll,nevtxcurr,nxevtcurr
      double precision,dimension(maxmes,nv) ::Z,P,X0,X2
      double precision,dimension(nv) ::Xprob
      double precision,dimension(nv,nv) ::Ut,Ut1
      double precision,dimension(:,:),allocatable ::VC,Corr
      double precision,dimension(npm) :: b,b1
      double precision,dimension(:),allocatable ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: eps,det,eta0,ytemp
      double precision,dimension(-1:ntrtot-3)::splaa
      double precision ::temp,aa,bb,aa1,bb1,cc1,dd1,som,betai
      double precision ::Y4,fevt
      double precision,dimension(ng) ::fi,pi,fi1
      double precision,dimension(ns,ng) ::ppiy,ppi
      double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor

      double precision,dimension(nxevt)::bevt,Xevt
      double precision,dimension(nbevt)::bevtint
      double precision::varexpsurv,surv_glob
      integer::sumnrisq,ii,ke
      double precision,dimension(maxval(nprisq))::brisq
      double precision,dimension(ng,nbevt)::risq,surv,surv0,survint


      eps=1.D-20
      ppiy=0.D0
      ppi=0.D0
      b1=0.d0
      do k=1,npm
         b1(k)=b(k)
      end do

      allocate(VC(maxmes,maxmes),Corr(maxmes,maxmes),Vi(maxmes*(maxmes+1)/2))
      VC=0.d0
      Vi=0.d0
      Corr=0.d0



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
               ppiy=-1.d0
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

         
         if (prior(i).ne.0) then
            pi=0.d0
            pi(prior(i))=1.d0
         else
            pi=1.d0   
            if(nprob.gt.0) then
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
            end if
            
            do g=1,ng
               pi(g) = pi(g)*pprior(i,g)
            end do
            
         end if
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
         ! Xevt X0 et X2 ok



         ! creation des donnees transformees



         ! creation de Y1
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
                  ppiy=-1.d0
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
                  ppiy=-1.d0
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
         !sumyig=0
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
                  nevtxcurr=0
                  do k=1,nv                  
                     if (idtdv(k).ne.1) then
                        if (idcom(k).eq.1) then  
                           if (idspecif(k).eq.2) then 
                              l=l+1
                              bevt(l)=b1(nprob+nrisqtot+nevtxcurr+g)
                              Xevt(l)=X(it+1,k)
                           else 
                              l=l+1
                              bevt(l)=b1(nprob+nrisqtot+nevtxcurr+1)
                              Xevt(l)=X(it+1,k)
                           end if

                        else 

                           if (idspecif((ke-1)*nv+k).eq.1) then   
                              l=l+1
                              bevt(l)=b1(nprob+nrisqtot+nevtxcurr+nxcurr(k)+ 1)
                              Xevt(l)=X(it+1,k)
                              nxcurr(k)=nxcurr(k)+1
                           end if
                           if (idspecif((ke-1)*nv+k).eq.2) then  
                              l=l+1
                              bevt(l)=b1(nprob+nrisqtot+nevtxcurr+nxcurr(k)+g)
                              Xevt(l)=X(it+1,k)
                              nxcurr(k)=nxcurr(k)+ng
                           end if
                        end if

                     else

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
                  ppiy=-1.d0
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
            fi1(g)=dexp(fi(g))



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
            fi(g)=fi(g)+fevt-surv_glob

            fi(g)=exp(fi(g))

            !ppi(i,g) = pi(g)*fi(g)
            !ppiy(i,g) = pi(g)*fi1(g)

         end do ! fin boucle classe



         !sumyig=DOT_PRODUCT(pi,fi1)
         !sumig=DOT_PRODUCT(pi,fi)

         do g=1,ng
            !ppiy(i,g)=ppiy(i,g)/sumyig
            !ppi(i,g)=ppi(i,g)/sumig
            ppiy(i,g)=pi(g)*fi1(g)/DOT_PRODUCT(pi,fi1)
            ppi(i,g)=pi(g)*fi(g)/DOT_PRODUCT(pi,fi)
         end do

         it = it + nmes(i)

      end do ! fin boucle sujet


456   continue

      deallocate(VC,Vi,Corr)

      return

      end subroutine postprob_comp







      subroutine splines_estime_comp(time,nsim,ke)
        use commun_comp
        implicit none

        integer::i,k,n,l,ke
        double precision::ht,htm,h2t,ht2,ht3,hht,h,hh,h2,h3,h4,h3m,h2n, &
        hn,hh2,hh3
        integer::nsim
        double precision,dimension(nsim)::time


        Tmm_est=0.d0
        Tmm1_est=0.d0
        Tmm2_est=0.d0
        Tmm3_est=0.d0
        Tim_est=0.d0
        Tim1_est=0.d0
        Tim2_est=0.d0
        Tim3_est=0.d0

        l=0


        zi(-2,ke)=zi(1,ke)
        zi(-1,ke)=zi(1,ke)
        zi(0,ke)=zi(1,ke)
        zi(nz(ke)+1,ke)=zi(nz(ke),ke)
        zi(nz(ke)+2,ke)=zi(nz(ke),ke)
        zi(nz(ke)+3,ke)=zi(nz(ke),ke)
        n=nz(ke)+2
!------------------- time ---------------------------
      Do i=1,nsim
        do k=2,n-2
           if ((time(i).ge.zi(k-1,ke)).and.  &
               time(i).lt.zi(k,ke)) then
               l=k-1
           end if
        end do

        if (time(i).eq.zi(n-2,ke)) then
           l=n-3
        end if

        ht = time(i)-zi(l,ke)
        htm = time(i)-zi(l-1,ke)
        h2t = time(i)-zi(l+2,ke)
        ht2 = zi(l+1,ke)-time(i)
        ht3 = zi(l+3,ke)-time(i)
        hht = time(i)-zi(l-2,ke)
        h = zi(l+1,ke)-zi(l,ke)
        hh = zi(l+1,ke)-zi(l-1,ke)
        h2 = zi(l+2,ke)-zi(l,ke)
        h3 = zi(l+3,ke)-zi(l,ke)
        h4 = zi(l+4,ke)-zi(l,ke)
        h3m = zi(l+3,ke)-zi(l-1,ke)
        h2n = zi(l+2,ke)-zi(l-1,ke)
        hn = zi(l+1,ke)-zi(l-2,ke)
        hh3 = zi(l+1,ke)-zi(l-3,ke)
        hh2 = zi(l+2,ke)-zi(l-2,ke)

        if (time(i).ne.zi(n-2,ke)) then

           Tmm3_est(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
           Tmm2_est(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                  +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))  &
                  +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmm1_est(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                  +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))   &
                  +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmm_est(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

        end if

        if (time(i).eq.zi(n-2,ke)) then

           Tmm3_est(i) = 0.d0
           Tmm2_est(i) = 0.d0
           Tmm1_est(i) = 0.d0
           Tmm_est(i) = 4.d0/h

        end if

        Tim3_est(i) = (0.25d0*(time(i)-zi(l-3,ke))*Tmm3_est(i)) &
                         +(0.25d0*hh2*Tmm2_est(i))        &
            +(0.25d0*h3m*Tmm1_est(i))+(0.25d0*h4*Tmm_est(i))
        Tim2_est(i) = (0.25d0*hht*Tmm2_est(i))  &
            +(h3m*Tmm1_est(i)*0.25d0)+(h4*Tmm_est(i)*0.25d0)
        Tim1_est(i) = (htm*Tmm1_est(i)*0.25d0)+(h4*Tmm_est(i)*0.25d0)
        Tim_est(i) = ht*Tmm_est(i)*0.25d0

       End Do

      end subroutine splines_estime_comp







      subroutine fct_risq_estime_comp(k,brisq,time,nsim,g,risq,surv)


        use commun_comp


      implicit none
! risq contient le risque instantane d'evenement
! surv contient le risq cumule d'evenement et non la fonction de survie directe


      integer ::i,g,ll,kk,ii,l,j,nsim,k
      double precision, dimension(nprisq(k))::brisq
      double precision,dimension(nsim*ng,nbevt):: risq,surv
      double precision::som
      double precision,dimension(nsim)::time


      l=0
      do i=1,nsim

! change CPL 20022015

         if (typrisq(k).eq.2.and.logspecif.eq.1) then
            surv(nsim*(g-1)+i,k)=brisq(1)*(time(i))**brisq(2)
            risq(nsim*(g-1)+i,k)=brisq(1)*brisq(2) &
                 *(time(i))**(brisq(2)-1)
         end if
         if (typrisq(k).eq.2.and.logspecif.eq.0) then
            surv(nsim*(g-1)+i,k)=(brisq(1)*time(i))**brisq(2)
            risq(nsim*(g-1)+i,k)=brisq(1)*brisq(2) &
                 *(brisq(1)*time(i))**(brisq(2)-1)
         end if
! end change CPL 20022015
         if (typrisq(k).eq.1) then
            do j=1,nz(k)-1
               som=0.d0
               do l=1,j-1
                  som=som+brisq(l)*(zi(l+1,k)-zi(l,k))
               end do
               if (time(i).ge.zi(j,k).and.time(i).le.zi(j+1,k)) then
                  surv(nsim*(g-1)+i,k)=som+brisq(j)*(time(i)-zi(j,k))
                  risq(nsim*(g-1)+i,k)=brisq(j)
               end if
            end do
         end if
               
         if (typrisq(k).eq.3) then
            !------------ survie et risq pour Tsurv ----------------
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

            surv(nsim*(g-1)+i,k)=som+brisq(ll)*Tim3_est(i)+brisq(ll+1)*Tim2_est(i) &
                 +brisq(ll+2)*Tim1_est(i)+brisq(ll+3)*Tim_est(i)
            risq(nsim*(g-1)+i,k)=brisq(ll)*Tmm3_est(i)+brisq(ll+1)*Tmm2_est(i)     &
                 +brisq(ll+2)*Tmm1_est(i)+brisq(ll+3)*Tmm_est(i)

         end if
      end do



    end subroutine fct_risq_estime_comp





     subroutine transfo_estimee_comp(b,npm,nsim,marker,transfY)

      use commun_comp

      implicit none

      integer::kk,nsim,npm,j,k,l
      double precision,dimension(nsim)::marker,transfY
      double precision,dimension(ntrtot)::splaa,Xspl
      double precision::aa1,bb1,dd1,aa,bb,betai,eps,pas,ytemp,cc1
      double precision, dimension(npm)::b,b1
      double precision ::ht,htm,ht2,ht3,hht,h,hh,h2,h3,h2n,hn
      double precision,dimension(nsim)::mmm,mmm1,mmm2,iim,iim1,iim2



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


        splaa(1)=b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)
        do kk=2,ntrtot
           splaa(kk)=b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+kk)**2
             end do

!           write(*,*)'ntr',(b1(nprob+nrisq
!     &       +nvarxevt+ncssg+ncg*ng+kk),kk=1,ntrtot)


      !calcul de H(y) pour remplacer call estim_splines..
        do j=1,nsim
! ou se trouve la valeur
         l=0

         do k = 2,ntrtot-2
               if ((marker(j).ge.zitr(k-1)).and.(marker(j).lt.zitr(k))) then
                  l=k-1
               end if
        end do

           if (marker(j).eq.zitr(ntrtot-2)) then
               l=ntrtot-3
            end if

!         if (l.lt.1.or.l.gt.ntrtot-1) then
!            write(*,*)'probleme estim splines',l
!            write(*,*)'j=',j,'test(j)',test(j)
!            stop
!         end if


               ht2 = zitr(l+1)-marker(j)
               htm= marker(j)-zitr(l-1)
               ht = marker(j)-zitr(l)
               ht3 = zitr(l+2)-marker(j)
               hht = marker(j)-zitr(l-2)
               h = zitr(l+1)-zitr(l)
               hh= zitr(l+1)-zitr(l-1)
               hn= zitr(l+1)-zitr(l-2)
               h2n=zitr(l+2)-zitr(l-1)
               h2= zitr(l+2)-zitr(l)
               h3= zitr(l+3)-zitr(l)

               if (marker(j).ne.zitr(ntrtot-2)) then
                  mmm2(j) = (3.d0*ht2*ht2)/(hh*h*hn)
                  mmm1(j) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
                  mmm(j)  = (3.d0*ht*ht)/(h3*h2*h)
               end if
               if (marker(j).eq.zitr(ntrtot-2)) then
                  mmm2(j) = 0.d0
                  mmm1(j) = 0.d0
                  mmm(j)  = 3.d0/h
               end if

               iim2(j)=hht*mmm2(j)/(3.d0)+ h2n*mmm1(j)/(3.d0) &
      +h3*mmm(j)/(3.d0)
               iim1(j)=htm*mmm1(j)/(3.d0)+h3*mmm(j)/(3.d0)

               iim(j)=ht*mmm(j)/(3.d0)

!-------- transformation :

            Xspl=0.d0
            Xspl(1)=1
            do k=2,l
               Xspl(k)=1
            end do
            Xspl(l+1)=iim2(j)
            Xspl(l+2)=iim1(j)
            Xspl(l+3)=iim(j)
            transfY(j)= dot_product(Xspl,splaa)
      end do
        !fin H(y)

        else if (idlink.eq.1) then

            aa1=exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1))/ &
             (1+exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)))
            bb1=exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2))/ &
             (1+exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2)))
            bb1=aa1*(1.d0-aa1)*bb1
            cc1=b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+3)
            dd1=abs(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+4))

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
                    transfY(j)=(marker(j)-b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1))/ &
                         abs(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2))
                 end do
        end if
        end subroutine transfo_estimee_comp








     subroutine residuals_comp(b1,npm,ppi,resid_m,pred_m_g,resid_ss, &
      pred_ss_g,pred_RE,Yobs)

        use commun_comp
!      use optim

      implicit none

      integer ::i,j,k,l,m,g,l2,m2,jj,npm,j1,j2,ii,ll
      integer ::ier,nmoins,it,kk
      double precision,dimension(maxmes,nv) ::Z,P,X0,X2
      double precision,dimension(nv) ::Xprob,err2
      double precision,dimension(nv,nv) ::Ut,Ut1
      double precision,dimension(:,:),allocatable ::VC,Corr,VC1
      double precision,dimension(maxmes,maxmes)::SigmaE,CovDev
      double precision,dimension(npm) ::b1
      double precision,dimension(:),allocatable ::Vi
      double precision,dimension(nv) :: b0,b2,bprob
      double precision :: eps,det
      double precision ::temp
      double precision,dimension(nv,maxmes)::Valea
      double precision,dimension(maxmes) :: mu,Y1,Y2,pred1,err1,tcor
      double precision,dimension(ng) :: pi
      double precision,dimension(nobs)::resid_m,resid_ss
      double precision,dimension(ns*nea)::pred_RE
      double precision,dimension(nobs*ng)::pred_m_g,pred_ss_g
      double precision,dimension(ns,ng) ::PPI
      double precision,dimension(nobs)::Yobs
      double precision,dimension(-1:(ntrtot-3))::splaa
      double precision::aa1,bb1,dd1,aa,bb,betai,ytemp,som,cc1



      allocate(VC(maxmes,maxmes),VC1(maxmes,maxmes),Corr(maxmes,maxmes),Vi(maxmes*(maxmes+1)/2))
      VC=0.d0
      VC1=0.d0
      Vi=0.d0
      Corr=0.d0
      SigmaE=0.d0
      CovDev=0.d0



      eps=1.D-20

!----------- rappel des parametres utilises ---------




!      write(*,*)'nvc',nvc,'nea',nea,'nwg',nwg,'nef',nef
! creation de Ut, decomposition de cholesky pour G

      Ut=0.d0
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

!         write(*,*)'Ut',Ut

! ----------- boucle sur les individus -------------

      resid_m=0.d0
      pred_m_g=0.d0
      resid_ss=0.d0
      pred_ss_g=0.d0
      pred_RE=0.d0
      Yobs=0.d0
      pred_RE=0.d0
      it=0
      kk=0
      do i=1,ns

! -------- creation de Vi = ZiGZi'+Ri+se*seIni ----------

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
        SigmaE=0.d0
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
                  if(idlink.ne.-1) SigmaE(j1,j2) = 1
                  if(idlink.eq.-1) SigmaE(j1,j2) = b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)**2
               end if
               if (ncor.eq.1) then 
                  Corr(j1,j2) = Corr(j1,j2)+ & 
                       b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor)**2 &
                          * min(tcor(j1),tcor(j2))
               else if (ncor.eq.2) then
                  Corr(j1,j2) = Corr(j1,j2)+ & 
                       b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor)**2 &
                          *exp(-b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+1)* &
                          abs(tcor(j1)-tcor(j2)))
               end if
            end do
         end do         
                 

!     creation de Y1

        Y1=0.d0
         if (idlink.eq.0) then  ! Linear link

            do j=1,nmes(i)
               !nmestot=nmestot+1
               Y1(j)=dble(Y(it+j)-b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)) &
                     /abs(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2))
               Yobs(it+j)=Y1(j)
            end do

         elseif (idlink.eq.1) then  ! Beta link


            aa1=exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1))/ &
             (1+exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)))
            bb1=exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2))/ &
             (1+exp(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2)))
            bb1=aa1*(1.d0-aa1)*bb1
            cc1=b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+3)
            dd1=abs(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+4))

            aa=aa1*aa1*(1-aa1)/bb1-aa1
            bb=aa*(1-aa1)/aa1


            do j=1,nmes(i)

               !nmestot=nmestot+1

               ytemp=(dble(Y(it+j))-minY+epsY)/(maxY-minY+2*epsY)
               Y1(j)=(betai(aa,bb,ytemp)-cc1)/dd1

               if (Y1(j).eq.999.d0) then
                  do k=1,nmes(i)
                     resid_m(it+k)=9999.d0
                     pred_m_g(it+k)=9999.d0
                     Yobs(it+k)=9999.d0
                     resid_ss(it+k)=9999.d0
                     pred_ss_g(it+k)=9999.d0
                  end do
                  do k=1,nea
                     pred_RE((i-1)*nea+k)=9999.d0
                  end do
                  goto 654
               else
                  Yobs(it+j)=Y1(j)
               end if

            end do

         elseif (idlink.eq.2) then ! Splines link


            bb=b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)

            do kk=2,ntrtot
               splaa(kk-3)=b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+kk)**2
            end do

            do j=1,nmes(i)

               !nmestot=nmestot+1

               ll=0
               if (Y(it+j).eq.zitr(ntrtot-2)) then
                  ll=ntrtot-3
               end if
               som=0.d0
               do kk = 2,ntrtot-2
                  if ((Y(it+j).ge.zitr(kk-1)).and.(Y(it+j).lt.zitr(kk))) then
                     ll=kk-1
                  end if
               end do
               if (ll.lt.1.or.ll.gt.ntrtot-3) then
               do k=1,nmes(i)
                  resid_m(it+k)=9999.d0
                  pred_m_g(it+k)=9999.d0
                  Yobs(it+k)=9999.d0
                  resid_ss(it+k)=9999.d0
                  pred_ss_g(it+k)=9999.d0
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

               Y1(j)=bb+som +splaa(ll-2)*im2(indiceY(it+j)) &
                 +splaa(ll-1)*im1(indiceY(it+j))+splaa(ll)*im(indiceY(it+j))

               Yobs(it+j)=Y1(j)


            end do

         else if (idlink.eq.-1) then
            do j=1,nmes(i)
               Y1(j)=dble(Y(it+j))
               Yobs(it+j)=Y1(j)
            end do

         end if






!     creation de P=Zi*Ut et V=P*P' que si non spec aux classes
!     debut du calcul de la vraisemblance
!     cas 1 : ng=1

         if (ng.eq.1) then


            Valea=0.d0
            VC=0.d0
            P=0.d0
            P=MATMUL(Z,Ut)
            Valea=MATMUL(Ut,transpose(P))
            VC=MATMUL(P,transpose(P))+Corr+SigmaE
            CovDev=MATMUL(P,transpose(P))+Corr

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
                  resid_m(it+j)=9999.d0
                  pred_m_g(it+j)=9999.d0

                  resid_ss(it+j)=9999.d0
                  pred_ss_g(it+j)=9999.d0
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
            m=0
            X0=0.d0
            do k=1,nv
               if (idg(k).ne.0) then
                  l=l+1
                  do j=1,nmes(i)
                     X0(j,l)=dble(X(it+j,k))
                  end do
                  if(k.eq.1.and.idlink.ne.-1) then
                     b0(l)=0.d0
                  else
                     m=m+1
                     b0(l)=b1(nprob+nrisqtot+nvarxevt+m)
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
               resid_m(it+j)=Y2(j)
               pred_m_g(it+j)=mu(j)

               resid_ss(it+j)=Y1(j)-pred1(j)
               pred_ss_g(it+j)=pred1(j)
            end do


            do k=1,nea
               pred_RE((i-1)*nea+k)=err2(k)
            end do



!     cas 2 :  ng>1  composantes
         else

!     transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
            if (prior(i).ne.0) then
                   pi=0.d0
               pi(prior(i))=1.d0
            else
               pi=1.d0   
               if(nprob.gt.0) then
                  Xprob=0.d0
                  !Xprob(1)=1.d0
                  l=0
                  do k=1,nv
                     if (idprob(k).eq.1) then
                        l=l+1
                        Xprob(l)=X(it+1,k)
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
               
               do g=1,ng
                  pi(g) = pi(g)*pprior(i,g)
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


            do g=1,ng
               nmoins=0
               l2=0
               m2=0
               b0=0.d0
               b2=0.d0
               do k=1,nv
                  if (idg(k).eq.1) then
                     if (k.eq.1.and.idlink.ne.-1) then
                        m2=m2+1
                        b0(m2)=0.d0
                     else
                        m2=m2+1
                        b0(m2)=b1(nprob+nrisqtot+nvarxevt+nmoins+1)
                        nmoins=nmoins+1
                     end if
                  else if (idg(k).eq.2) then
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

               mu=0.d0
               mu=matmul(X0,b0)+matmul(X2,b2)
               VC=0.d0
               P=0.d0
               Ut1=Ut
               if (nwg.ne.0) then
                  if (g.eq.ng) then
                     Ut1=Ut
                  else
                     Ut1=Ut*b1(nprob+nrisqtot+nvarxevt+nef+nvc+g)
                  end if
               end if
               P=0.d0
               Valea=0.d0
               P=MATMUL(Z,Ut1)
               Valea=MATMUL(Ut1,transpose(P))
               VC=MATMUL(P,transpose(P))+Corr+SigmaE
               CovDEv=MATMUL(P,transpose(P))+Corr


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
                     resid_m(it+j)=9999.d0
                     resid_ss(it+j)=9999.d0
                     do l=1,ng
                        pred_ss_g((l-1)*nobs+it+j)=9999.d0
                        pred_m_g((l-1)*nobs+it+j)=9999.d0
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
                  pred_m_g((g-1)*nobs+it+j)=mu(j)
                  pred_ss_g((g-1)*nobs+it+j)=pred1(j)

                  resid_ss(it+j)=resid_ss(it+j)   &
                      +ppi(i,g)*(Y1(j)-pred1(j))
                  resid_m(it+j)=resid_m(it+j)+pi(g)*(Y2(j))

               end do
               do k=1,nea
                  pred_RE((i-1)*nea+k)=pred_RE((i-1)*nea+k)+ppi(i,g)*err2(k)
               end do

            end do

         end if


         it=it+nmes(i)
      end do

 654  continue

      deallocate(VC,VC1,Vi,Corr)
      return

      end subroutine residuals_comp










!------------------------------------------------------------
!                   test independance Y et T
!------------------------------------------------------------


      subroutine scoretest_comp(b,npm,statglob,statevt)

        use commun_comp
!        use optim


        IMPLICIT NONE
        integer ::i,j,k,l,m,g,l2,m2,jj,npm,ll,ii !,id,jd
        integer ::ier,nmoins,kk,j1,j2,q,iqq
        integer::ke,sumnrisq,nxevtcurr,nevtxcurr
        double precision::statglob
        double precision,dimension(nbevt)::statevt
        double precision,dimension(ns,nea*nbevt)::stat
        double precision,dimension(nea*nbevt)::statscore
        double precision,dimension(nea*nbevt,nea*nbevt)::varscore
        double precision,dimension(nea*nbevt*(nea*nbevt+1)/2)::vars
        double precision,dimension(nea)::RE
        double precision,dimension(nea,maxmes)::BZT
        double precision,dimension(maxmes,nv) ::X00,X2
        double precision,dimension(maxmes,nea) ::Z,P
        double precision,dimension(nvarprob) ::Xprob,bprob 
        double precision,dimension(nea,nea) ::Ut,Ut1
        double precision,dimension(maxmes,maxmes) ::VC,Corr
        double precision,dimension(npm) :: b,b1
        double precision,dimension(maxmes*(maxmes+1)/2) ::Vi
        double precision,dimension(nv) :: b0,b2
        double precision,dimension(nxevt)::Xevt,bevt
        double precision,dimension(nbevt)::bevtint,risqcum
        double precision,dimension(maxval(nrisq))::brisq

        double precision :: eps,det,som,temp,eta0,expo,surv_glob,fevt,varexpsurv
        double precision ::Y4,Li,ytemp !,beta_densite
        !double precision,dimension(ng)::expo
        double precision,dimension(maxmes) :: mu,Y1,Y2,Y3,tcor
        double precision,dimension(ng) :: pi
        double precision,dimension(-1:ntrtot-3)::splaa
        double precision::aa1,bb1,dd1,aa,bb,betai,cc1
        double precision,dimension(ng,nbevt)::risq,surv,surv0,survint

        nmes_curr=0
        eps=1.d-7


        b1=0.d0
        do k=1,npm
           b1(k)=b(k)
        end do


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

        stat=0.d0
        statscore=0.d0

        do i=1,ns


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
                    surv0(g,ke)=surv0(g,ke)*exp(b1(nprob+sumnrisq+nprisq(ke)+g))
                 end if
              end do
              sumnrisq = sumnrisq + nrisq(ke)
           end do





           ! -------- creation de Vi = ZiGZi'+se*seIni ----------
           ! creation de Zi

           Z=0.d0
           l=0
           do k=1,nv
              if (idea(k).eq.1) then
                 l=l+1
                 do j=1,nmes(i)
                    Z(j,l)=dble(X(nmes_curr+j,k))
                 end do
              end if
           end do

           !matrice Corr=Ri+s2*I

           Corr=0.d0
           tcor=0.d0
           if (ncor.gt.0) then
              do k=1,nv
                 if (idcor(k).eq.1) then
                    do j=1,nmes(i)
                       tcor(j) = X(nmes_curr+j,k)
                    end do
                 end if
              end do
              do j1=1,nmes(i)
                 do j2=1,nmes(i)
                    if (ncor.eq.1) then
                       Corr(j1,j2) = Corr(j1,j2)+b1(nprob+nrisqtot+nvarxevt &
                            +nef+nvc+nwg+ncor)**2*min(tcor(j1),tcor(j2))
                    else if (ncor.eq.2) then
                       Corr(j1,j2) = Corr(j1,j2)+b1(nprob+nrisqtot+nvarxevt &
                            +nef+nvc+nwg+ncor)**2*exp(-b1(nprob+nrisqtot+nvarxevt &
                            +nef+nvc+nwg+1)*abs(tcor(j1)-tcor(j2)))
                    end if
                    if(j1.eq.j2) then
                       if(idlink.eq.-1) then 
                          Corr(j1,j2) = Corr(j1,j2)+b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)**2
                       else
                          Corr(j1,j2) = Corr(j1,j2) + 1 !variance de l'erreur
                       end if
                    end if
                 end do
              end do
           end if

         if(ncor.eq.0) then
            do j1=1,nmes(i)
               if(idlink.ne.-1) Corr(j1,j1) = 1 ! variance de l'erreur
               if(idlink.eq.-1) Corr(j1,j1) = b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+1)**2
            end do
         end if



           ! creation de Y1
           Y1=0.d0
           splaa=0.d0

           if (idlink.eq.0) then  ! Linear link

              do j=1,nmes(i)
                 Y1(j)=(dble(Y(nmes_curr+j))-b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg &
                      +ncor+1))/abs(b1(nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+2))
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

                 ytemp=(dble(Y(nmes_curr+j))-minY+epsY)/(maxY-minY+2*epsY)
                 Y1(j)=(betai(aa,bb,ytemp)-cc1)/dd1


                 if (Y1(j).eq.999.d0) then
                    statglob=9999.d0
                    !print*,"pb betai"
                    goto 654
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

                 if (Y(nmes_curr+j).eq.zitr(ntrtot-2)) then
                    ll=ntrtot-3
                 end if

                 som=0.d0
                 do kk = 2,ntrtot-2
                    if ((Y(nmes_curr+j).ge.zitr(kk-1)).and. &
                         (Y(nmes_curr+j).lt.zitr(kk))) then
                       ll=kk-1
                    end if
                 end do

                 if (ll.lt.1.or.ll.gt.ntrtot-3) then
                    statglob=9999.d0
                    !print*,"pb ll"
                    goto 654
                 end if
                 if (ll.gt.1) then
                    do ii=2,ll
                       som=som+splaa(ii-3)
                    end do
                 end if



                 Y1(j)=eta0+som +splaa(ll-2)*im2(indiceY(nmes_curr+j)) &
                      +splaa(ll-1)*im1(indiceY(nmes_curr+j))&
                      + splaa(ll)*im(indiceY(nmes_curr+j))


                 !                write(*,*)'Y',Y1(j),j,jacobien
              end do
           else if (idlink.eq.-1) then
              do j=1,nmes(i)
                 Y1(j)=dble(Y(nmes_curr+j))
              end do
           end if ! fin transfos



           ! creation de P=Zi*Ut et V=P*P' que si non spec aux classes

           if (nwg.eq.0.OR.NG.EQ.1) then
              P=0.d0
              P=MATMUL(Z,Ut)
              VC=0.d0
              VC=MATMUL(P,transpose(P))+Corr

              BZT=0.d0
              BZT=matmul(Ut,transpose(P))


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
                 !print*,"pb inv vi"
                 statglob=9999.d0
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



           ! contribution individuelle a la vraisemblance




           if (prior(i).ne.0) then
              pi=0.d0
              pi(prior(i))=1.d0
           else
              pi=1.d0
              if(ng.gt.1) then
                 pi=1.d0   
                 if(nprob.gt.0) then
                    ! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
                    Xprob=0.d0
                    !Xprob(1)=1
                    l=0
                    do k=1,nv
                       if (idprob(k).eq.1) then
                          l=l+1
                          Xprob(l)=X(nmes_curr+1,k)
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
                       !write(*,*)"pi(g)",pi(g)
                    end do
                    
                    pi(ng)=1/(1+temp)
                    
                    do g=1,ng-1
                       pi(g)=pi(g)*pi(ng)
                    end do
                 end if
                 
                 do g=1,ng
                    pi(g) = pi(g)*pprior(i,g)
                 end do
                
              end if
              
           end if

           ! creation des vecteurs de variables explicatives
           l=0
           m=0
           q=0
           X00=0.d0
           X2=0.d0
           do k=1,nv
              if (idg(k).eq.2) then
                 l=l+1
                 do j=1,nmes(i)
                    X2(j,l)=dble(X(nmes_curr+j,k))
                 end do
              else if (idg(k).eq.1) then
                 m=m+1
                 do j=1,nmes(i)
                    X00(j,m)=dble(X(nmes_curr+j,k))
                 end do
              end if
           end do

           b2=0.d0
           b0=0.d0
           Li=0.d0
           do g=1,ng

              if(nxevt.ne.0) then
                 bevt=0.d0
                 bevtint=0.d0
                 Xevt=0.d0
                 l=0
                 nxcurr=0 ! for each covariate, number of previous effects before ke
                  do ke=1,nbevt      
                     nevtxcurr=0 ! number of previous effects in b1 before k
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




              ! variance covariance si spec aux classes :

              Ut1=Ut
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


                 BZT=0.d0
                 BZT=matmul(Ut1,transpose(P))

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
                    !print*,"pb inv vi"
                    statglob=9999.d0
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
              Y4=DOT_PRODUCT(Y2,Y3)

              RE=0.d0
              RE=matmul(BZT,y3)

              expo=0.d0
              risqcum=0.d0
              surv_glob=0.d0
              nxevtcurr=0
              fevt=0.d0

              do ke=1,nbevt

                 varexpsurv=0.d0
                 if (nxevt.ne.0) then 
                    varexpsurv=DOT_PRODUCT(Xevt((nxevtcurr+1):(nxevtcurr+nxevtspec(ke)))&
                         ,bevt((nxevtcurr+1):(nxevtcurr+nxevtspec(ke))))
                 end if

                 risqcum(ke) = -(survint(g,ke)+exp(bevtint(ke))*(surv(g,ke)- &
                      survint(g,ke)))*exp(varexpsurv)

                 Surv_glob=surv_glob+ survint(g,ke)*exp(varexpsurv)+ &
                      exp(bevtint(ke)+varexpsurv)*(surv(g,ke)-survint(g,ke))

                 if (Devt(i).eq.ke) then
                    fevt=risq(g,ke)*exp(varexpsurv)                     
                    if (ind_survint(i).eq.1) then
                       fevt=fevt*exp(bevtint(ke))
                    end if
                    risqcum(ke)=1+risqcum(ke)
                 end if
                 nxevtcurr=nxevtcurr+nxevtspec(ke)

              end do

               if (Devt(i).ne.0) then
                  expo=expo+pi(g)*fevt*exp(&
                       +(-det-Y4)/2.d0-surv_glob)
               else
                  expo=expo+pi(g)*exp((-det-Y4)/2.d0-surv_glob)
               end if
               
               Li = Li + expo
               
               do ke=1,nbevt

                  do j=1,nea
                     stat(i,(ke-1)*nea+j)=stat(i,(ke-1)*nea+j)+expo*risqcum(ke)*RE(j)
                  end do
                  
               end do

            end do ! fin boucle classe


           !stats de test
           do ke=1,nbevt

              do j=1,nea
                 stat(i,(ke-1)*nea+j)=stat(i,(ke-1)*nea+j)/Li
                 statscore((ke-1)*nea+j)=statscore((ke-1)*nea+j)+stat(i,(ke-1)*nea+j)
              end do

           end do
           nmes_curr=nmes_curr+nmes(i)

        end do ! fin boucle i


        !calcul var(u)=sum(u_i*u_i')-u*u'/N
        varscore=0.d0
        do j=1,nea*nbevt
           do k=1,nea*nbevt
              do i=1,ns
                 varscore(j,k)=varscore(j,k)+stat(i,j)*stat(i,k)
              end do
              varscore(j,k)=varscore(j,k)-statscore(j)*statscore(k)/dble(ns)
           end do
        end do

        !inverser var(u)
        vars=0.d0
        jj=0
        do k=1,nea*nbevt
           do j=1,nea*nbevt
              jj=j+k*(k-1)/2
              vars(jj)=varscore(j,k)
           end do
        end do

        ier=0
        det=0.d0

        call dsinv(vars,nea*nbevt,eps,ier,det)
        if(ier.eq.-1) then
           !print*,"pb inv vars",eps,ier,det
           statglob=9999.d0
           goto 654
        end if


        !stat du test global u'*var(u)**-1*u
        statglob=0.d0
        do j=1,nea*nbevt
           do k=1,nea*nbevt
              if(k.ge.j) then
                 iqq=k*(k-1)/2+j
              else
                 iqq=j*(j-1)/2+k
              end if
              statglob=statglob+statscore(j)*vars(iqq)*statscore(k)
           end do
        end do

   
        !test pour chaque evt
        statevt=0.d0
        do ke=1,nbevt
           do j=1,nea
              do k=1,nea
                 if(k.ge.j) then
                    iqq=((ke-1)*nea+k)*((ke-1)*nea+k-1)/2+((ke-1)*nea+j)
                 else
                    iqq=((ke-1)*nea+j)*((ke-1)*nea+j-1)/2+((ke-1)*nea+k)
                 end if

                 statevt(ke)=statevt(ke)+statscore((ke-1)*nea+j)*vars(iqq)*statscore((ke-1)*nea+k)

              end do
           end do
        end do



654     continue


      end subroutine scoretest_comp




      
  subroutine loglikjointlcmm(Y0,X0,Prior0,pprior0,Tentr0,Tevt0,Devt0,ind_survint0 &
       ,idprob0,idea0,idg0,idcor0,idcom0,idspecif0,idtdv0,idlink0 &
       ,epsY0,nbzitr0,zitr0,uniqueY0,nvalSPL0,indiceY0 &
       ,typrisq0,risqcom0,nz0,zi0 &
       ,ns0,ng0,nv0,nobs0,nmes0,nbevt0,nea0,nwg0,ncor0,idiag0,idtrunc0,logspecif0 &
       ,npm0,b0,ppi0,ppitest0,resid_m &
       ,resid_ss,pred_m_g,pred_ss_g,pred_RE &
       ,time,risq_est,risqcum_est,marker,transfY,nsim,Yobs,statglob,statevt &
       ,fix0,nfix0,bfix0,estim0,loglik)


      use commun_comp

      IMPLICIT NONE

      !Declaration des variables en entree
      integer,intent(in)::ns0,ng0,nv0,nobs0,nea0,npm0,nsim,nfix0
      integer,intent(in)::nwg0,ncor0,idiag0,idtrunc0,logspecif0
      double precision,dimension(nobs0),intent(in)::Y0
      integer,dimension(nobs0),intent(in)::indiceY0
      double precision,dimension(nobs0*nv0),intent(in)::X0
      double precision, dimension(ns0),intent(in)::Tentr0,Tevt0
      integer, dimension(ns0),intent(in)::ind_survint0,nmes0,Devt0,prior0
      double precision, dimension(ns0*ng0), intent(in) :: pprior0
      integer,dimension(nv0),intent(in)::idprob0,idea0,idg0,idcor0 
      integer,intent(in)::idlink0,nbzitr0,nvalSPL0,nbevt0
      integer, dimension(nv0),intent(in)::idcom0,idtdv0
      integer, dimension(nv0*nbevt0),intent(in)::idspecif0
      double precision,dimension(nbzitr0),intent(in)::zitr0
      double precision,dimension(nvalSPL0),intent(in)::uniqueY0 
      integer,dimension(nbevt0),intent(in)::typrisq0,risqcom0,nz0
      double precision,dimension(maxval(nz0),nbevt0),intent(in)::zi0    
      double precision,intent(in)::epsY0
      integer,dimension(npm0+nfix0),intent(in)::fix0
      double precision,dimension(nfix0),intent(in)::bfix0
      integer,intent(in)::estim0 

      !Declaration des variable en entree et sortie
      double precision,dimension(npm0), intent(inout) :: b0

      !Declaration des variables en sortie
      double precision,intent(out)::loglik
      double precision,dimension(ns0*ng0),intent(out)::ppi0,ppitest0
      double precision,dimension(nobs0),intent(out)::resid_m,resid_ss,Yobs
      double precision,dimension(nobs0*ng0),intent(out)::pred_m_g
      double precision,dimension(nobs0*ng0),intent(out)::pred_ss_g
      double precision,dimension(ns0*nea0),intent(out)::pred_RE
      double precision,dimension(nsim),intent(out)::marker,transfY,time
      double precision,dimension(nsim*ng0,nbevt0),intent(out)::risq_est,risqcum_est
      double precision,intent(out)::statglob
      double precision,dimension(nbevt0),intent(out)::statevt

      !Variables locales
      integer::jtemp,i,g,j,ier,k,ktemp,ig,it,ke,sumnrisq,dejaspl
      integer::nbfix,npmtot0,k2,id,jd
      double precision::thi,thj
      double precision,dimension(ns0,ng0)::PPI,ppiy
      double precision,dimension(npm0+nfix0)::btot
      double precision,external::vrais_comp


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
      statglob=0.d0
      statevt=0.d0

      ! nmbre max de mesures pour 1 sujet
      maxmes=0
      do i=1,ns0
         if (nmes0(i).gt.maxmes) then
            maxmes=nmes0(i)
         end if
      end do


      ! alloc pour partie modele mixte
      minY=zitr0(1)
      maxY=zitr0(nbzitr0)

      rangeY=0
      !if (Ydiscret.eq.1) rangeY=maxY-minY


      epsY=epsY0
      idlink=idlink0
      nvalSPL=nvalSPL0
      if (idlink.eq.-1) ntrtot=1 !sans trasnfo
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


      allocate(Y(nobs0),idprob(nv0),X(nobs0,nv0) &
      ,idea(nv0),idg(nv0),idcor(nv0),nmes(ns0),prior(ns0) &
      ,idcom(nv0),idspecif(nv0*nbevt0),idtdv(nv0),pprior(ns0,ng0))

      !  alloc pour partie survie

      allocate(Tsurv0(ns0),Tsurv(ns0),Tsurvint(ns0),ind_survint(ns0),Devt(ns0))
      allocate(risqcom(nbevt0),typrisq(nbevt0),nz(nbevt0),nprisq(nbevt0),nrisq(nbevt0), &
           nxevtspec(nbevt0),nevtparx(nv0),nxcurr(nv0))





      nbevt=nbevt0    
      typrisq=typrisq0
      risqcom=risqcom0
      idtrunc=idtrunc0
      Tsurv0=Tentr0   
      Tsurv=Tevt0    
      devt=devt0    
      ind_survint=ind_survint0
      logspecif=logspecif0 
! Tsurvint defini plus loin
      Tsurvint=Tsurv

      ! zi : contient noeuds pour hazard (ou min et max si Weibull)
      if(any(typrisq.eq.3)) then
         allocate(zi(-2:maxval(nz0)+3,nbevt))
      else
         allocate(zi(maxval(nz0),nbevt))
      end if



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
      prior=0
      pprior=0.d0
      nmes=nmes0
      Y=0.d0
      X=0.d0
      idprob=0 !classmb
      idea=0   !random
      idg=0    !fixed
      idcor=0  !cor
      idcom=0 !survival
      idtdv=0 !timedepvar
      idspecif=0  !survival
      idspecif(1:nv*nbevt)=idspecif0(1:nv*nbevt)

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
               nmes(i)=nmes0(i)
               prior(i)=prior0(i)
               do j=1,nmes(i)
                  jtemp=jtemp+1
                  Y(jtemp)=Y0(jtemp)
               end do
            end if

            do j=1,nmes(i)
               ktemp=ktemp+1
               it=it+1
               X(it,k)=X0(ktemp)
            end do
            
            do g=1,ng                 
               pprior(i,g)=pprior0((i-1)*ng+g)
            end do
         end do
      end do

! definition Tsurvint 
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


      ! nvarxevt = nombre total de coef pour survie (sans prm hazard)
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


      if((ng.eq.1.and.ncg.gt.0).or.(ng.eq.1.and.nprob.gt.0)) then
         go to 1236
      end if


      if (idiag.eq.1) then
         nvc=nea
      else if(idiag.eq.0) then
         nvc=nea*(nea+1)/2
      end if


      nef=ncssg+ncg*ng
      if(idlink.ne.-1) nef=nef-1
      npmtot=nprob+nrisqtot+nvarxevt+nef+nvc+nwg+ncor+ntrtot

      if (npmtot.ne.npmtot0) then
         goto 1236
      end if


      if (nwg.gt.0) then
         do i=1,nwg
            btot(nprob+nrisqtot+nvarxevt+nef+nvc+i)=abs(btot(nprob+nrisqtot+nvarxevt+nef+nvc+i))
         end do
      end if

      ! creer base de splines si transfo splines
      if (idlink.eq.2) then 
         call design_splines_comp(ier)
         if (ier.eq.-1) then          
            go to 1236
         end if
      end if

      ! creer base de splines si au moins un hazard splines
      dejaspl=0
      if(any(typrisq.eq.3)) then
            allocate(Tmm(ns),Tmm1(ns),Tmm2(ns),Tmm3(ns),Tim(ns)         &
                 ,Tim1(ns),Tim2(ns),Tim3(ns),Tmm0(ns),Tmm01(ns),Tmm02(ns)  &
                 ,Tmm03(ns),Tim0(ns),Tim01(ns),Tim02(ns),Tim03(ns),          &
                 Tmmt(ns),Tmmt1(ns),Tmmt2(ns),Tmmt3(ns),Timt(ns),Timt1(ns) &
                 ,Timt2(ns),Timt3(ns))

            do ke=1,nbevt
               if(typrisq(ke).eq.3 .and. dejaspl.eq.0) then
                  call splines(ke)
                  dejaspl=1
               end if
            end do

      end if

      if (estim0.eq.1) then
         
         IF (npm0.eq.1) then
            go to 1589
         else
            
            id=0
            jd=0
            thi=0.d0
            thj=0.d0
            
            loglik=vrais_comp(b0,npm0,id,thi,jd,thj)
            
         end if
         
      else            

         !  injecter le b estime dans btot
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
         
         
         ! calculs post-estimation
         

      ! proba d'appartenance aux classes
         ppi=0.d0
         ppiy=0.d0
         if(ng>1) then
            call postprob_comp(btot,npmtot,ppi,ppiy)
         end if
!print*,"fin postprob"
      ! residus et predictions des effets aleatoires
         call residuals_comp(btot,npmtot,ppi,resid_m,pred_m_g,resid_ss,pred_ss_g,pred_RE,Yobs)

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
         if(nea.gt.0) then
            call scoretest_comp(btot,npmtot,statglob,statevt)
         else
            statglob=9999.d0
            statevt=9999.d0
         end if


      ! grille de valeur pour les risques de base

        allocate(brisq_est(maxval(nprisq)))

        dejaspl=0
        if(any(typrisq.eq.3)) then
           allocate(Tmm_est(nsim),Tmm1_est(nsim),Tmm2_est(nsim),Tmm3_est(nsim), &
                Tim_est(nsim),Tim1_est(nsim),Tim2_est(nsim),Tim3_est(nsim))

           do ke=1,nbevt
              if(typrisq(ke).eq.3 .and. dejaspl.eq.0) then
                 call splines_estime_comp(time,nsim,ke)
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
              
              call fct_risq_estime_comp(ke,brisq_est,time,nsim,g,risq_est,risqcum_est)

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

!print*,"fin risq"

        ! grille de valeurs pour la transfo
        call transfo_estimee_comp(btot,npmtot,nsim,marker,transfY)

     end if

     

 1589 continue



!--------------- liberation de la memoire ------------------------


      if (any(typrisq.eq.3)) then
         deallocate(Tmm,Tmm1,Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,    &
              Tmm01,Tmm02,Tmm03,Tim0,Tim01,Tim02,Tim03,Tmmt,Tmmt1,     &
              Tmmt2,Tmmt3,Timt,Timt1,Timt2,Timt3)
      endif


 1236 continue

      deallocate(Y,X,idprob,idea,idg,idcor,nmes,Tsurv0,Tsurv,Tsurvint &
     ,ind_survint,zi,devt,prior,pprior,zitr,mm,mm1,mm2,im,im1,im2 &
     ,indiceY,uniqueY,risqcom,typrisq,nz,nprisq,nrisq,idcom,idspecif,idtdv &
     ,nxevtspec,nevtparx,nxcurr)


      deallocate(bfix,fix)

      return


end subroutine loglikjointlcmm
