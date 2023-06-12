

    !!! signification des variables de commun_mpj !!!

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
!!$  idlink : type de transformation (0=linear,1=beta,2=splines)
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


     module commun_mpj

      implicit none
      integer,save::nbK,ns,ng,nv2,ncssg,ncg,  &
      nvarprob,maxmes,nobs,nvarxevt,nySPL,ntotvalSPL, &
      nbevt,logspecif,idtrunc,nvdepsurv,nrisqtot,nxevt,npmtot
      integer,save::ntrtot,nprob,nmescur
      integer,dimension(:),allocatable,save::nef,ncontr,nvc,nw,ncor,nerr, &
           nalea,ntr,ny,nea,idiag,nv,contraint
      integer,dimension(:),allocatable,save::risqcom,typrisq,nz,nprisq,nrisq,nxevtspec,nevtparx,nxcurr
      double precision,dimension(:),allocatable,save::Y,uniqueY,minY,maxY,epsY
      integer,dimension(:),allocatable,save::indiceY
      double precision,dimension(:,:),allocatable,save ::X,Xns,zi,zitr
      double precision,dimension(:),allocatable,save::Tsurv0,Tsurv,Tsurvint
      integer,dimension(:),allocatable,save::Devt,ind_survint
      integer,dimension(:),allocatable,save ::idea,idg,idcor,idcontr
      integer,dimension(:),allocatable,save::idcom,idspecif,idtdv,idprob,idlink
      integer,dimension(:,:),allocatable,save::nmes,nmesparK
      integer,dimension(:),allocatable,save::prior,nvalSPL
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
    end module commun_mpj





! ===================================================================
! ===================  vraisemblance sujet i  =======================
! ===================================================================

double precision function vrais_mpj_i(b,npm,id,thi,jd,thj,i)
  
!  use parameters
  use commun_mpj
!  use optim

  implicit none

  integer::i,id,jd,npm
  double precision::thi,thj
  double precision,dimension(npm)::b

  integer::k,l,m,sumnrisq,ke,g,j,yk,nevtxcurr,nxevtcurr,sumMesYk,sumntr,numSPL
  integer::sumny,tmp,kk,ll,ii,j1,j2,jj,ier,nmoins,q,sumparK,sumnv
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
  double precision,dimension(maxval(nea),maxval(nea))::Ut
  double precision,dimension(maxmes,maxval(nea))::Z,P
  double precision,dimension(maxmes,maxmes)::Corr,VC
  double precision,dimension(maxmes)::tcor,mu,Y2,Y3
  double precision,dimension(nbK*maxmes)::Y1
  double precision,dimension(maxmes*(maxmes+1)/2)::Vi
  double precision,dimension(maxmes,(maxval(ncontr)+sum(idcontr)))::X01
  double precision,dimension(maxval(ncontr)+sum(idcontr))::b01
  double precision,dimension(maxmes,maxval(nv))::X00,X2
  double precision,dimension(maxval(nv))::b0,b2

  vrais_mpj_i=0.d0
  !print*,"vrais pour i=",i

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
!              if(i.eq.1) then
!print*,"tr lin prm 1 :",b1(nprob+nrisqtot+nvarxevt&
!                   +tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+ &
!                   sumntr+1)
!print*,"tr lin prm 2 :",b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+ &
!     nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr+2)
!end if
              jacobien = jacobien - log(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ &
                   ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr+2))
           end do

        else if (idlink(sumny+yk).eq.1) then  ! Beta link
!if(i.eq.1) then
!print*,"tr beta prm 1:",b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+&
 !               ncor(k)+nerr(k)+nalea(k)+sumntr+1)
!print*,"tr beta prm 2:",b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+&
!                ncor(k)+nerr(k)+nalea(k)+sumntr+2)
!print*,"tr beta prm 3:",b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+&
!                ncor(k)+nerr(k)+nalea(k)+sumntr+3)
!print*,"tr beta prm 4:",b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+&
!     ncor(k)+nerr(k)+nalea(k)+sumntr+4)
!end if
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
                 vrais_mpj_i=-1.d9
                 !print*,"-1.d9 Y1=999"
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
!if(i<=5 .and. id==0 .and. jd==0) print*,"ll=",ll,"im=",im(1:10),"im1=",im1(1:10),"im2=",im2(1:10)
           
              if (ll.lt.1.or.ll.gt.ntr(sumny+yk)-3) then          
                 vrais_mpj_i=-1.d9
                 !print*,"-1.d9 ll<1 ou ll>ntrtot-3",Y(nmescur+sumMesYk+sumparK+j)
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
        else if(idlink(sumny+yk).eq.-1) then
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
     sumparK=0
     sumnv=0
     do k=1,nbK
        
        sumntr=0

        !! variance des effets aleatoires
        Ut=0.d0
        if (nea(k)>0) then
           
           if(contraint(k).eq.2) then
              Ut=0.d0
              Ut(1,1)=1.d0
              if (nea(k)>1) then 

                 If (idiag(k).eq.1) then
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
        end if
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
       
        sumMesYk = 0
        do yk=1,ny(k)
           do j1=1,nmes(i,sumny+yk)
              if(contraint(k).eq.1) then
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
           sumMesYk = sumMesYk + nmes(i,sumny+yk)
           sumntr = sumntr + ntr(sumny+yk) ! sert a incrementer tmp
        end do


        if(nmesparK(i,k).eq.0) goto 100

        
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
           vrais_mpj_i=-1.d9
           !print*,"-1.d9 dsinv"
           !print*,"b=",b
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

        ! variables en effet fixe
        b0=0.d0
        b01=0.d0
        l=0
        ll=0
        m=0
        X00=0.d0
        X01=0.d0
        do kk=1,nv(k)
           if (idg(sumnv+kk).ne.0) then
              l=l+1
              do j=1,nmesparK(i,k)
                 X00(j,l)=dble(X(nmescur+sumparK+j,sumnv+kk))  
              end do
                 if (kk.eq.1.and.contraint(k).ne.0) then
                    b0(l)=0.d0 ! intercept fixe a 0
                 else
                    ll=ll+1
                    b0(l)=b1(nprob+nrisqtot+nvarxevt+tmp+ll)
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
        
        ! calcul de la vrais_Y
        mu=0.d0
        Y2=0.d0
        Y3=0.d0
        Y4=0.d0
        mu=matmul(X00,b0)+matmul(X01,b01)
        do j=1,nmesparK(i,k)
           Y2(j) = Y1(sumparK+j)-mu(j)
        end do
       ! Y2=Y1-mu
        Y3=matmul(VC,Y2)
        Y4=DOT_PRODUCT(Y2,Y3) 

      !  if(i.lt.5) then
      !     print*,"Y4=",Y4
      !     print*,"b0=",b0, " b01=",b01
      !     print*,"X00=",X00
      !     print*,"Z=",Z
      !     print*,"sumnv=",sumnv
      !     print*,"ncurr=",nmescur+sumparK+sumMesYk
      !     print*,"nv(k)=" ,nv(k)
      !     print*,"mu=",mu
      !     print*,"VC=",VC
      !  end if
        vrais_Y = vrais_Y-nmesparK(i,k)*dlog(dble(2*3.14159265))-det-Y4
        !if(k.eq.1) vrais_Y1 = vrais_Y1 -nmesparK(i,k)*dlog(dble(2*3.14159265))-det-Y4
        !if(k.eq.2) vrais_Y2 = vrais_Y2 -nmesparK(i,k)*dlog(dble(2*3.14159265))-det-Y4

100 continue
        
        sumny = sumny + ny(k)
        tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr
        sumparK = sumparK + nmesparK(i,k)
        sumnv = sumnv + nv(k)
     end do ! fin boucle nbK

     ! !! log-vraisemblance totale
     vrais_mpj_i = vrais_surv + vrais_Y/2.d0 + jacobien

     ! entree retardee
     if(idtrunc.eq.1) then 
        vrais_mpj_i = vrais_mpj_i - entretard
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
     vrais_mpj_i=0.d0
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

              if(contraint(k).eq.2) then
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

           sumMesYk=0
           do yk=1,ny(k)
              do j1=1,nmes(i,sumny+yk)
                 if(contraint(k).eq.1) then
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
              sumntr = sumntr + ntr(sumny+yk) ! sert a incrementer tmp
           end do

           if(nmesparK(i,k).eq.0) goto 200

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
              vrais_mpj_i=-1.d9
              !print*,"-1.d9 dsinv"
              !print*,"i=",i, "nmespark=",nmesparK(i,k)
              !print*,"Vi=",Vi
              !print*,"b=",b
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
                 if (kk.eq.1.and.contraint(k).ne.0) then
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
                 if (kk.eq.1.and.contraint(k).ne.0) then
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

200        continue
           
           sumny = sumny + ny(k)
           tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr
           sumparK=sumparK+nmesparK(i,k)
           sumnv = sumnv +nv(k)
        end do ! fin boucle K

        ! sommer sur les classes
        expo = expo + pi(g)*vrais_surv*vrais_Y
        
     end do ! fin boucle classes
     
     ! log-vraisemblance totale
     vrais_mpj_i = -sum(nmesparK(i,:))/2.d0*dlog(dble(2*3.14159265)) + log(expo) + jacobien
     ! if(i.le.3) then
     !    print*,"log(expo)=",log(expo)
     !    print*,"jac=",jacobien
     !    print*,"nmes=",sum(nmesparK(i,:))
     !    print*,"vrais_mpj_i=",vrais_mpj_i
     ! endif
     
     ! entree retardee
     if(idtrunc.eq.1) then 
        vrais_mpj_i = vrais_mpj_i - log(retard)
     end if
    ! print*,"apres idtrunc"
  end if

  
654 continue

  return
  
end function vrais_mpj_i






double precision function vrais_mpj(b,m,id,thi,jd,thj)


  use commun_mpj,only:ns,nmes,nmescur

  implicit none

  integer::m,i,id,jd
  double precision::thi,thj,vrais_mpj_i,temp
  double precision,dimension(m)::b

  
  nmescur=0
  vrais_mpj=0.d0
  do i=1,ns
     !call intpr("dans vrais tot i=",-1,i,1)
     !print*,vrais_mpj_i(b,m,id,thi,jd,thj,i) 
     temp=vrais_mpj_i(b,m,id,thi,jd,thj,i)  
     !call dblepr("vrais=",-1,temp,1) 
     !if(i.eq.1)write(*,*)"vrais= ",temp        
     vrais_mpj = vrais_mpj + temp
     !print*,"i=",i," vrais_mpj=",vrais_mpj
     !if (temp.ne.-1.d9 .or. temp/temp.ne.1) then 
     if (abs(temp+1.d9).lt.1.d-6 .or. temp.ne.temp) then 
        !write(*,*)"i=",i,"vrais= ",temp," id=",id," jd=",jd
        vrais_mpj = -1.d9
        ! if(verbose==1) write(*,*)"i=",i,"vrais= ",temp
        goto 541
     end if
     nmescur = nmescur + sum(nmes(i,:))
    ! print*,"nmescur=",nmescur
  end do
541 continue
    ! print*,"vrais_mpj=",vrais_mpj
  return

end function vrais_mpj

          



! =============================================
! subroutine de creation de design matrix
! =============================================

      subroutine design_splines_mpj (ier)

      use commun_mpj

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
        
   end if

end do


765     continue

end subroutine design_splines_mpj






subroutine splines_mpj(k)
  use commun_mpj
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

     if (abs(Tsurv(i)-zi(n-2,k)).lt.1.d-6) then
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

     if (abs(Tsurv(i)-zi(n-2,k)).gt.1.d-6) then

        Tmm3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
        Tmm2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
             +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))  &
             +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
        Tmm1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
             +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))   &
             +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
        Tmm(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

     end if

     if (abs(Tsurv(i)-zi(n-2,k)).lt.1.d-6) then

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

        if (abs(Tsurv0(i)-zi(n-2,k)).lt.1.d-6) then
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

        if (abs(Tsurv0(i)-zi(nz(k)-2,k)).gt.1.d-6) then

           Tmm03(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))

           Tmm02(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn))  &
                +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))   &
                +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmm01(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))    &
                +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmm0(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

        end if

        if (abs(Tsurv0(i)-zi(n-2,k)).lt.1.d-6) then

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

        if (abs(Tsurvint(i)-zi(nz(k)-2,k)).lt.1.d-6) then
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

        if (abs(Tsurvint(i)-zi(nz(k)-2,k)).gt.1.d-6) then

           Tmmt3(i) = ((4.d0*ht2*ht2*ht2)/(h*hh*hn*hh3))
           Tmmt2(i) = ((4.d0*hht*ht2*ht2)/(hh2*hh*h*hn)) &
                +((-4.d0*h2t*htm*ht2)/(hh2*h2n*hh*h))     &
                +((4.d0*h2t*h2t*ht)/(hh2*h2*h*h2n))
           Tmmt1(i) = (4.d0*(htm*htm*ht2)/(h3m*h2n*hh*h)) &
                +((-4.d0*htm*ht*h2t)/(h3m*h2*h*h2n))       &
                +((4.d0*ht3*ht*ht)/(h3m*h3*h2*h))
           Tmmt(i) = 4.d0*(ht*ht*ht)/(h4*h3*h2*h)

        end if

        if (abs(Tsurvint(i)-zi(nz(k)-2,k)).lt.1.d-6) then

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

end subroutine splines_mpj



subroutine fct_risq_mpj_i(i,k,brisq,g,risq,surv,surv0,survint)

  use commun_mpj

  implicit none

  integer::i,k,g
  double precision,dimension(nprisq(k))::brisq
  double precision,dimension(ng,nbevt)::risq,surv,surv0,survint

  integer::j,l,ll,kk,ii
  double precision::som

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
     if (abs(Tsurv(i)-zi(nz(k),k)).lt.1.d-6) then
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
        if (abs(Tsurv0(i)-zi(nz(k),k)).lt.1.d-6) then
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
        if (abs(Tsurvint(i)-zi(nz(k),k)).lt.1.d-6) then
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


end subroutine fct_risq_mpj_i



subroutine postprob_mpj(b,npm,ppi,ppiy)


  use commun_mpj
!  use optim

  implicit none

  integer::npm
  double precision,dimension(npm)::b
  double precision,dimension(ns,ng) ::ppiy,ppi

  integer::i,k,l,m,sumnrisq,ke,g,j,yk,nevtxcurr,nxevtcurr,sumMesYk,sumntr,numSPL
  integer::sumny,tmp,kk,ll,ii,j1,j2,jj,ier,nmoins,q,sumparK,sumnv,it
  double precision::eps,temp,fevt,varexpsurv,surv_glob
  double precision,dimension(npmtot)::b1
  double precision,dimension(ng)::pi
  double precision,dimension(nvarprob)::Xprob,bprob 
  double precision,dimension(maxval(nrisq))::brisq
  double precision,dimension(ng,nbevt)::risq,surv,surv0,survint
  double precision,dimension(nxevt)::Xevt,bevt
  double precision,dimension(nbevt)::bevtint
  double precision,dimension(-1:maxval(ntr)-3)::splaa
  double precision::aa1,bb1,dd1,aa,bb,betai,cc1,ytemp,som,eta0,vrais_Y,det,Y4
  double precision,dimension(maxval(nea),maxval(nea))::Ut
  double precision,dimension(maxmes,maxval(nea))::Z,P
  double precision,dimension(:,:),allocatable::Corr,VC
  double precision,dimension(maxmes)::tcor,mu,Y2,Y3
  double precision,dimension(nbK*maxmes)::Y1
  double precision,dimension(:),allocatable::Vi
  double precision,dimension(maxmes,(maxval(ncontr)+sum(idcontr)))::X01
  double precision,dimension(maxval(ncontr)+sum(idcontr))::b01
  double precision,dimension(maxmes,maxval(nv))::X00,X2
  double precision,dimension(maxval(nv))::b0,b2
  double precision,dimension(ng)::fi,fi1


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

  it=0
  do i=1,ns

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

     !! transformation des Y
     Y1=0.d0

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
                 Y1(sumparK+sumMesYk+j)=(dble(Y(it+sumMesYk+sumparK+j))-b1(nprob+nrisqtot+nvarxevt&
                      +tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+ &
                      sumntr+1)) &
                      /abs(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+ &
                      nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr+2))

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

                 ytemp=(dble(Y(it+sumMesYk+sumparK+j))-minY(sumny+yk)+epsY(sumny+yk)) &
                      /(maxY(sumny+yk)-minY(sumny+yk)+2*epsY(sumny+yk))
                 Y1(sumMesYk+sumparK+j)=(betai(aa,bb,ytemp)-cc1)/dd1

                 if (abs(Y1(sumMesYk+sumparK+j)-999.d0).lt.1.d-6) then
                    !vrais_mpj_i=-1.d9
                    !print*,"-1.d9 Y1=999"
                    goto 456
                 end if

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

              do j=1,nmes(i,sumny+yk)
                 ll=0

                 if (abs(Y(it+sumMesYk+sumparK+j)-zitr(ntr(sumny+yk)-2,numSPL)).lt.1.d-6) then
                    ll=ntr(sumny+yk)-3
                 end if

                 som=0.d0
                 do kk = 2,ntr(sumny+yk)-2
                    if ((Y(it+sumMesYk+sumparK+j).ge.zitr(kk-1,numSPL)).and. &
                         (Y(it+sumMesYk+sumparK+j).lt.zitr(kk,numSPL))) then
                       ll=kk-1
                    end if
                 end do

                 if (ll.lt.1.or.ll.gt.ntr(sumny+yk)-3) then          
                    !vrais_mpj_i=-1.d9
                    !print*,"-1.d9 ll<1 ou ll>ntrtot-3"
                    goto 456
                 end if
                 if (ll.gt.1) then
                    do ii=2,ll
                       som=som+splaa(ii-3)
                    end do
                 end if

                 Y1(sumMesYk+sumparK+j)=eta0+som +splaa(ll-2)*im2(indiceY(it+sumMesYk+sumparK+j)) &
                      +splaa(ll-1)*im1(indiceY(it+sumMesYk+sumparK+j))&
                      +splaa(ll)*im(indiceY(it+sumMesYk+sumparK+j))

              end do
           else if(idlink(sumny+yk).eq.-1) then
              do j=1,nmes(i,sumny+yk)
                 Y1(sumMesYk+sumparK+j)=Y(it+sumMesYk+sumparK+j)
              end do
           end if
           sumntr = sumntr + ntr(sumny+yk)
           sumMesYk = sumMesYk + nmes(i,sumny+yk)
        end do
        sumny = sumny + ny(k)
        tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr
        sumparK = sumparK + nmesparK(i,k)
     end do

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
       ! pig ok
     end if

     ! calcul des proba a posteriori de chaque classe
     fi=0.d0
     fi1=0.d0
     do g=1,ng

        if(nbevt.gt.0) then
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
                 fevt = log(risq(g,ke)) + varexpsurv                     
                 if (ind_survint(i).eq.1) then
                    fevt = fevt + bevtint(ke)
                 end if
              end if

              nxevtcurr=nxevtcurr+nxevtspec(ke)
           end do

           fi(g) = fevt-surv_glob

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

              if(contraint(k).eq.2) then
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
              if(nw(g).ne.0.and.g.lt.ng) then
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
                    Z(j,l)=dble(X(it+sumparK+j,sumnv+kk))
                 end do
              end if
           end do

           !print*,"Z ok"
           !! matrice Corr (cor+err+alea)
           Corr=0.d0
           tcor=0.d0
           if(ncor(k).gt.0) then

              do kk=1,nv(k)
                 if (idcor(sumnv+kk).eq.1) then
                    do j=1,nmesparK(i,k)
                       tcor(j) = X(it+sumparK+j,sumnv+kk)
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

           sumMesYk=0
           do yk=1,ny(k)
              do j1=1,nmes(i,sumny+yk)
                 if(contraint(k).eq.1) then
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
              sumntr = sumntr + ntr(sumny+yk) ! sert a incrementer tmp
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
              !vrais_mpj_i=-1.d9
              !print*,"-1.d9 dsinv"
              ! print*,"i=",i, "nmespark=",nmesparK(i,k)
              ! print*,"Vi=",Vi
              ! print*,"b=",b
              ! print*,"bfix=",bfix
              ! print*,"id=",id,"thi=",thi
              ! print*,"jd=",jd,"thj=",thj
              ! print*,"fix=",fix
              goto 456
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
                    X00(j,m)=dble(X(it+sumparK+j,sumnv+kk))
                 end do
                 if (kk.eq.1.and.contraint(k).ne.0) then
                    b0(m)=0.d0
                 else
                    b0(m)=b1(nprob+nrisqtot+nvarxevt+tmp+nmoins+1)
                    nmoins=nmoins+1
                 end if

                 ! avec mixture
              else if (idg(sumnv+kk).eq.2) then
                 l=l+1
                 do j=1,nmesparK(i,k)
                    X2(j,l)=dble(X(it+sumparK+j,sumnv+kk))
                 end do
                 if (kk.eq.1.and.contraint(k).ne.0) then
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
                       X01(sumMesYk+j,(q-1)*ny(k)+yk) = dble(X(it+sumparK+sumMesYk+j,sumnv+kk))
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

           ! calcul de la vrais pour Yik
           mu=0.d0
           y2=0.d0
           y3=0.d0
           y4=0.d0
           mu=matmul(X00,b0)+matmul(X2,b2)+matmul(X01,b01)

           do j=1,nmesparK(i,k)
              Y2(j) = Y1(sumparK+j)-mu(j)
           end do

           Y3=Matmul(VC,Y2)
           Y4=0.d0
           Y4=DOT_PRODUCT(Y2,Y3)
          

           fi(g) = fi(g) + (-nmesparK(i,k)*log(dble(2*3.14159265)) -det -Y4)/2.d0 

           fi1(g) = fi1(g) + (-nmesparK(i,k)*log(dble(2*3.14159265)) -det -Y4)/2.d0

           sumny = sumny + ny(k)
           tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr
           sumparK=sumparK+nmesparK(i,k)
           sumnv = sumnv +nv(k)
        end do ! fin boucle K

        fi(g) = exp(fi(g))
        fi1(g) = exp(fi1(g))
        
     end do ! fin boucle classe

     do g=1,ng
        ppiy(i,g)=pi(g)*fi1(g)/DOT_PRODUCT(pi,fi1)
        ppi(i,g)=pi(g)*fi(g)/DOT_PRODUCT(pi,fi)
     end do

     it = it + sum(nmes(i,:))

  end do !  fin boucle sujet

456 continue

  deallocate(VC,Vi,Corr)

  return

end subroutine postprob_mpj



subroutine residuals_mpj(b1,npm,ppi,resid_m,pred_m_g,resid_ss, &
     pred_ss_g,pred_RE,pred_RE_Y,Yobs)

  use commun_mpj
!  use optim

  implicit none

  integer ::i,j,k,l,m,g,jj,npm,j1,j2,ii,ll,yk,kk,q,iea,ialea
  integer ::ier,nmoins,it,sumea,sumalea,sumntr,sumparK,sumny,numSPL,tmp,sumMesYk,sumnv
  double precision,dimension(maxmes,maxval(nv)) ::X0,X2
  double precision,dimension(maxmes,maxval(nea)) ::Z,P
  double precision,dimension(maxval(nv)) ::err2
  double precision,dimension(nvarprob) ::Xprob, bprob
  double precision,dimension(maxval(nea),maxval(nea)) ::Ut
  double precision,dimension(maxmes,maxmes) ::VC,Corr
  double precision,dimension(maxmes,maxmes)::Valea,SigmaE,CovDev
  double precision,dimension(maxval(nea),maxmes)::covUY
  double precision,dimension(npm) ::b1
  double precision,dimension(:),allocatable ::Vi
  double precision,dimension(maxval(nv)) :: b0,b2
  double precision ::eps,det,temp
  double precision,dimension(nbK*maxmes) :: Y1
  double precision,dimension(maxmes) :: mu,Y2,pred1,err1,tcor
  double precision,dimension(ng) :: pi
  double precision,dimension(nobs)::resid_m,resid_ss
  double precision,dimension(ns*sum(nea(:)))::pred_RE
  double precision,dimension(ns*sum(nalea(:)))::pred_RE_Y
  double precision,dimension(nobs*ng)::pred_m_g,pred_ss_g
  double precision,dimension(ns,ng) ::PPI
  double precision,dimension(nobs)::Yobs
  double precision,dimension(-1:(ntrtot-3))::splaa
  double precision::aa1,bb1,dd1,aa,bb,betai,ytemp,som,cc1,eta0
  double precision,dimension(maxval(nalea),maxmes)::err3
  double precision,dimension(maxval(nalea))::err4 
  double precision,dimension(maxmes,(maxval(ncontr)+sum(idcontr)))::X01
  double precision,dimension(maxval(ncontr)+sum(idcontr))::b01

  allocate(Vi(maxmes*(maxmes+1)/2))

  VC=0.d0
  Vi=0.d0
  Corr=0.d0
  SigmaE=0.d0
  CovDev=0.d0
  CovUY=0.d0

  eps=1.D-20
  resid_m=0.d0
  pred_m_g=0.d0
  resid_ss=0.d0
  pred_ss_g=0.d0
  pred_RE=0.d0
  Yobs=0.d0
  pred_RE=0.d0
  it=0
  iea=0
  ialea=0
  kk=0

  do i=1,ns

     !! transformation des Y
     Y1=0.d0
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
                 Y1(sumparK+sumMesYk+j)=(dble(Y(it+sumMesYk+sumparK+j))-b1(nprob+nrisqtot+nvarxevt&
                      +tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+ &
                      sumntr+1)) &
                      /abs(b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+ &
                      nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr+2))
                 Yobs(it+sumparK+sumMesYk+j) = Y1(sumparK+sumMesYk+j)
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

                 ytemp=(dble(Y(it+sumMesYk+sumparK+j))-minY(sumny+yk)+epsY(sumny+yk)) &
                      /(maxY(sumny+yk)-minY(sumny+yk)+2*epsY(sumny+yk))
                 Y1(sumMesYk+sumparK+j)=(betai(aa,bb,ytemp)-cc1)/dd1
                 Yobs(it+sumparK+sumMesYk+j) = Y1(sumparK+sumMesYk+j)

                 if (abs(Y1(sumMesYk+sumparK+j)-999.d0).lt.1.d-6) then
                    !vrais_mpj_i=-1.d9
                    !print*,"-1.d9 Y1=999"
                    goto 456
                 end if

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

              do j=1,nmes(i,sumny+yk)
                 ll=0

                 if (abs(Y(it+sumMesYk+sumparK+j)-zitr(ntr(sumny+yk)-2,numSPL)).lt.1.d-6) then
                    ll=ntr(sumny+yk)-3
                 end if

                 som=0.d0
                 do kk = 2,ntr(sumny+yk)-2
                    if ((Y(it+sumMesYk+sumparK+j).ge.zitr(kk-1,numSPL)).and. &
                         (Y(it+sumMesYk+sumparK+j).lt.zitr(kk,numSPL))) then
                       ll=kk-1
                    end if
                 end do

                 if (ll.lt.1.or.ll.gt.ntr(sumny+yk)-3) then          
                    !vrais_mpj_i=-1.d9
                    !print*,"-1.d9 ll<1 ou ll>ntrtot-3"
                    goto 456
                 end if
                 if (ll.gt.1) then
                    do ii=2,ll
                       som=som+splaa(ii-3)
                    end do
                 end if

                 Y1(sumMesYk+sumparK+j)=eta0+som +splaa(ll-2)*im2(indiceY(it+sumMesYk+sumparK+j)) &
                      +splaa(ll-1)*im1(indiceY(it+sumMesYk+sumparK+j))&
                      +splaa(ll)*im(indiceY(it+sumMesYk+sumparK+j))
                 Yobs(it+sumparK+sumMesYk+j) = Y1(sumparK+sumMesYk+j)

              end do
           else if(idlink(sumny+yk).eq.-1) then
              do j=1,nmes(i,sumny+yk)
                 Y1(sumMesYk+sumparK+j)=Y(it+sumMesYk+sumparK+j)
                 Yobs(it+sumparK+sumMesYk+j) = Y1(sumparK+sumMesYk+j)
              end do
           end if
           sumntr = sumntr + ntr(sumny+yk)
           sumMesYk = sumMesYk + nmes(i,sumny+yk)
        end do
        sumny = sumny + ny(k)
        tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr
        sumparK = sumparK + nmesparK(i,k)
     end do

     sumny=0
     sumnv=0
     sumparK=0
     tmp=0
     sumea=0
     sumalea=0
     
     !! debut calcul des residus
     if(ng.eq.1) then

        do k=1,nbK
           sumntr=0

           !! variance des effets aleatoires
           Ut=0.d0
           if (nea(k)>0) then

              if(contraint(k).eq.2) then
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
           end if

           !! variables en effet aleatoire
           Z=0.d0
           l=0
           do kk=1,nv(k)
              if (idea(sumnv+kk).eq.1) then
                 l=l+1
                 do j=1,nmesparK(i,k)
                    Z(j,l)=dble(X(it+sumparK+j,sumnv+kk))
                 end do
              end if
           end do

           !! matrice Corr
           Corr=0.d0
           tcor=0.d0
           if(ncor(k).gt.0) then

              do kk=1,nv(k)
                 if (idcor(sumnv+kk).eq.1) then
                    do j=1,nmesparK(i,k)
                       tcor(j) = X(it+sumparK+j,sumnv+kk)
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


           sumMesYk = 0
           Valea=0.d0
           SigmaE=0.d0
           do yk=1,ny(k)
              do j1=1,nmes(i,sumny+yk)
                 if(contraint(k).eq.1) then
                    SigmaE(sumMesYk+j1,sumMesYk+j1) = 1
                 else
                    SigmaE(sumMesYk+j1,sumMesYk+j1) =   b1(nprob+nrisqtot+nvarxevt+tmp+&
                         nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+yk)**2 !variance de l'erreur yk
                 end if
                 if (nalea(k).eq.ny(k)) then ! intercept aleatoire de yk
                    do j2=1,nmes(i,sumny+yk)
                       Valea(sumMesYk+j1,sumMesYk+j2) = b1(nprob+nrisqtot+nvarxevt+tmp+&
                            nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+yk)**2
                    end do
                 end if
              end do
              sumMesYk = sumMesYk + nmes(i,sumny+yk)
              sumntr = sumntr + ntr(sumny+yk) ! sert a incrementer tmp
           end do

           covUY=0.d0
           VC=0.d0
           P=0.d0
           CovDev=0.d0

           P = MATMUL(Z,Ut)
           covUY = MATMUL(Ut,transpose(P))
           VC = MATMUL(P,transpose(P))+Corr+Valea+SigmaE
           CovDev=MATMUL(P,transpose(P))+Corr+Valea

           !     Vi en vecteur
           j=0
           do jj=1,nmesparK(i,k)
              do kk=jj,nmesparK(i,k)
                 j=jj+kk*(kk-1)/2
                 Vi(j)=VC(jj,kk)
              end do
           end do

           CALL dsinv(Vi,nmesparK(i,k),eps,ier,det)
           if (ier.eq.-1) then
              do j=1,nmesparK(i,k)
                 resid_m(it+sumparK+j)=9999.d0
                 pred_m_g(it+sumparK+j)=9999.d0
                 resid_ss(it+sumparK+j)=9999.d0
                 pred_ss_g(it+sumparK+j)=9999.d0
              end do
              do kk=1,nea(k)
                 pred_RE(iea+sumea+kk)=9999.d0
              end do
              do kk=1,nalea(k)
                 pred_RE_Y(ialea+sumalea+kk)=9999.d0
              end do
              goto 456
           end if

           !     retransformation du vecteur Vi en matrice :
           VC=0.d0
           do jj=1,nmesparK(i,k)
              do kk=1,nmesparK(i,k)
                 if (kk.ge.jj) then
                    VC(jj,kk)=Vi(jj+kk*(kk-1)/2)
                 else
                    VC(jj,kk)=Vi(kk+jj*(jj-1)/2)
                 end if
              end do
           end do


           ! variables en effet fixe
           b0=0.d0
           b01=0.d0
           l=0
           ll=0
           m=0
           X0=0.d0
           X01=0.d0
           do kk=1,nv(k)
              if (idg(sumnv+kk).ne.0) then
                 l=l+1
                 do j=1,nmesparK(i,k)
                    X0(j,l)=dble(X(it+sumparK+j,sumnv+kk))  
                 end do
                 if (kk.eq.1.and.contraint(k).ne.0) then
                    b0(l)=0.d0 ! intercept fixe a 0
                 else
                    ll=ll+1
                    b0(l)=b1(nprob+nrisqtot+nvarxevt+tmp+ll)
                 end if
              end if

              !contrast : 
              if (idcontr(sumnv+kk).ne.0) then
                 m=m+1
                 sumMesYk=0
                 do yk=1,ny(k)
                    ! creation matrice design des contrastes: X01
                    do j=1,nmes(i,sumny+yk)
                       X01(sumMesYk+j,(m-1)*ny(k)+yk) = dble(X(it+sumparK+sumMesYk+j,sumnv+kk))
                    end do
                    sumMesYk=sumMesYk+nmes(i,sumny+yk)
                    ! creation vecteur parms des contrastes: b01
                    if (yk<ny(k)) then
                       b01((m-1)*ny(k)+yk)=b1(nprob+nrisqtot+nvarxevt+tmp+nef(k)+(m-1)*(ny(k)-1)+yk)
                    else
                       b01((m-1)*ny+ny) =-sum(b1(nprob+nrisqtot+nvarxevt+tmp+ &
                            nef(k)+(m-1)*(ny(k)-1)+1 &
                            :nprob+nrisqtot+nvarxevt+tmp+nef(k)+(m-1)*(ny(k)-1)+ny(k)-1))
                    end if
                 end do
              end if
           end do


           mu=matmul(X0,b0)+matmul(X01,b01)

           do j=1,nmesparK(i,k)
              Y2(j)=Y1(sumparK+j)-mu(j)
           end do

           err1=0.d0
           err1=MATMUL(VC,Y2)
           err2=0.d0
           err2=MATMUL(covUY,err1)

           pred1=0.d0
           pred1=MATMUL(CovDev,err1)

           ! residus
           do j=1,nmesparK(i,k)
              resid_m(it+sumparK+j)=Y2(j)
              pred_m_g(it+sumparK+j)=mu(j)
              pred_ss_g(it+sumparK+j)=mu(j)+pred1(j)
              resid_ss(it+sumparK+j)=Y2(j)-pred1(j)
           end do

           !pred_RE
           do kk=1,nea(k)
              pred_RE(iea+sumea+kk)=err2(kk)
           end do

           !pred_RE_Y 
           err3=0.d0
           sumMesYk=0
           if (nalea(k).eq.ny(k)) then
              do yk=1,ny(k)
                 err3(yk,:) = Valea(sumMesYk+nmes(i,sumny+yk),:)
                 sumMesYk=sumMesYk+nmes(i,sumny+yk)
              end do
           end if
           err4=0.d0
           err4= matmul(err3,err1)
           do kk=1,nalea(k)
              pred_RE_Y(ialea+sumalea+kk) = err4(kk)
           end do


           tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr
           sumea = sumea + nea(k)
           sumalea = sumalea + nalea(k)
           sumparK = sumparK + nmesparK(i,k)
           sumnv = sumnv + nv(k)
           sumny = sumny+ ny(k)
        end do ! fin boucle k

     else ! ng>1

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
        ! pig ok

        ! boucle sur les classes
        do g=1,ng

           sumny=0
           tmp=0
           sumparK=0
           sumnv=0
           sumalea=0
           sumea=0
           do k=1,nbK
              !     print*,"k=",k
              sumntr=0

              !! variance des effets aleatoires
              Ut=0.d0
              if (nea(k)>0) then

                 if(contraint(k).eq.2) then
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
                 if(nw(g).ne.0.and.g.lt.ng) then
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
                       Z(j,l)=dble(X(it+sumparK+j,sumnv+kk))
                    end do
                 end if
              end do

              !print*,"Z ok"
              !! matrice Corr (cor+err+alea)
              Corr=0.d0
              tcor=0.d0
              if(ncor(k).gt.0) then

                 do kk=1,nv(k)
                    if (idcor(sumnv+kk).eq.1) then
                       do j=1,nmesparK(i,k)
                          tcor(j) = X(it+sumparK+j,sumnv+kk)
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

              sumMesYk=0
              SigmaE=0.d0
              Valea=0.d0
              do yk=1,ny(k)
                 do j1=1,nmes(i,sumny+yk)
                    if(contraint(k).eq.1) then
                       SigmaE(sumMesYk+j1,sumMesYk+j1) = 1
                    else
                       SigmaE(sumMesYk+j1,sumMesYk+j1) = b1(nprob+nrisqtot+nvarxevt+tmp+&
                            nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+yk)**2 !variance de l'erreur yk
                    end if
                    if (nalea(k).eq.ny(k)) then ! intercept aleatoire de yk
                       do j2=1,nmes(i,sumny+yk)
                          Valea(sumMesYk+j1,sumMesYk+j2) = b1(nprob+nrisqtot+nvarxevt+&
                               tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+yk)**2
                       end do
                    end if
                 end do
                 !   print*,"Corr",yk," ok"
                 sumMesYk = sumMesYk + nmes(i,sumny+yk)
                 sumntr = sumntr + ntr(sumny+yk) ! sert a incrementer tmp
              end do


              ! matrice de variance
              covUY=0.d0
              VC=0.d0
              P=0.d0
              CovDev=0.d0

              P = MATMUL(Z,Ut)
              covUY = MATMUL(Ut,transpose(P))
              VC = MATMUL(P,transpose(P))+Corr+Valea+SigmaE
              CovDev=MATMUL(P,transpose(P))+Corr+Valea


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
                 !print*,"-1.d9 dsinv i=",i,"g=",g
                 do j=1,nmesparK(i,k)
                    resid_m(it+sumparK+j)=9999.d0
                    pred_m_g(nobs*(g-1)+it+sumparK+j)=9999.d0
                    resid_ss(it+sumparK+j)=9999.d0
                    pred_ss_g(nobs*(g-1)+it+sumparK+j)=9999.d0
                 end do
                 do kk=1,nea(k)
                    pred_RE(iea+sumea+kk)=9999.d0
                 end do
                 do kk=1,nalea(k)
                    pred_RE_Y(ialea+sumalea+kk)=9999.d0
                 end do
                 goto 456
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


              ! variables en effet fixe
              l=0
              m=0
              q=0
              b2=0.d0
              b0=0.d0
              nmoins=0
              X0=0.d0
              X2=0.d0
              X01=0.d0
              b01=0.d0
              do kk=1,nv(k)

                 ! sans mixture
                 if (idg(sumnv+kk).eq.1) then
                    m=m+1
                    do j=1,nmesparK(i,k)
                       X0(j,m)=dble(X(it+sumparK+j,sumnv+kk))
                    end do
                    if (kk.eq.1.and.contraint(k).ne.0) then
                       b0(m)=0.d0
                    else
                       b0(m)=b1(nprob+nrisqtot+nvarxevt+tmp+nmoins+1)
                       nmoins=nmoins+1
                    end if

                    ! avec mixture
                 else if (idg(sumnv+kk).eq.2) then
                    l=l+1
                    do j=1,nmesparK(i,k)
                       X2(j,l)=dble(X(it+sumparK+j,sumnv+kk))
                    end do
                    if (kk.eq.1.and.contraint(k).ne.0) then
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
                          X01(sumMesYk+j,(q-1)*ny(k)+yk) = dble(X(it+sumparK+sumMesYk+j,sumnv+kk))
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


              mu=matmul(X0,b0)+matmul(X01,b01)+matmul(X2,b2)

              do j=1,nmesparK(i,k)
                 Y2(j)=Y1(sumparK+j)-mu(j)
              end do

              err1=0.d0
              err1=MATMUL(VC,Y2)
              err2=0.d0
              err2=MATMUL(covUY,err1)           
              pred1=0.d0
              pred1=MATMUL(CovDev,err1)

              ! residus
              do j=1,nmesparK(i,k)
                 pred_m_g(nobs*(g-1)+it+sumparK+j)=mu(j)
                 pred_ss_g(nobs*(g-1)+it+sumparK+j)=mu(j)+pred1(j)

                 resid_m(it+sumparK+j)=resid_m(it+sumparK+j)+pi(g)*(Y2(j))
                 resid_ss(it+sumparK+j)=resid_ss(it+sumparK+j)+ppi(i,g)*(Y2(j)-pred1(j))
              end do

              !pred_RE
              do kk=1,nea(k)
                 pred_RE(iea+sumea+kk)=pred_RE(iea+sumea+kk)+ppi(i,g)*err2(kk)
              end do

              !pred_RE_Y 
              err3=0.d0
              if (nalea(k).eq.ny(k)) then
                 do yk=1,ny(k)
                    err3(yk,:) = Valea(nmesparK(i,k),:)
                 end do
              end if
              err4=0.d0
              err4= matmul(err3,err1)
              do kk=1,nalea(k)
                 pred_RE_Y(ialea+sumalea+kk) = pred_RE_Y(ialea+sumalea+kk)+ppi(i,g)*err4(kk)
              end do


              tmp = tmp + nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr
              sumea = sumea + nea(k)
              sumalea = sumalea + nalea(k)
              sumparK = sumparK + nmesparK(i,k)
              sumnv = sumnv + nv(k)
              sumny = sumny + ny(k)
           end do ! fin boucle k

        end do ! fin boucle classe

     end if

456  continue

     it = it + sum(nmesparK(i,:))
     iea = iea + sum(nea(:))
     ialea = ialea + sum(nalea(:))

  end do ! fin boucle sujet

  deallocate(Vi)

end subroutine residuals_mpj




subroutine splines_estime_mpj(time,nsim,ke)
  use commun_mpj
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

end subroutine splines_estime_mpj





subroutine fct_risq_estime_mpj(k,brisq,time,nsim,g,risq,surv)


  use commun_mpj


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

end subroutine fct_risq_estime_mpj





subroutine transfo_estimee_mpj(b,npm,nsim,marker,transfY)

  use commun_mpj

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
end subroutine transfo_estimee_mpj







subroutine loglikmpjlcmm(K0,ny0,nbevt0,ng0,ns0,Y0,nobs0,X0,nv0,Xns0,nv20, &
     Prior0,Tentr0,Tevt0,Devt0,ind_survint0, &
     idnv0,idnv20,idspecif0,idlink0, &
     epsY0,nbzitr0,zitr0,uniqueY0,nvalSPL0,indiceY0, &
     typrisq0,risqcom0,nz0,zi0, &
     nmes0,nea0,nw0,ncor0,nalea0,idiag0,idtrunc0,logspecif0, &
     npm0,b0,ppi0,ppitest0,resid_m, &
     resid_ss,pred_m_g,pred_ss_g,pred_RE,pred_RE_Y, &
     time,risq_est,risqcum_est,marker,transfY,nsim,Yobs,statscoretest, &
     fix0,contrainte0,nfix0,bfix0,estim0,loglik)


  use commun_mpj

  IMPLICIT NONE

  !Declaration des variables en entree
  integer,intent(in)::K0,nbevt0,ns0,ng0,nv20,nobs0,npm0,nsim
  integer,dimension(K0)::ny0,nv0,nea0,nw0,ncor0,nalea0,idiag0,contrainte0
  integer,dimension(ns0,sum(ny0(:)))::nmes0   
  integer,intent(in)::idtrunc0,logspecif0,estim0,nfix0
  double precision,dimension(nfix0),intent(in)::bfix0
  double precision,dimension(nobs0),intent(in)::Y0
  integer,dimension(nobs0),intent(in)::indiceY0
  double precision,dimension(nobs0*sum(nv0(:))),intent(in)::X0
  double precision,dimension(ns0*nv20),intent(in)::Xns0
  double precision, dimension(ns0),intent(in)::Tentr0,Tevt0
  integer, dimension(ns0),intent(in)::prior0,Devt0,ind_survint0
  integer,dimension(4*sum(nv0(:))),intent(in)::idnv0
  integer,dimension(3*nv20),intent(in)::idnv20
  integer, dimension(nv20*nbevt0),intent(in)::idspecif0
  integer, dimension(sum(ny0(:))),intent(in)::idlink0,nbzitr0,nvalSPL0
  double precision, dimension(sum(ny0(:))),intent(in)::epsY0
  double precision,dimension(maxval(nbzitr0),sum(ny0(:))),intent(in)::zitr0
  double precision,dimension(sum(nvalSPL0(:))),intent(in)::uniqueY0 
  integer,dimension(nbevt0),intent(in)::typrisq0,risqcom0,nz0
  double precision,dimension(maxval(nz0),nbevt0),intent(in)::zi0
  integer,dimension(npm0+nfix0),intent(in)::fix0

  !Declaration des variable en entree et sortie
  double precision,dimension(npm0), intent(inout) :: b0

  !Declaration des variables en sortie
  double precision,intent(out)::loglik
  double precision,dimension(ns0*ng0),intent(out)::ppi0,ppitest0
  double precision,dimension(nobs0),intent(out)::resid_m,resid_ss,Yobs
  double precision,dimension(nobs0*ng0),intent(out)::pred_m_g
  double precision,dimension(nobs0*ng0),intent(out)::pred_ss_g
  double precision,dimension(ns0*sum(nea0(:))),intent(out)::pred_RE
  double precision,dimension(ns0*sum(nalea0(:))),intent(out)::pred_RE_Y
  double precision,dimension(nsim*sum(ny0)),intent(out)::marker,transfY
  double precision,dimension(nsim),intent(out)::time
  double precision,dimension(nsim*ng0,nbevt0),intent(out)::risq_est,risqcum_est
  double precision,dimension(1+nbevt0),intent(out)::statscoretest

  !Variables locales
  integer::jtemp,i,g,j,ier,k,ktemp,ig,ke,sumnrisq,dejaspl,nbfix,nytot
  integer::mi,k1,l,neftot,nvctot,ncontrtot,nwtot,ncortot,nerrtot,naleatot
  integer::sumny,m,jj,ll,npmtot0,k2,id,jd
  double precision::thi,thj
  double precision,dimension(ns0,ng0)::PPI,ppiy
  double precision,dimension(npm0+nfix0)::btot
  double precision,external::vrais_mpj


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

  ! alloc pour partie modele mixte
  nbK=K0
  nytot = sum(ny0(:))
  allocate(minY(nytot),maxY(nytot),idlink(nytot),ntr(nytot),epsY(nytot))
  allocate(nef(nbK),ncontr(nbK),nea(nbK),nvc(nbK),nw(nbK),ncor(nbK),nerr(nbK),nalea(nbK),ny(nbK),idiag(nbK),nv(nbK),&
       contraint(nbK))

  contraint=contrainte0
  ny=ny0
  epsY=epsY0
  nySPL=0 
  do k=1,nytot
     idlink(k)=idlink0(k)
     minY(k)=zitr0(1,k)
     maxY(k)=zitr0(nbzitr0(k),k)
     if (idlink(k).eq.2) then
        nySPL=nySPL+1
     end if
  end do
  if(nySPL>0) then 
     allocate(nvalSPL(nySPL))
     nvalSPL(1:nySPL)=nvalSPL0(1:nySPL)
  else
     allocate(nvalSPL(1))
     nvalSPL(1) = 0
  end if

  if(all(idlink/=2)) then
     ntotvalSPL=0
     allocate(zitr(1,1))
     allocate(mm(1),mm1(1),mm2(1),im(1),im1(1),im2(1))
     mm(1)=0.d0
     mm1(1)=0.d0
     mm2(1)=0.d0
     im(1)=0.d0
     im1(1)=0.d0
     im2(1)=0.d0
  else
     ntotvalSPL=sum(nvalSPL(:))
     allocate(zitr(-1:(maxval(nbzitr0)+2),nySPL))
     allocate(mm(ntotvalSPL),mm1(ntotvalSPL),mm2(ntotvalSPL),im(ntotvalSPL),im1(ntotvalSPL),im2(ntotvalSPL))
  end if


  zitr=0.d0  
  k1=0
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
  end do

  allocate(Y(nobs0),X(nobs0,sum(nv0(:))),Xns(ns0,nv20),prior(ns0),  &
       idea(sum(nv0(:))),idg(sum(nv0(:))),idcontr(sum(nv0(:))),&
       idcor(sum(nv0(:))), nmes(ns0,nytot),&
       idprob(nv20),idcom(nv20),idspecif(nv20*nbevt0),idtdv(nv20), &
       uniqueY(ntotvalSPL),indiceY(nobs0))

  ns=ns0
  ng=ng0
  nv=nv0
  nv2=nv20
  nobs=nobs0
  nbevt=nbevt0
  idtrunc=idtrunc0

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

  if (ntotvalSPL.gt.0) uniqueY(1:ntotvalSPL)=uniqueY0(1:ntotvalSPL)
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

     typrisq=typrisq0
     risqcom=risqcom0
     idtrunc=idtrunc0
     Tsurv0=Tentr0   
     Tsurv=Tevt0    
     devt=devt0    
     ind_survint=ind_survint0
     logspecif=logspecif0
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
           go to 1236
        end if

        nef(k)=ncssg+ng*ncg
        ncontr(k)=ncontr(k)+idcontr(ll+l)*(ny(k)-1)
        if(contraint(k).ne.0) nef(k)=nef(k)-1
        nea(k)=nea(k)+idea(ll+l)
        if(idiag(k).eq.1) then
           nvc(k)=nea(k)
        else
           nvc(k)=nea(k)*(nea(k)+1)/2
        end if
        if(contraint(k).eq.2) nvc(k)=nvc(k)-1
        if(contraint(k).ne.1) nerr(k)=ny(k)
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

  npmtot0 = npm0+nfix0

  if (npmtot.ne.npmtot0) then
     !print*,"npmtot",npmtot,"npmtot0",npmtot0
     !print*,nprob,nrisqtot,nvarxevt,nef,ncontr,nvc,nw,ncor,nerr,nalea,ntrtot
     goto 1236
  end if



  ! creer base de splines si transfo splines
  if (any(idlink.eq.2)) then 
     call design_splines_mpj(ier)
     if (ier.eq.-1) then          
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
  bfix(1:nbfix)=bfix0(1:nbfix)


  ! changer prm de B par prm a estimer
  ! tmp=0
  ! sumny=0
  ! sumntr=0
  ! do k=1,nbK

  !    if (idiag(k).eq.1) then
  !       do j=1,nvc(k)
  !          btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j)=dsqrt(abs(btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j)))
  !       end do
  !    else

  !       if(contraint.eq.2) then
  !          mvc(1)=1.d0
  !          do j=1,nvc(k)
  !             mvc(1+j)=btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j)
  !          end do
  !       else
  !          do j=1,nvc(k)
  !             mvc(j)=btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+j)
  !          end do
  !       end if

  !       eps=1.d-20
  !       call DMFSD(mvc,nea(k),eps,ier)

  !       if(contraint.eq.2) then
  !          do j=1,nvc(k)
  !             btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+j)=mvc(1+j)
  !          end do
  !       else
  !          do j=1,nvc(k)
  !             btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+j)=mvc(j)
  !          end do
  !       end if
  !    end if

  !    if (nw(k).gt.0) then
  !       do j=1,nw(k)
  !          btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+j)= &
  !               abs(btot(nprob+nrisqtot+nvarxevt+tmp+nef(k)+ncontr(k)+nvc(k)+j))
  !       end do
  !    end if

  !    sumntr=0
  !    do l=1,ny(k)
  !       sumntr = sumntr+ntr(sumny+l)
  !    end do

  !    tmp = tmp+nef(k)+ncontr(k)+nvc(k)+nw(k)+ncor(k)+nerr(k)+nalea(k)+sumntr
  !    sumny = sumny+ny(k)
  ! end do

  if (estim0.eq.1) then
     
     IF (npm0.eq.1) then
        go to 1589
     else
        
        id=0
        jd=0
        thi=0.d0
        thj=0.d0
       
        loglik=vrais_mpj(b0,npm0,id,thi,jd,thj)
        
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
     

     ! proba d'appartenance aux classes
     ppi=0.d0
     ppiy=0.d0
     if(ng>1) then
        call postprob_mpj(btot,npmtot,ppi,ppiy)
     end if

     ! residus et predictions des effets aleatoires
     call residuals_mpj(btot,npmtot,ppi,resid_m,pred_m_g,resid_ss,pred_ss_g,pred_RE,pred_RE_Y,Yobs)

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
        call transfo_estimee_mpj(btot,npmtot,nsim,marker,transfY)
     end if

     !print*,"fin transfo"



  end if



  !--------------- liberation de la memoire ------------------------


  if(nbevt.gt.0) then
     if (any(typrisq.eq.3)) then
        deallocate(Tmm,Tmm1,Tmm2,Tmm3,Tim,Tim1,Tim2,Tim3,Tmm0,    &
             Tmm01,Tmm02,Tmm03,Tim0,Tim01,Tim02,Tim03,Tmmt,Tmmt1,     &
             Tmmt2,Tmmt3,Timt,Timt1,Timt2,Timt3)
     endif
  end if

1589 continue
  
  deallocate(bfix,fix)

1236 continue

  deallocate(nmesparK,minY,maxY,idlink,ntr,epsY,nef,ncontr,nea,nvc,nw,ncor,nerr,&
       nalea,ny,idiag,nvalSPL,nv,contraint)

  deallocate(Y,X,Xns,idprob,idea,idg,idcontr,idcor,nmes,prior,zitr,mm,mm1,mm2,im,im1,&
       im2,indiceY,uniqueY,idcom,idspecif,idtdv)

  if(nbevt.gt.0) then
     deallocate(Tsurv0,Tsurv,Tsurvint,ind_survint,zi,devt,risqcom,typrisq, &
          nz,nprisq,nrisq,nxevtspec,nevtparx,nxcurr)
  else
     deallocate(nrisq)
  end if


  return


end subroutine loglikmpjlcmm


