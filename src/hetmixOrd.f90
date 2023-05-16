
!==========================================================
!
!      Latent class mixed model for curvilinear ordinal outcome
!           usingF a threshold transformation
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
!                                      21/06/2010
!===========================================================
! - Version fortran 90
!






! minY et maxY donnent le range du Y
! ide0 est de taille maxY-minY et dit si parm estime ou pas.
! Dans B, toutes les valeurs, en vérifiant que des 0 pour les non estimes



! ------------------------------------------------
!  Residus non faits au 23/06/2011




!----------------------------------------------------------
!
!- Module COMMUN avec les donnees dynamiques
!
!----------------------------------------------------------




module communo

  implicit none
  integer,save ::ns,ng,nv,idiag,ncssg,nvc,nea,ncg,nwg,npmtot2 &
       ,nprob,nvarprob,maxmes,nobs,ntrtot,nrisq,nvarxevt,nef
  double precision,dimension(:),allocatable,save::Y
  double precision,dimension(:,:),allocatable,save ::X
  integer,dimension(:),allocatable,save ::idea,idg,idprob,ide
  integer,dimension(:),allocatable,save :: nmes,prior
  double precision,dimension(:,:),allocatable,save ::pprior
  integer,save :: minY,maxY
  integer,dimension(:),allocatable,save::fix
  double precision,dimension(:),allocatable,save::bfix

end module communo



module donnees_indiv

  implicit none
  double precision,dimension(:,:),allocatable::Ut1
  double precision,dimension(:),allocatable::mu
  double precision,dimension(:,:),allocatable::Z
  integer::numpat,nmescur
  integer,parameter ::nf=1
  double precision,dimension(:),allocatable::b1
end module donnees_indiv



!-----------------------------------------------------------
!                        FUNCPA
!------------------------------------------------------------


double precision function funcpao(b,npm,id,thi,jd,thj)

  use communo
  !use optim
  use donnees_indiv



  IMPLICIT NONE
  integer ::i,j,k,l,m,g,l2,m2,id,jd,npm,it,kk
  integer ::nmoins
  double precision,dimension(maxmes,nv) ::X00,X2
  double precision,dimension(nv) ::Xprob
  double precision,dimension(nea,nea) ::Ut
  double precision,dimension(nea) ::Xea
  double precision,dimension(npm) :: b
  double precision,dimension(nv) :: b0,b2,bprob
  double precision :: vrais,thi,thj
  double precision :: expo,temp,vraisobs
  double precision,dimension(ng) :: pi


  allocate(Ut1(nea,nea),mu(maxmes),Z(maxmes,nea),b1(npmtot2))


  m=0
  l=0
  b1=0.d0
  do k=1,nef+nvc+nwg
     if (fix(k).eq.0) then 
        l=l+1
        b1(k)=b(l)
        if (id.eq.l) then
           b1(k)=b1(k)+thi
        end if
        if (jd.eq.l) then
           b1(k)=b1(k)+thj
        end if
     else
        m=m+1
        b1(k)=bfix(m)
     end if
  end do


  ! changement Cecile
  kk=nef+nwg+nvc
  do k=1,maxY-minY
     if (ide(k).eq.1) then
        kk=kk+1
        if (fix(kk).eq.0) then
           l=l+1
           b1(nef+nwg+nvc+k)=b(l)
           if (id.eq.l) then
              b1(nef+nwg+nvc+k)=b1(nef+nwg+nvc+k)+thi
           end if
           if (jd.eq.l) then
              b1(nef+nwg+nvc+k)=b1(nef+nwg+nvc+k)+thj
           end if
        else 
           m=m+1
           b1(nef+nwg+nvc+k)=bfix(m)
        end if
     else
        b1(nef+nwg+nvc+k)=0.d0
     end if
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
           Ut(j,k)=b1(nef+k+j*(+j-1)/2)
        end do
     end do
  end if

  Ut1=Ut


  ! ----------- boucle sur les individus -------------
  it=0
  vrais=0.d0
  do i=1,ns

     numpat=i
     nmescur=it

     ! -------- creation de Vi = ZiGZi'+se*seIni ----------
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




     ! cas 1 : ng=1

     if (ng.eq.1) then


        b0=0.d0
        l=0
        X00=0.d0
        do k=1,nv
           if (idg(k).ne.0) then
              l=l+1
              do j=1,nmes(i)
                 X00(j,l)=dble(X(it+j,k))
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
        mu=matmul(X00,b0)



        if (nea.gt.0) then 
           temp=vraisobs()      
        else
           Xea=0.d0
           call vraistot(nea,Xea,nf,temp)
        end if


        if(temp.lt.1.d-300) then
           !               write(*,*)'temp a 1.d-300',i
           temp=1.d-300
        end if

        vrais=vrais+log(temp)

        !            if(thi.eq.0.and.thj.eq.0) then 
        !               write(*,*)'i',i,temp,log(temp),vrais
        !            end if

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
           if (nwg.ne.0.and.g.ne.ng) then
              Ut1=Ut*abs(b1(nef+nvc+g))
           end if

           mu=0.d0
           mu=matmul(X00,b0)+matmul(X2,b2)




           if (nea.gt.0) then 
              temp=vraisobs()      
           else
              Xea=0.d0
              call vraistot(nea,Xea,nf,temp)
           end if

           if (temp.lt.1.d-300) then

              !               write(*,*)'temp a 1.d-300',i
              temp=1.d-300
           end if
           expo = expo+pi(g)*temp


        end do
        vrais=vrais+log(expo)

     end if

     it=it+nmes(i)
  end do

  ! FIN BOUCLE SUJET

  funcpao=vrais


  deallocate(Ut1,mu,Z,b1)

  !            if(thi.eq.0.and.thj.eq.0) then 
  !               write(*,*)'funcpao',funcpao
  !            end if




  return

end function funcpao



!-----------------------------------------------------------
!                        FUNCPIO
!------------------------------------------------------------


double precision function funcpio(b,npm,id,thi,jd,thj,i)

  use communo
  !use optim
  use donnees_indiv



  IMPLICIT NONE
  integer ::i,j,k,l,m,g,l2,m2,id,jd,npm,kk
  integer ::nmoins
  double precision,dimension(maxmes,nv) ::X00,X2
  double precision,dimension(nv) ::Xprob
  double precision,dimension(nea,nea) ::Ut
  double precision,dimension(nea) ::Xea
  double precision,dimension(npm) :: b
  double precision,dimension(nv) :: b0,b2,bprob
  double precision :: vrais,thi,thj
  double precision :: expo,temp,vraisobs
  double precision,dimension(ng) :: pi

  allocate(Ut1(nea,nea),mu(maxmes),Z(maxmes,nea),b1(npmtot2))


  m=0
  l=0
  b1=0.d0
  do k=1,nef+nvc+nwg
     if (fix(k).eq.0) then 
        l=l+1
        b1(k)=b(l)
        if (id.eq.l) then
           b1(k)=b1(k)+thi
        end if
        if (jd.eq.l) then
           b1(k)=b1(k)+thj
        end if
     else
        m=m+1
        b1(k)=bfix(m)
     end if
  end do


  ! changement Cecile
  kk=nef+nwg+nvc
  do k=1,maxY-minY
     if (ide(k).eq.1) then
        kk=kk+1
        if (fix(kk).eq.0) then
           l=l+1
           b1(nef+nwg+nvc+k)=b(l)
           if (id.eq.l) then
              b1(nef+nwg+nvc+k)=b1(nef+nwg+nvc+k)+thi
           end if
           if (jd.eq.l) then
              b1(nef+nwg+nvc+k)=b1(nef+nwg+nvc+k)+thj
           end if
        else 
           m=m+1
           b1(nef+nwg+nvc+k)=bfix(m)
        end if
     else
        b1(nef+nwg+nvc+k)=0.d0
     end if
  end do




  !      l=nef+nwg+nvc
  !      do k=1,maxY-minY
  !         if (ide(k).eq.1) then
  !            l=l+1
  !            b1(nef+nwg+nvc+k)=b(l)
  !            if (id.eq.l) then
  !               b1(nef+nwg+nvc+k)=b1(nef+nwg+nvc+k)+thi
  !            end if
  !            if (jd.eq.l) then
  !               b1(nef+nwg+nvc+k)=b1(nef+nwg+nvc+k)+thj
  !            end if
  !         else
  !            b1(nef+nwg+nvc+k)=0.d0
  !         end if
  !      end do





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

  vrais=0.d0
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
           ! idg ne 0 pour l'iintercept forcement donc on met le parm a 0
           if (k.eq.1) then
              b0(l)=0.d0
           else
              b0(l)=b1(nprob+l-1)
           end if
        end if
     end do


     mu=0.d0
     mu=matmul(X00,b0)



     if (nea.gt.0) then 
        temp=vraisobs()      
     else
        Xea=0.d0
        call vraistot(nea,Xea,nf,temp)
     end if


     if(temp.lt.1.d-300) then
        !               write(*,*)'temp a 1.d-300',i
        temp=1.d-300
     end if

     vrais=vrais+log(temp)


     !            if(thi.eq.0.and.thj.eq.0) then 
     !               write(*,*)'i',i,temp,log(temp),vrais
     !            end if

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

        do g=1,ng
           pi(g) = pi(g)*pprior(i,g)
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
           temp=vraisobs()      
        else
           Xea=0.d0
           call vraistot(nea,Xea,nf,temp)
        end if

        if (temp.lt.1.d-300) then

           !               write(*,*)'temp a 1.d-300',i
           temp=1.d-300
        end if
        expo = expo+pi(g)*temp


     end do
     vrais=vrais+log(expo)

  end if

  funcpio=vrais


  deallocate(Ut1,mu,Z,b1)

  return

end function funcpio








! ================================================================
!
!
!     Integration numerique : calcul de funcpa pour 1 sujet
!     dans vraisobs()
!
! =================================================================

double precision function vraisobs()

!  use lois_normales
  use donnees_indiv
  use communo


  implicit none
  integer::j
  external ::vraistot
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
     !        if (nea.gt.1) then

     MINPTS=30
     MAXPTS=500               !  100 !500
     EPSABS=1.d-100
     EPSREL=1.d-100
     RESTAR=0
     call  hrmsym( NDIM2, NF2, MINPTS, MAXPTS, vraistot, EPSABS, &
          EPSREL, RESTAR, RESULT, ABSERR2, NEVAL, IFAIL, WORK)

     if (result(1).le. 1.d-300) then
        result(1)=1.d-300
     end if
     vraisobs=result(1)

     return

  else

     !           if (nea.eq.1) then 

     npg=30
     ! on definit les points
     call gausshermite(gauss,npg)


     ! boucle pour faire l'integration multiple
     do j=1,npg
        Xea(1)=gauss(1,j)
        call vraistot(nea,Xea,nf2,funvls)
        result(1)=result(1)+funvls*gauss(2,j)
     end do

     if (result(1).le. 1.d-300) then
        result(1)=1.d-300
     end if

     vraisobs=result(1)
     return

  end if




end function vraisobs






! ===============================================================
!
!     Vraisemblance sachant les effets aleatoires
! ===============================================================

subroutine vraistot(NDIM2,Xea,nf2,FUNVLS)

  use donnees_indiv
  use communo
  !      use lois_normales

  implicit none
  double precision ::vraisind,gamma1,gamma0,alnorm
  double precision :: funvls,sup,inf
  integer::ndim2,nf2,i,k,j,ind
  logical ::upper
  double precision,dimension(NDIM2)::Xea
  double precision,dimension(NDIM2)::ui
  double precision, dimension(maxmes)::mu1


  !     write(*,*)'entree vraistot',NDIM2,Xea,NF2,FUNVLS,numpat,nea,nf

  i=numpat
  nf2=nf2



  mu1=0.d0
  ui=0.d0
  if (nea.ge.1) then 
     ui=MATMUL(Ut1,Xea)
     mu1=mu+MATMUL(Z,ui)
  else              
     mu1=mu          
  end if

  vraisind=1.d0

  !      write(*,*)'i',i,(ui(j),j=1,nea)
  !      write(*,*)'mu1',(mu1(j),j=1,nmes(i))
  !      write(*,*)'mu1',mu1

  upper=.false.


  do j=1,nmes(i)
     ind=0
     if (Y(nmescur+j).eq.minY) then
        gamma0=(b1(nef+nvc+nwg+1)-mu1(j)  )
        vraisind=vraisind*(alnorm(gamma0,upper))
        ind=1
     else
        sup=b1(nef+nvc+nwg+1)
        inf=sup
        do k=1,maxY-minY-1
           sup=sup+ b1(nef+nvc+nwg+k+1)*b1(nef+nvc+nwg+k+1)
           if (Y(nmescur+j).eq.dble(minY+k)) then
              gamma1=(sup-mu1(j))
              gamma0=(inf-mu1(j))
              vraisind=vraisind*(alnorm(gamma1,upper) &
                   -alnorm(gamma0,upper))
              ind=1
           end if
           inf=sup
        end do
        if (Y(nmescur+j).eq.maxY) then
           gamma0=(sup-mu1(j))
           vraisind=vraisind*(1.d0-alnorm(gamma0,upper))
           ind=1
        end if

     end if


     !         write(*,*)'Y',i,j,Y(nmescur+j),gamma0,gamma1

     !         if (ind.eq.0) then
     !            write(*,*)'problem with modality for subject',i &
     !           ,'and occasion',j
     !            write(*,*)'Y=',Y(nmescur+j),'while ordinal outcome in'  &
     !                ,minY, maxY
     !            stop
     !         endif

  end do



  FUNVLS=vraisind

end subroutine vraistot



!!!!!!!!!!!!!!!!!!!!!!!!  fin vraistot !!!!!!
































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

subroutine postprobo(b,npm,PPI)
  use communo
!  use optim
  use donnees_indiv
  implicit none

  integer ::i,j,k,l,m,g,l2,m2,it,npm,nmoins
  double precision,dimension(maxmes,nv) ::X0,X2
  double precision,dimension(nv) ::Xprob
  double precision,dimension(nea,nea) ::Ut    
  double precision,dimension(nea) ::Xea
  double precision,dimension(npm) :: b
  double precision,dimension(nv) :: b0,b2,bprob
  double precision :: temp,f,vraisobs
  double precision,dimension(ng) ::fi,pi
  double precision,dimension(ns,ng) ::PPI




  allocate(Ut1(nea,nea),mu(maxmes),Z(maxmes,nea),b1(npmtot2))


  PPI=0.D0

  b1=0.d0
  do k=1,nef+nvc+nwg
     b1(k)=b(k)
  end do

  l=nef+nwg+nvc
  do k=1,maxY-minY
     if (ide(k).eq.1) then
        l=l+1
        b1(nef+nwg+nvc+k)=b(l)
     else
        b1(nef+nwg+nvc+k)=0.d0
     end if
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

  Ut1=Ut

  ! ----------- boucle sur les individus -------------
  it=0
  do i=1,ns




     numpat=i
     nmescur=it

     ! -------- creation de Vi = ZiGZi'+se*seIni ----------

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
        if (nwg.ne.0.and.g.ne.ng) then
           Ut1=Ut*abs(b1(nef+nvc+g))
        end if


        mu=0.d0
        mu=matmul(X0,b0)+matmul(X2,b2)

        if (nea.gt.0) then 
           temp=vraisobs() 
        else
           Xea=0.d0
           call vraistot(nea,Xea,nf,temp)
        end if
        if(temp.lt.1.d-300) then
           temp=1.d-300
        end if

        fi(g)=temp

     end do
     f=DOT_PRODUCT(pi,fi)
     do g=1,ng
        PPI(i,g)=pi(g)*fi(g)/f
     end do

     it=it+nmes(i)

  end do



  deallocate(Ut1,mu,Z,b1)
  return

end subroutine postprobo






!---------------------------------------------------------------

!residuals non fait

!---------------------------------------------------------------






!---------------------------------------------------------------

!transfo estimee

!---------------------------------------------------------------

subroutine transfo_estimee_ord(b,npm,marker,transfY)

  use communo

  implicit none

  double precision,dimension(2*(maxY-minY+1))::marker,transfY
  integer::npm,k,l
  double precision, dimension(npm)::b
  double precision, dimension(npmtot2)::b1
  double precision::inf,sup

  b1=0.d0
  do k=1,nef+nvc+nwg
     b1(k)=b(k)
  end do

  l=nef+nwg+nvc
  do k=1,maxY-minY
     if (ide(k).eq.1) then
        l=l+1
        b1(nef+nwg+nvc+k)=b(l)
     else
        b1(nef+nwg+nvc+k)=0.d0
     end if
  end do

  !       write(*,*)'infos',minY,maxY,nsim,npm
  !       write(*,*)'b',(b1(j),j=1,npm)

  marker=0.d0
  transfY=0.d0

  l=nef+nwg+nvc
  marker(1)=minY
  marker(2)=minY
  transfY(1)=-1.d10
  transfY(2)=b1(l+1)
  sup=b1(l+1)
  if (maxY-minY.gt.1) then
     do k=1,maxY-minY-1
        marker(2*k+1)=minY+k
        marker(2*k+2)=minY+k
        inf=sup
        sup=sup+b1(l+1+k)*b1(l+1+k)
        transfY(2*k+1)=inf
        transfY(2*k+2)=sup
     end do
  end if
  marker(2*(maxY-minY)+1)=maxY
  transfY(2*(maxY-minY)+1)=sup
  marker(2*(maxY-minY)+2)=maxY
  transfY(2*(maxY-minY)+2)=1.d10




end subroutine transfo_estimee_ord






!===========================================================
!      DERIV pour UACV 
!===========================================================


subroutine computUACVo(b,m,rlindiv,vopt,UACV)

  ! Calcul du gradient et de la matrice remplacant la hessienne
  ! par Fisher scoring empirique
  use communo
  use donnees_indiv,only:nmescur

  IMPLICIT NONE

  double precision::funcpio, &
       rldiscret,UACV,trace,th0,thn,th
  integer::m,i,k,id,j
  double precision,dimension(m,1)::Uscore
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
  DO i=1,ns
     Uscore=0.d0
     rlindiv(i)=funcpio(b,m,id,th0,id,th0,i)
     rldiscret=rldiscret+rlindiv(i)
     do k=1,m
        th=DMAX1(1.d-6, 1.d-4 * DABS(b(k)))
        thn=-1.D0*th
        Uscore(k,1)=-(funcpio(b,m,k,th,id,th0,i) &
             - funcpio(b,m,k,thn,id,th0,i))/(2.d0*th)
     END DO

     !           write(*,*)'Uscore',Uscore
     !           write(*,*)'Uscore2',Uscore2

     J_cond=J_cond+MATMUL(Uscore,transpose(Uscore))
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

  !        write(*,*)'MAT',MAT
  !       write(*,*)'J_cond',J_cond
  !      write(*,*)'H_1',H_1

  return

end subroutine computUACVo




subroutine logliklcmmord(Y0,X0,Prior0,pprior0,idprob0,idea0,idg0,ns0,ng0,nv0,nobs0, &
     nea0,nmes0,idiag0,nwg0,npm0,b0,ppi0,resid_m, &
     resid_ss,pred_m_g,pred_ss_g,pred_RE, &
     minY0,maxY0,ide0,marker,transfY,UACV,rlindiv,V,fix0,nfix0,bfix0,estim0,loglik)


  use communo

  IMPLICIT NONE

  !Declaration des variables en entree
  integer,intent(in)::nv0,nfix0,estim0
  integer,intent(in)::minY0,maxY0
  integer,dimension(maxY0-minY0),intent(in)::ide0

  integer, intent(in)::ns0,ng0,nobs0,idiag0,nwg0,npm0,nea0
  integer, dimension(nv0),intent(in)::idea0,idg0,idprob0
  integer, dimension(ns0),intent(in)::nmes0,prior0
  double precision, dimension(ns0*ng0), intent(in) :: pprior0
  double precision,dimension(nobs0),intent(in)::Y0
  double precision,dimension(nobs0*nv0),intent(in)::X0
  integer,dimension(npm0+nfix0),intent(in)::fix0
  double precision,dimension(nfix0),intent(in)::bfix0
  !Declaration des variable en entree et sortie
  double precision, dimension(npm0), intent(in) :: b0
  double precision,dimension(npm0*(npm0+3)/2)::V ! pour computeUACV
  !Declaration des variables en sortie
  double precision,intent(out)::loglik
  double precision,dimension(ns0*ng0),intent(out)::ppi0
  double precision,dimension(nobs0),intent(out)::resid_m,resid_ss
  double precision,dimension(nobs0*ng0),intent(out)::pred_m_g
  double precision,dimension(nobs0*ng0),intent(out)::pred_ss_g
  double precision,dimension(ns0*nea0),intent(out)::pred_RE
  double precision,dimension(2*(maxY0-minY0+1)),intent(out)::marker,transfY
  double precision,dimension(ns0),intent(out)::rlindiv
  !Variables locales
  integer::jtemp,i,g,j,ij,k,ktemp,ig,nmestot,it,npmtot,nbfix
  integer::k2,id,jd,npmtot0
  double precision::thi,thj
  double precision::UACV
  double precision,dimension(ns0,ng0)::PPI
  double precision,dimension(npm0+nfix0)::btot
  double precision,external::funcpao





  ! sorties initialisees

  ppi0=0.d0
  UACV=0.d0
  pred_ss_g=0.d0
  pred_m_g=0.d0
  pred_RE=0.d0
  marker=0.d0
  transfY=0.d0
  resid_m=0.d0
  resid_ss=0.d0
  rlindiv=0.d0
  loglik=0.d0


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

  minY=minY0
  maxY=maxY0

  allocate(Y(nobs0),idprob(nv0),X(nobs0,nv0) &
       ,idea(nv0),idg(nv0),nmes(ns0),prior(ns0),pprior(ns0,ng0),ide(maxY-minY))

  ! enrigstrement pour les modules
  ns=ns0
  ng=ng0
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
  nmes=0
  Y=0.d0
  X=0.d0
  idprob=0
  idea=0
  idg=0
  nmestot=0
  ktemp=0
  do k=1,nv
     idprob(k)=idprob0(k)
     idea(k)=idea0(k)
     idg(k)=idg0(k)

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
        
        do g=1,ng                 
           pprior(i,g)=pprior0((i-1)*ng+g)
        end do
     end do
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
  nprob=0 !ng-1
  nvarprob=min(ng-1,1)
  do k=1,nv
     if (idg(k).eq.1) then
        ncssg=ncssg+1      ! nb var. sans melange
     else if (idg(k).eq.2) then
        ncg=ncg+1      ! nb var. sans melange
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




  ntrtot=0
  ide=0
  do k=1,maxY0-minY0
     if(ide0(k).eq.1) then
        ntrtot=ntrtot+1
        ide(k)=1
     end if
  end do

  !     write(*,*)'ide',(ide(k),k=1,maxY-minY)
  !     write(*,*)'maxY,minY',maxY,minY,ntrtot


  nef=nprob+ncssg+ncg*ng+nvarxevt+nrisq-1
  npmtot=nef+nvc+nwg+ntrtot
  npmtot2=nef+nvc+nwg+maxY-minY


  !    if (nwg.gt.0) then
  !       do i=1,nwg
  !          btot(nef+nvc+i)=abs(btot(nef+nvc+i))
  !       end do
  !    end if



  ! lancement de l'optimisation

  IF (estim0.eq.1) then
     id=0
     jd=0
     thi=0.d0
     thj=0.d0

     loglik=funcpao(b0,npm0,id,thi,jd,thj)

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


     call transfo_estimee_ord(btot,npmtot,marker,transfY)


     ! probas posteriori

     !      write(*,*)'avant postprob'

     call computUACVo(b0,npm0,rlindiv,v,UACV)

     if (ng.gt.1) then
        call postprobo(btot,npmtot,PPI)
     end if

     ig=0
     ij=0
     do i=1,ns
        do g=1,ng0
           ig=ig+1
           ppi0(ig)=PPI(i,g)
        end do
     end do


  end if


1589 continue

  deallocate(Y,X,idprob,idea,idg,nmes,prior,pprior,ide)
  deallocate(fix,bfix)

  return
end subroutine logliklcmmord
