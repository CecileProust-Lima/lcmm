!====================================================================


! Integration points and weights with MC or GH method


!===================================================================


subroutine Integ(X0, idprob, idea, idg, idcor, ng, ncor, nv, maxmes, idiag, nwg, &
     npm, b1, nsim, methInteg, points, weights)

  IMPLICIT NONE

  ! in input
  integer, intent(in) :: ng, ncor, nv, maxmes, idiag, nwg, npm, nsim, methInteg
  double precision, dimension(maxmes * nv), intent(in) :: X0
  integer, dimension(nv), intent(in) :: idprob, idea, idg, idcor
  double precision, dimension(npm), intent(in) :: b1

  ! for computation
  integer :: j, k, l, m, g, l2, m2, jj, npm2, j1, j2
  integer :: ier, nmoins, kk, nvarprob, ncg, ncssg, nea, nef, nvc, nprob
  double precision, dimension(nv,nv) :: Ut, Ut1
  double precision, dimension(:,:), allocatable :: VC, Z, P, R, X00, X2
  double precision, dimension(:), allocatable :: Vi, ysim, usim, mu, tcor
  double precision, dimension(nv) :: b0, b2
  double precision :: eps, SX, x22
  double precision,dimension(2,51) :: gauss

  ! for output
  double precision, dimension(maxmes * ng * nsim), intent(out) :: points
  double precision, dimension(nsim), intent(out) :: weights

  
  ! GetRNGstate
  call getrand()


  allocate(ysim(maxmes), usim(maxmes), mu(maxmes), tcor(maxmes), Vi(maxmes*(maxmes+1)/2), &
       VC(maxmes,maxmes), Z(maxmes,nv), P(maxmes,nv), R(maxmes,maxmes), X00(maxmes,nv), X2(maxmes,nv))


  eps=1.d-20


  ! creation des parametres

  nea = 0
  ncg = 0
  ncssg = 0
  nprob = 0 
  nvarprob = 0
  do k = 1,nv
     if (idg(k).eq.1) then
        ncssg = ncssg + 1      ! nb var. sans melange
     else if (idg(k).eq.2) then
        ncg = ncg + 1      ! nb var. dans melange
     end if
     nea = nea + idea(k)
     nprob = nprob + (idprob(k)) * (ng-1)
     nvarprob = nvarprob + idprob(k)
  end do

  if((ng.eq.1.and.ncg.gt.0).or.(ng.eq.1.and.nprob.gt.0)) then
     points = 9999.d0
     weights = 9999.d0
     go to 654
  end if


  !  nb effets fixes = nb effets fixes sans melange
  !                  + ng fois le nb de var dans melange

  nvc = nea
  if(idiag.eq.0) then
     nvc = nea * (nea + 1) / 2
  end if

  nef = ncssg + ncg * ng
  npm2 = nprob + nef + nvc + nwg + ncor + 1



  if (npm.ne.npm2) then
     points = 9999.d0
     weights = 9999.d0
     goto 654
  end if


!========== debut des calculs =======================
!   Variance of random-effects


  Ut = 0.d0
  If (idiag.eq.1) then
     do j = 1,nea
        do k = 1,nea
           if (j.eq.k) then
              Ut(j,k) = b1(nprob + nef + j)
           else
              Ut(j,k) = 0.d0
           end if
        end do
     end do
  end if

  If (idiag.eq.0) then
     do j = 1,nea
        do k = 1,j
           Ut(j,k) = b1(nprob + nef + k + j * (j-1) / 2)
        end do
     end do
  end if
  Ut1 = Ut

  points = 0.d0
  weights = 0.d0


! -------- creation de Vi = ZiGZi'+se*seIni ----------
! creation de Zi

  Z = 0.d0
  l = 0
  do k = 1,nv
     if (idea(k).eq.1) then
        l = l + 1
        do j = 1,maxmes
           Z(j,l) = dble(X0(maxmes * (k - 1) + j))
        end do
     end if
  end do
  
  ! creation de Ri
  R = 0.d0
  tcor = 0.d0
  if (ncor.gt.0) then
     do k = 1,nv
        if (idcor(k).eq.1) then
           do j = 1,maxmes
              tcor(j) = X0(maxmes * (k - 1) + j)
           end do
        end if
     end do
  end if
  do j1 = 1,maxmes
     do j2 = 1,maxmes
        if (ncor.eq.1) then 
           R(j1,j2) = b1(nprob + nef + nvc + nwg + ncor) * b1(nprob + nef + nvc + nwg + ncor) * &
                min(tcor(j1), tcor(j2))
        else if (ncor.eq.2) then
           R(j1,j2) = b1(nprob + nef + nvc + nwg + ncor) * b1(nprob + nef + nvc + nwg + ncor) * &
                exp(-b1(nprob + nef + nvc + nwg + 1) * abs(tcor(j1) - tcor(j2)))
               end if
            end do
         end do     



! creation de P=Zi*Ut et V=P*P' que si non spec aux classes

  if (nwg.eq.0.OR.NG.EQ.1) then
     P = 0.d0
     P = MATMUL(Z, Ut)
     VC = 0.d0
     VC = MATMUL(P, transpose(P)) + R
     ! j'ajoute l'erreur apres


     ! Vi en vecteur
     jj = 0
     Vi = 0.d0
     do j = 1,maxmes
        VC(j,j) = VC(j,j) + b1(npm) * b1(npm) ! erreur de mesure
        do k = j,maxmes
           jj = j + k * (k - 1) / 2
           Vi(jj) = VC(j,k)
        end do
     end do

     if (methInteg.eq.1) then 
        CALL DMFSD(Vi, maxmes, EPS, IER)
        if (ier.eq.-1) then
           points = 9999.d0
           weights = 9999.d0
           goto 654
        end if
        
        VC = 0.d0
        do j = 1,maxmes
           do k = 1,j
              VC(j,k) = Vi(k + j * (j - 1) / 2)
           end do
        end do
     end if
  end if
     

  ! calcul des points et poids

  if (ng.eq.1) then ! cas 1 : ng=1

     b0 = 0.d0
     l = 0
     X00 = 0.d0
     do k = 1,nv
        if (idg(k).ne.0) then
           l = l + 1
           do j = 1,maxmes
              X00(j,l) = dble(X0(maxmes * (k - 1) + j))
           end do
           b0(l) = b1(nprob + l)
        end if
     end do
     
     mu=0.d0
     mu=matmul(X00,b0)

     
     if (methInteg.eq.1) then
        ! i.e. methode de MC a faire
        do l = 1,nsim
           usim = 0.d0
           ysim = 0.d0
           do m = 1,maxmes
              SX = 1.d0
              call bgos(SX, 0, usim(m), x22, 0.d0)
           end do
           
           ysim = mu + MATMUL(VC, usim)
           do j = 1,maxmes
              points((l - 1) * maxmes + j) = ysim(j)
           end do

           weights(l) = 1.d0 / dble(nsim)

        end do
        
     else

        call gausshermite(gauss,nsim)
        
        do l = 1,nsim
           do j = 1,maxmes
              points((l - 1) * maxmes + j) = mu(j) + sqrt(VC(j,j)) * gauss(1, l)
           end do
           
           weights(l) = gauss(2, l)
           
        end do
     end if

     
  else ! cas 2 :  ng>1  composantes

     ! creation des vecteurs de variables explicatives
     l = 0
     m = 0
     X00 = 0.d0
     X2 = 0.d0
     do k = 1,nv
        if (idg(k).eq.2) then
           l = l + 1
           do j = 1,maxmes
              X2(j, l) = dble(X0(maxmes * (k - 1) + j))
           end do
        else if (idg(k).eq.1) then
           m = m + 1
           do j = 1,maxmes
              X00(j, m) = dble(X0(maxmes * (k - 1) + j))
           end do
        end if
     end do

     b2 = 0.d0
     b0 = 0.d0
     do g = 1,ng
        nmoins = 0
        l2 = 0
        m2 = 0
        do k = 1,nv
           if (idg(k).eq.1) then
              m2 = m2 + 1
              b0(m2) = b1(nprob + nmoins + 1)
              nmoins = nmoins + 1
           else if (idg(k).eq.2) then
              l2 = l2 + 1
              b2(l2) = b1(nprob + nmoins + g)
              nmoins = nmoins + ng
           end if
        end do

!        write(*,*)'b2',b2
!        write(*,*)'b0',b0

        ! variance covariance si spec aux classes :

        Ut1 = Ut
        if (nwg.ne.0) then
           Ut1 = 0.d0
           if (g.eq.ng) then
              Ut1 = Ut
           else
              Ut1 = Ut * abs(b1(nprob + nef + nvc + g))
           end if

           P = 0.d0
           P = MATMUL(Z, Ut1)
           VC = 0.d0
           VC = MATMUL(P, transpose(P)) + R

           ! Vi en vecteur
           Vi = 0.d0
           jj = 0
           do j = 1,maxmes
              VC(j,j) = VC(j,j) + b1(npm) * b1(npm)
              do k = j,maxmes
                 jj = j + k * (k - 1) / 2
                 Vi(jj) = VC(j,k)
              end do
           end do

           ! cholesky de Vi
           if (methInteg.eq.1) then 
              CALL DMFSD(Vi, maxmes, EPS, IER)
              if (ier.eq.-1) then
                 points = 9999.d0
                 weights = 9999.d0
                 goto 654
              end if
              
              VC = 0.d0
              do j = 1,maxmes
                 do k = 1,j
                    VC(j,k) = Vi(k + j * (j - 1) / 2)
                 end do
              end do
           end if
        end if

        mu = 0.d0
        mu = matmul(X00, b0) + matmul(X2, b2)
        
        if (methInteg.eq.1) then

           ! i.e. methode de MC a faire
           do l = 1,nsim
              usim = 0.d0
              ysim = 0.d0
              do m = 1,maxmes
                 SX = 1.d0
                 call bgos(SX, 0, usim(m), x22, 0.d0)
              end do
              
              ysim = mu + MATMUL(VC, usim)
              
              do j = 1,maxmes
                 points((g - 1) * maxmes * nsim + (l - 1) * maxmes + j) = ysim(j)
              end do
              
              weights(l) = 1.d0 / dble(nsim)
           end do
           

        else

           call gausshermite(gauss,nsim)

           do l = 1,nsim
              do j = 1,maxmes
                 points((g - 1) * maxmes * nsim + (l - 1) * maxmes + j) =  mu(j) + sqrt(VC(j,j)) * gauss(1, l)
              end do

              weights(l) = gauss(2, l)
           end do
           
        end if
           
     end do ! fin g

  end if
  
654 continue
     
  deallocate(Z,P,R,X00,X2,ysim,usim,mu,tcor,VC,Vi)
  
  ! PutRNGstate
  call putrand()
  
  
  return
  
end subroutine Integ
















