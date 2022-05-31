

subroutine predictmult(X0,idprob,idea,idg,idcor,idcontr &
     ,ng,ncor,nalea,nv,ny,maxmes,idiag,nwg,npm,b1,epsY,idlink &
     ,nbzitr,zitr0,modalite,nbmod,nsim,methInteg,Ymarg)  

!  use optim
  IMPLICIT NONE

  ! in input
  integer,intent(in)::ng,ncor,nv,maxmes,idiag,nwg,npm,nsim,methInteg,nalea,ny
  double precision,dimension(maxmes*nv),intent(in) ::X0
  integer,dimension(nv),intent(in)::idprob,idea,idg,idcor,idcontr
  integer,dimension(ny),intent(in)::idlink,nbzitr,nbmod
  double precision,dimension(npm),intent(in)::b1
  double precision,dimension(ny),intent(in) ::epsY
  double precision,dimension(maxval(nbzitr),ny),intent(in)::zitr0
  integer, dimension(sum(nbmod)),intent(in)::modalite

  ! for computation
  integer ::j,k,l,m,g,l2,m2,jj,npm2,j1,j2,yk,sumntrtot,q,sumnbmod
  integer ::ier,nmoins,kk,nrisq,nvarxevt,niter,nvarprob,ncg,ncssg,nea,nef,nvc,nprob,ncontr
!  double precision,dimension(nv) ::Xprob,  bprob
!  double precision::temp
!  double precision,dimension(ng) :: pi
  double precision,dimension(nv,nv) ::Ut,Ut1

  double precision,dimension(:,:),allocatable::VC,Z,P,Corr,X00,X2,X01,Sigma
  double precision,dimension(:),allocatable::Vi,ysim,usim,mu,tcor,usim2,asim,wsim
  double precision,dimension(maxmes*ny)::wi
  double precision,dimension(nv) :: b0,b2,uu
  double precision,dimension(nv*ny) :: b01 
  double precision :: eps,ytemp2,x22,ai
  double precision ::ytemp,diff,SX,beta
  double precision::aa1,bb1,dd1,aa,bb,cc1
  double precision,dimension(:),allocatable::zitr,splaa
  double precision::xinbta,INV_ISPLINES,beta_ln
  double precision,dimension(2,51)::gauss
  double precision,dimension(ny)::minY,maxY
  integer,dimension(ny)::ntrtot
  double precision::alnorm

  ! for output
  double precision,dimension(maxmes*ny,ng),intent(out) ::Ymarg

!     print*,"dans predictMult"

  !============= recup des places de parametres

  allocate(ysim(maxmes*ny),usim(maxmes*ny),mu(maxmes*ny),tcor(maxmes),Vi(maxmes*ny*(maxmes*ny+1)/2), &
       VC(maxmes*ny,maxmes*ny),Z(maxmes*ny,nv),P(maxmes*ny,nv),Corr(maxmes*ny,maxmes*ny), &
       X00(maxmes*ny,nv),X2(maxmes*ny,nv),X01(maxmes*ny,nv*ny), usim2(nv), Sigma(maxmes,maxmes), &
       asim(ny), wsim(maxmes*ny))

!       print*,"apres allocate1"
  ! en prevision de l'extension au conjoint
  nrisq=0
  nvarxevt=0
  ! fin en prevision

 do yk=1,ny 
 
  minY(yk)=zitr0(1,yk)
  maxY(yk)=zitr0(nbzitr(yk),yk)
  
!  if(verbose==1) print*,"min=",minY(yk),"max=",maxY(yk)

  if (idlink(yk).eq.0) ntrtot(yk)=2
  if (idlink(yk).eq.1) ntrtot(yk)=4
  if (idlink(yk).eq.2) ntrtot(yk)=nbzitr(yk)+2
  if (idlink(yk).eq.3) ntrtot(yk)=nbmod(yk)-1
  end do
  
!  if(verbose==1) print*,"ntrtot=",ntrtot
  
  if(all(idlink/=2)) then
         allocate(zitr(1),splaa(1))
  else
         allocate(zitr(-1:maxval(ntrtot)),splaa(-1:(maxval(ntrtot)-3)))
  end if
  
!  print*,"apres allocate2"
  zitr=0.d0
  splaa=0.d0
  ysim=0.d0
  usim=0.d0
  X2=0.d0
  X01=0.d0

  eps=1.d-20


!  write(*,*)'entree'
!  write(*,*)idprob
!  write(*,*)idea
!  write(*,*)idg
!  write(*,*)ng,nv,maxmes,idiag,nwg,npm
!  write(*,*)b1
!  write(*,*)'epsY',epsY,idlink,nbzitr
!  write(*,*)'zitr',zitr0
!  write(*,*)nsim,methInteg
!  write(*,*)'revu'


! if (verbose==1 .and. methInteg==1) print*,"Monte Carlo method"
! if (verbose==1 .and. methInteg==0) print*,"Gauss Hermite method"

  ! creation des parametres

  nea=0
  ncg=0
  ncssg=0
  nprob=0
  ncontr=0
  nvarprob=0
  do k=1,nv
     if (idg(k).eq.1) then 
        ncssg=ncssg+1      ! nb var. sans melange
     else if (idg(k).eq.2) then
        ncg=ncg+1      ! nb var. sans melange
     end if
     nea=nea+idea(k)
     nprob=nprob+(idprob(k))*(ng-1)
     nvarprob=nvarprob+idprob(k)
     ncontr=ncontr+idcontr(k)*(ny-1)
  end do

  if((ng.eq.1.and.ncg.gt.0).or.(ng.eq.1.and.nprob.gt.0)) then
     ymarg=9999.d0
!     if(verbose==1) print*,"pb : ng=",ng,"ncg =",ncg,"nprob=",nprob
     go to 654
  end if


  !  nb effets fixes = nb effets fixes sans melange
  !                  + ng fois le nb de var dans melange

  nvc=nea-1
  if(idiag.eq.0) then
     nvc=nea*(nea+1)/2-1
  end if

  nef=nprob+ncssg+ncg*ng+nvarxevt+nrisq-1+ncontr
  npm2=nef+nvc+nwg+ncor+ny+nalea+sum(ntrtot)

!if(verbose==1) then
! print*,"calcules : "
! print*,"nv=",nv
! print*,"idea=",idea,"idcontr=",idcontr
! print*,"ncontr=",ncontr,"nvc=",nvc,"nea=",nea
! print*,nprob,ncssg,ncg*ng,nvarxevt,nrisq,ncontr,nvc,nwg,ncor,ny,nalea,sum(ntrtot)
!end if

  if (npm.ne.npm2) then
     ymarg=9999.d0
!     if (verbose==1) print *,"pb npm : ","en entree :",npm,"calculé:",npm2
     goto 654
  end if


!========== debut des calculs =======================
!   Variance of random-effects


      Ut=0.d0
      Ut(1,1)=1.d0
      if (nea>1) then 
      
      If (idiag.eq.1) then
         do j=2,nea
            do k=2,nea
               if (j.eq.k) then
                  Ut(j,k)=b1(nef+j-1)
               else
                  Ut(j,k)=0.d0
               end if
            end do
         end do
      end if

      If (idiag.eq.0) then
         do j=2,nea
            do k=1,j
                 Ut(j,k)=b1(nef+k-1+j*(j-1)/2)
            end do
         end do
      end if
      
      end if
  Ut1=Ut

!  print*,"Ut ok"

  ymarg=0.d0


! -------- creation de Vi = ZiGZi'+Corr ----------
! creation de Zi

  Z=0.d0
  l=0
  do k=1,nv
     if (idea(k).eq.1) then
        l=l+1
        do yk=1,ny
          do j=1,maxmes
             Z(maxmes*(yk-1)+j,l)=dble(X0(maxmes*(k-1)+j))
          end do
        end do
     end if
  end do
!    print*,"Z ok"
  ! creation de Corr (contient erreurs de mesure, variance effet aleatoire sur test et AR/BM)  
        Corr=0.d0
        tcor=0.d0
        if (ncor.gt.0) then
           do k=1,nv
              if (idcor(k).eq.1) then
                 do j=1,maxmes
                    tcor(j) = X0(maxmes*(k-1)+j)
                 end do
              end if
           end do
         end if
        do yk=1,ny
         do j1=1,maxmes
            do j2=1,maxmes      
               if (ncor.eq.1) then 
                  Corr(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2) = b1(nef+nvc+nwg+ncor)*b1(nef+nvc+nwg+ncor)*&
                  min(tcor(j1),tcor(j2))
               else if (ncor.eq.2) then
                  Corr(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2) = b1(nef+nvc+nwg+ncor)*b1(nef+nvc+nwg+ncor)*&
                  exp(-b1(nef+nvc+nwg+1)*abs(tcor(j1)-tcor(j2)))
               end if

               if(all(idlink.ne.3)) then
                  if(nalea.eq.ny) Corr(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2) = Corr(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2) +&
                       b1(nef+nvc+nwg+ncor+ny+yk)**2 !variance de alpha_k
                  
                  if(j1.eq.j2) Corr(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2) = Corr(maxmes*(yk-1)+j1,maxmes*(yk-1)+j2)+&
                       b1(nef+nvc+nwg+ncor+yk) **2 ! variance de l'erreur
               end if
            end do
         end do
      end do


!      print*,"Corr ok"

        ! creation de P=Zi*Ut et V=P*P' que si non spec aux classes
           

      if(all(idlink.ne.3)) then
         
         if (nwg.eq.0.OR.NG.EQ.1) then
            P=0.d0
            P=MATMUL(Z,Ut)
            VC=0.d0
            VC=MATMUL(P,transpose(P))+Corr
            !       print*,"VC ok"
            ! Vi en vecteur
            
            jj=0
            Vi=0.d0
            do j=1,maxmes*ny
               do k=j,maxmes*ny
                  jj=j+k*(k-1)/2
                  Vi(jj)=VC(j,k)
               end do
            end do
            
            if (methInteg.eq.1) then 
               CALL DMFSD(Vi,maxmes*ny,EPS,IER)
               if (ier.eq.-1) then
                  ymarg=9999.d0
                  !           if (verbose==1) print*,"pb cholesky"
                  goto 654
               end if
               
               VC=0.d0
               do j=1,maxmes*ny
                  do k=1,j
                     VC(j,k)=Vi(k+j*(j-1)/2)
                  end do
               end do
            end if
         end if

      else

         if(any(idlink.eq.3)) then

            VC=0.d0
            do yk=1,ny
               
               Sigma=0.d0
               do j1=1,maxmes
                  do j2=1,maxmes 
                     if(nalea.eq.ny) Sigma(j1,j2) = b1(nef+nvc+nwg+ncor+ny+yk)**2 !variance de alpha_k
                     
                     if(j1.eq.j2) Sigma(j1,j2) = Sigma(j1,j2)+&
                          b1(nef+nvc+nwg+ncor+yk) **2 ! variance de l'erreur
                  end do
               end do
               
               jj=0
               Vi=0.d0
               do j=1,maxmes
                  do k=j,maxmes
                     jj=j+k*(k-1)/2
                     Vi(jj)=Sigma(j,k)
                  end do
               end do
               
               CALL DMFSD(Vi,maxmes,EPS,IER)
               if (ier.eq.-1) then
                  ymarg=9999.d0
                  goto 654
               end if
               
               VC=0.d0
               do j=1,maxmes
                  do k=1,j
                     VC((yk-1)*maxmes+j,(yk-1)*maxmes+k)=Vi(k+j*(j-1)/2)
                  end do
               end do
               
            end do
            
         end if
         
      end if
      

  ! cas 1 : ng=1

  if (ng.eq.1) then

     b0=0.d0
     l=0
     X00=0.d0
     b01=0.d0
     X01=0.d0
     m=0
     do k=1,nv
        if (idg(k).ne.0) then
           l=l+1
           do yk=1,ny
           do j=1,maxmes
              X00(maxmes*(yk-1)+j,l)=dble(X0(maxmes*(k-1)+j))
           end do
           end do
           ! idg ne 0 pour l'intercept forcement donc on met le parm a 0
           if (k.eq.1) then
              b0(l)=0.d0
           else
              b0(l)=b1(nprob+l-1)
           end if
        end if
        
         if (idcontr(k).ne.0) then
           m=m+1
           do yk=1,ny
! creation matrice design des contrastes: X01
                do j=1,maxmes
                  X01(maxmes*(yk-1)+j,(m-1)*ny+yk) = dble(X0(maxmes*(k-1)+j))
                end do
! creation vecteur parms des contrastes: b01
                if (yk<ny) THEN
                   b01((m-1)*ny+yk)=b1(nef-ncontr+(m-1)*(ny-1)+yk)
                else
                  b01((m-1)*ny+ny) =-sum(b1(nef-ncontr+(m-1)*(ny-1)+1 &
                             :nef-ncontr+(m-1)*(ny-1)+ny-1))
                end if
            end do          
         end if
        
     end do

     mu=0.d0
     mu=matmul(X00,b0)+matmul(X01,b01)
!        print*,"mu ok"
!    if (verbose==1) print*,"mu=",mu
    
     if(methInteg.eq.1) then !MC
        !         print*,"debut MC"
        do l=1,nsim
           
           if(any(idlink.ne.0)) then
              
              if(any(idlink.eq.3)) then

                 !! simuler les effets aleatoires
                 usim2=0.d0
                 ysim=0.d0
                 do m=1,nea
                    SX=1.d0
                    call bgos(SX,0,usim2(m),x22,0.d0)
                 end do
                 uu = 0.d0
                 uu = matmul(Ut,usim2)
                 ysim=matmul(X00,b0)+matmul(X01,b01)+matmul(Z,uu)

                 !! simuler le BM ou AR
                 if(ncor.gt.0) then
                    wsim=0.d0
                    x22=0.d0
                    SX=1.d0
                    do j=1,maxmes*ny
                       call bgos(SX,0,wsim(j),x22,0.d0)
                    end do
                    wi=0.d0
                    wi=matmul(Corr,wsim)
                    ysim = ysim + wi
                 end if
                 
            else
               
              usim=0.d0
              ysim=0.d0
              do m=1,maxmes*ny
                 SX=1.d0
                 call bgos(SX,0,usim(m),x22,0.d0)
              end do
              ysim=mu+MATMUL(VC,usim)

           end if
          end if

          sumntrtot=0
          sumnbmod = 0
          do yk=1,ny
             
             !          print*,"boucle sur idlink, yk=",yk, "nsim=",l
             if (idlink(yk).eq.0 .and. l.eq.1) then  ! Linear link
                
                aa = b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)
                bb = b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2)
                do j=1,maxmes
                   Ymarg(maxmes*(yk-1)+j,1) = mu(maxmes*(yk-1)+j)*bb+aa
                end do
                
                !sumntrtot = sumntrtot +ntrtot(yk) 
                
             else if (idlink(yk).eq.3) then

                 !! simuler l'EA specifique au test
                 ai=0.d0
                 if(nalea.gt.0) then
                    call bgos(SX,0,asim(yk),x22,0.d0)
                    ai = b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)*asim(yk)
                 end if
                 
                aa = b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)
                do j=1,maxmes
                   if(nalea.gt.0) ysim((yk-1)*maxmes+j) = ysim((yk-1)*maxmes+j)+ai
                   ytemp = modalite(sumnbmod + nbmod(yk)) - &
                        alnorm((aa-ysim((yk-1)*maxmes+j))/abs(b1(nef+nvc+nwg+ncor+yk)),.false.)
                   Ymarg(maxmes*(yk-1)+j,1) = Ymarg(maxmes*(yk-1)+j,1) + ytemp / dble(nsim)
                end do

                if(nbmod(yk).gt.2) then
                   do jj=2,nbmod(yk)-1
                      aa = aa + b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+jj)**2
                      do j=1,maxmes
                         ytemp = (modalite(sumnbmod + jj) - modalite(sumnbmod + jj+1))* &
                              alnorm((aa-ysim((yk-1)*maxmes+j))/abs(b1(nef+nvc+nwg+ncor+yk)),.false.)
                         Ymarg(maxmes*(yk-1)+j,1) = Ymarg(maxmes*(yk-1)+j,1) + ytemp / dble(nsim)                             
                      end do
                   end do
                end if
                
                sumnbmod = sumnbmod + nbmod(yk)
                
             else if (idlink(yk).eq.1) then  ! Beta link
                
                aa1=exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1))/ &
                     (1+exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)))
                bb1=exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2))/ &
                     (1+exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2)))
                bb1=aa1*(1.d0-aa1)*bb1
                
                cc1=abs(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+3))
                
                dd1=abs(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+4))
                
                aa=aa1*aa1*(1-aa1)/bb1-aa1
                bb=aa*(1-aa1)/aa1
                beta=beta_ln(aa,bb)

                !! simuler une N(0,Sigma)
                if(any(idlink.eq.3)) then
                   usim=0.d0
                   do m=1,maxmes
                      SX=1.d0
                      call bgos(SX,0,usim(m),x22,0.d0)
                   end do

                   Sigma=0.d0
                   do j1=1,maxmes
                      do j2=1,j1
                         Sigma(j1,j2) = VC((yk-1)*maxmes+j1,(yk-1)*maxmes+j2)
                      end do
                   end do

                   Vi=0.d0
                   Vi = matmul(Sigma,usim)
                   
                   do j=1,maxmes
                      ysim(maxmes*(yk-1)+j) = ysim(maxmes*(yk-1)+j) + Vi(j)
                   end do
                end if
                   
                
                do j=1,maxmes
                   ytemp=ysim(maxmes*(yk-1)+j)*dd1+cc1
                   if (ytemp.lt.0) then
                      ytemp=0.d0
                   end if
                 if (ytemp.gt.1) then
                    ytemp=1.d0
                 end if
                 ymarg(maxmes*(yk-1)+j,1)=ymarg(maxmes*(yk-1)+j,1)+xinbta(aa,bb,beta,ytemp,ier)/dble(nsim)
                 if (ier.ne.0.or.ymarg(maxmes*(yk-1)+j,1).eq.9999.d0) then
                    ymarg(maxmes*(yk-1)+j,1)=9999.d0
                    !                    if (verbose==1) print*,"pb beta"
                 end if
              end do
              
              !sumntrtot = sumntrtot +ntrtot(yk) 
              
           else if (idlink(yk).eq.2) then ! Splines link
              
              zitr=0.d0
              zitr(1:nbzitr(yk))=zitr0(1:nbzitr(yk),yk)
              zitr(-1)=zitr(1)
              zitr(0)=zitr(1)
              zitr(ntrtot(yk)-1)=zitr(ntrtot(yk)-2)
              zitr(ntrtot(yk))=zitr(ntrtot(yk)-1)
              
              !       if(verbose==1 .and. l==1) print*,"yk=",yk,"zitr=",zitr
              
              bb=b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)
              do kk=2,ntrtot(yk)
                 splaa(kk-3)=b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+kk)**2
              end do
              
              !! simuler une N(0,Sigma)
              if(any(idlink.eq.3)) then
                 usim=0.d0
                 do m=1,maxmes
                    SX=1.d0
                    call bgos(SX,0,usim(m),x22,0.d0)
                 end do
                   
                 Sigma=0.d0
                 do j1=1,maxmes
                    do j2=1,j1
                       Sigma(j1,j2) = VC((yk-1)*maxmes+j1,(yk-1)*maxmes+j2)
                    end do
                 end do

                 Vi=0.d0
                 Vi = matmul(Sigma,usim)
                 
                 do j=1,maxmes
                    ysim(maxmes*(yk-1)+j) = ysim(maxmes*(yk-1)+j) + Vi(j)
                 end do
              end if
                 
              do j=1,maxmes
                 niter=0
                 diff=0.d0
                 ier=0
                 ytemp=INV_ISPLINES(ysim(maxmes*(yk-1)+j),splaa,bb,nbzitr(yk),zitr,ier,niter,diff)
                 if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3).or.ymarg(maxmes*(yk-1)+j,1).eq.9999.d0) then
                    ymarg(maxmes*(yk-1)+j,1)=9999.d0
                    !                       if (verbose==1) print*,"pb inversion splines"
                 else
                    ymarg(maxmes*(yk-1)+j,1)=ymarg(maxmes*(yk-1)+j,1)+ytemp/dble(nsim)
                 end if
                 
              end do
              
              !sumntrtot = sumntrtot+ntrtot(yk)     
              
           end if
           sumntrtot = sumntrtot+ntrtot(yk) 
        end do
     end do
     ! fin MC
     
  else !GH
     !      print*,"debut GH"
     call gausshermite(gauss,nsim) ! on a les nsim poids et points d'integration dans gauss
 
   sumntrtot=0
   do yk=1,ny
!   print*,"boucle pour idlink, yk=",yk
     if (idlink(yk).eq.0) then  ! Linear link
         aa = b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)
         bb = b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2)
         do j=1,maxmes
           Ymarg(maxmes*(yk-1)+j,1) = mu(maxmes*(yk-1)+j)*bb+aa
         end do
         sumntrtot = sumntrtot +ntrtot(yk) 
         
     else if (idlink(yk).eq.1) then  ! Beta link
     
        aa1=exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1))/ &
             (1+exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)))
        bb1=exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2))/ &
             (1+exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2)))
        bb1=aa1*(1.d0-aa1)*bb1

        cc1=abs(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+3))

        dd1=abs(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+4))

        aa=aa1*aa1*(1-aa1)/bb1-aa1
        bb=aa*(1-aa1)/aa1
        beta=beta_ln(aa,bb)     
      
     
      do j=1,maxmes
        do l=1,nsim 
                ytemp=mu(maxmes*(yk-1)+j)+sqrt(VC(maxmes*(yk-1)+j,maxmes*(yk-1)+j))*gauss(1,l)
                ytemp=ytemp*dd1+cc1

                 if (ytemp.lt.0) then
                    ytemp=0.d0
                 end if
                 if (ytemp.gt.1) then
                    ytemp=1.d0
                 end if
                 ymarg(maxmes*(yk-1)+j,1)=ymarg(maxmes*(yk-1)+j,1)+xinbta(aa,bb,beta,ytemp,ier)*gauss(2,l)
                if (ier.ne.0.or.ymarg(maxmes*(yk-1)+j,1).eq.9999.d0) then
                    ymarg(maxmes*(yk-1)+j,1)=9999.d0
!                    if (verbose==1) print*,"pb beta"
                 end if
        end do
      end do
      sumntrtot = sumntrtot +ntrtot(yk) 
  
   else if (idlink(yk).eq.2) then ! Splines link
       k=0
        zitr=0.d0
        zitr(1:nbzitr(yk))=zitr0(1:nbzitr(yk),yk)
        zitr(-1)=zitr(1)
        zitr(0)=zitr(1)
        zitr(ntrtot(yk)-1)=zitr(ntrtot(yk)-2)
        zitr(ntrtot(yk))=zitr(ntrtot(yk)-1)

        bb=b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)
!        if(verbose==1) print*,"bb=",bb,"indice",nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1
        splaa=0.d0
        do kk=2,ntrtot(yk)
           splaa(kk-3)=b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+kk)**2
        end do 

      ! if(verbose==1) print*,"yk=",yk,"zitr=",zitr,"splaa=",splaa
       
       
      do j=1,maxmes
        do l=1,nsim 
         ier=0
         niter=0
         diff=0.d0
         ytemp = mu(maxmes*(yk-1)+j)+sqrt(VC(maxmes*(yk-1)+j,maxmes*(yk-1)+j))*gauss(1,l) 
         ytemp2=INV_ISPLINES(ytemp,splaa,bb,nbzitr(yk),zitr,ier,niter,diff)
         
        if ((ier.eq.3).or.(ier.ne.1 .and. diff.gt.1.d-3).or.ymarg(maxmes*(yk-1)+j,1).eq.9999.d0) then
                ymarg(maxmes*(yk-1)+j,1)=9999.d0
   !             if (verbose==1 .and. k==0) print*,"pb splines :","ier=",ier,"diff=",diff,"ymarg=",&
   !             ymarg(maxmes*(yk-1)+j,1),"ytemp=",ytemp,"VC=",VC(maxmes*(yk-1)+j,maxmes*(yk-1)+j)
 !               k=k+1
!print*,"ymarg",ymarg(maxmes*(yk-1)+j,1)
        else
                ymarg(maxmes*(yk-1)+j,1)=ymarg(maxmes*(yk-1)+j,1)+ytemp2*gauss(2,l)
        end if
 
        end do
      end do 
      
      !     if(verbose==1.and. k>0) print*,k,"fois pb splines"
       
       sumntrtot = sumntrtot + ntrtot(yk)
    end if
  end do

end if
  ! fin if methInteg
  !fin cas ng=1

 else ! cas ng>1
     ! transformation des  pig=exp(Xbg)/(1+somme(Xbk,k=1,G-1))
!     Xprob=0.d0
!     Xprob(1)=1
!     l=0
!     do k=1,nv
!        if (idprob(k).eq.1) then
!           l=l+1
!           Xprob(1+l)=X0(maxmes*(k-1)+1)
!        end if
!     end do
!     pi=0.d0
!     temp=0.d0
!     Do g=1,ng-1
!        bprob=0.d0
!        do k=1,nvarprob
!           bprob(k)=b1((k-1)*(ng-1)+g)
!        end do

!        temp=temp+exp(DOT_PRODUCT(bprob,Xprob))

!        pi(g)=exp(DOT_PRODUCT(bprob,Xprob))

!     end do

!     pi(ng)=1/(1+temp)

!     do g=1,ng-1
!        pi(g)=pi(g)*pi(ng)
!     end do

! creation des vecteurs de variables explicatives
            l=0
            m=0
            q=0
            X00=0.d0
            X2=0.d0
            X01=0.d0
            b01=0.d0
            do k=1,nv
               if (idg(k).eq.2) then
                  l=l+1
                  do yk=1,ny
                  do j=1,maxmes
                        X2(maxmes*(yk-1)+j,l)=dble(X0(maxmes*(k-1)+j))
                  end do
                  end do
               else if (idg(k).eq.1) then
                  m=m+1
                  do yk=1,ny
                  do j=1,maxmes
                     X00(maxmes*(yk-1)+j,m)=dble(X0(maxmes*(k-1)+j))
                  end do
                  end do
               end if
               
               !contrast : 
               if (idcontr(k).ne.0) then
                  q=q+1
                  do yk=1,ny
! creation matrice design des contrastes: X01
                     do j=1,maxmes
                        X01(maxmes*(yk-1)+j,(q-1)*ny+yk) = dble(X0(maxmes*(k-1)+j))
                     end do
! creation vecteur parms des contrastes: b01
                     if (yk<ny) THEN
                        b01((q-1)*ny+yk)=b1(nef-ncontr+(q-1)*(ny-1)+yk)
                     else
                        b01((q-1)*ny+ny) =-sum(b1((nef-ncontr+(q-1)*(ny-1)+1) &
                             :(nef-ncontr+(q-1)*(ny-1)+ny-1)))
                     end if
                  end do                  
               end if          
            end do
!             print*,"contrast ok"
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

!           print*,"mixture ok"

        ! variance covariance si spec aux classes :

               Ut1=Ut
        if (nwg.ne.0) then
           Ut1=0.d0
           if (g.eq.ng) then
              Ut1=Ut
           else
              Ut1=Ut*abs(b1(nef+nvc+g))
           end if

           if(all(idlink.ne.3)) then 
           P=0.d0
           P=MATMUL(Z,Ut1)
           VC=0.d0
           VC=MATMUL(P,transpose(P))+Corr
!             print*,"VC ok"
           ! Vi en vecteur
           Vi=0.d0
           jj=0
           do j=1,maxmes*ny
              do k=j,maxmes*ny
                 jj=j+k*(k-1)/2
                 Vi(jj)=VC(j,k)
              end do
           end do

           if (methInteg.eq.1) then 
              CALL DMFSD(Vi,maxmes*ny,EPS,IER)
              if (ier.eq.-1) then
                 ymarg=9999.d0
     !            if (verbose==1) print*,"pb cholesky"
                 goto 654
              end if
              
              VC=0.d0
              do j=1,maxmes*ny
                 do k=1,j
                    VC(j,k)=Vi(k+j*(j-1)/2)
                 end do
              end do
           end if
        end if
     end if

        mu=0.d0
        mu=matmul(X00,b0)+matmul(X2,b2)+matmul(X01,b01)
!              print*,"mu ok"
        !if (verbose==1) print*,"mu=",mu        
!        write(*,*)'mu'

    if(methInteg.eq.1) then !MC
!    print*,"debut MC"
      do l=1,nsim
         
         if(any(idlink.ne.0)) then

            if(any(idlink.eq.3)) then

               !! simuler les effets aleatoires
               usim2=0.d0
               ysim=0.d0
               do m=1,nea
                  SX=1.d0
                  call bgos(SX,0,usim2(m),x22,0.d0)
               end do
               uu = 0.d0
               uu = matmul(Ut1,usim2)
               ysim = mu+matmul(Z,uu)
               
               !! simuler le BM ou AR
               if(ncor.gt.0) then
                  wsim=0.d0
                  x22=0.d0
                  SX=1.d0
                  do j=1,maxmes*ny
                     call bgos(SX,0,wsim(j),x22,0.d0)
                  end do
                  wi=0.d0
                  wi=matmul(Corr,wsim)
                  ysim = ysim + wi
               end if

            else
               
               usim=0.d0
               ysim=0.d0
               do m=1,maxmes*ny
                  SX=1.d0
                  call bgos(SX,0,usim(m),x22,0.d0)
               end do
               ysim=mu+MATMUL(VC,usim)
               
            end if
         end if

!    if (verbose==1) print*,"ysim=",ysim
      sumntrtot=0
      sumnbmod = 0
      do yk=1,ny
!          print*,"boucle idlink, yk et nsim=",yk,nsim 
         if (idlink(yk).eq.0 .and. l.eq.1) then  ! Linear link
            
            aa = b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)
            bb = b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2)
            do j=1,maxmes
               Ymarg(maxmes*(yk-1)+j,g) = mu(maxmes*(yk-1)+j)*bb+aa
            end do     
            
         else if (idlink(yk).eq.3) then ! ordinal

           !! simuler l'EA specifique au test
           ai=0.d0
           if(nalea.gt.0) then
              call bgos(SX,0,asim(yk),x22,0.d0)
              ai = b1(nprob+nef+ncontr+nvc+nwg+ncor+ny+yk)*asim(yk)
           end if
           
           aa = b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)
           do j=1,maxmes
              if(nalea.gt.0) ysim((yk-1)*maxmes+j) = ysim((yk-1)*maxmes+j)+ai
              ytemp = modalite(sumnbmod + nbmod(yk)) - &
                   alnorm((aa-ysim((yk-1)*maxmes+j))/abs(b1(nef+nvc+nwg+ncor+yk)),.false.)
              Ymarg(maxmes*(yk-1)+j,g) = Ymarg(maxmes*(yk-1)+j,g) + ytemp / dble(nsim)
           end do
           
           if(nbmod(yk).gt.2) then
              do jj=2,nbmod(yk)-1
                 aa = aa + b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+jj)**2
                 do j=1,maxmes
                    ytemp = (modalite(sumnbmod + jj) - modalite(sumnbmod + jj+1))* &
                         alnorm((aa-ysim((yk-1)*maxmes+j))/abs(b1(nef+nvc+nwg+ncor+yk)),.false.)
                    Ymarg(maxmes*(yk-1)+j,g) = Ymarg(maxmes*(yk-1)+j,g) + ytemp / dble(nsim)                             
                 end do
              end do
           end if
           
           sumnbmod = sumnbmod + nbmod(yk)
                
        else if (idlink(yk).eq.1) then  ! Beta link
              
           aa1=exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1))/ &
                   (1+exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)))
              bb1=exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2))/ &
                   (1+exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2)))
              bb1=aa1*(1.d0-aa1)*bb1
              
              cc1=abs(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+3))
              
              dd1=abs(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+4))
              
              aa=aa1*aa1*(1-aa1)/bb1-aa1
              bb=aa*(1-aa1)/aa1
              beta=beta_ln(aa,bb)
              
              !! simuler une N(0,Sigma)
              if(any(idlink.eq.3)) then
                 usim=0.d0
                 do m=1,maxmes
                    SX=1.d0
                    call bgos(SX,0,usim(m),x22,0.d0)
                 end do
                 
                 Sigma=0.d0
                 do j1=1,maxmes
                    do j2=1,j1
                       Sigma(j1,j2) = VC((yk-1)*maxmes+j1,(yk-1)*maxmes+j2)
                    end do
                 end do
                 
                 Vi=0.d0
                 Vi = matmul(Sigma,usim)
                 
                 do j=1,maxmes
                    ysim(maxmes*(yk-1)+j) = ysim(maxmes*(yk-1)+j) + Vi(j)
                 end do
              end if
              
              do j=1,maxmes
                 ytemp=ysim(maxmes*(yk-1)+j)*dd1+cc1
                 if (ytemp.lt.0) then
                    ytemp=0.d0
                 end if
                 if (ytemp.gt.1) then
                    ytemp=1.d0
                 end if
                 ymarg(maxmes*(yk-1)+j,g)=ymarg(maxmes*(yk-1)+j,g)+xinbta(aa,bb,beta,ytemp,ier)/dble(nsim)
                 if (ier.ne.0 .or. ymarg(maxmes*(yk-1)+j,g).eq.9999.d0) then
                    ymarg(maxmes*(yk-1)+j,g)=9999.d0
!                    if (verbose==1) print*,"pb beta"
                 end if
              end do

      else if (idlink(yk).eq.2) then ! Splines link
         
        zitr=0.d0
        zitr(1:nbzitr(yk))=zitr0(1:nbzitr(yk),yk)
        zitr(-1)=zitr(1)
        zitr(0)=zitr(1)
        zitr(ntrtot(yk)-1)=zitr(ntrtot(yk)-2)
        zitr(ntrtot(yk))=zitr(ntrtot(yk)-1)

        bb=b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)
        do kk=2,ntrtot(yk)
           splaa(kk-3)=b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+kk)**2
        end do
  
        !! simuler une N(0,Sigma)
        if(any(idlink.eq.3)) then
           usim=0.d0
           do m=1,maxmes
              SX=1.d0
              call bgos(SX,0,usim(m),x22,0.d0)
           end do
           
           Sigma=0.d0
           do j1=1,maxmes
              do j2=1,j1
                 Sigma(j1,j2) = VC((yk-1)*maxmes+j1,(yk-1)*maxmes+j2)
              end do
           end do
           
           Vi=0.d0
           Vi = matmul(Sigma,usim)
           
           do j=1,maxmes
              ysim(maxmes*(yk-1)+j) = ysim(maxmes*(yk-1)+j) + Vi(j)
           end do
        end if
        
        do j=1,maxmes
           niter=0
           diff=0.d0
           ier=0
           ytemp=INV_ISPLINES(ysim(maxmes*(yk-1)+j),splaa,bb,nbzitr(yk),zitr,ier,niter,diff)
           if ((ier.eq.3).or.(ier.ne.1 .and. diff.gt.1.d-3).or.ymarg(maxmes*(yk-1)+j,g).eq.9999.d0) then
              ymarg(maxmes*(yk-1)+j,g)=9999.d0
              !                       if (verbose==1) print*,"pb splines"
           else
              ymarg(maxmes*(yk-1)+j,g)=ymarg(maxmes*(yk-1)+j,g)+ytemp/dble(nsim)
           end if
           
        end do
        
     end if
  ! fin if idlink
        
     sumntrtot = sumntrtot+ntrtot(yk)     
  
  end do
end do
  
  else !GH
!       print*,"debut GH"
    call gausshermite(gauss,nsim)

   sumntrtot=0
   do yk=1,ny
!   print*,"boucle idlink, yk=",yk
     if (idlink(yk).eq.0) then  ! Linear link
         aa = b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)
         bb = b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2)
         do j=1,maxmes
           Ymarg(maxmes*(yk-1)+j,g) = mu(maxmes*(yk-1)+j)*bb+aa
         end do
         sumntrtot = sumntrtot +ntrtot(yk) 
         
     else if (idlink(yk).eq.1) then  ! Beta link
     
        aa1=exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1))/ &
             (1+exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)))
        bb1=exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2))/ &
             (1+exp(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+2)))
        bb1=aa1*(1.d0-aa1)*bb1

        cc1=abs(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+3))

        dd1=abs(b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+4))

        aa=aa1*aa1*(1-aa1)/bb1-aa1
        bb=aa*(1-aa1)/aa1
        beta=beta_ln(aa,bb)     
      
     
      do j=1,maxmes
        do l=1,nsim 
                ytemp=mu(maxmes*(yk-1)+j)+sqrt(VC(maxmes*(yk-1)+j,maxmes*(yk-1)+j))*gauss(1,l)
                ytemp=ytemp*dd1+cc1

                 if (ytemp.lt.0) then
                    ytemp=0.d0
                 end if
                 if (ytemp.gt.1) then
                    ytemp=1.d0
                 end if
                 ymarg(maxmes*(yk-1)+j,g)=ymarg(maxmes*(yk-1)+j,g)+xinbta(aa,bb,beta,ytemp,ier)*gauss(2,l)
                if (ier.ne.0 .or. ymarg(maxmes*(yk-1)+j,g).eq.9999.d0) then
                    ymarg(maxmes*(yk-1)+j,g)=9999.d0
!                    if (verbose==1) print*,"pb beta"
                 end if
        end do
      end do
      sumntrtot = sumntrtot +ntrtot(yk) 
  
   else if (idlink(yk).eq.2) then ! Splines link
     
        zitr=0.d0
        zitr(1:nbzitr(yk))=zitr0(1:nbzitr(yk),yk)
        zitr(-1)=zitr(1)
        zitr(0)=zitr(1)
        zitr(ntrtot(yk)-1)=zitr(ntrtot(yk)-2)
        zitr(ntrtot(yk))=zitr(ntrtot(yk)-1)

        bb=b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+1)
        do kk=2,ntrtot(yk)
           splaa(kk-3)=b1(nef+nvc+nwg+ncor+ny+nalea+sumntrtot+kk)**2
        end do 

      do j=1,maxmes
        do l=1,nsim 
         ier=0
         niter=0
         diff=0.d0
         ytemp = mu(maxmes*(yk-1)+j)+sqrt(VC(maxmes*(yk-1)+j,maxmes*(yk-1)+j))*gauss(1,l)
         ytemp2=INV_ISPLINES(ytemp,splaa,bb,nbzitr(yk),zitr,ier,niter,diff)
        if ((ier.eq.3).or.(ier.ne.1.and.diff.gt.1.d-3).or.ymarg(maxmes*(yk-1)+j,g).eq.9999.d0) then
                ymarg(maxmes*(yk-1)+j,g)=9999.d0
!                if (verbose==1) print*,"pb splines"
        else
                ymarg(maxmes*(yk-1)+j,g)=ymarg(maxmes*(yk-1)+j,g)+ytemp2*gauss(2,l)
        end if  

        end do
      end do 
      sumntrtot = sumntrtot +ntrtot(yk) 
  
    end if
  end do
  
 end if
 !fin if methInteg
 
 end do 
 ! fin do ng
 
 end if
 ! fin if ng
!     print*,"fin boucles"
 
do yk=1,ny

  if (idlink(yk).eq.1) then
   do g=1,ng
     do j=1,maxmes
        if (ymarg(maxmes*(yk-1)+j,g).ne.9999.d0) then
        ymarg(maxmes*(yk-1)+j,g)=ymarg(maxmes*(yk-1)+j,g)*(maxY(yk)-minY(yk)+2*epsY(yk))+minY(yk)-epsY(yk)
       if (ymarg(maxmes*(yk-1)+j,g).lt.minY(yk)) ymarg(maxmes*(yk-1)+j,g)=minY(yk)
        if (ymarg(maxmes*(yk-1)+j,g).gt.maxY(yk)) ymarg(maxmes*(yk-1)+j,g)=maxY(yk)
        end if
     end do
   end do
 
 else
   do g=1,ng
      do j=1,maxmes
        if (ymarg(maxmes*(yk-1)+j,g).ne.9999.d0) then
        ymarg(maxmes*(yk-1)+j,g)=ymarg(maxmes*(yk-1)+j,g)
        end if
     end do
   end do

 end if
end do

654 continue
!        print*,"avant delocate"
  deallocate(zitr,splaa)
  deallocate(Z,P,Corr,X00,X2,X01,ysim,usim,mu,tcor,VC,Vi,usim2,Sigma,wsim,asim)


  return

end subroutine predictmult

