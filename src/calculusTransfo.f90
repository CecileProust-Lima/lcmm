

subroutine calculustransfo(b,npm,ny,idlink,ntrtot,imoins,zitr0,maxnbzitr,marker,nsim,minY,maxY,epsY,ide,dimide,transfY)

  implicit none

  integer::kk,nsim,npm,j,k,sumntrtot,yk,imoins,ny,maxnbzitr,l
  double precision,dimension(ny*nsim)::marker,transfY
  integer, dimension(ny)::idlink,ntrtot,dimide
  integer,dimension(maxval(dimide),ny)::ide
  double precision,dimension(maxval(ntrtot))::splaa
  double precision::aa1,bb1,dd1,aa,bb,betai,eps,ytemp,cc1
  double precision, dimension(npm)::b
  double precision,dimension(ny)::minY,maxY,epsY
  double precision,dimension(maxnbzitr,ny)::zitr0
  double precision,dimension(:),allocatable::zitr
  double precision,dimension(maxval(ntrtot))::Xspl
  double precision,dimension(nsim)::mmm,mmm1,mmm2,iim,iim1,iim2
  double precision ::ht,htm,ht2,ht3,hht,h,hh,h2,h3,h2n,hn,inf,sup
  double precision, dimension(npm-minval(ntrtot)+maxval(dimide))::b1

  eps=1.D-20

  !      write(*,*)'infos',minY,maxY,nsim,npm,ntrtot
  !       write(*,*)'b',(b(j),j=1,npm)

  transfY=0.d0

  !      write(*,*)(marker(j),j=1,nsim)

  sumntrtot = 0

  if (any(idlink.eq.2)) then
     allocate(zitr(-1:maxval(ntrtot)))
  else
     allocate(zitr(1))
  end if

  zitr=0.d0
  mmm=0.d0
  mmm1=0.d0
  mmm2=0.d0
  iim=0.d0
  iim1=0.d0
  iim2=0.d0
  Xspl=0.d0

  b1=0.d0
  inf=0.d0
  sup=0.d0

  do yk=1,ny

     if (idlink(yk).eq.2) then

        zitr=0.d0
        mmm=0.d0
        mmm1=0.d0
        mmm2=0.d0
        iim=0.d0
        iim1=0.d0
        iim2=0.d0
        Xspl=0.d0
        splaa=0.d0

        zitr(1:(ntrtot(yk)-2))=zitr0(1:(ntrtot(yk)-2),yk)
        zitr(-1)=zitr(1)
        zitr(0)=zitr(1)
        zitr(ntrtot(yk)-1)=zitr(ntrtot(yk)-2)
        zitr(ntrtot(yk))=zitr(ntrtot(yk)-1)


        splaa(1)=b(imoins+sumntrtot+1)
        do kk=2,ntrtot(yk)
           splaa(kk)=b(imoins+sumntrtot+kk)*b(imoins+sumntrtot+kk)
        end do

        !calcul de H(y)
        do j=1,nsim
           ! ou se trouve la valeur
           l=0

           do k = 2,ntrtot(yk)-2
              if ((marker((yk-1)*nsim+j).ge.zitr(k-1)).and.(marker((yk-1)*nsim+j).lt.zitr(k))) then
                 l=k-1
              end if
           end do

           if (marker((yk-1)*nsim+j).eq.zitr(ntrtot(yk)-2)) then
              l=ntrtot(yk)-3
           end if

           !         if (l.lt.1.or.l.gt.ntrtot-1) then
           !            write(*,*)'probleme estim splines',l
           !            write(*,*)'j=',j,'test(j)',test(j)
           !            stop
           !         end if


           ht2 = zitr(l+1)-marker((yk-1)*nsim+j)
           htm= marker((yk-1)*nsim+j)-zitr(l-1)
           ht = marker((yk-1)*nsim+j)-zitr(l)
           ht3 = zitr(l+2)-marker((yk-1)*nsim+j)
           hht = marker((yk-1)*nsim+j)-zitr(l-2)
           h = zitr(l+1)-zitr(l)
           hh= zitr(l+1)-zitr(l-1)
           hn= zitr(l+1)-zitr(l-2)
           h2n=zitr(l+2)-zitr(l-1)
           h2= zitr(l+2)-zitr(l)
           h3= zitr(l+3)-zitr(l)

           if (marker((yk-1)*nsim+j).ne.zitr(ntrtot(yk)-2)) then
              mmm2(j) = (3.d0*ht2*ht2)/(hh*h*hn)
              mmm1(j) = (3.d0*htm*ht2)/(h2n*hh*h)+(3.d0*ht*ht3)/(h2*h*h2n)
              mmm(j)  = (3.d0*ht*ht)/(h3*h2*h)
           end if
           if (marker((yk-1)*nsim+j).eq.zitr(ntrtot(yk)-2)) then
              mmm2(j) = 0.d0
              mmm1(j) = 0.d0
              mmm(j)  = 3.d0/h
           end if

           iim2(j)=hht*mmm2(j)/(3.d0)+ h2n*mmm1(j)/(3.d0) &
                +h3*mmm(j)/(3.d0)
           iim1(j)=htm*mmm1(j)/(3.d0)+h3*mmm(j)/(3.d0)

           iim(j)=ht*mmm(j)/(3.d0)

           !-------- transformation  :

           Xspl=0.d0
           Xspl(1)=1
           do k=2,l
              Xspl(k)=1
           end do
           Xspl(l+1)=iim2(j)
           Xspl(l+2)=iim1(j)
           Xspl(l+3)=iim(j)
           transfY((yk-1)*nsim+j)= dot_product(Xspl,splaa)
        end do
        !fin H(y)




     else if (idlink(yk).eq.1) then

        aa1=exp(b(imoins+sumntrtot+1))/ &
             (1+exp(b(imoins+sumntrtot+1)))
        bb1=exp(b(imoins+sumntrtot+2))/ &
             (1+exp(b(imoins+sumntrtot+2)))
        bb1=aa1*(1.d0-aa1)*bb1
        cc1=b(imoins+sumntrtot+3)
        dd1=abs(b(imoins+sumntrtot+4))

        aa=aa1*aa1*(1-aa1)/bb1-aa1
        bb=aa*(1-aa1)/aa1

        do j=1,nsim
           ytemp=(marker((yk-1)*nsim+j)-minY(yk)+epsY(yk))/(maxY(yk)-minY(yk)+2*epsY(yk))
           transfY((yk-1)*nsim+j)=(betai(aa,bb,ytemp)-cc1)/dd1
           !                if (transfY((yk-1)*nsim+j).eq.999.d0) then
           !                    write(*,*)'problem'
           !                end if

        end do


     else if (idlink(yk).eq.0) then

        do j=1,nsim
           transfY((yk-1)*nsim+j)=(marker((yk-1)*nsim+j)-b(imoins+sumntrtot+1)) &
                /abs(b(imoins+sumntrtot+2))
        end do

     else if (idlink(yk).eq.3) then

        ! On suppose ici qu'on a dans marker les memes valeurs que dans estimlink (ie toutes les modalites du min au max)

        b1=0.d0
        do k=1,imoins
           b1(k)=b(k)
        end do

        l=imoins+sumntrtot
        do k=1,INT(maxY(yk))-INT(minY(yk))
           if (ide(k,yk).eq.1) then
              l=l+1
              b1(imoins+sumntrtot+k)=b(l)
           else
              b1(imoins+sumntrtot+k)=0.d0
           end if
        end do

        !       write(*,*)'infos',minY,maxY,nsim,npm
        !       write(*,*)'b',(b1(j),j=1,npm)


        l=imoins+sumntrtot
        transfY((yk-1)*nsim+1)=-1.d10
        transfY((yk-1)*nsim+2)=b1(l+1)
        sup=b1(l+1)
        if (maxY(yk)-minY(yk).gt.1) then
           do k=1,INT(maxY(yk))-INT(minY(yk))-1
              inf=sup
              sup=sup+b1(l+1+k)*b1(l+1+k)
              transfY((yk-1)*nsim+2*k+1)=inf
              transfY((yk-1)*nsim+2*k+2)=sup
           end do
        end if
        transfY((yk-1)*nsim+2*(dimide(yk))+1)=sup
        transfY((yk-1)*nsim+2*(dimide(yk))+2)=1.d10
        if(nsim.gt.2*dimide(yk)+2) then
           do k = 2*dimide(yk)+3, nsim
              transfY((yk-1)*nsim+k) = 1.d10
           end do
        end if

     end if
     sumntrtot = sumntrtot + ntrtot(yk)
  end do
  !     write(*,*)(marker(j),j=1,ny*nsim)
  !     write(*,*)(transfY(j),j=1,ny*nsim)

  deallocate(zitr)

end subroutine calculustransfo
