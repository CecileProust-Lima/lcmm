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


!===================================================================
!
!
      !  MODULE LOIS_NORMALES
!
!
!===================================================================





        module lois_normalesc

        contains


!===================================================================
!
!
!  SUbroutine de Genz
!
!
!===================================================================


        subroutine hrmsym(ndim,nf,minpts,maxpts,funsub,epsabs,epsrel,&
        restar,result,abserr,neval,ifail,work)


!!!!BEGIN PROLOGUE HRMSYM
!
!***AUTHOR
!            Alan Genz,
!            Department of Mathematics
!            Washington State University
!            Pullman, WA 99164-3113, USA
!            Email: alangenz@wsu.edu
!
!      Reference: Genz, A., and Keister, B. (1996), Fully Symmetric
!      Interpolatory Rules for Multiple Integrals over Infinite
!      Regions with Gaussian Weight, J. Comp. Appl. Math., 71, 299-309.
!
!***KEYWORDS automatic multidimensional integrator,
!            n-dimensional region (-infin, infin)^n
!           Gaussian weight
!***PURPOSE  The routine calculates an approximation to a given
!            vector of definite integrals
!
!      infin     infin
!     I    ...  I     w(X)(F ,F ,...,F  ) DX(NDIM)...DX(2)DX(1),
!     -infin    -infin      1  2      NF
!
!       where F = F (X ,X ,...,X    ), I = 1,2,...,NF,
!              I   I  1  2      NDIM
!
!       w(X) = EXP(-( X(1)**2 + ... + X(NDIM)**2 )/2)/SQRT(2*PI)**NDIM.
!
!***DESCRIPTION Computation of integrals over infinite regions with
!               Gaussian weight.
!
!   ON ENTRY
!
!     NDIM   Integer, number of variables, with 1 < NDIM <= 1000.
!     NF     Integer, number of components of the integral.
!     MINPTS Integer, minimum number of function evaluations.
!     MAXPTS Integer, maximum number of function evaluations.
!     FUNSUB Externally declared subroutine for computing
!            all components of the integrand at the given
!            evaluation point.
!            It must have parameters (NDIM,X,NF,FUNVLS)
!            Input parameters:
!              NDIM   Integer that defines the dimension of the
!                     integral.
!              X      Real array of dimension NDIM
!                     that defines the evaluation point.
!              NF Integer that defines the number of components of I.
!            Output parameter:
!              FUNVLS Real array of dimension NF
!                     that defines NF components of the integrand.
!
!     EPSABS Real, requested absolute accuracy.
!     EPSREL Real, requested relative accuracy.
!     RESTAR Integer.
!            If RESTAR = 0, this is the first attempt to compute
!            the integral.
!            If RESTAR = 1, then we restart a previous attempt.
!     WORK   Real array of working storage, must have dimension at
!            least 2*NF + (NF+1)*NUMSMS. NUMSMS is the number
!            of NDIM partitions of the integers 0,...,RULE, after
!            RULE+1 approximations to the integrals have been computed.
!            The required value for NUMSMS will never exceed 10000.
!
!   ON RETURN
!
!     RESULT Real array of dimension NF.
!            Approximations to all components of the integral.
!     ABSERR Real array of dimension NF.
!            Estimates of absolute accuracies.
!     NEVAL  Integer, number of function evaluations used.
!     IFAIL  Integer.
!            IFAIL = 0 for normal exit, when
!              ABSERR(K) <= MAX( EPSABS, ABS(RESULT(K))*EPSREL ) for
!              all K, 0 < K <= NF, with <= MAXPTS function values.
!            IFAIL = 1 if MAXPTS was too small to obtain the required
!              accuracy. In this case values of RESULT are returned
!              with estimated absolute accuracies ABSERR.
!     WORK   Real array of working storage, that contains information
!            that might be needed for subsequent calls of HRMSYM. This
!            array should not be modified.
!
!***ROUTINES CALLED
!***END PROLOGUE HRMSYM
!
!   Global variables.
!
        implicit none

        external::funsub
        integer::ndim,nf,minpts,maxpts,restar,neval,ifail
        double precision::epsabs,epsrel
        double precision,dimension(*)::result,abserr,work
!
!   Local variables.
!
        integer::i,intcls,numsms
        integer,save::rule,mnrule
        integer,parameter::maxrul=25
        double precision::differ,rntcls
        double precision,save::cltotl
        double precision,parameter::epower=1.5

        ifail=1
        if (ndim .le. 4) then
!
!     Call product Gauss-Hermite Rule
!
                call hermit(ndim,nf,minpts,maxpts,funsub,epsabs,&
                 epsrel,restar,result,abserr,neval,ifail,work)

        else

                if (restar.eq.0) then

                        mnrule = -1
                        rule = 0
                        cltotl = 0

                        do i=1,nf
                                work(i) = 0
                                work(i+nf) = 0
                        end do

                end if

                neval = 0

                do while (neval .le. maxpts .and. rule .le. maxrul .and. (ifail.gt.0 .or. neval.lt.minpts))

                        call hrmtrl(ndim,nf,funsub,mnrule,rule,result,intcls,&
                        work(2*nf+1),work(3*nf+1),numsms)

                        rntcls = real(intcls)**epower
                        cltotl = cltotl + rntcls

                        do i=1,nf
                                differ = (result(i) - work(i))/cltotl
                                work(i) = work(i) + rntcls*differ
                                work(nf+i) = (cltotl - rntcls)*(work(nf+i)/cltotl + rntcls*differ**2)
                        end do

                        ifail = 0

                        do i= 1,nf

                                if (rule .gt. 0) then
                                        abserr(i) = sqrt(work(nf+i))
                                else
                                        abserr(i) = abs(result(i))
                                end if

                                if (abserr(i).gt.max(epsabs,epsrel*abs(result(i)))) then
                                        ifail = 1
                                end if
                        end do

                        neval = neval + intcls
                        rule = rule + 1
                end do
        end if
!
!***END HRMSYM
!
        end subroutine hrmsym







        subroutine hrmtrl(s,n,f,minord,maxord,intval,intcls,work,fulsms,numsms)
!  MULTIDIMENSIONAL FULLY SYMMETRIC RULE INTEGRATION SUBROUTINE
!
!   THIS SUBROUTINE COMPUTES A SEQUENCE OF FULLY SYMMETRIC RULE
!   APPROXIMATIONS TO A FULLY SYMMETRIC MULTIPLE INTEGRAL. WRITTEN BY
!     Alan Genz
!     Department of Mathematics
!     Washington State University
!     Pullman, Washington 99164-3113  USA
!     Telephone: 509-335-2131
!     Electronic Mail: alangenz@wsu.edu
!     Fax: 509-335-1188
!
!**************  PARAMETERS FOR HRMTRL  ********************************
!*****INPUT PARAMETERS
!  S       INTEGER number of variables, with 0 < S <= 1000.
!  N       Integer, number of components of the integral.
!  F       EXTERNALly declared user defined integrand subroutine.
!          it must have parameters (S,X,N,F), where X is a REAL S-array
!          and F is a real N-array.
!  MINORD  INTEGER minimum order parameter.  On entry MINORD specifies
!          the current highest order approximation to the integral,
!          available in the array INTVAL.  For the first call
!          MINORD should be set to -1. Otherwise a previous call is
!          assumed that computed INTVAL. On exit MINORD is set to MAXORD.
!  MAXORD  INTEGER maximum order parameter, must be greater than MINORD
!          and not exceed 25. The subroutine computes approximations of
!          polynomial degree 2*MAXORD+1.
!******OUTPUT PARAMETERS
!  INTVAL  REAL array of length N.  Upon successful exit
!          INTVAL(1),..., INTVAL(N) are approximations to the components
!          of the integral.  These are all approximations of polynomial
!          degree 2*MAXORD+1.
!  INTCLS  INTEGER number of F values needed for INTVAL
!  WORK    REAL working storage array.
!  FULSMS  REAL working storage array with dimension (N+1,*). On exit
!          FULSMS(I,J) contains the fully symmetric basic rule sum
!          indexed by the jth s-partition of the integers 0,...,MAXORD,
!          for the Ith component of the integrand.
!          FULSMS(N+1,J) contains number of points for the fully
!          symmetric basic rule sum indexed by the Jth S-partition of
!          the integers 0,...,MAXORD.
!  NUMSMS  INTEGER number of S-partitions of the integers 0,...,MAXORD.
!***********************************************************************
        implicit none

        external::f
        integer::s,n,minord,maxord,numsms,d,i,modofm,intcls,prtcnt
        integer,parameter::maxdim=1000,maxrdm=25
        double precision::fulwgt
        double precision,dimension(maxdim)::x
        double precision,dimension(n)::intval
        double precision,dimension(n+1,*)::fulsms
        double precision,dimension(*)::work
        integer,dimension(maxdim)::m,k

        d= minord + 1
        intcls = 0


        if (d.eq.0) then
                do i=1,n
                        intval(i)=0.d0
                end do
        end if
!
!***  Begin loop for each D
!      for each D find all distinct partitions M with |M| <= D
!
        do while (d.le.min(maxord,maxrdm))

                prtcnt=0
                call nxpart(prtcnt,s,m,modofm)

                do while (modofm.le.d)
!
!***  Calculate updated weight for partition M alnd
!***     fully symmetric sums ( when necessary )
!
                        fulwgt = fweight(s,x,m,k,modofm,d)

                        if (d.eq.modofm) then
                        do i=1,n
                                fulsms(i,prtcnt)=0.d0
                        end do
                        fulsms(n+1,prtcnt)=0.d0
                        end if

                        if (fulsms(n+1,prtcnt).eq.0 .and.fulwgt .ne.0) then
                                call fulsmh(s,m,n,f,fulsms(1,prtcnt),x,work)
                                intcls=intcls+INT(fulsms(n+1,prtcnt))
                        end if

                        do i=1,n
                                intval(i)=intval(i)+fulwgt*fulsms(i,prtcnt)
                        end do

                        call nxpart(prtcnt,s,m,modofm)

                end do
!
!***  End loop for each D
!
                d=d+1
        end do


        minord = maxord
        numsms =prtcnt - 1

        end subroutine hrmtrl






!##############################################################################################






!##########################################################################################

        double precision function fweight(s,intrps,m,k,modofm,d)

        implicit none
!
!***  Subroutine to update weight for partition m
!

        integer::s,i,l,d,modofm
        integer,parameter::nzrmax=17
        integer,parameter::maxrdm=25
        integer,dimension(s)::m,k
        double precision::momnkn
        double precision,dimension(s)::intrps
        double precision,dimension(0:maxrdm,0:maxrdm),save::momprd
        double precision,dimension(0:nzrmax)::a,g


!
!     Modified moments and generators for 1 + 2 + 6 + 10 + 16 = 35 point
!        degree 51 rule, with degree 1, 5 and 19 imbedded rules.
!

        DATA A /  2*1D0, 0D0, 6D0,-0.48378475125832451D2, 3*0D0, 34020D0, &
        -0.98606453173677489D6, 5*0D0, 0.12912054173706603D13, &
        -0.11268664521456168D15, 0.29248520348796280D16 /
        DATA G( 0),G( 1)/0D0                  , 0.17320508075688773D1/
        DATA G( 2),G( 3)/0.41849560176727319D1, 0.74109534999454084D0/
        DATA G( 4),G( 5)/0.28612795760570581D1, 0.63633944943363700D1/
        DATA G( 6),G( 7)/0.12304236340273060D1, 0.51870160399136561D1/
        DATA G( 8),G( 9)/0.25960831150492022D1, 0.32053337944991945D1/
        DATA G(10),G(11)/0.90169397898903025D1, 0.24899229757996061D0/
        DATA G(12),G(13)/0.79807717985905609D1, 0.22336260616769417D1/
        DATA G(14),G(15)/0.71221067008046167D1, 0.36353185190372782D1/
        DATA G(16),G(17)/0.56981777684881096D1, 0.47364330859522971D1/
        DATA MOMPRD(0,0) / 0D0 /


        if (momprd(0,0).eq. 0) then
!
!***  Calculate moments
!
                do l = 0,maxrdm
                        do i= 0,maxrdm
                        momprd(l,i) = 0
                        end do
                end do

                momprd(0,0) = a(0)

                do l = 0,nzrmax
                        momnkn = 1
                        do i = 1,nzrmax
                                if (i .le. l) then
                                        momnkn = momnkn*(g(l)**2 - g(i-1)**2)
                                else
                                        momnkn = momnkn*(g(l)**2 - g(i)**2)
                                end if
                                if (i .ge. l) momprd(l,i) = a(i)/momnkn
                        end do
                end do
        end if
!
!     Determine Updated Weight Contribution
!
        do i = 2,s
                intrps(i) = 0
                k(i) = m(i)
        end do

        k(1) = d - modofm + m(1)
 10   intrps(1) = momprd(m(1),k(1))

        do i = 2,s
                intrps(i) = intrps(I) + momprd(m(i),k(i))*intrps(i-1)
                intrps(i-1) = 0
                k(1) = k(1) - 1
                k(i) = k(i) + 1
                if (k(1) .ge. m(1)) goto 10
                k(1) = k(1) + k(i) - m(i)
                k(i) = m(i)
        end do

        fweight = intrps(s)

        end function fweight



        !###########################################################################################








        subroutine fulsmh(s,m,n,f,fulsms,x,funval)

        implicit none

!
!***  To compute fully symmetric basic rule sums
!
        integer::s,ix,lx,i,l,mi,ml,il,n,sumcls
        integer,parameter::mx=17
        integer,dimension(*)::m
        double precision::intwgt
        double precision,dimension(*)::x,fulsms,funval
        double precision,dimension(0:mx),save::g
        external:: f

!
!        Generators for 1 + 2 + 6 + 10 + 16 = 35 point
!        degree 51 rule, with degree 1, 5, 15 and 29 imbedded rules.
!

        DATA G( 0),G( 1)/0D0                  , 0.17320508075688773D1/
        DATA G( 2),G( 3)/0.41849560176727319D1, 0.74109534999454084D0/
        DATA G( 4),G( 5)/0.28612795760570581D1, 0.63633944943363700D1/
        DATA G( 6),G( 7)/0.12304236340273060D1, 0.51870160399136561D1/
        DATA G( 8),G( 9)/0.25960831150492022D1, 0.32053337944991945D1/
        DATA G(10),G(11)/0.90169397898903025D1, 0.24899229757996061D0/
        DATA G(12),G(13)/0.79807717985905609D1, 0.22336260616769417D1/
        DATA G(14),G(15)/0.71221067008046167D1, 0.36353185190372782D1/
        DATA G(16),G(17)/0.56981777684881096D1, 0.47364330859522971D1/

! ai ajoute la valeur initiale suivante:
        lx=0

        intwgt=1
        do i=1,s
                if (m(i).ne.0) intwgt=intwgt/2
        end do

        sumcls=0

        do i=1,n
                fulsms(i)=0
        end do

!
!*******  Compute centrally symmetric sum for permutation M
!
 10   do i = 1, s
         x(i) = -g(m(i))
      end do

 20   sumcls = sumcls + 1

        call f(s,x,n,funval)

        do i = 1,n
                fulsms(i) = fulsms(i) + intwgt*funval(i)
        end do

        do i = 1, s
                x(i) = -x(i)
                if (x(i) .gt. 0 ) go to 20
        end do
!
!*******  END Integration loop for M
!
!
!*******  Find next distinct permutation of M and loop back
!          to compute next centrally symmetric sum
!
        do i = 2, s
                if (m(i-1) .gt. m(i)) then

                        mi = m(i)
                        ix = i - 1

                        if (i .gt. 2) then
                                do l = 1,ix/2
                                        ml = m(l)
                                        il = i - l
                                        m(l) = m(il)
                                        m(il) = ml
                                        if ( ml .le. mi ) ix = ix - 1
                                        if (m(l) .gt. mi) lx = l
                                end do
                                if (m(ix) .le. mi) ix = lx
                        end if

                        m(i) = m(ix)
                        m(ix) = mi
                        goto 10
                end if
        end do
!
!***  END Loop for permutations of M and associated sums
!
!
!*** Restore original order to M.
!
        do i= 1,s/2
                mi = m(i)
                m(I) = m(s-i+1)
                m(s-i+1) = mi
        end do

        fulsms(n+1) = sumcls

        end subroutine fulsmh


!######################################################################################""

        subroutine nxpart(prtcnt,s,m,modofm)

        implicit none
!
!***  SUBROUTINE TO COMPUTE THE NEXT S PARTITION
!
        integer::s, prtcnt,modofm,i,msum,l
        integer,dimension(s)::m


        if (prtcnt .eq. 0) then

                do i = 1,s
                        m(i) = 0
                end do

                prtcnt = 1
                modofm = 0
        else
                prtcnt = prtcnt + 1
                msum = m(1)

                do i= 2,s
                        msum = msum + m(i)
                        if (m(1) .le. m(I) + 1) then
                                m(i) = 0
                        else
                                m(1) = msum - (i-1)*(m(i)+1)
                                do l= 2,i
                                        m(l) = m(i) + 1
                                end do
                                return
                        end if
                end do

                m(1) = msum + 1
                modofm = m(1)
        end if

        end subroutine nxpart

!###########################################################################################

        subroutine mltrul(ndim,numfun,funsub,np,point,weight,intval,funs,x,ic)

        implicit none
!
!     Computes product integration rule
!
        external::funsub
        integer::ndim,np,numfun,ici,i
        double precision::wtprod
        double precision,dimension(*)::point,weight,intval,x,funs,ic

        do i = 1,ndim
                ic(i) = 1.d0
        end do

        do i = 1,numfun
                intval(i) = 0.d0
        end do

 10   wtprod = 1.d0

        do i= 1,ndim
                ici = INT(ic(i))
                x(i) = point(ici)
                wtprod = wtprod*weight(ici)
        end do

        call funsub(ndim,x,numfun,funs)

        do i = 1,numfun
                intval(i) = intval(i) + wtprod*funs(i)
        end do

        do i = 1,ndim
                ic(i) = ic(i) + 1.d0
                if (ic(i) .le. np) goto 10
                ic(i) = 1.d0
        end do

        end subroutine mltrul

!###########################################################################"""

        subroutine hermit(ndim,numfun,minpts,maxpts,funsub,espabs,epsrel,restar, &
        result,abserr,neval,ifail,work)

        implicit none

!***BEGIN PROLOGUE HERMIT
!***CATEGORY NO. H2B1A1
!***AUTHOR
!
!     Alan Genz
!     Department of Mathematics
!     Washington State University
!     Pullman, Washington 99164-3113  USA
!     Telephone: 509-335-2131
!     Electronic Mail: alangenz@wsu.edu
!     Fax: 509-335-1188
!
!***KEYWORDS automatic multidimensional integrator, Gaussian weight
!            n-dimensional region (-infin, infin)^n
!***PURPOSE  The routine calculates an approximation to a vector of
!            definite integrals using Gauss-Hermite product rules.
!
!      infin     infin
!     I    ...  I       w(X) (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
!     -infin    -infin         1  2      NUMFUN
!
!       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN,
!              I   I  1  2      NDIM
!
!       w(X) = EXP(-( X(1)**2 + ... + X(NDIM)**2 )/2)/SQRT(2*PI)**NDIM.
!
!***DESCRIPTION Computation of integrals over infinite regions with
!               Gaussian weight using Gauss-Hermite product rules.
!
!   ON ENTRY
!
!     NDIM   Integer, number of variables.
!     NUMFUN Integer, number of components of the integral.
!     MINPTS Integer, minimum number of function evaluations.
!     MAXPTS Integer, maximum number of function evaluations.
!     FUNSUB Externally declared subroutine for computing components of
!            the integrand at the given evaluation point.
!            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
!            Input parameters:
!              NDIM   Integer, number of dimensions for the integral.
!              X      Real NDIM-array for the evaluation point.
!              NUMFUN Integer, number of components of I.
!            Output parameter:
!              FUNVLS Real NUMFUN-array for integrand components.
!
!     EPSABS Real, requested absolute accuracy.
!     EPSREL Real, requested relative accuracy.
!     RESTAR Integer.
!            If RESTAR = 0, this is the first attempt to compute
!            the integral.
!            If RESTAR = 1, then we restart a previous attempt.
!     WORK   Real array of working storage, must have dimensions at
!            least 2*NUMFUN+2*NDIM
!
!   ON RETURN
!
!     RESULT Real NUMFUN-array of approximations to all components
!            of the integral.
!     ABSERR Real NUMFUN-array of estimates of absolute accuracies.
!     NEVAL  Integer, number of function evaluations used.
!     IFAIL  Integer.
!            IFAIL = 0 for normal exit, when
!              ABSERR(K) <= MAX( EPSABS, ABS(RESULT(K))*EPSREL ) for
!              all K, 0 < K <= NUMFUN, with <= MAXPTS function values.
!            IFAIL = 1 if MAXPTS was too small to obtain the required
!              accuracy. In this case values of RESULT are returned
!              with estimated absolute accuracies ABSERR.
!
!***ROUTINES CALLED MLTRUL
!***END PROLOGUE HERMIT
!
!   Global variables.
!
        external::funsub
        integer::ndim,numfun,minpts,maxpts,restar,neval,ifail,i
        integer,save::rule
        integer,parameter::maxrul=50
        double precision::espabs,epsrel
        double precision,dimension(*)::result,abserr,work
        double precision,dimension(maxrul)::point,weight
        double precision,dimension(25,50),save::w,t
!
!   Local variables.
!

!
!     Gauss Hermite Weights and Points, N = 1,50
!
        DATA   W(1, 1), T(1, 1) / 1D0, 0D0 /
        DATA ( W(I, 2), T(I, 2), I = 1, 1) /0.5000000000000001D+00, 0.1000000000000000D+01/
        DATA ( W(I, 3), T(I, 3), I = 1, 2) / &
        0.1666666666666667D+00, 0.1732050807568877D+01, &
        0.6666666666666664D+00, 0.1107367643833737D-15/
        DATA ( W(I, 4), T(I, 4), I = 1, 2) / &
        0.4587585476806855D-01, 0.2334414218338977D+01, &
        0.4541241452319317D+00, 0.7419637843027258D+00/
        DATA ( W(I, 5), T(I, 5), I = 1, 3) /&
        0.1125741132772071D-01, 0.2856970013872805D+01,&
        0.2220759220056126D+00, 0.1355626179974265D+01,&
        0.5333333333333342D+00, 0.9386691848789097D-16/
        DATA ( W(I, 6), T(I, 6), I = 1, 3) / &
        0.2555784402056243D-02, 0.3324257433552119D+01,&
        0.8861574604191447D-01, 0.1889175877753710D+01,&
        0.4088284695560291D+00, 0.6167065901925933D+00/
        DATA ( W(I, 7), T(I, 7), I = 1, 4) / &
        0.5482688559722184D-03, 0.3750439717725742D+01,&
        0.3075712396758645D-01, 0.2366759410734542D+01,&
        0.2401231786050126D+00, 0.1154405394739968D+01,&
        0.4571428571428575D+00, 0.2669848554723344D-16/
        DATA ( W(I, 8), T(I, 8), I = 1, 4) /&
        0.1126145383753679D-03, 0.4144547186125893D+01,&
        0.9635220120788268D-02, 0.2802485861287542D+01,&
        0.1172399076617590D+00, 0.1636519042435109D+01,&
        0.3730122576790775D+00, 0.5390798113513754D+00/
        DATA ( W(I, 9), T(I, 9), I = 1, 5) /&
        0.2234584400774664D-04, 0.4512745863399781D+01,&
        0.2789141321231769D-02, 0.3205429002856470D+01,&
        0.4991640676521780D-01, 0.2076847978677829D+01,&
        0.2440975028949394D+00, 0.1023255663789133D+01,&
        0.4063492063492066D+00, 0.0000000000000000D+00/
        DATA ( W(I,10), T(I,10), I = 1, 5) / &
        0.4310652630718282D-05, 0.4859462828332311D+01,&
        0.7580709343122187D-03, 0.3581823483551927D+01,&
        0.1911158050077027D-01, 0.2484325841638954D+01,&
        0.1354837029802680D+00, 0.1465989094391158D+01,&
        0.3446423349320194D+00, 0.4849357075154977D+00/
        DATA ( W(I,11), T(I,11), I = 1, 6) /&
        0.8121849790214922D-06, 0.5188001224374871D+01,&
        0.1956719302712241D-03, 0.3936166607129977D+01,&
        0.6720285235537304D-02, 0.2865123160643646D+01,&
        0.6613874607105794D-01, 0.1876035020154847D+01,&
        0.2422402998739701D+00, 0.9288689973810635D+00,&
        0.3694083694083690D+00, 0.0000000000000000D+00/
        DATA ( W(I,12), T(I,12), I = 1, 6) /&
        0.1499927167637166D-06, 0.5500901704467746D+01,&
        0.4837184922590630D-04, 0.4271825847932281D+01,&
        0.2203380687533207D-02, 0.3223709828770096D+01,&
        0.2911668791236414D-01, 0.2259464451000800D+01,&
        0.1469670480453302D+00, 0.1340375197151617D+01,&
        0.3216643615128298D+00, 0.4444030019441390D+00/
        DATA ( W(I,13), T(I,13), I = 1, 7) /&
        0.2722627642805887D-07, 0.5800167252386502D+01,&
        0.1152659652733391D-04, 0.4591398448936520D+01,&
        0.6812363504429268D-03, 0.3563444380281636D+01,&
        0.1177056050599653D-01, 0.2620689973432215D+01,&
        0.7916895586044999D-01, 0.1725418379588239D+01,&
        0.2378715229641365D+00, 0.8566794935194499D+00,&
        0.3409923409923412D+00, 0.2011511664336819D-15/
        DATA ( W(I,14), T(I,14), I = 1, 7) /&
        0.4868161257748367D-08, 0.6087409546901291D+01,&
        0.2660991344067620D-05, 0.4896936397345567D+01,&
        0.2003395537607445D-03, 0.3886924575059772D+01,&
        0.4428919106947401D-02, 0.2963036579838668D+01,&
        0.3865010882425336D-01, 0.2088344745701943D+01,&
        0.1540833398425136D+00, 0.1242688955485464D+01,&
        0.3026346268130198D+00, 0.4125904579546022D+00/
        DATA ( W(I,15), T(I,15), I = 1, 8) /&
        0.8589649899633300D-09, 0.6363947888829836D+01,&
        0.5975419597920602D-06, 0.5190093591304780D+01,&
        0.5642146405189029D-04, 0.4196207711269018D+01,&
        0.1567357503549958D-02, 0.3289082424398766D+01,&
        0.1736577449213763D-01, 0.2432436827009758D+01,&
        0.8941779539984458D-01, 0.1606710069028730D+01,&
        0.2324622936097322D+00, 0.7991290683245483D+00,&
        0.3182595182595181D+00, 0.0000000000000000D+00/
        DATA ( W(I,16), T(I,16), I = 1, 8) /&
        0.1497814723161838D-09, 0.6630878198393126D+01,&
        0.1309473216286842D-06, 0.5472225705949343D+01,&
        0.1530003216248727D-04, 0.4492955302520013D+01,&
        0.5259849265739089D-03, 0.3600873624171548D+01,&
        0.7266937601184742D-02, 0.2760245047630703D+01,&
        0.4728475235401395D-01, 0.1951980345716333D+01,&
        0.1583383727509496D+00, 0.1163829100554964D+01,&
        0.2865685212380120D+00, 0.3867606045005573D+00/
        DATA ( W(I,17), T(I,17), I = 1, 9) /&
        0.2584314919374932D-10, 0.6889122439895331D+01,&
        0.2808016117930569D-07, 0.5744460078659410D+01,&
        0.4012679447979839D-05, 0.4778531589629983D+01,&
        0.1684914315513387D-03, 0.3900065717198010D+01,&
        0.2858946062284621D-02, 0.3073797175328194D+01,&
        0.2308665702571097D-01, 0.2281019440252989D+01,&
        0.9740637116272111D-01, 0.1509883307796740D+01,&
        0.2267063084689769D+00, 0.7518426007038956D+00,&
        0.2995383701266057D+00, 0.0000000000000000D+00/
        DATA ( W(I,18), T(I,18), I = 1, 9) /&
        0.4416588769358736D-11, 0.7139464849146476D+01,&
        0.5905488478836554D-08, 0.6007745911359599D+01,&
        0.1021552397636983D-05, 0.5054072685442739D+01,&
        0.5179896144116204D-04, 0.4188020231629400D+01,&
        0.1065484796291652D-02, 0.3374736535778089D+01,&
        0.1051651775194131D-01, 0.2595833688911239D+01,&
        0.5489663248022256D-01, 0.1839779921508646D+01,&
        0.1606853038935128D+00, 0.1098395518091501D+01,&
        0.2727832346542882D+00, 0.3652457555076979D+00/
        DATA ( W(I,19), T(I,19), I = 1,10) /&
        0.7482830054057162D-12, 0.7382579024030434D+01,&
        0.1220370848447449D-08, 0.6262891156513252D+01,&
        0.2532220032092866D-06, 0.5320536377336039D+01,&
        0.1535114595466674D-04, 0.4465872626831029D+01,&
        0.3785021094142701D-03, 0.3664416547450636D+01,&
        0.4507235420342067D-02, 0.2898051276515753D+01,&
        0.2866669103011841D-01, 0.2155502761316934D+01,&
        0.1036036572761442D+00, 0.1428876676078373D+01,&
        0.2209417121991433D+00, 0.7120850440423796D+00,&
        0.2837731927515210D+00, 0.4118522463420039D-15/
        DATA ( W(I,20), T(I,20), I = 1,10) /&
        0.1257800672437914D-12, 0.7619048541679760D+01,&
        0.2482062362315163D-09, 0.6510590157013660D+01,&
        0.6127490259983006D-07, 0.5578738805893195D+01,&
        0.4402121090230841D-05, 0.4734581334046057D+01,&
        0.1288262799619300D-03, 0.3943967350657311D+01,&
        0.1830103131080496D-02, 0.3189014816553389D+01,&
        0.1399783744710099D-01, 0.2458663611172367D+01,&
        0.6150637206397690D-01, 0.1745247320814126D+01,&
        0.1617393339840001D+00, 0.1042945348802752D+01,&
        0.2607930634495551D+00, 0.3469641570813557D+00/
        DATA ( W(I,21), T(I,21), I = 1,11) /&
        0.2098991219565665D-13, 0.7849382895113822D+01,&
        0.4975368604121770D-10, 0.6751444718717456D+01,&
        0.1450661284493093D-07, 0.5829382007304472D+01,&
        0.1225354836148259D-05, 0.4994963944782024D+01,&
        0.4219234742551696D-04, 0.4214343981688420D+01,&
        0.7080477954815349D-03, 0.3469846690475375D+01,&
        0.6439697051408779D-02, 0.2750592981052372D+01,&
        0.3395272978654278D-01, 0.2049102468257161D+01,&
        0.1083922856264195D+00, 0.1359765823211230D+01,&
        0.2153337156950595D+00, 0.6780456924406435D+00,&
        0.2702601835728773D+00, 0.0000000000000000D+00/
        DATA ( W(I,22), T(I,22), I = 1,11) /&
        0.3479460647877136D-14, 0.8074029984021710D+01,&
        0.9841378982346056D-11, 0.6985980424018808D+01,&
        0.3366514159458310D-08, 0.6073074951122888D+01,&
        0.3319853749814059D-06, 0.5247724433714421D+01,&
        0.1334597712680954D-04, 0.4476361977310866D+01,&
        0.2622833032559635D-03, 0.3741496350266517D+01,&
        0.2808761047577212D-02, 0.3032404227831676D+01,&
        0.1756907288080571D-01, 0.2341759996287707D+01,&
        0.6719631142889003D-01, 0.1664124839117906D+01,&
        0.1619062934136754D+00, 0.9951624222712152D+00,&
        0.2502435965869353D+00, 0.3311793157152742D+00/
        DATA ( W(I,23), T(I,23), I = 1,12) /&
        0.5732383167802038D-15, 0.8293386027417354D+01,&
        0.1922935311567786D-11, 0.7214659435051866D+01,&
        0.7670888862399855D-09, 0.6310349854448401D+01,&
        0.8775062483861979D-07, 0.5493473986471793D+01,&
        0.4089977244992140D-05, 0.4730724197451473D+01,&
        0.9340818609031275D-04, 0.4004775321733304D+01,&
        0.1167628637497855D-02, 0.3305040021752963D+01,&
        0.8579678391465647D-02, 0.2624323634059181D+01,&
        0.3886718370348111D-01, 0.1957327552933424D+01,&
        0.1120733826026210D+00, 0.1299876468303978D+01,&
        0.2099596695775429D+00, 0.6484711535344957D+00,&
        0.2585097408088385D+00, 0.0000000000000000D+00/
        DATA ( W(I,24), T(I,24), I = 1,12) /&
        0.9390193689041782D-16, 0.8507803519195264D+01,&
        0.3714974152762395D-12, 0.7437890666021664D+01,&
        0.1718664927964866D-09, 0.6541675005098631D+01,&
        0.2267461673480609D-07, 0.5732747175251204D+01,&
        0.1217659745442582D-05, 0.4978041374639117D+01,&
        0.3209500565274598D-04, 0.4260383605019904D+01,&
        0.4647187187793975D-03, 0.3569306764073560D+01,&
        0.3976608929181313D-02, 0.2897728643223314D+01,&
        0.2112634440896754D-01, 0.2240467851691752D+01,&
        0.7206936401717838D-01, 0.1593480429816420D+01,&
        0.1614595128670001D+00, 0.9534219229321088D+00,&
        0.2408701155466405D+00, 0.3173700966294525D+00/
        DATA ( W(I,25), T(I,25), I = 1,13) /&
        0.1530038997998690D-16, 0.8717597678399592D+01,&
        0.7102103037003980D-13, 0.7656037955393078D+01,&
        0.3791150000477161D-10, 0.6767464963809719D+01,&
        0.5738023868899356D-08, 0.5966014690606704D+01,&
        0.3530152560245470D-06, 0.5218848093644280D+01,&
        0.1067219490520254D-04, 0.4508929922967284D+01,&
        0.1777669069265268D-03, 0.3825900569972490D+01,&
        0.1757850405263803D-02, 0.3162775679388193D+01,&
        0.1085675599146230D-01, 0.2514473303952205D+01,&
        0.4337997016764489D-01, 0.1877058369947839D+01,&
        0.1148809243039517D+00, 0.1247311975616789D+01,&
        0.2048510256503405D+00, 0.6224622791860757D+00,&
        0.2481693511764858D+00, 0.0000000000000000D+00/
        DATA ( W(I,26), T(I,26), I = 1,13) /&
        0.2480694260393664D-17, 0.8923051727828243D+01,&
        0.1344547649663596D-13, 0.7869426697637738D+01,&
        0.8242809443163844D-11, 0.6988088770623415D+01,&
        0.1424293237988014D-08, 0.6193693483796317D+01,&
        0.9986755573314568D-07, 0.5453615383857833D+01,&
        0.3443413612308114D-05, 0.4750947483085378D+01,&
        0.6557558694333818D-04, 0.4075427214412228D+01,&
        0.7442025763604303D-03, 0.3420156373999979D+01,&
        0.5302198015682246D-02, 0.2780138499509748D+01,&
        0.2459766565712125D-01, 0.2151530090121648D+01,&
        0.7622953220630281D-01, 0.1531215708695402D+01,&
        0.1605865456137948D+00, 0.9165450413386282D+00,&
        0.2324707356300776D+00, 0.3051559707592978D+00/
        DATA ( W(I,27), T(I,27), I = 1,14) /&
        0.4003364766550257D-18, 0.9124421250672931D+01,&
        0.2522363250873417D-14, 0.8078349274534165D+01,&
        0.1768236511219616D-11, 0.7203876611910644D+01,&
        0.3472604702845840D-09, 0.6416154934562174D+01,&
        0.2761933918447925D-07, 0.5682760761629052D+01,&
        0.1080581533683832D-05, 0.4986906410679802D+01,&
        0.2339557671566820D-04, 0.4318417671936682D+01,&
        0.3028398259361645D-03, 0.3670472986492407D+01,&
        0.2471872445961970D-02, 0.3038150251871036D+01,&
        0.1321102584046355D-01, 0.2417683983162542D+01,&
        0.4748957556274264D-01, 0.1806045213138672D+01,&
        0.1169962651750102D+00, 0.1200683354549981D+01,&
        0.2000149701605136D+00, 0.5993548807899852D+00,&
        0.2389778937255043D+00, 0.0000000000000000D+00/
        DATA ( W(I,28), T(I,28), I = 1,14) /&
        0.6432547438801930D-19, 0.9321937814408766D+01,&
        0.4691765569500354D-15, 0.8283069540861421D+01,&
        0.3745901035176660D-12, 0.7415125286176065D+01,&
        0.8326609843882241D-10, 0.6633731493950435D+01,&
        0.7479362584613589D-08, 0.5906656325824994D+01,&
        0.3304864449926482D-06, 0.5217223673447450D+01,&
        0.8093584057145153D-05, 0.4555340384596974D+01,&
        0.1188285381401779D-03, 0.3914253725963635D+01,&
        0.1104305927857598D-02, 0.3289106970171833D+01,&
        0.6752459709030160D-02, 0.2676201879526944D+01,&
        0.2793578476788097D-01, 0.2072582674144621D+01,&
        0.7977336601159966D-01, 0.1475781736957922D+01,&
        0.1594181936613094D+00, 0.8836525629929802D+00,&
        0.2248886297506769D+00, 0.2942517144887133D+00/
        DATA ( W(I,29), T(I,29), I = 1,15) /&
        0.1029341808721942D-19, 0.9515812006947357D+01,&
        0.8657491667957282D-16, 0.8483826557846555D+01,&
        0.7842840425658472D-13, 0.7622102722480985D+01,&
        0.1965709944734762D-10, 0.6846722135707994D+01,&
        0.1986123546067053D-08, 0.6125635348243716D+01,&
        0.9868968543560684D-07, 0.5442271089178505D+01,&
        0.2721127828058099D-05, 0.4786611062352805D+01,&
        0.4508394026980976D-04, 0.4151964855100983D+01,&
        0.4743663738893483D-03, 0.3533533770990993D+01,&
        0.3297972210833669D-02, 0.2927678154322886D+01,&
        0.1559400577786720D-01, 0.2331504884065565D+01,&
        0.5121083528871912D-01, 0.1742616232610662D+01,&
        0.1185603192669036D+00, 0.1158946149400189D+01,&
        0.1954459569679010D+00, 0.5786461780331649D+00,&
        0.2307372767004873D+00, 0.0000000000000000D+00/
        DATA ( W(I,30), T(I,30), I = 1,15) /&
        0.1640807008117853D-20, 0.9706235997359524D+01,&
        0.1585560944966296D-16, 0.8680837722732207D+01,&
        0.1624080129972436D-13, 0.7825051744352813D+01,&
        0.4573425871326147D-11, 0.7055396866960296D+01,&
        0.5178459467189710D-09, 0.6339997686869597D+01,&
        0.2882175154047618D-07, 0.5662381850082873D+01,&
        0.8909088868621158D-06, 0.5012600596486518D+01,&
        0.1657998163067346D-04, 0.4384020365898051D+01,&
        0.1965129439848249D-03, 0.3771894423159236D+01,&
        0.1544707339866097D-02, 0.3172634639420402D+01,&
        0.8295747557723240D-02, 0.2583402100229274D+01,&
        0.3111177018350134D-01, 0.2001858612956431D+01,&
        0.8278683671562172D-01, 0.1426005658374115D+01,&
        0.1580469532090208D+00, 0.8540733517109733D+00,&
        0.2179999718155776D+00, 0.2844387607362094D+00/
        DATA ( W(I,31), T(I,31), I = 1,16) /&
        0.2605973854893011D-21, 0.9893385708986649D+01,&
        0.2883352367857899D-17, 0.8874301409488794D+01,&
        0.3328468324148409D-14, 0.8024193227361653D+01,&
        0.1049603362311349D-11, 0.7260000488890867D+01,&
        0.1327251483589731D-09, 0.6550014268765684D+01,&
        0.8243931619119761D-08, 0.5877855885986261D+01,&
        0.2845610088162858D-06, 0.5233641511712708D+01,&
        0.5923202317686233D-05, 0.4610789797323995D+01,&
        0.7871624069602249D-04, 0.4004600901491224D+01,&
        0.6960312713792868D-03, 0.3411532415843158D+01,&
        0.4221717767270697D-02, 0.2828792768157509D+01,&
        0.1796787584344161D-01, 0.2254095000754410D+01,&
        0.5456725889447496D-01, 0.1685497905069052D+01,&
        0.1196831096958545D+00, 0.1121297374047009D+01,&
        0.1911320047746435D+00, 0.5599475878410030D+00,&
        0.2232941387424060D+00, 0.9191380905810332D-16/
        DATA ( W(I,32), T(I,32), I = 1,16) /&
        0.4124607489018384D-22, 0.1007742267422945D+02,&
        0.5208449591960853D-18, 0.9064399210702408D+01,&
        0.6755290223670036D-15, 0.8219728765382246D+01,&
        0.2378064855777808D-12, 0.7460755754121516D+01,&
        0.3347501239801238D-10, 0.6755930830540704D+01,&
        0.2312518412074224D-08, 0.6088964309076983D+01,&
        0.8881290713105934D-07, 0.5450033273623426D+01,&
        0.2059622103953437D-05, 0.4832604613244488D+01,&
        0.3055980306089618D-04, 0.4232021109995410D+01,&
        0.3025570258170642D-03, 0.3644781249880835D+01,&
        0.2062051051307883D-02, 0.3068135169013122D+01,&
        0.9903461702320572D-02, 0.2499840415187396D+01,&
        0.3410984772609194D-01, 0.1938004905925718D+01,&
        0.8534480827208071D-01, 0.1380980199272144D+01,&
        0.1565389937575984D+00, 0.8272849037797656D+00,&
        0.2117055698804795D+00, 0.2755464192302757D+00/
        DATA ( W(I,33), T(I,33), I = 1,17) /&
        0.6506889970402893D-23, 0.1025849562613868D+02,&
        0.9349155921250728D-19, 0.9251297851734609D+01,&
        0.1358447599822066D-15, 0.8411842935668213D+01,&
        0.5323023871225154D-13, 0.7657866034784416D+01,&
        0.8316046910555386D-11, 0.6957971061087896D+01,&
        0.6369261724538402D-09, 0.6295953125159368D+01,&
        0.2712480030928315D-07, 0.5662046690089215D+01,&
        0.6982937054010693D-06, 0.5049763451908822D+01,&
        0.1152282979948381D-04, 0.4454485185983318D+01,&
        0.1271924628565943D-03, 0.3872747224621657D+01,&
        0.9695342064240463D-03, 0.3301836743259241D+01,&
        0.5227605765843533D-02, 0.2739550282026145D+01,&
        0.2030404475704196D-01, 0.2184038256077331D+01,&
        0.5758631173578715D-01, 0.1633699795932273D+01,&
        0.1204510521056565D+00, 0.1087107916669903D+01,&
        0.1870581852279859D+00, 0.5429533766656418D+00,&
        0.2165276496896071D+00, 0.0000000000000000D+00/
        DATA ( W(I,34), T(I,34), I = 1,17) /&
        0.1023327129805438D-23, 0.1043674187100505D+02,&
        0.1668144375578546D-19, 0.9435150833760799D+01,&
        0.2708054709262291D-16, 0.8600705233983431D+01,&
        0.1177937790622538D-13, 0.7851517590699072D+01,&
        0.2036657679770825D-11, 0.7156339259134842D+01,&
        0.1724305566745258D-09, 0.6499046353927969D+01,&
        0.8117409040122166D-08, 0.5869927588596779D+01,&
        0.2312034264322868D-06, 0.5262536481734633D+01,&
        0.4227725748387717D-05, 0.4672290691400979D+01,&
        0.5182712643366873D-04, 0.4095758971954162D+01,&
        0.4399649667746255D-03, 0.3530261634074169D+01,&
        0.2650824238310194D-02, 0.2973629650303989D+01,&
        0.1155073894677711D-01, 0.2424050509756231D+01,&
        0.3692346295804431D-01, 0.1879964366420186D+01,&
        0.8751168135862399D-01, 0.1339990625522619D+01,&
        0.1549420965148600D+00, 0.8028734607837125D+00,&
        0.2059249366691133D+00, 0.2674391839515571D+00/
        DATA ( W(I,35), T(I,35), I = 1,18) /&
        0.1604619191790137D-24, 0.1061228847764259D+02,&
        0.2959542628709020D-20, 0.9616099851106188D+01,&
        0.5354094198066113D-17, 0.8786471736571588D+01,&
        0.2578597004420350D-14, 0.8041881508402966D+01,&
        0.4921158318497817D-12, 0.7351222593955577D+01,&
        0.4592895709255784D-10, 0.6698448669188526D+01,&
        0.2383123333705684D-08, 0.6073899909835108D+01,&
        0.7486460178518943D-07, 0.5471169042499284D+01,&
        0.1511924078161164D-05, 0.4885706922446520D+01,&
        0.2050985370117561D-04, 0.4314112806712937D+01,&
        0.1931447146139865D-03, 0.3753736849678044D+01,&
        0.1294831077498348D-02, 0.3202440651527510D+01,&
        0.6300195959720374D-02, 0.2658444295466304D+01,&
        0.2258156121393648D-01, 0.2120223604836576D+01,&
        0.6029658086774051D-01, 0.1586437891088980D+01,&
        0.1209324521970301D+00, 0.1055876792225099D+01,&
        0.1832085621911518D+00, 0.5274192342262778D+00,&
        0.2103411454127607D+00, 0.2220535009031490D-16/
        DATA ( W(I,36), T(I,36), I = 1,18) /&
        0.2509037634634927D-25, 0.1078525331238753D+02,&
        0.5222366200862934D-21, 0.9794276019583023D+01,&
        0.1050294193474738D-17, 0.8969286534562617D+01,&
        0.5587113978649920D-15, 0.8229115367471579D+01,&
        0.1174029255105150D-12, 0.7542793039211395D+01,&
        0.1204744586548017D-10, 0.6894347646173911D+01,&
        0.6871083387212228D-09, 0.6274168326809511D+01,&
        0.2373822777365743D-07, 0.5675884710106672D+01,&
        0.5278315192800947D-06, 0.5094978513857614D+01,&
        0.7896978086723001D-05, 0.4528076990600175D+01,&
        0.8220025562410442D-04, 0.3972557341929988D+01,&
        0.6107548335511020D-03, 0.3426308595129129D+01,&
        0.3304134538435289D-02, 0.2887579695004719D+01,&
        0.1321657401560215D-01, 0.2354877715992540D+01,&
        0.3955236976655980D-01, 0.1826896577986744D+01,&
        0.8934249750438648D-01, 0.1302464954480165D+01,&
        0.1532910133997116D+00, 0.7805064920524665D+00,&
        0.2005920064390222D+00, 0.2600079252490002D+00/
        DATA ( W(I,37), T(I,37), I = 1,19) /&
        0.3912701900272739D-26, 0.1095574594356165D+02,&
        0.9167997950194993D-22, 0.9969800945691102D+01,&
        0.2045033189096993D-18, 0.9149282977696753D+01,&
        0.1198838092763842D-15, 0.8413364679488879D+01,&
        0.2767202859289950D-13, 0.7731209035776908D+01,&
        0.3114550633754001D-11, 0.7086915685017345D+01,&
        0.1947529536776704D-09, 0.6470920475390299D+01,&
        0.7379429831202629D-08, 0.5876887892594201D+01,&
        0.1801392729842884D-06, 0.5300328473991422D+01,&
        0.2963204697508371D-05, 0.4737895299703953D+01,&
        0.3397941877127311D-04, 0.4186990225627911D+01,&
        0.2788064134714879D-03, 0.3645526993515297D+01,&
        0.1670452623119233D-02, 0.3111780274991816D+01,&
        0.7424836460756516D-02, 0.2584284718210775D+01,&
        0.2478561700046269D-01, 0.2061764482625703D+01,&
        0.6272612686101582D-01, 0.1543082026656558D+01,&
        0.1211816464812620D+00, 0.1027199336691835D+01,&
        0.1795672590244481D+00, 0.5131472568284307D+00,&
        0.2046562495907946D+00, 0.0000000000000000D+00/
        DATA ( W(I,38), T(I,38), I = 1,19) /&
        0.6086019568424894D-27, 0.1112386843494987D+02,&
        0.1601586834974089D-22, 0.1014278766112967D+02,&
        0.3953752210235847D-19, 0.9326584757395571D+01,&
        0.2548653789376282D-16, 0.8594764136386107D+01,&
        0.6447906524854284D-14, 0.7916616928480361D+01,&
        0.7941720700519787D-12, 0.7276311665555140D+01,&
        0.5431469647397207D-10, 0.6664328864414711D+01,&
        0.2251484789660455D-08, 0.6074366042059564D+01,&
        0.6017570329946914D-07, 0.5501960756811018D+01,&
        0.1085199304765951D-05, 0.4943790029802265D+01,&
        0.1366659848918500D-04, 0.4397278309932714D+01,&
        0.1234207566124386D-03, 0.3860361738163258D+01,&
        0.8159914125966156D-03, 0.3331338060544006D+01,&
        0.4014342501780207D-02, 0.2808766413901917D+01,&
        0.1488362222069118D-01, 0.2291397563624214D+01,&
        0.4200051997416258D-01, 0.1778123425136563D+01,&
        0.9088415555258827D-01, 0.1267939114780749D+01,&
        0.1516111462601078D+00, 0.7599132481367392D+00,&
        0.1956519870413640D+00, 0.2531636353359348D+00/
        DATA ( W(I,39), T(I,39), I = 1,20) /&
        0.9443344575063092D-28, 0.1128971604447632D+02,&
        0.2784787505225604D-23, 0.1031334144275699D+02,&
        0.7592450542206611D-20, 0.9501306853782864D+01,&
        0.5370701458462819D-17, 0.8773438698067656D+01,&
        0.1486129587733307D-14, 0.8099152213152410D+01,&
        0.1998726680623669D-12, 0.7462682377866057D+01,&
        0.1491720010448932D-10, 0.6854552519790964D+01,&
        0.6748646364787737D-09, 0.6268491549685291D+01,&
        0.1969872920599359D-07, 0.5700062454271859D+01,&
        0.3884118283713278D-06, 0.5145964544094186D+01,&
        0.5356584901373566D-05, 0.4603643074288992D+01,&
        0.5307742112417270D-04, 0.4071054595947238D+01,&
        0.3859381697690359D-03, 0.3546517669768644D+01,&
        0.2093837438880663D-02, 0.3028613314092798D+01,&
        0.8588083029362219D-02, 0.2516115854339318D+01,&
        0.2690624148395815D-01, 0.2007943065265498D+01,&
        0.6490157458316884D-01, 0.1503118900228415D+01,&
        0.1212421280631244D+00, 0.1000744572356018D+01,&
        0.1761190277014499D+00, 0.4999751896072681D+00,&
        0.1994086534474408D+00, 0.1536204102353874D-15/
        DATA ( W(I,40), T(I,40), I = 1,20) /&
        0.1461839873869467D-28, 0.1145337784154873D+02,&
        0.4820467940200524D-24, 0.1048156053467427D+02,&
        0.1448609431551587D-20, 0.9673556366934033D+01,&
        0.1122275206827074D-17, 0.8949504543855559D+01,&
        0.3389853443248306D-15, 0.8278940623659475D+01,&
        0.4968088529197761D-13, 0.7646163764541459D+01,&
        0.4037638581695192D-11, 0.7041738406453829D+01,&
        0.1989118526027766D-09, 0.6459423377583766D+01,&
        0.6325897188548972D-08, 0.5894805675372016D+01,&
        0.1360342421574886D-06, 0.5344605445720084D+01,&
        0.2048897436081474D-05, 0.4806287192093873D+01,&
        0.2221177143247582D-04, 0.4277826156362752D+01,&
        0.1770729287992397D-03, 0.3757559776168985D+01,&
        0.1055879016901825D-02, 0.3244088732999869D+01,&
        0.4773544881823334D-02, 0.2736208340465433D+01,&
        0.1653784414256937D-01, 0.2232859218634873D+01,&
        0.4427455520227679D-01, 0.1733090590631720D+01,&
        0.9217657917006089D-01, 0.1236032004799159D+01,&
        0.1499211117635710D+00, 0.7408707252859313D+00,&
        0.1910590096619904D+00, 0.2468328960227240D+00/
        DATA ( W(I,41), T(I,41), I = 1,21) /&
        0.2257863956583089D-29, 0.1161493725433746D+02,&
        0.8308558938782992D-25, 0.1064753678631932D+02,&
        0.2746891228522292D-21, 0.9843433249157988D+01,&
        0.2326384145587187D-18, 0.9123069907984480D+01,&
        0.7655982291966812D-16, 0.8456099083269388D+01,&
        0.1220334874202772D-13, 0.7826882004053867D+01,&
        0.1077818394935909D-11, 0.7226022663732790D+01,&
        0.5769853428092003D-10, 0.6647308470747191D+01,&
        0.1994794756757339D-08, 0.6086349164878472D+01,&
        0.4667347708107243D-07, 0.5539884440458126D+01,&
        0.7658186077982435D-06, 0.5005396683404125D+01,&
        0.9058608622433030D-05, 0.4480878331594004D+01,&
        0.7894719319504627D-04, 0.3964684028033266D+01,&
        0.5158014443431912D-03, 0.3455432217780992D+01,&
        0.2561642428649777D-02, 0.2951937016381193D+01,&
        0.9777902738208298D-02, 0.2453159345907049D+01,&
        0.2893721174793441D-01, 0.1958170711977291D+01,&
        0.6684765935446599D-01, 0.1466125457295967D+01,&
        0.1211489170115104D+00, 0.9762387671800500D+00,&
        0.1728495310506020D+00, 0.4877685693194347D+00,&
        0.1945450277536008D+00, 0.2585532684499631D-15/
        DATA ( W(I,42), T(I,42), I = 1,21) /&
        0.3479841758734498D-30, 0.1177447255645880D+02,&
        0.1426197845863333D-25, 0.1081135621818894D+02,&
        0.5178070329449428D-22, 0.1001103095231325D+02,&
        0.4785541849652557D-19, 0.9294235815925036D+01,&
        0.1712817711028008D-16, 0.8630736540442662D+01,&
        0.2963878417982294D-14, 0.8004954459331870D+01,&
        0.2839400066530831D-12, 0.7407531683161432D+01,&
        0.1648408975918176D-10, 0.6832282984221202D+01,&
        0.6182418956905496D-09, 0.6274839704458262D+01,&
        0.1570405641380752D-07, 0.5731959941924930D+01,&
        0.2800444926136691D-06, 0.5201142761234950D+01,&
        0.3605310211304410D-05, 0.4680396489703305D+01,&
        0.3425737772910867D-04, 0.4168091525525806D+01,&
        0.2445266890869028D-03, 0.3662862441045986D+01,&
        0.1329885905572134D-02, 0.3163540283549010D+01,&
        0.5573924561218487D-02, 0.2669104116603810D+01,&
        0.1816809011555160D-01, 0.2178645208762502D+01,&
        0.4638274179115778D-01, 0.1691339732834912D+01,&
        0.9325376860459540D-01, 0.1206427277827926D+01,&
        0.1482345543992958D+00, 0.7231933449391704D+00,&
        0.1867743488620199D+00, 0.2409545348665914D+00/
        DATA ( W(I,43), T(I,43), I = 1,22) /&
        0.5352075224079761D-31, 0.1193205730105631D+02,&
        0.2438509785586889D-26, 0.1097309952495776D+02,&
        0.9705961174387400D-23, 0.1017643700187913D+02,&
        0.9772188868027905D-20, 0.9463096735522324D+01,&
        0.3797453851791020D-17, 0.8802954705722227D+01,&
        0.7121219095950012D-15, 0.8180490511448660D+01,&
        0.7386428797509359D-13, 0.7586383052570858D+01,&
        0.4641647576106650D-11, 0.7014473354182730D+01,&
        0.1884798249393598D-09, 0.6460413330895460D+01,&
        0.5186705460818634D-08, 0.5920978461471519D+01,&
        0.1003006675541783D-06, 0.5393683423125381D+01,&
        0.1402088827949991D-05, 0.4876551284706595D+01,&
        0.1448871182961238D-04, 0.4367966934806672D+01,&
        0.1126824016383807D-03, 0.3866579655692208D+01,&
        0.6691683916969021D-03, 0.3371235817148052D+01,&
        0.3070012395787847D-02, 0.2880930777450376D+01,&
        0.1098373316809114D-01, 0.2394773428842528D+01,&
        0.3087516724563331D-01, 0.1911959274951712D+01,&
        0.6858704943109335D-01, 0.1431749364451692D+01,&
        0.1209303683096624D+00, 0.9534532756297385D+00,&
        0.1697454597838771D+00, 0.4764148861084949D+00,&
        0.1900207247825863D+00, 0.0000000000000000D+00/
        DATA ( W(I,44), T(I,44), I = 1,22) /&
        0.8215250893174055D-32, 0.1208776070905845D+02,&
        0.4153634860809548D-27, 0.1113284252424393D+02,&
        0.1809481070701256D-23, 0.1033973350764328D+02,&
        0.1981515599992759D-20, 0.9629741154666469D+01,&
        0.8346676981526894D-18, 0.8972848703628870D+01,&
        0.1693427110357361D-15, 0.8353592294951005D+01,&
        0.1898518041004563D-13, 0.7762686386163553D+01,&
        0.1289064815963923D-11, 0.7193997236488432D+01,&
        0.5656565204102218D-10, 0.6643196399727595D+01,&
        0.1683045805404500D-08, 0.6107075817029749D+01,&
        0.3522080139402352D-07, 0.5583164829645622D+01,&
        0.5334152162808288D-06, 0.5069500234589782D+01,&
        0.5980506497871435D-05, 0.4564480302199131D+01,&
        0.5055024343502252D-04, 0.4066767790495557D+01,&
        0.3269029726329101D-03, 0.3575222999425195D+01,&
        0.1636879726834301D-02, 0.3088855994142003D+01,&
        0.6407990218344719D-02, 0.2606791462987426D+01,&
        0.1976568083382755D-01, 0.2128242119185467D+01,&
        0.4833422729538417D-01, 0.1652487989479395D+01,&
        0.9414471679496565D-01, 0.1178859803328855D+01,&
        0.1465614441714095D+00, 0.7067252381063451D+00,&
        0.1827650568597313D+00, 0.2354771181719222D+00/
        DATA ( W(I,45), T(I,45), I = 1,23) /&
        0.1258601266761721D-32, 0.1224164801738294D+02,&
        0.7049452438028993D-28, 0.1129065655801674D+02,&
        0.3355898215846036D-24, 0.1050099761933589D+02,&
        0.3990921764567262D-21, 0.9794252095357150D+01,&
        0.1819415421781025D-18, 0.9140507651218687D+01,&
        0.3987388394988911D-16, 0.8524355348620958D+01,&
        0.4823870286948576D-14, 0.7936544056928108D+01,&
        0.3532976789571211D-12, 0.7370964332170609D+01,&
        0.1672373578472200D-10, 0.6823306517591415D+01,&
        0.5370107061670978D-09, 0.6290378188903940D+01,&
        0.1213731971819947D-07, 0.5769722503692224D+01,&
        0.1987368663904533D-06, 0.5259389088450044D+01,&
        0.2412167676023244D-05, 0.4757788618871960D+01,&
        0.2210680714129846D-04, 0.4263596248773437D+01,&
        0.1552889103199588D-03, 0.3775684997351202D+01,&
        0.8463565255385183D-03, 0.3293078264001083D+01,&
        0.3614802525284821D-02, 0.2814914962782011D+01,&
        0.1219644981036976D-01, 0.2340423205074869D+01,&
        0.3271890357088906D-01, 0.1868899890086888D+01,&
        0.7014033220973423D-01, 0.1399694432579061D+01,&
        0.1206095596750258D+00, 0.9321954002025563D+00,&
        0.1667945553649214D+00, 0.4658191781783683D+00,&
        0.1857980420096402D+00, 0.0000000000000000D+00/
        DATA ( W(I,46), T(I,46), I = 1,23) /&
        0.1924664627046275D-33, 0.1239378079201923D+02,&
        0.1192246272132870D-28, 0.1144660885260523D+02,&
        0.6192859698073460D-25, 0.1066030193428052D+02,&
        0.7986127260533585D-22, 0.9956707572498630D+01,&
        0.3934545090204791D-19, 0.9306015173119082D+01,&
        0.9300353507700858D-17, 0.8692869193236842D+01,&
        0.1212244415430154D-14, 0.8108051845054231D+01,&
        0.9561363845435019D-13, 0.7545477116064648D+01,&
        0.4874112963243041D-11, 0.7000853362533224D+01,&
        0.1686106290702903D-09, 0.6471003045526279D+01,&
        0.4108222480852658D-08, 0.5953482378238384D+01,&
        0.7258481863429432D-07, 0.5446353016216582D+01,&
        0.9517578349775918D-06, 0.4948037176754664D+01,&
        0.9436598821520223D-05, 0.4457221460031489D+01,&
        0.7183252612530643D-04, 0.3972790546235975D+01,&
        0.4250585899081316D-03, 0.3493779976317873D+01,&
        0.1975260370663049D-02, 0.3019341533379349D+01,&
        0.7268729480993403D-02, 0.2548717171165183D+01,&
        0.2132401414735712D-01, 0.2081218863808104D+01,&
        0.5013852830109192D-01, 0.1616212619554741D+01,&
        0.9487419222577377D-01, 0.1153105444663487D+01,&
        0.1449090160374568D+00, 0.6913343908693306D+00,&
        0.1790029030973513D+00, 0.2303570440689780D+00/
        DATA ( W(I,47), T(I,47), I = 1,24) /&
        0.2937991918931770D-34, 0.1254421721021957D+02,&
        0.2009633017503842D-29, 0.1160076284239976D+02,&
        0.1137324853690828D-25, 0.1081771486308368D+02,&
        0.1588168143815981D-22, 0.1011718100450289D+02,&
        0.8443830093861547D-20, 0.9469449861345275D+01,&
        0.2149644881290608D-17, 0.8859217846077781D+01,&
        0.3014376026444397D-15, 0.8277299513807538D+01,&
        0.2556494006562142D-13, 0.7717631482085797D+01,&
        0.1401228890589135D-11, 0.7175939408289193D+01,&
        0.5213279817493255D-10, 0.6649059958212412D+01,&
        0.1366926136038253D-08, 0.6134561715294918D+01,&
        0.2601160240819878D-07, 0.5630517648304252D+01,&
        0.3677447240564392D-06, 0.5135360748767347D+01,&
        0.3936411635330397D-05, 0.4647788224870149D+01,&
        0.3239911617951789D-04, 0.4166695488690111D+01,&
        0.2076565749155176D-03, 0.3691129181403769D+01,&
        0.1047278482271871D-02, 0.3220252777152719D+01,&
        0.4191754574650512D-02, 0.2753320728634848D+01,&
        0.1340827900837288D-01, 0.2289658541588672D+01,&
        0.3446881631958411D-01, 0.1828647032406020D+01,&
        0.7152609366400937D-01, 0.1369709566513547D+01,&
        0.1202053639683183D+00, 0.9123014234783601D+00,&
        0.1639855806134269D+00, 0.4559006613049348D+00,&
        0.1818448921796476D+00, 0.0000000000000000D+00/
        DATA ( W(I,48), T(I,48), I = 1,24) / &
        0.4477155473876067D-35, 0.1269301231543958D+02,&
        0.3376456143135106D-30, 0.1175317846161662D+02,&
        0.2079058971418782D-26, 0.1097330095851353D+02,&
        0.3139476644799137D-23, 0.1027574158173351D+02,&
        0.1798854916237414D-20, 0.9630885686949055D+01,&
        0.4925465109698832D-18, 0.9023480280417658D+01,&
        0.7419993598065164D-16, 0.8444371322523445D+01,&
        0.6756677274657091D-14, 0.7887517316500372D+01,&
        0.3975805958471478D-12, 0.7348660565908911D+01,&
        0.1588360981246298D-10, 0.6824651320745980D+01,&
        0.4474287153436318D-09, 0.6313069914905501D+01,&
        0.9154055746313672D-08, 0.5811999987702533D+01,&
        0.1392791689559962D-06, 0.5319884620419087D+01,&
        0.1606393695542611D-05, 0.4835430885878957D+01,&
        0.1426609232424675D-04, 0.4357544108790293D+01,&
        0.9881804917611217D-04, 0.3885281117248394D+01,&
        0.5395865846271863D-03, 0.3417816048684776D+01,&
        0.2343082111760394D-02, 0.2954414688880120D+01,&
        0.8149696855351926D-02, 0.2494414743040025D+01,&
        0.2283821210009132D-01, 0.2037210303259033D+01,&
        0.5180518354956036D-01, 0.1582239319674964D+01,&
        0.9546340056142952D-01, 0.1128973231509649D+01,&
        0.1432824569931762D+00, 0.6769081422063483D+00,&
        0.1754635418118663D+00, 0.2255570729801638D+00/
        DATA ( W(I,49), T(I,49), I = 1,25) /&
        0.6811389123116583D-36, 0.1284021824817448D+02,&
        0.5655218285268972D-31, 0.1190391240788757D+02,&
        0.3783656387961836D-27, 0.1112712121198839D+02,&
        0.6170433961110509D-24, 0.1043245459795100D+02,&
        0.3805275010242669D-21, 0.9790392369513130D+01,&
        0.1119149962220974D-18, 0.9185730837014976D+01,&
        0.1808768646225968D-16, 0.8609346484896321D+01,&
        0.1765991900170034D-14, 0.8055219008757504D+01,&
        0.1113983833798787D-12, 0.7519106753949869D+01,&
        0.4771647040124784D-11, 0.6997872987021261D+01,&
        0.1441769692730330D-09, 0.6489109229832087D+01,&
        0.3166143293671931D-08, 0.5990909213118121D+01,&
        0.5175289227690812D-07, 0.5501725495034997D+01,&
        0.6419581525894787D-06, 0.5020274351287227D+01,&
        0.6139401305210235D-05, 0.4545470293840825D+01,&
        0.4586356164456607D-04, 0.4076379533242183D+01,&
        0.2705406889774859D-03, 0.3612185969516418D+01,&
        0.1271499268760875D-02, 0.3152165704938866D+01,&
        0.4796631268346652D-02, 0.2695667488564915D+01,&
        0.1461268197482865D-01, 0.2242097365824569D+01,&
        0.3612646737149015D-01, 0.1790906348772901D+01,&
        0.7276104720165333D-01, 0.1341580271678923D+01,&
        0.1197332834080540D+00, 0.8936312256129486D+00,&
        0.1613082628631473D+00, 0.4465901175404579D+00,&
        0.1781337719310835D+00, 0.0000000000000000D+00/
        DATA ( W(I,50), T(I,50), I = 1,25) /&
        0.1034607500576990D-36, 0.1298588445541555D+02,&
        0.9443414659584510D-32, 0.1205301838092448D+02,&
        0.6856280758924735D-28, 0.1127923332148262D+02,&
        0.1206044550761014D-24, 0.1058738174919177D+02,&
        0.7995094477915292D-22, 0.9948035709637500D+01,&
        0.2522482807168144D-19, 0.9346039593575728D+01,&
        0.4368171816201588D-17, 0.8772299579514598D+01,&
        0.4566698246800344D-15, 0.8220815907982127D+01,&
        0.3083828687005300D-13, 0.7687362406712500D+01,&
        0.1414228936126661D-11, 0.7168814837853899D+01,&
        0.4576636712310442D-10, 0.6662775399018720D+01,&
        0.1077060789389039D-08, 0.6167347388659921D+01,&
        0.1888225976835208D-07, 0.5680992291033284D+01,&
        0.2514609880838772D-06, 0.5202434993399912D+01,&
        0.2584937658949391D-05, 0.4730598550228594D+01,&
        0.2078485175734569D-04, 0.4264557843038109D+01,&
        0.1321726328668984D-03, 0.3803505741742012D+01,&
        0.6708280619787080D-03, 0.3346727774732429D+01,&
        0.2738160896935348D-02, 0.2893582727707738D+01,&
        0.9045054154849623D-02, 0.2443487452654017D+01,&
        0.2430481286424306D-01, 0.1995904709795124D+01,&
        0.5334352453170102D-01, 0.1550333214338771D+01,&
        0.9593054035810168D-01, 0.1106299289397183D+01,&
        0.1416854132499443D+00, 0.6633496795082918D+00,&
        0.1721258519924433D+00, 0.2210451816445435D+00/
!
!
!***END HERMIT
!


      if (restar .eq. 0) rule = 1
      neval= 0

 10   if (neval + rule**ndim .le. maxpts .and. rule .lt. maxrul) then
                do i = 1,rule/2
                        point(i) = -t(i,rule)
                        weight(i) =  w(i,rule)
                        point(rule-i+1) = t(i,rule)
                        weight(rule-i+1) = w(i,rule)
                end do

                if (mod(rule, 2) .eq. 1) then
                        point(rule/2 + 1 ) = 0
                        weight(rule/2 + 1 ) = w(rule/2 + 1,rule)
                end if

                call mltrul(ndim,numfun,funsub,rule,point,weight,result,work,&
                work(numfun+1),work(numfun+ndim+1))

                neval = neval + rule**ndim
                ifail = 0

                do i = 1,numfun

                        if (rule .gt. 1) then
                                abserr(i) = abs(result(i) - work(2*ndim+numfun+i))
                        else
                                abserr(i) = abs(result(i))
                        end if

                        work(2*ndim+numfun+i) = result(i)
                        if (abserr(i) .gt. max(espabs,epsrel*abs(result(i)))) then
                                 ifail = 1
                        end if
                end do

                rule = rule + 1
                if (ifail .gt. 0 .or. neval .lt. minpts) goto 10
        end if

      end subroutine hermit



        end module lois_normalesc





!===========================================================
!
      !SUBROUTINES
!







!===========================================================
!
      !SUBROUTINES
!
!===========================================================


!-----------------------------------------------------------
!                    SUBROUTINE PRINCIPALE: HETMIXCONT
!------------------------------------------------------------



      subroutine hetmixcont(Y0,X0,Prior0,idprob0,idea0,idg0,idcor0,ns0,ng0,nv0,nobs0 &
          ,nea0,nmes0,idiag0,nwg0,ncor0,npmtot0,btot,Vopt,vrais,ni,istop,gconv,ppi0,resid_m &
          ,resid_ss,pred_m_g,pred_ss_g,pred_RE,convB,convL,convG,maxiter0 &
          ,epsY0,idlink0,nbzitr0,zitr0,marker,transfY,nsim0,Yobs,Ydiscret,vraisdiscret,UACV,rlindiv,pbH0,fix0)

      use parameters
      use communc
      use optim

      IMPLICIT NONE

        !Declaration des variables en entree
      integer,intent(in)::nv0,maxiter0,Ydiscret
      double precision,intent(in)::epsY0
      integer,intent(in)::idlink0,nbzitr0
      double precision,dimension(nbzitr0),intent(in)::zitr0

      integer, intent(in)::ns0,ng0,nobs0,idiag0,nwg0,npmtot0,nea0,nsim0,ncor0
      integer, dimension(nv0),intent(in)::idea0,idg0,idprob0,idcor0
      integer, dimension(ns0),intent(in)::nmes0,prior0
      double precision,dimension(nobs0),intent(in)::Y0
      double precision,dimension(nobs0*nv0),intent(in)::X0
      double precision,intent(in)::convB,convL,convG
      integer,dimension(npmtot0),intent(in)::pbH0,fix0
          
        !Declaration des variable en entree et sortie
      double precision, dimension(npmtot0), intent(inout) :: btot
        !Declaration des variables en sortie
      double precision,intent(out)::vrais,vraisdiscret,UACV
      double precision,dimension(3),intent(out)::gconv
      double precision,dimension(ns0*ng0),intent(out)::ppi0
      double precision,dimension(ns0),intent(out)::rlindiv
      double precision,dimension(nobs0),intent(out)::resid_m,resid_ss,Yobs
      double precision,dimension(nobs0*ng0),intent(out)::pred_m_g
      double precision,dimension(nobs0*ng0),intent(out)::pred_ss_g
      double precision,dimension(ns0*nea0),intent(out)::pred_RE
      double precision,dimension(nsim0),intent(out)::marker,transfY
      double precision,dimension(npmtot0*(npmtot0+1)/2),intent(out)::Vopt
      integer, intent(out)::ni,istop
        !Variables locales
      integer::jtemp,i,g,j,ij,npm,ier,k,ktemp,ig,nmestot,it,id,nbfix
      double precision::eps,ca,cb,dd,thi
      double precision,dimension(ns0,ng0)::PPI
      double precision,dimension(npmtot0)::mvc,b
      double precision,dimension(npmtot0*(npmtot0+3)/2)::V
      double precision,external::funcpac

         !print*,"debut fortran"
! sorties initialisees

       ppi0=0.d0
       Vopt=0.d0
       gconv=0.d0
       pred_ss_g=0.d0
       pred_m_g=0.d0
       pred_RE=0.d0
       marker=0.d0
       transfY=0.d0
       resid_m=0.d0
       resid_ss=0.d0
       vrais=0.d0
       ni=0
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

      epsa=convB
      epsb=convL
      epsd=convG
      maxiter=maxiter0

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

      eps=1.d-20
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
         istop=12
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

              !print*,"npm=",npm

!       write(*,*)'idlink',idlink
!       write(*,*)'idea',idea
!       write(*,*)'NVC',nvc



      if (idiag.eq.1) then
         DO j=1,nvc
            btot(nef+j)=dsqrt(abs(btot(nef+j)))
         END DO
      end if

! si idiag=0, on met dans le vecteur des parms, les parms
! de la transformee de Cholesky

      if (idiag.eq.0) then

         DO j=1,nvc
            mvc(j)=btot(nef+j)
         END DO
            !print*,"avant dmfsd"
         CALL dmfsd(mvc,nea,EPS,IER)
         !print*,"apres dmfsd"
         DO j=1,nvc
            btot(nef+j)=mvc(j)
         END DO
      end if
      if (nwg.gt.0) then
         do i=1,nwg
            btot(nef+nvc+i)=abs(btot(nef+nvc+i))
         end do
      end if


      if (idlink.eq.2) then
         call design_splines(ier)
         if (ier.eq.-1) then
            istop=9
            go to 1589
         end if
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
      k=0
      do j=1,npmtot
        if(fix0(j).eq.0) then
           k=k+1
           pbH(k)=pbH0(j)
        end if
     end do



! lancement de l'optimisation

      IF (npm.eq.1) then
         istop=10
         go to 1589
      else
         ca=0.d0
         cb=0.d0
         dd=0.d0
         
!         write(*,*) "avant appel marq98 :"
!         write(*,*) "nef=",nef
!         write(*,*) "nprob+ncssg+ncg*ng+nvarxevt+nrisq-1=",nprob,"+",ncssg,"+",ncg,"*",ng,"+",nvarxevt,"+",nrisq,"-",1
!         write(*,*) "nvc=",nvc
!         write(*,*) "nwg=",nwg
!         write(*,*) "ntrtot=",ntrtot
 !        write(*,*) "ncor=",ncor
!         write(*,*) "npm=",npm
!         write(*,*) "npmtot=",npmtot
!         write(*,*) "idg0=",idg0
!         write(*,*) "idg=",idg   
!write(*,*) "B=",b,"bfix=",bfix  
!write(*,*) "Btot=",btot

         
         call marq98(b,npm,ni,V,vrais,ier,istop,ca,cb,dd,funcpac)

!         write(*,*)
!         write(*,*)'    FIN OPTIMISATION  ..... '
!         write(*,*)'istop',istop,'vrais',vrais,'b=',b

         gconv=0.d0
         gconv(1)=ca
         gconv(2)=cb
         gconv(3)=dd
         vopt(1:(npm*(npm+1)/2))=V(1:(npm*(npm+1)/2))

!  injecter le b estime dans btot
         k=0
         do j=1,npmtot
            if(fix0(j).eq.0) then
               k=k+1
               btot(j)=b(k)
            end if
         end do


         if (istop.eq.1.or.istop.eq.2.or.istop.eq.3) then  
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


!         else
!            ig=0
!            ij=0
!            do i=1,ns
!               do g=1,ng0
!                  ig=ig+1
!                  ppi0(ig)=0.d0
!               end do
!            end do

            if (Ydiscret.eq.1.and.ncor.eq.0) then
               id=0
               thi=0.d0
               call vrais_discret(b,npm,id,thi,id,thi,vraisdiscret)
 !              write(*,*)'vrais_discret',vraisdiscret

!               call vrais_cont(b,npm,id,thi,id,thi,vraisdiscret)
!               write(*,*)'vrais_cont',vraisdiscret

!               write(*,*)'avant UACV'
               call computUACV(b,npm,rlindiv,vopt,UACV)
!               write(*,*)'arpes UACV'
            end if

         end if

      end if


!      write(*,*)'avant deallocate'

 1589 continue

      deallocate(Y,X,idprob,idea,idg,idcor,nmes,prior)

      deallocate(zitr,mm,mm1,mm2,im,im1,im2)

      deallocate(pbH,fix,bfix)


!      write(*,*)'fin'
      return
    end subroutine hetmixcont



!-----------------------------------------------------------
!                        FUNCPA
!------------------------------------------------------------


      double precision function funcpac(b,npm,id,thi,jd,thj)

      use communc
      use optim
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
      use optim
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

      use optim

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






! =================================================================
!              Densite d'une beta
!=================================================================

      double precision Function beta_densite(X,a,b)

      implicit none

      double precision :: beta,a,b,gammln,X

      beta=exp(gammln(a+b)-gammln(a)-gammln(b))

      beta_densite=((X)**(a-1))*((1-X)**(b-1))*beta

      return

      end Function beta_densite


! Calcul de gamma (gammln(a))

      double precision Function gammln(xx)

! retourne la valeur ln(gamma(xx)) pour xx>0
      implicit none

      integer::j
      double precision:: ser,stp,tmp,x,y,xx
      double precision,dimension(6)::cof

      save cof,stp

      data cof,stp/76.18009172947146d0,-86.50532032941677d0, &
      24.01409824083091d0, -1.231739572450155d0,.1208650973866179d-2, &
      -.5395239384953d-5,2.5066282746310005d0/

      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
         y=y+1.d0
         ser=ser+cof(j)/y
      end do
      gammln=tmp+log(stp*ser/x)

      return

      end Function gammln






! =================================================================
!              CDF incomplete d'une beta
!=================================================================


      double precision Function betai(a,b,x)

      implicit none

      double precision :: a,b,x,bt,betaCF,gammln,temp

      if(x.lt.0.d0.or.x.gt.1.d0) then
            betai=999.d0
            return
      end if
      if (x.eq.0..or.X.eq.1.) then
         bt=0.
      else
         bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
      end if

      if (x.lt.(a+1.)/(a+B+2.)) then
         temp=betaCF(a,b,x)
         if (temp.eq.999.d0) then
            betai=999.d0
            return
         end if
         betai=bt*temp/A
         return
      else
        temp=betaCF(b,a,1.-x)
        if (temp.eq.999.d0) then
            betai=999.d0
            return
         end if
         betai=1.-bt*temp/b
         return
      end if

      end Function betai

! betaCF est utilisï¿½ par betai

      double precision Function betaCF(a,b,x)

      implicit none

      integer ::m,m2
      integer,parameter ::maxit=100
      double precision ::a,b,x,aa,c,d,del,h,qab,qam,qap
      double precision,parameter::eps=3.e-7,fpmin=1.e-30


      qab=a+b
      qap=a+1
      qam=a-1
      c=1.       ! first step of Lentz's method
      d=1.-qab*x/qap
      if (abs(d).lt.fpmin) d=fpmin
      d=1./d
      h=d
      do m=1,maxit
         m2=2*m
         aa=m*(b-m)*x/((qam+m2)*(a+m2))
         d=1.+aa*d        ! one step (the even one) of the recurrence
         if (abs(d).lt.fpmin) d=fpmin
         c=1.+aa/c
         if (abs(c).lt.fpmin) c=fpmin
         d=1./d
         h=h*d*c
         aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
         d=1.+aa*d        ! next step of the recurrence (the odd one)
         if (abs(d).lt.fpmin) d=fpmin
         c=1.+aa/c
         if (abs(c).lt.fpmin) c=fpmin
         d=1./d
         del=d*c
         h=h*del
         if (abs(del-1.).lt.eps) goto 1
      end do

!  pause 'a or b too big, or maxit too small in betaCF'
            betaCF=999.d0
!            write(*,*)'problem'
            return

 1    betaCF=h
      return
      end Function betaCF




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
        use optim
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

        use lois_normalesc
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
              call gaussher(gauss,npg)

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
      use optim
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
