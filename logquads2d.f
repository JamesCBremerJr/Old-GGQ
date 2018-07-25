c       Construct two sets of  tensor product quadrature rules, one for 
c       discretizing integral operators with kernels of the form
c
c             log(|x-y|) \psi(y) + \phi(y)
c
c       with \psi and \phi smooth, and a second set for integral operators
c       with kernels of the form
c
c            \Grad_x log|x-y| \cdot F(y) | phi(y)
c
        implicit double precision (a-h,o-z)
        double precision, allocatable :: xsdisc(:),ysdisc(:),whtsdisc(:)
        double precision, allocatable :: xsbdy(:),ysbdy(:),whtsbdy(:)
        double precision, allocatable :: xsall(:,:),ysall(:,:)
        double precision, allocatable :: xslege(:),whtslege(:)

        double precision, allocatable :: whtsall(:,:)
        double precision, allocatable :: xs(:),ys(:),whts(:)
        integer, allocatable          :: nquadall(:)
        double precision              :: amatr(2,2)
c
        eps       = 1.0d-20
        epsadap   = 1.0d-20
        erraccept = 1.0d-20
        ifcheck   = 0
c
        ndegree = 8           ! degree of the discretization polynomials
        npoly   = 12          ! degree of polynomials to integrate
        nmax    = 999         ! maximum possible singular quadrature size
!
        nlege   = 12

        allocate(xslege(nlege),whtslege(nlege))
        call legequad(nlege,xslege,whtslege)
c     
c       Fetch the discretization quadrature rule for the interior
c
        call discquad_info(ndisc)
        allocate(xsdisc(ndisc),ysdisc(ndisc),whtsdisc(ndisc))
        call discquad(ndisc,xsdisc,ysdisc,whtsdisc)

        call prinf("ndisc    = *",ndisc,1)
        call prin2("xsdisc   = *",xsdisc,ndisc)
        call prin2("ysdisc   = *",ysdisc,ndisc)
        call prin2("whtsdisc = *",whtsdisc,ndisc)
c
c       Construct the boundary quadrature rule
c
        nquadbdy = nlege*4
        allocate(xsbdy(nquadbdy),ysbdy(nquadbdy),whtsbdy(nquadbdy))

        idx = 0

        do i=1,nlege
        idx          = idx + 1
        xsbdy(idx)   = -1.0d0
        ysbdy(idx)   = xslege(i)
        whtsbdy(idx) = whtslege(i)
        end do

        do i=1,nlege
        idx          = idx + 1
        xsbdy(idx)   = xslege(i)
        ysbdy(idx)   = 1.0d0
        whtsbdy(idx) = whtslege(i)
        end do

        do i=1,nlege
        idx = idx + 1
        xsbdy(idx)   = 1.0d0
        ysbdy(idx)   = xslege(i)
        whtsbdy(idx) = whtslege(i)
        end do
c
        do i=1,nlege
        idx = idx + 1
        xsbdy(idx)   = xslege(i)
        ysbdy(idx)   = -1.0d0
        whtsbdy(idx) = whtslege(i)
        end do
c
        call prinf("nquadbdy = *",nquadbdy,1)
        call prin2("xsbdy    = *",xsbdy,nquadbdy)
        call prin2("ysbdy    = *",ysbdy,nquadbdy)
        call prin2("whtsbdy  = *",whtsbdy,nquadbdy)
c
c       Write out the discretization quadrature rule 
c
        iw = 1001
        open(iw,FILE='quads.f90')
c
 0100 format(A)
 0150 format("double precision :: ",A,"(",I3.3,",",I3.3,")")
 0160 format("double precision :: ",A,"(",I3.3,")")
 0175 format("integer          :: ",A,"(",I3.3,")")
 0300 format(A,I3.3)
 0400 format("xs(",I3.3,")   = ",D44.36)
 0500 format("ys(",I3.3,")   = ",D44.36)
 0600 format("whts(",I3.3,") = ",D44.36)
c
        write(iw,0100) "subroutine discquad(nquad,xs,ys,whts)"
        write(iw,0100) "implicit double precision (a-h,o-z)"
        write(iw,0100) "integer :: nquad"
        write(iw,0160) "xs",ndisc
        write(iw,0160) "ys",ndisc
        write(iw,0160) "whts",ndisc
        write(iw,0300) "nquad     = ",ndisc

        do i=1,ndisc
        write(iw,0400) i,xsdisc(i)
        end do
c
        do i=1,ndisc
        write(iw,0500) i,ysdisc(i)
        end do
c
        do i=1,ndisc
        write(iw,0600) i,whtsdisc(i)
        end do
c
        write(iw,0100) "end subroutine"
        write(iw,*)    ""
        write(iw,*)    ""
c
        write(iw,0100) "subroutine bdyquad(nquad,xs,ys,whts)"
        write(iw,0100) "implicit double precision (a-h,o-z)"
        
        write(iw,0160) "xs",nquadbdy
        write(iw,0160) "ys",nquadbdy
        write(iw,0160) "whts",nquadbdy
        write(iw,0100) "integer          :: nquad"

        write(iw,0300) "nquad     = ",nquadbdy

        do i=1,nquadbdy
        write(iw,0400) i,xsbdy(i)
        end do
c
        do i=1,nquadbdy
        write(iw,0500) i,ysbdy(i)
        end do
c
        do i=1,nquadbdy
        write(iw,0600) i,whtsbdy(i)
        end do
c
        write(iw,0100) "end subroutine"
        write(iw,*)    ""
        write(iw,*)    ""

        call flush(iw)

        close(iw)
c
        allocate(xsall(nmax,ndisc+nquadbdy),ysall(nmax,ndisc+nquadbdy))
        allocate(whtsall(nmax,ndisc+nquadbdy))
        allocate(nquadall(ndisc+nquadbdy))
c
 1100 format("nquadsall(",I3.3,")    = ",I3.3)
 1200 format("xsall(",I3.3,",",I3.3,")   = ",D44.36)
 1300 format("ysall(",I3.3,",",I3.3,")   = ",D44.36)
 1400 format("whtsall(",I3.3,",",I3.3,") = ",D44.36)
 2100 format(A," = ",I3.3)
c
c       Build the first set of singular quadrature rules (for log 
c       singularities)
c
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(xs,ys,whts,i,x1,x2,ier,nquad,
!$OMP!    t1,t2,errmax)
        allocate(xs(10000),ys(10000),whts(10000))
!$OMP   DO
        do i=1,ndisc
        x1    = xsdisc(i)
        x2    = ysdisc(i)
        call prini("idisc = ",i)
        call elapsed(t1)
        call diagquad(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
        call elapsed(t2)
        call prin2("time = *",t2-t1,1)

        if (nquad .gt. nmax) then
        print *,"QUADRATURE SIZED EXCEEDED NMAX"
        stop
        endif

        xsall(1:nquad,i)   = xs(1:nquad)
        ysall(1:nquad,i)   = ys(1:nquad)
        whtsall(1:nquad,i) = whts(1:nquad)
        nquadall(i)        = nquad
        end do
!$OMP   END DO
!$OMP   END PARALLEL
c
c
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(xs,ys,whts,i,x1,x2,ier,nquad,
!$OMP!    t1,t2)
        allocate(xs(10000),ys(10000),whts(10000))
!$OMP   DO
        do i=1,nlege
        x1    = xsbdy(i)
        x2    = ysbdy(i)
        call elapsed(t1)
        call diagquad(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
        call elapsed(t2)
        call prin2("time = *",t2-t1,1)
        if (nquad .gt. nmax) then
        print *,"QUADRATURE SIZED EXCEEDED NMAX"
        stop
        endif

        xsall(1:nquad,i+ndisc)   = xs(1:nquad)
        ysall(1:nquad,i+ndisc)   = ys(1:nquad)
        whtsall(1:nquad,i+ndisc) = whts(1:nquad)
        nquadall(i+ndisc)        = nquad
        end do
!$OMP   END DO
!$OMP   END PARALLEL
c
c       Rotate the boundary quadratures to construct the others
c
        amatr(1,1) = 0
        amatr(1,2) = 1
        amatr(2,1) = -1
        amatr(2,2) = 0

        do i=1,nlege
        nq = nquadall(ndisc+i)
        nquadall(ndisc+i+nlege) = nq

        do j=1,nq
        x = xsall(j,i+ndisc)
        y = ysall(j,i+ndisc)
        wht = whtsall(j,i+ndisc)

        xsall(j,i+ndisc+nlege)   = amatr(1,1)*x + amatr(1,2)*y
        ysall(j,i+ndisc+nlege)   = amatr(2,1)*x + amatr(2,2)*y
        whtsall(j,i+ndisc+nlege) = wht
        end do
        end do
c
        amatr(1,1) = -1
        amatr(1,2) = 0
        amatr(2,1) = 0
        amatr(2,2) = 1
c
        do i=1,nlege
        nq = nquadall(ndisc+i)
        nquadall(ndisc+i+2*nlege) = nq
        do j=1,nq
        x = xsall(j,i+ndisc)
        y = ysall(j,i+ndisc)
        wht = whtsall(j,i+ndisc)
        xsall(j,i+ndisc+2*nlege)   = amatr(1,1)*x + amatr(1,2)*y
        ysall(j,i+ndisc+2*nlege)   = amatr(2,1)*x + amatr(2,2)*y
        whtsall(j,i+ndisc+2*nlege) = wht
        end do
        end do
c
c
        amatr(1,1) = 0
        amatr(1,2) = 1
        amatr(2,1) = 1
        amatr(2,2) = 0

        do i=1,nlege
        nq = nquadall(ndisc+i)
        nquadall(ndisc+i+3*nlege) = nq
        do j=1,nq
        x   = xsall(j,i+ndisc)
        y   = ysall(j,i+ndisc)
        wht = whtsall(j,i+ndisc)

        xsall(j,i+ndisc+3*nlege)   = amatr(1,1)*x + amatr(1,2)*y
        ysall(j,i+ndisc+3*nlege)   = amatr(2,1)*x + amatr(2,2)*y
        whtsall(j,i+ndisc+3*nlege) = wht
        end do
        end do
c
c       Check the quadrature rules
c
        if (ifcheck .eq. 1) then
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(i,x1,x2,nq,errmax)
!$OMP   DO
        do i=1,ndisc
        x1    = xsdisc(i)
        x2    = ysdisc(i)
        nq    = nquadall(i)
        call diagquad_check(epsadap,x1,x2,npoly,nq,
     -    xsall(:,i),ysall(:,i),whtsall(:,i),errmax)

        if (errmax .gt. erraccept) then
        call prin2("x1 = *",x1,1)
        call prin2("x2 = *",x2,1)
        call prin2("errmax = *",errmax,1)
        call prina("diagquad error too high*")
        endif
        end do
!$OMP END DO
!$OMP END PARALLEL
c

c
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(i,x1,x2,nq,errmax)
!$OMP   DO
        do i=1,nquadbdy
        x1    = xsbdy(i)
        x2    = ysbdy(i)
        nq = nquadall(ndisc+i)
        call diagquad_check(epsadap,x1,x2,npoly,nq,
     -    xsall(:,i+ndisc),ysall(:,i+ndisc),whtsall(:,i+ndisc),errmax)
c
        if (errmax .gt. erraccept) then
        call prin2("x1 = *",x1,1)
        call prin2("x2 = *",x2,1)
        call prin2("errmax = *",errmax,1)
        call prina("diagquad bdy error too high*")
        endif
        end do
!$OMP END DO
!$OMP END PARALLEL
c
        endif

c
c       Write the first set of singular rules to the disc
c
        iw = 1001
        open(iw,FILE='quads.f90',STATUS='OLD',ACCESS='append')

        write(iw,0100) "subroutine singquads(nquadsall,xsall," //
     -    "ysall,whtsall)"
        write(iw,0100) "implicit double precision (a-h,o-z)"
        write(iw,0150) "xsall",nmax,ndisc+nquadbdy
        write(iw,0150) "ysall",nmax,ndisc+nquadbdy
        write(iw,0150) "whtsall",nmax,ndisc+nquadbdy
        write(iw,0175) "nquadsall",ndisc+nquadbdy
c

        do i=1,ndisc+nquadbdy
        nq = nquadall(i)
        write(iw,1100) i,nq
        do j=1,nq
        write(iw,1200) j,i,xsall(j,i)
        write(iw,1300) j,i,ysall(j,i)
        write(iw,1400) j,i,whtsall(j,i)
        end do
        end do
c
        write(iw,0100) "end subroutine"
        write(iw,*)    ""
        write(iw,*)    ""

        close(iw)
        call flush(iw)
c
c       Construct the second set of singular rules
c
        nn = ndisc
c
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(xs,ys,whts,i,x1,x2,ier,nquad)
        allocate(xs(10000),ys(10000),whts(10000))
!$OMP   DO
        do i=1,ndisc
        x1    = xsdisc(i)
        x2    = ysdisc(i)
        call diagquad2(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
        if (nquad .gt. nmax) then
        print *,"QUADRATURE SIZED EXCEEDED NMAX"
        stop
        endif

        xsall(1:nquad,i)   = xs(1:nquad)
        ysall(1:nquad,i)   = ys(1:nquad)
        whtsall(1:nquad,i) = whts(1:nquad)
        nquadall(i)        = nquad
        end do
!$OMP   END DO
!$OMP   END PARALLEL

c
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(xs,ys,whts,i,x1,x2,ier,nquad,
!$OMP!    t1,t2)
        allocate(xs(10000),ys(10000),whts(10000))
!$OMP   DO
        do i=1,nlege
        x1    = xsbdy(i)
        x2    = ysbdy(i)
        call elapsed(t1)
        call diagquad2(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
        call elapsed(t2)
        call prin2("time = *",t2-t1,1)
        if (nquad .gt. nmax) then
        print *,"QUADRATURE SIZED EXCEEDED NMAX"
        stop
        endif
c
        xsall(1:nquad,i+ndisc)   = xs(1:nquad)
        ysall(1:nquad,i+ndisc)   = ys(1:nquad)
        whtsall(1:nquad,i+ndisc) = whts(1:nquad)
        nquadall(i+ndisc)        = nquad
        end do
!$OMP   END DO
!$OMP   END PARALLEL
c
c       Rotate the boundary quadratures to construct the others
c
        amatr(1,1) = 0
        amatr(1,2) = 1
        amatr(2,1) = -1
        amatr(2,2) = 0

        do i=1,nlege
        nq = nquadall(ndisc+i)
        nquadall(ndisc+i+nlege) = nq

        do j=1,nq
        x = xsall(j,i+ndisc)
        y = ysall(j,i+ndisc)
        wht = whtsall(j,i+ndisc)

        xsall(j,i+ndisc+nlege)   = amatr(1,1)*x + amatr(1,2)*y
        ysall(j,i+ndisc+nlege)   = amatr(2,1)*x + amatr(2,2)*y
        whtsall(j,i+ndisc+nlege) = wht
        end do
        end do
c
        amatr(1,1) = -1
        amatr(1,2) = 0
        amatr(2,1) = 0
        amatr(2,2) = 1
c
        do i=1,nlege
        nq = nquadall(ndisc+i)
        nquadall(ndisc+i+2*nlege) = nq
        do j=1,nq
        x = xsall(j,i+ndisc)
        y = ysall(j,i+ndisc)
        wht = whtsall(j,i+ndisc)
        xsall(j,i+ndisc+2*nlege)   = amatr(1,1)*x + amatr(1,2)*y
        ysall(j,i+ndisc+2*nlege)   = amatr(2,1)*x + amatr(2,2)*y
        whtsall(j,i+ndisc+2*nlege) = wht
        end do
        end do
c
c
        amatr(1,1) = 0
        amatr(1,2) = 1
        amatr(2,1) = 1
        amatr(2,2) = 0

        do i=1,nlege
        nq = nquadall(ndisc+i)
        nquadall(ndisc+i+3*nlege) = nq
        do j=1,nq
        x   = xsall(j,i+ndisc)
        y   = ysall(j,i+ndisc)
        wht = whtsall(j,i+ndisc)

        xsall(j,i+ndisc+3*nlege)   = amatr(1,1)*x + amatr(1,2)*y
        ysall(j,i+ndisc+3*nlege)   = amatr(2,1)*x + amatr(2,2)*y
        whtsall(j,i+ndisc+3*nlege) = wht
        end do
        end do
c
c      Check the second set of diagonal quadratures
c
        if (ifcheck .eq. 1) then
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(i,x1,x2,nq,errmax)
!$OMP   DO
        do i=1,ndisc
        x1    = xsdisc(i)
        x2    = ysdisc(i)
        nq    = nquadall(i)
        call diagquad2_check(epsadap,x1,x2,npoly,nq,
     -    xsall(:,i),ysall(:,i),whtsall(:,i),errmax)
c
        if (errmax .gt. erraccept) then
        call prin2("x1 = *",x1,1)
        call prin2("x2 = *",x2,1)
        call prin2("errmax = *",errmax,1)
        call prina("diagquad2 error too high*")
        endif

        end do
!$OMP   END DO
!$OMP END PARALLEL

c
c
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(i,x1,x2,nq,errmax)
!$OMP   DO

        do i=1,nquadbdy
        x1    = xsbdy(i)
        x2    = ysbdy(i)
        nq = nquadall(ndisc+i)
        call diagquad2_check(epsadap,x1,x2,npoly,nq,
     -    xsall(:,i+ndisc),ysall(:,i+ndisc),whtsall(:,i+ndisc),errmax)
c
        if (errmax .gt. erraccept) then
        call prin2("x1 = *",x1,1)
        call prin2("x2 = *",x2,1)
        call prin2("errmax = *",errmax,1)
        call prina("diagquad2 bdy error too high*")
        endif

        end do
!$OMP   END DO
!$OMP END PARALLEL
c
        endif
c
c       Write the second set of singular rules to the disc
c
        iw = 1001
        open(iw,FILE='quads.f90',STATUS='OLD',ACCESS='append')

        write(iw,0100) "subroutine singquads2(nquadsall,xsall," //
     -    "ysall,whtsall)"
        write(iw,0100) "implicit double precision (a-h,o-z)"
        write(iw,0150) "xsall",nmax,ndisc+nquadbdy
        write(iw,0150) "ysall",nmax,ndisc+nquadbdy
        write(iw,0150) "whtsall",nmax,ndisc+nquadbdy
        write(iw,0175) "nquadsall",ndisc+nquadbdy
c

        do i=1,ndisc+nquadbdy
        nq = nquadall(i)
        write(iw,1100) i,nq
        do j=1,nq
        write(iw,1200) j,i,xsall(j,i)
        write(iw,1300) j,i,ysall(j,i)
        write(iw,1400) j,i,whtsall(j,i)
        end do
        end do
c
        write(iw,0100) "end subroutine"
        write(iw,*)    ""
        write(iw,*)    ""

        close(iw)
        call flush(iw)
c
        iw = 1001
        open(iw,FILE='quads.f90',STATUS='OLD',ACCESS='append')

        write(iw,0100) "subroutine discquad_info(ndegree," //
     -    "nquad,nlege,naux,nmax)"
        write(iw,0100) "implicit double precision (a-h,o-z)"
        write(iw,2100) "ndegree",ndegree
        write(iw,2100) "nquad",ndisc
        write(iw,2100) "nlege",nlege
        write(iw,2100) "naux",npoly
        write(iw,2100) "nmax",nmax
        write(iw,0100) "end subroutine"
        write(iw,*)    ""
        write(iw,*)    ""
        close(iw)
c
        end program



        subroutine diagquad(ier,eps,npoly,x1,x2,
     -    nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension                     :: ab0(2,100)
        double precision, allocatable :: disc(:),coefs(:,:)
        double precision, allocatable :: xs0(:),whts0(:),rints(:)
        double precision, allocatable :: xs1(:),whts1(:),ab(:,:)
        double precision, allocatable :: xs2(:),whts2(:),errs(:)
        double precision              :: xs(1),ys(1),whts(1)
        double precision              :: verts(2,3)
        integer, allocatable          :: idxs(:)
        external funuser1,funuser2,funeval,funquad1,funquad2
c
c       Build a tensor product quadrature on [-1,1] x [-1,1] for integrals of
c       the form
c
c           1     1
c       \int  \int  (  log( (x1-y1)^2 + (x2-y2)^2 )  P_i1(y1) P_j1(y2)  +
c           -1    1
c
c                 P_i2(y1) P_j2(y2) ) dy1 dy2
c
c       where (x1,x2) is a  user-specified point in [-1,1] x [-1,1], P_k
c       denotes the Legendre polynomial of degree k and
c
c           0 <= i1+j1 < = npoly,   0 <= i2+j2 <= npoly.
c
c       Input parameters:
c         npoly - the degree of the polynomials 
c         (x1,x2) - the target node
c
c       Output parameters:
c         ier - an error return code
c             ier   =  0 indicates successful execution
c             ier   =  
c         (nquad,xs,ys,whts) - the output quadrature rule
c 
        ier = 0

c
c       First build a quadrature for integrating the functions F_ij defined via
c
c                                  1
c         F_ij(y2) = P_j2(y2) \int   log ( (x1-y1)^2 + (x2-y2)^2 ) P_j(y1) dy1
c                                 -1
c
c       and the functions P_j(y2),  where 0<= i1+j1, j2 <= npoly and P_n denotes the 
c       Legendre polynomial of degree n.
c
        call prin2("in diagquad, x1 = *",x1,1)
        call prin2("in diagquad, x2 = *",x2,1)
        call prinf("in diagquad, npoly = *",npoly,1)
c
        ldisc = 100 000 000
        allocate(disc(ldisc))
        allocate(xs1(1000),whts1(1000),rints(1000),idxs(1000))
        allocate(xs2(1000),whts2(1000))
c
        kdisc = 30
        nquad = 0
 0100 continue
        a      = -1.0d0-x2
        b      =  1.0d0-x2
        nfuns1 = (npoly+1)*(npoly+2)/2 + npoly+1
c
        if (a .eq. 0 .OR. b .eq. 0) then
        nints0    = 1
        ab0(1,1)  = a
        ab0(2,1)  = b
        else
        nints0    = 2
        ab0(1,1)  = a
        ab0(2,1)  = 0
c
        ab0(1,2)  = 0
        ab0(2,2)  = b
        endif

        kdisc     = kdisc+1
        ifadap    = 1
c
        call prinf("nfuns1 = *",nfuns1,1)
        call legedisc(ier,ifadap,nints0,ab0,kdisc,eps,nfuns1,
     -    funuser1,npoly,x1,x2,eps,disc,ldisc,ncoefs,lkeep)
        if (ier .ne. 0) then
        call prinf("after legedisc, ier = *",ier,1)
        stop
        endif
c
c       Gram-Schmidt the input functions.
c
        allocate(coefs(ncoefs,nfuns1))
c
        call legegs(eps,disc,nfuns1,funuser1,npoly,x1,x2,eps,
     1    krank,coefs)
        call prinf("after legegs, krank=*",krank,1)
c
c
c       Fetch the initial, oversampled, quadrature rule.
c
        allocate(xs0(ncoefs),whts0(ncoefs))
        call legedisc_quad(disc,nquad0,xs0,whts0)
c
c       Compute the Chebyshev quadrature.
c
        call chebquad(krank,funeval,disc,coefs,krank,par4,
     1    nquad0,xs0,whts0,nquad1,xs1,whts1,rints)
c
c
c       Build the Gaussian quadrature.
c
        ifaccept = 0
        ngoal    = 0
c
        call gaussquad(eps,krank,rints,funeval,disc,coefs,krank,
     1    par4,nquad1,xs1,whts1,a,b,ngoal,ifaccept)
c
        deallocate(coefs,xs0,whts0)
c
        if (ier .ne. 0) goto 0100
c        ngoal = (krank+1)/2+1
c        if (nquad1 .gt. ngoal) goto 0100
        do i=1,nquad1
        if (whts1(i) .le. 0) goto 0100
        if (xs1(i)   .lt. a .OR. xs1(i) .gt. b) goto 0100
        end do

c
        do i=1,nquad1
        idxs(i) = i
        end do
        call quicksort(nquad1,xs1,idxs)
        whts1(1:nquad1)  = whts1(idxs(1:nquad1))
        xs1(1:nquad1) = xs1(1:nquad1)
        call prin2("xs1    = *",xs1,nquad1)
        call prin2("whts1  = *",whts1,nquad1)
c
c       For each point y2 in the preceeding quadrature rule, build
c       a quadrature for the functions
c
c             1
c        \int    ( log ( (x1-y1)^2 + (x2-y2)^2 ) P_i(y1)  ) dy1 
c            -1
c
c       where 0<= i <= npoly.
c
        do idisc=1,nquad1
        y2    = xs1(idisc)
        kdisc = 30
 0200 continue
        a      = -1.0d0-x1
        b      =  1.0d0-x1
c
        if (a .eq. 0 .OR. b .eq. 0) then
        nints0    = 1
        ab0(1,1)  = a
        ab0(2,1)  = b
        else
        nints0    = 2
        ab0(1,1)  = a
        ab0(2,1)  = 0
        ab0(1,2)  = 0
        ab0(2,2)  = b
        endif

        kdisc     = kdisc+1
        nfuns2    = (npoly+1)*2
        ifadap    = 1
c
        call legedisc(ier,ifadap,nints0,ab0,kdisc,eps,nfuns2,
     -    funuser2,npoly,x1,x2,y2,disc,ldisc,ncoefs,lkeep)
c
c       Gram-Schmidt the input functions.
c
        allocate(coefs(ncoefs,nfuns2))
c
        call prinf("before legegs, nfuns2 = *",nfuns2,1)
        call legegs(eps,disc,nfuns2,funuser2,npoly,x1,x2,y2,
     1    krank,coefs)
        call prinf("after legegs, krank=*",krank,1)
c
c
c       Fetch the initial, oversampled, quadrature rule.
c
        allocate(xs0(ncoefs),whts0(ncoefs))
        call legedisc_quad(disc,nquad0,xs0,whts0)
c
c       Compute the Chebyshev quadrature.
c
        call chebquad(krank,funeval,disc,coefs,krank,par4,
     1    nquad0,xs0,whts0,nquad2,xs2,whts2,rints)
c
c       Build the Gaussian quadrature.
c
        ifaccept = 0
        ngoal    = 0
c
        call gaussquad(eps,krank,rints,funeval,disc,coefs,krank,
     1    par4,nquad2,xs2,whts2,a,b,ngoal,ifaccept)
        deallocate(coefs,xs0,whts0)
c
        if (ier .ne. 0) goto 0200
c        ngoal = (krank+1)/2+1
c        if (nquad2 .gt. ngoal) goto 0100
        do i=1,nquad2
        if (whts2(i) .le. 0) goto 0200
        if (xs2(i)   .lt. a .OR. xs2(i) .gt. b) goto 0200
        end do
c
        do i=1,nquad2
        idxs(i) = i
        end do
        call quicksort(nquad2,xs2,idxs)
        whts2(1:nquad2)  = whts2(idxs(1:nquad2))

        call prin2("xs2    = *",xs2,nquad2)
        call prin2("whts2  = *",whts2,nquad2)
c
        do i=1,nquad2
        nquad       = nquad+1
        xs(nquad)   = xs2(i)       + x1
        ys(nquad)   = xs1(idisc)   + x2
        whts(nquad) = whts1(idisc)*whts2(i)
        end do
c
        end do
c
        call prina("*")
        call prinf("final nquad = *",nquad,1)
        call prin2("final xs    = *",xs,nquad)
        call prin2("final ys    = *",ys,nquad)
        call prin2("final whts  = *",whts,nquad)

        end



        subroutine diagquad_check(epsadap,x1,x2,npoly,nquad,xs,ys,whts,
     -    errmax)
        implicit double precision (a-h,o-z)
        dimension xs(nquad),ys(nquad),whts(nquad)
        dimension verts(2,3)
        double precision, allocatable :: errs(:)
        external funquad1,funquad2
c
c       Thouroughly test one of the diagonal quadrature rules
c
        allocate(errs( (npoly+1)*(npoly+2) ) )
        idx = 0

        do nn1=0,npoly
        do i1=0,nn1
        j1 = nn1-i1
c
        verts(1,1) = -1
        verts(2,1) = -1
        verts(1,2) =  1
        verts(2,2) = -1
        verts(1,3) =  1
        verts(2,3) =  1
c
        call adaptri(ier,epsadap,verts,funquad1,x1,x2,i1,j1,
     1    val1,nquad01)
c
        verts(1,1) = -1
        verts(2,1) = -1
        verts(1,2) = -1
        verts(2,2) =  1
        verts(1,3) =  1
        verts(2,3) =  1
c
        call adaptri(ier,epsadap,verts,funquad1,x1,x2,i1,j1,
     1    val2,nquad02)
c
        val0  = val1+val2
c
        val  = 0
        do i=1,nquad
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
        call funquad1(x,y,x1,x2,i1,j1,val00)
        val = val + wht*val00
        end do
c
        idx = idx + 1
        errs(idx) = abs(val-val0)
c
        end do
        end do
c
        do nn1=0,npoly
        do i1=0,nn1
        j1 = nn1-i1
c
        verts(1,1) = -1
        verts(2,1) = -1
        verts(1,2) =  1
        verts(2,2) = -1
        verts(1,3) =  1
        verts(2,3) =  1
c
        call adaptri(ier,epsadap,verts,funquad2,x1,x2,i1,j1,
     1    val1,nquad01)
c
        verts(1,1) = -1
        verts(2,1) = -1
        verts(1,2) = -1
        verts(2,2) =  1
        verts(1,3) =  1
        verts(2,3) =  1
c
        call adaptri(ier,epsadap,verts,funquad2,x1,x2,i1,j1,
     1    val2,nquad02)
c
        val0  = val1+val2
c
        val  = 0
        do i=1,nquad
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
        call funquad2(x,y,x1,x2,i1,j1,val00)
        val = val + wht*val00
        end do
c
        idx = idx + 1
        errs(idx) = abs(val-val0)

c
        end do
        end do
c
        errmax = maxval(errs)
        call prin2("x1 = *",x1,1)
        call prin2("x2 = *",x2,1)
        call prin2("errs = *",errs,idx)
        call prin2("errmax = *",errmax,1)
        call prina("*")
c
        end subroutine



        subroutine funeval(t,disc,coefs,nfuns,par4,vals,ders,ifders)
        implicit double precision (a-h,o-z)
        dimension vals(1),ders(1)
        call leged_eval(disc,nfuns,coefs,t,vals,ders,ifders)
        end


        subroutine funquad0(y1,val,x1,x2,y2,ideg,par5,par6,par7,
     -    par8)
        implicit double precision (a-h,o-z)
        call lege0(ideg,x1+y1,pol,der)
        dd  = log( (y1)**2 + (y2)**2)
        val = dd*pol
        end subroutine


        subroutine funuser1(y2,vals,npoly,x1,x2,eps00)
        implicit double precision (a-h,o-z)
        double precision vals(1),pols(0:100)
        external funquad0
c
        a    = -1.0d0-x1
        b    =  1.0d0-x1
        m    = 30
        ideg = 0
c     
        eps = eps00/20
c
        call lege(npoly+1,y2+x2,pols)
c
        idx  = 0 
        do ndeg=0,npoly
        do ideg1=0,ndeg
        ideg2 = ndeg-ideg1
        call adapgauss(ier,a,b,eps,m,val,funquad0,x1,x2,y2,
     1    ideg1,par5,par6,par7,par8)
        idx       = idx + 1
        vals(idx) = val*pols(ideg2)
        end do
        end do

        do ideg=0,npoly
        idx = idx + 1
        vals(idx) = pols(ideg)
        end do

        return
        end subroutine



        subroutine funuser2(y1,vals,npoly,x1,x2,y2)
        implicit double precision (a-h,o-z)
        dimension ns(1),ys2(1),pols(0:100),vals(1)
        external funquad1
c
c
        call lege(npoly+1,y1+x1,pols)
c
        dd = log((y1)**2 + (y2)**2)
        idx   = 0
c
        do i=0,npoly
        idx       = idx + 1
        vals(idx) = pols(i)*dd
        end do
c
        do i=0,npoly
        idx       = idx + 1
        vals(idx) = pols(i)
        end do
c
        end subroutine


        subroutine funquad1(y1,y2,x1,x2,i1,j1,val)
        implicit double precision (a-h,o-z)
        dimension pols1(0:100)
        dimension pols2(0:100)
        call lege(i1,y1,pols1)
        call lege(j1,y2,pols2)
        val = log( (x1-y1)**2 + (x2-y2)**2) * pols1(i1)*pols2(j1)
        end subroutine


        subroutine funquad2(y1,y2,x1,x2,i1,j1,val)
        implicit double precision (a-h,o-z)
        dimension pols1(0:100)
        dimension pols2(0:100)
        call lege(i1,y1,pols1)
        call lege(j1,y2,pols2)
        val = pols1(i1)*pols2(j1)
        end subroutine



        subroutine diagquad2(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
c
c       Build a tensor product quadrature on [-1,1] x [-1,1] for integrals of
c       the form
c
c           1     1
c       \int  \int  (  (x1-y1) / ( (x1-y1)^2 + (x2-y2)^2 )  P_i1(y1) P_j1(y2)  +
c           -1    1
c         
c                      (x2-y2) / ( (x1-y1)^2 + (x2-y2)^2 )  P_i2(y1) P_j2(y2)  +
c
c
c                 P_i3(y1) P_j3(y2) ) dy1 dy2
c
c       where (x1,x2) is a  user-specified point in (-1,1) x (-1,1), P_k
c       denotes the Legendre polynomial of degree k and
c
c           0 <= i1+j1 < = npoly,  0 <= i2+j2 <= npoly,  0 <= i3+j3 <= npoly.
c
c
c       Input parameters:
c         npoly - the degree of the polynomials 
c         (x1,x2) - the specified target node
c
c       Output parameters:
c         ier - an error return code
c             ier   =  0 indicates successful execution
c             ier   =  
c         (nquad,xs,ys,whts) - the output quadrature rule
c 
        dimension thetas(4),rect(4),ab0(2,100),verts(2,3)
        double precision, allocatable :: disc(:),coefs(:,:)
        double precision, allocatable :: xslege(:),whtslege(:)
        double precision, allocatable :: xs0(:),whts0(:)
        double precision, allocatable :: xs1(:),whts1(:)
        double precision, allocatable :: rints(:),errs(:)
        integer, allocatable          :: idxs(:)
        external funuser3,funeval,funquad3,funquad4,funquad5
c
        allocate(xs1(1000),whts1(1000),rints(1000),idxs(1000))

        ier = 0
        pi  = acos(-1.0d0)
c
        nlege = npoly/2+1
        allocate(xslege(nlege),whtslege(nlege))
        call legequad(nlege,xslege,whtslege)
c
        call prin2("in diagquad2 , x1 = *",x1,1)
        call prin2("in diagquad2 , x2 = *",x2,1)

c
c       Determine the coordinates
c
        a   = -1-x1
        b   =  1-x1
c
        c   = -1-x2
        d   =  1-x2
c
        rect(1)   = a
        rect(2)   = c
        rect(3)   = b
        rect(4)   = d
c
        thetas(1) = atan2(d,b)
        thetas(2) = atan2(d,a)
        thetas(3) = atan2(c,a)
        thetas(4) = atan2(c,b)
c
        do i=1,4
        if (thetas(i) .lt. 0) thetas(i) = thetas(i) + 2*pi
        end do
c
        ifadap   = 1
        kdisc    = 30
        nints0   = 5
        ab0(1,1) = 0.0d0
        ab0(2,1) = thetas(1)
        ab0(1,2) = thetas(1)
        ab0(2,2) = thetas(2)
        ab0(1,3) = thetas(2)
        ab0(2,3) = thetas(3)
        ab0(1,4) = thetas(3)
        ab0(2,4) = thetas(4)
        ab0(1,5) = thetas(4)
        ab0(2,5) = 2*pi
c
        ldisc = 100 000 000
        allocate(disc(ldisc))
c
        nfuns3 = (npoly+2)*(npoly+3)/2*3
        call prinf("nfuns3 = *",nfuns3,1)
 0100 continue
        kdisc = kdisc+1
c
        call legedisc(ier,ifadap,nints0,ab0,kdisc,eps,nfuns3,
     -    funuser3,npoly,thetas,rect,par4,disc,ldisc,ncoefs,lkeep)
        if (ier .ne. 0) then
        call prinf("after legedisc, ier = *",ier,1)
        stop
        endif
c
        call prinf("after legedisc, ncoefs = *",ncoefs,1)
c
c       Gram-Schmidt the input functions.
c
        allocate(coefs(ncoefs,nfuns3))
        call legegs(eps,disc,nfuns3,funuser3,npoly,thetas,rect,par4,
     1    krank,coefs)
        call prinf("after legegs, krank=*",krank,1)
c
c       Fetch the initial, oversampled, quadrature rule.
c
        allocate(xs0(ncoefs),whts0(ncoefs))
        call legedisc_quad(disc,nquad0,xs0,whts0)
c
c       Compute the Chebyshev quadrature.
c
        call chebquad(krank,funeval,disc,coefs,krank,par4,
     1    nquad0,xs0,whts0,nquad1,xs1,whts1,rints)
c
        ifaccept = 0
        ngoal    = 0
        a        = 0
        b        = 2*pi
c
        call gaussquad(eps,krank,rints,funeval,disc,coefs,krank,
     1    par4,nquad1,xs1,whts1,a,b,ngoal,ifaccept)
        deallocate(coefs,xs0,whts0)
c
        if (ier .ne. 0) goto 0100
c        ngoal = (krank+1)/2+1
c        if (nquad1 .gt. ngoal) goto 0100
        do i=1,nquad1
        if (whts1(i) .le. 0) goto 0100
        if (xs1(i) .lt. a .OR. xs1(i) .gt. b) goto 0100
        end do

c
c       Build the tensor product quadrature
c
        nquad = 0

        do i=1,nquad1

        t    = xs1(i)
        twht = whts1(i)

        r1   = 0
        call rfun(rect,thetas,t,r2)

        do j=1,nlege
        r    =   xslege(j)*(r2-r1)/2 + (r2-r1)/2
        rwht = whtslege(j)*(r2-r1)/2
c
        x    = r*cos(t)
        y    = r*sin(t)
        wht  = twht*rwht*r
c
        nquad = nquad+1
        xs(nquad)   = x      + x1
        ys(nquad)   = y      + x2
        whts(nquad) = wht
        end do
        end do
c
        call prinf("in diagquad2, final nquad = *",nquad,1)
        call prin2("in diagquad2, final xs    = *",xs,nquad)
        call prin2("in diagquad2, final ys    = *",ys,nquad)
        call prin2("in diagquad2, final whts  = *",whts,nquad)

         end

        subroutine diagquad2_check(epsadap,x1,x2,npoly,
     -    nquad,xs,ys,whts,errmax)
        implicit double precision (a-h,o-z)
        double precision :: xs(nquad),ys(nquad),whts(nquad)
        double precision, allocatable :: errs(:)
        dimension verts(2,3)
c     
c       Test the newly minted quadrature rule thoroughly via comparison
c       with results obtained with adaptive quadrature
c
        nn = (npoly+1)*(npoly+2)/2*3

        allocate(errs(nn))
c
        idx = 0

        do nn1=0,npoly
        do i1=0,nn1
        j1 = nn1-i1
c
        verts(1,1) = -1-x1
        verts(2,1) = -1-x2
        verts(1,2) =  1-x1
        verts(2,2) = -1-x2
        verts(1,3) =  1-x1
        verts(2,3) =  1-x2
c
        call adaptri(ier,epsadap,verts,funquad3,x1,x2,i1,j1,
     1    val1,nquad01)
        if (ier .ne. 0) then
           print *,"!"
           stop
        endif
c
        verts(1,1) = -1-x1
        verts(2,1) = -1-x2
        verts(1,2) = -1-x1
        verts(2,2) =  1-x2
        verts(1,3) =  1-x1
        verts(2,3) =  1-x2
c
        call adaptri(ier,epsadap,verts,funquad3,x1,x2,i1,j1,
     1    val2,nquad02)
        if (ier .ne. 0) then
           print *,"!"
           stop
        endif
c
        val0  = val1+val2
c
        val  = 0
        do i=1,nquad
        x   = xs(i)-x1
        y   = ys(i)-x2
        wht = whts(i)
        call funquad3(x,y,x1,x2,i1,j1,val00)
        val = val + wht*val00
        end do
c
c
        idx = idx + 1
        errs(idx) = abs(val-val0)
c
        end do
        end do
c
        do nn1=0,npoly
        do i1=0,nn1
        j1 = nn1-i1
c
        verts(1,1) = -1-x1
        verts(2,1) = -1-x2
        verts(1,2) =  1-x1
        verts(2,2) = -1-x2
        verts(1,3) =  1-x1
        verts(2,3) =  1-x2
c
        call adaptri(ier,epsadap,verts,funquad4,x1,x2,i1,j1,
     1    val1,nquad01)
        if (ier .ne. 0) then
           print *,"!"
           stop
        endif
c
        verts(1,1) = -1-x1
        verts(2,1) = -1-x2
        verts(1,2) = -1-x1
        verts(2,2) =  1-x2
        verts(1,3) =  1-x1
        verts(2,3) =  1-x2
c
        call adaptri(ier,epsadap,verts,funquad4,x1,x2,i1,j1,
     1    val2,nquad02)
        if (ier .ne. 0) then
           print *,"!"
           stop
        endif

        val0 = val1+val2
        val  = 0
        do i=1,nquad
        x   = xs(i)-x1
        y   = ys(i)-x2
        wht = whts(i)
        call funquad4(x,y,x1,x2,i1,j1,val00)
        val = val + wht*val00
        end do
c
        idx = idx + 1
        errs(idx) = abs(val-val0)
c
        end do
        end do
c
        do nn1=0,npoly
        do i1=0,nn1
        j1 = nn1-i1
c
        verts(1,1) = -1
        verts(2,1) = -1
        verts(1,2) =  1
        verts(2,2) = -1
        verts(1,3) =  1
        verts(2,3) =  1
c
        call adaptri(ier,epsadap,verts,funquad5,x1,x2,i1,j1,
     1    val1,nquad01)
        if (ier .ne. 0) then
           print *,"!"
           stop
        endif

c
        verts(1,1) = -1
        verts(2,1) = -1
        verts(1,2) = -1
        verts(2,2) =  1
        verts(1,3) =  1
        verts(2,3) =  1
c
        call adaptri(ier,epsadap,verts,funquad5,x1,x2,i1,j1,
     1    val2,nquad02)
        if (ier .ne. 0) then
           print *,"!"
           stop
        endif

c
        val0  = val1+val2
c
        val  = 0
        do i=1,nquad
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
        call funquad5(x,y,x1,x2,i1,j1,val00)
        val = val + wht*val00
        end do
c
        idx = idx + 1
        errs(idx) = abs(val-val0)
c
        end do
        end do
c
        errmax = maxval(errs)

        call prin2("x1 = *",x1,1)
        call prin2("x2 = *",x1,1)
        call prin2("errs = *",errs,idx)
        call prin2("errmax = *",errmax,1)
cc        call prina("*")
c
c$$$        if (errmax .gt. 1.0d-16) then
c$$$        call prina("MAXIMUM ERROR EXCEEDED*")
c$$$        stop
c$$$        endif
        
c
        end


        subroutine funuser3(theta,vals,npoly,thetas,rect,par4)
        implicit double precision (a-h,o-z)
        dimension thetas(4),rect(4),vals(1)
c
        call rfun(rect,thetas,theta,r)

        idx = 0
c
c$$$        do i  = 0,npoly+1
c$$$        do j1 = 0,npoly
c$$$        do j2 = 0,npoly
c$$$        idx       = idx+1
c$$$
c$$$        if (j2 .eq. 0) then
c$$$        vals(idx) =  r**(i+1)*cos(j1*theta)
c$$$        else
c$$$        vals(idx) =  r**(i+1)*cos(j1*theta) * sin(j2*theta)
c$$$        endif
c$$$
c$$$        end do
c$$$        end do
c$$$        end do
c
        do i  = 0,npoly+1
        do j1 = 0,i
        j2 = i-j1
        idx       = idx+1
        vals(idx) = r**(i+1)*cos(theta)**j1*sin(theta)**j2
        end do
        end do
c
        do i  = 0,npoly+1
        do j1 = 0,i
        j2 = i-j1
        idx       = idx+1
        vals(idx) =  r**(i+1)*cos(theta)**j1 *sin(theta)**j2*cos(theta)
        end do
        end do
c
        do i  = 0,npoly+1
        do j1 = 0,i
        j2 = i-j1
        idx       = idx+1
        vals(idx) =  r**(i+1)*cos(theta)**j1*sin(theta)**j2*sin(theta)
        end do
        end do
c$$$c
c$$$        do i  = 0,npoly+1
c$$$        do j1 = 0,i
c$$$        j2 = i-j1
c$$$        idx       = idx+1
c$$$        vals(idx) = log**(r)*r**(i+1)*cos(theta)**j1*sin(theta)**j2
c$$$        end do
c$$$        end do

        end subroutine


        subroutine rfun(rect,thetas,theta,r)
        implicit double precision (a-h,o-z)
        dimension rect(4),thetas(4)

        a  = rect(1)
        c  = rect(2)
        b  = rect(3)
        d  = rect(4)
c
        if (theta .gt. thetas(4)) then
        r = b/cos(theta)
        elseif (theta .gt. thetas(3)) then
        r = c/sin(theta)
        elseif (theta .gt. thetas(2)) then
        r = a/cos(theta)
        elseif (theta .gt. thetas(1)) then
        r = d/sin(theta)
        else
        r = b/cos(theta)
        endif

        end subroutine



        subroutine funquad3(y1,y2,x1,x2,i1,j1,val)
        implicit double precision (a-h,o-z)
        dimension pols1(0:100)
        dimension pols2(0:100)
        call lege(i1,y1+x1,pols1)
        call lege(j1,y2+x2,pols2)
        dd = (y1) / ( (y1)**2 + (y2)**2 )
        val = dd * pols1(i1)*pols2(j1)

        end subroutine

        subroutine funquad4(y1,y2,x1,x2,i1,j1,val)
        implicit double precision (a-h,o-z)
        dimension pols1(0:100)
        dimension pols2(0:100)
        call lege(i1,y1+x1,pols1)
        call lege(j1,y2+x2,pols2)
        dd = (y2)/( (y1)**2 + (y2)**2 ) 
        val = dd * pols1(i1)*pols2(j1)
        end subroutine

        subroutine funquad5(y1,y2,x1,x2,i1,j1,val)
        implicit double precision (a-h,o-z)
        dimension pols1(0:100)
        dimension pols2(0:100)
        call lege(i1,y1,pols1)
        call lege(j1,y2,pols2)
        val = pols1(i1)*pols2(j1)
        end subroutine





cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine discquad_info(nquad0)
        implicit double precision (a-h,o-z)
        nquad0 = 52
        end subroutine


        subroutine discquad(nquad0,xs0,ys0,whts0)
        implicit double precision (a-h,o-z)
        double precision :: xs0(1),ys0(1),whts0(1)
c
c       Return to the user the discretization quadrature, which integrates
c       polynomials of order 2*ndegree on the rectangle [-1,1] x [-1,1].
c
c       The quadrature is accurate to around 30 digits.  
c
c       It is intended for use in discretizing integral operators by
!       representing their solutions as 8th order polynomials.
c
c                      Output Parameters:
c
c   nquad0 - this will be set to the length of the quadrature rule
c   (xs0,ys0) - these user-supplied arrays will contain the coordinates
c       of the quadrature nodes
c   whts0 - this user-supplied array will contain the (positive)
c       weights
c   
c
        dimension xs(52),ys(52),whts(52)        
        data nquad / 52 /
        data xs /
     -    0.595709534481334984968258155046D-01,
     -   -0.227986154105160583330241093646D+00,
     -    0.600996666306949497920716364961D+00,
     -    0.899352126918021066010099992948D+00,
     -   -0.986669755940225846627651472078D+00,
     -   -0.976248877535581195712571988672D+00,
     -    0.898954253364050910482878685520D+00,
     -   -0.455218061112388245497059635150D+00,
     -    0.697952302565951703906302170236D+00,
     -   -0.219470041816650633047291737158D+00,
     -    0.912266187505508504203200336835D+00,
     -    0.429388818960401874327302776813D+00,
     -    0.994553557340949754016655779132D+00,
     -   -0.911321637435155034819363454824D+00,
     -    0.889712226690185728221707579384D-01,
     -   -0.523844342751902321214790804564D+00,
     -    0.375695403756973699029595174365D+00,
     -   -0.844122099320538929776986898317D+00,
     -    0.736910885966299129758821858716D+00,
     -    0.509514079403359354303515994049D+00,
     -   -0.585837006729577258462052907303D+00,
     -   -0.785616781982602733758263303620D+00,
     -    0.140979513073424654237818794913D+00,
     -    0.987984838658164341063113321777D+00,
     -   -0.195401398340828032823825326342D+00,
     -   -0.718782510989325749368730792766D+00,
     -   -0.779348966307823564954955519465D+00,
     -    0.821037982719170895230912745308D-01,
     -   -0.184477574015358186976301823433D+00,
     -   -0.663824545977195158726413435793D+00,
     -   -0.913728884106098714298670668792D+00,
     -   -0.481346007008979756881408394808D+00,
     -    0.992107383124282573388293241092D+00,
     -    0.834813791806285356988213069303D+00,
     -    0.542517445741263930874392304827D+00,
     -   -0.507084824878077205157045501174D+00,
     -    0.178941947129298424684946933618D+00,
     -   -0.261453349505757331177424442658D+00,
     -    0.449160199288593304212458032116D+00,
     -   -0.182024135261913941285463584569D+00,
     -    0.379205496709569284373702003473D+00,
     -   -0.852529301294830953587057444579D+00,
     -   -0.727186145424036973111518124570D+00,
     -    0.906180038488740883301800198178D+00,
     -    0.973876692040707422155676165243D+00,
     -    0.693844090276462422396386972459D+00,
     -    0.860991271011021653190944230385D-01,
     -   -0.984548557629027586455257508164D+00,
     -   -0.965964031370258865804594092170D+00,
     -    0.912286945536719071206844083055D+00,
     -   -0.977694909717452212066829917651D+00,
     -    0.750120757616906900586825958791D+00/
        data ys /
     -   -0.689535805405646919002663483239D-01,
     -   -0.989124771806003910804562214942D+00,
     -   -0.367484934832922375648928777847D-01,
     -   -0.194606775085394016069830134723D+00,
     -   -0.959593163032965917037033349720D+00,
     -   -0.534001349721656143485264536171D+00,
     -   -0.974016067873542729464141877075D+00,
     -   -0.462102355939142299831437813532D+00,
     -   -0.450028960327200829995732484926D+00,
     -   -0.199801816686291300380418304252D+00,
     -   -0.643910307022118264717441968850D+00,
     -   -0.710266710943780059707379144216D+00,
     -   -0.321478062065767456353028059268D+00,
     -   -0.819021078614365693815229943608D+00,
     -   -0.495615135142315027438918718158D+00,
     -   -0.890191073145983906569560124040D+00,
     -   -0.251942935648320595893507362467D+00,
     -   -0.294678273131199118407137457468D+00,
     -   -0.856000358390583292812483365901D+00,
     -   -0.976531607234978933225320269321D+00,
     -   -0.940607666511545512571131747397D-01,
     -   -0.972013121342990592532014728216D+00,
     -   -0.908096326078309022235027565223D+00,
     -   -0.866030995664358036032006704901D+00,
     -   -0.742668088922330857582418078126D+00,
     -   -0.655633170883756139778169025070D+00,
     -    0.973969808448651180343754153805D+00,
     -    0.645333624384553718422789531263D-01,
     -    0.755111198862645749642122903570D+00,
     -    0.115468415287513987832298165919D+00,
     -    0.834591614484737089954174021260D+00,
     -    0.500964564503117494560330076117D+00,
     -    0.829729472514972049926736637682D+00,
     -    0.800522211339293868449208740686D-01,
     -    0.976112501669398133394581987984D+00,
     -    0.901379456863444628860640125050D+00,
     -    0.905977077045154644349975643589D+00,
     -    0.233674599046493686446443488904D+00,
     -    0.700590488387249751957313916403D+00,
     -    0.990541069619248078291670313177D+00,
     -    0.244076584425381010981673100673D+00,
     -    0.341352596492516425000960652095D+00,
     -    0.690614614929954706792134900088D+00,
     -    0.612241567065143478736817713980D+00,
     -    0.263855673266835337180277831011D+00,
     -    0.422268124277051047495238141553D+00,
     -    0.505619733358839230041293961490D+00,
     -    0.963689161403230562389735519594D+00,
     -    0.169186684439590326367049014821D-01,
     -    0.967327907048547942946479874605D+00,
     -    0.563190791429714154246918389041D+00,
     -    0.846785895025508492765999972283D+00/
        data whts /
     -    0.615281232546091849148905325633D-01,
     -    0.277574626539027077602847404324D-01,
     -    0.101949032789261702489909901828D+00,
     -    0.585430473028673511136389238768D-01,
     -    0.691647504897177509278309989318D-02,
     -    0.324967851871293063058988061809D-01,
     -    0.168470376004779304953726084377D-01,
     -    0.138375481028201034591794543967D+00,
     -    0.124460719402234548575317861368D+00,
     -    0.140465739804564229261602351541D+00,
     -    0.624492654785096935096637054296D-01,
     -    0.130249656491254347325207965620D+00,
     -    0.194033305322927622277072922220D-01,
     -    0.452970623713953478283506021525D-01,
     -    0.158601159010172776345470093530D+00,
     -    0.775826928453944343783396858536D-01,
     -    0.134641139527149578718183516271D+00,
     -    0.103673659218933042754113921637D+00,
     -    0.698236200830464292195244898939D-01,
     -    0.326324623889909371977603863151D-01,
     -    0.992431040086694358451654820131D-01,
     -    0.252737399713167042642801677928D-01,
     -    0.864315360650566644873925189463D-01,
     -    0.130204450455155208971844302674D-01,
     -    0.130497153584554189292248159508D+00,
     -    0.996819736787061016704764376112D-01,
     -    0.240016578722406330163980310464D-01,
     -    0.944960494793990777133013530024D-01,
     -    0.128093144418908326009820112500D+00,
     -    0.981047262199623563955626454904D-01,
     -    0.407817547706184255877150730520D-01,
     -    0.133427778894442783875949882847D+00,
     -    0.134842281610511311194429524676D-01,
     -    0.769880431627056043786047056251D-01,
     -    0.324813939129360386868179338272D-01,
     -    0.744213980433165021978994373745D-01,
     -    0.872784921380304707404297849913D-01,
     -    0.157081503430786621736640640203D+00,
     -    0.129648782098094093287614782707D+00,
     -    0.270547203072481796739620160488D-01,
     -    0.140424549898729083566087059419D+00,
     -    0.937646853218780675474811397482D-01,
     -    0.906948487106129705168995019889D-01,
     -    0.654220462216669162057596782103D-01,
     -    0.381843876612410217797338594854D-01,
     -    0.122109504310161140514977577169D+00,
     -    0.161076101681514160857336401800D+00,
     -    0.650443643251314782561366838004D-02,
     -    0.483388732616093902972680343073D-01,
     -    0.186030248524897467185471164651D-01,
     -    0.298772843894898525921809435685D-01,
     -    0.698146799751765205953934431233D-01/
!
        nquad0 = nquad
        do i=1,nquad0
        xs0(i)   = xs(i)
        ys0(i)   = ys(i)
        whts0(i) = whts(i)
        end do
        end


c$$$        subroutine discquad_info(nquad0)
c$$$        implicit double precision (a-h,o-z)
c$$$        nquad0 = 78
c$$$        end subroutine
c$$$
c$$$
c$$$        subroutine discquad(nquad0,xs0,ys0,whts0)
c$$$        implicit double precision (a-h,o-z)
c$$$        dimension xs(78),ys(78),whts(78)
c$$$        dimension xs0(1),ys0(1),whts0(1)
c$$$c
c$$$c       Return to the user a 77-point quadrature which integrates
c$$$c       polynomials of order 20 on the rectangle [-1,1] x [-1,1].
c$$$c
c$$$c       The quadrature is accurate to around 30 digits.  
c$$$c
c$$$c       It is intended for use in discretizing integral operators by
c$$$!       representing their solutions as 10th order polynomials.
c$$$c
c$$$c                      Output Parameters:
c$$$c
c$$$c   nquad0 - this will be set to the length of the quadrature rule
c$$$c   (xs0,ys0) - these user-supplied arrays will contain the coordinates
c$$$c       of the quadrature nodes
c$$$c   whts0 - this user-supplied array will contain the (positive)
c$$$c       weights
c$$$c   
c$$$        data xs /
c$$$     -    0.926555975905755513915971020626D+00,
c$$$     -   -0.913024526647190478369028753619D+00,
c$$$     -   -0.988383918890483819383402716217D+00,
c$$$     -   -0.897909422586895589416728775740D+00,
c$$$     -    0.918250528282333966739170174725D+00,
c$$$     -    0.993143912399823407705012257604D+00,
c$$$     -    0.601563907173806964601862595067D+00,
c$$$     -    0.591759998035315642428834249077D+00,
c$$$     -    0.989273837398036769501509838829D+00,
c$$$     -    0.884733440486114138685366797227D+00,
c$$$     -    0.987599205146504368430067331523D+00,
c$$$     -   -0.540417259542395087683575775109D+00,
c$$$     -   -0.778721260506808284854511539355D+00,
c$$$     -   -0.889304384214424968764066350354D+00,
c$$$     -   -0.839511846147327511251256435860D+00,
c$$$     -   -0.612685594580635033592251142141D+00,
c$$$     -   -0.941505382885246589850443105423D+00,
c$$$     -    0.984048738732380799043733621615D+00,
c$$$     -   -0.150010524775347637804222872163D+00,
c$$$     -   -0.877201282592437277489172004856D+00,
c$$$     -   -0.836917200379954706480569355748D+00,
c$$$     -   -0.945990308809321678878553206313D+00,
c$$$     -   -0.456193598944486462387801958428D+00,
c$$$     -    0.574676524396181193725809306587D+00,
c$$$     -   -0.446114543415753682566439186327D+00,
c$$$     -   -0.410309846957920086422436798277D+00,
c$$$     -   -0.642916721460959665305071493064D+00,
c$$$     -    0.807284370961669538482761738836D+00,
c$$$     -    0.221915152687594509474002166572D+00,
c$$$     -   -0.653296443512545709106091297676D+00,
c$$$     -   -0.818353338238848651312936934649D+00,
c$$$     -    0.147870626481583843302945488476D+00,
c$$$     -    0.123099056887539677578056155333D-01,
c$$$     -    0.797550310031424215456924128415D+00,
c$$$     -    0.919586820989622809256771707413D+00,
c$$$     -   -0.180461794924694458508021614405D+00,
c$$$     -   -0.296338789612619901800456842131D+00,
c$$$     -   -0.199396171849393871743703895950D+00,
c$$$     -   -0.986659079083316701998132121168D+00,
c$$$     -   -0.642046340295263110253769590402D+00,
c$$$     -    0.643200129865208274234781541213D+00,
c$$$     -    0.984615106335576606425482055215D+00,
c$$$     -    0.977416813702975436606890120051D+00,
c$$$     -    0.344272670575593763050503230385D+00,
c$$$     -    0.333447895749259321026740268278D+00,
c$$$     -    0.454657277404616748499685536652D+00,
c$$$     -    0.754189994881143297810597417724D-01,
c$$$     -   -0.253606236270413835788494038084D+00,
c$$$     -    0.947241867923324739557767719148D-01,
c$$$     -   -0.800249996718763461100192874774D+00,
c$$$     -    0.930962764468901076407877765860D+00,
c$$$     -   -0.981942622386643710915337719703D+00,
c$$$     -   -0.688121808834076768852645297433D+00,
c$$$     -   -0.461031209491973810695035990513D+00,
c$$$     -   -0.149213663175234621638448027073D+00,
c$$$     -    0.924011183305361011036081377753D+00,
c$$$     -    0.398529190358020520979892730693D+00,
c$$$     -   -0.428407343854623858635750309209D+00,
c$$$     -   -0.240172369916866889617323319703D+00,
c$$$     -    0.276893339778031360449601997239D+00,
c$$$     -    0.524397383009391746127601919110D+00,
c$$$     -    0.932996045706880613226821264688D+00,
c$$$     -    0.665912752284963279216215882303D-01,
c$$$     -    0.414344270105987592196560959318D+00,
c$$$     -    0.811236955485436692199089514073D+00,
c$$$     -   -0.996366768246143206855615424925D+00,
c$$$     -   -0.935245785900152428541161254202D+00,
c$$$     -   -0.953404675232094786009673726234D+00,
c$$$     -    0.731134516240156292380218685231D+00,
c$$$     -    0.791098009763760864638204388590D+00,
c$$$     -   -0.991652442640997640154153133121D+00,
c$$$     -    0.169860983241556642305791642017D+00,
c$$$     -   -0.679450346661891347463572565658D+00,
c$$$     -    0.628826776673926873172011341350D+00,
c$$$     -    0.795431813190396868786694476380D+00,
c$$$     -   -0.279211619091827389292393778096D+00,
c$$$     -   -0.359790108177762311921762039076D-01,
c$$$     -    0.640016632659688312578046148144D+00/
c$$$        data ys /
c$$$     -   -0.993648204713861939238672212463D+00,
c$$$     -    0.987130947373865906151551613803D+00,
c$$$     -   -0.602417236771653182297198661648D-01,
c$$$     -   -0.990078266619169016546422409097D+00,
c$$$     -    0.992170806540107267359623443082D+00,
c$$$     -    0.649923047829817170968616730613D+00,
c$$$     -   -0.988992904425822674075772433170D+00,
c$$$     -    0.988741360723661080483120087239D+00,
c$$$     -   -0.661885609912338212297293041116D+00,
c$$$     -   -0.589648677073402235776469461567D-01,
c$$$     -   -0.196936485442616018184760091396D+00,
c$$$     -   -0.122559070264894807019501360355D+00,
c$$$     -    0.662032776392395025277887792183D-01,
c$$$     -    0.354671122862718363330850358394D+00,
c$$$     -    0.524180112287403829686349341006D+00,
c$$$     -   -0.974131307081767833213995019516D+00,
c$$$     -    0.749589811011490065755033408537D+00,
c$$$     -   -0.932279757497518808467932279191D+00,
c$$$     -   -0.196546387804985927687699786736D+00,
c$$$     -   -0.158610129929298288536740967721D+00,
c$$$     -   -0.570090901257351846150866303548D+00,
c$$$     -    0.221950836855751208462952484345D+00,
c$$$     -    0.538549475095439271953139876694D+00,
c$$$     -    0.496023794868792503828407698271D-02,
c$$$     -   -0.525189588437537864756702545208D+00,
c$$$     -    0.638128797305795367434545215062D-01,
c$$$     -    0.973174302216364122273736045956D+00,
c$$$     -    0.636291109246586358535297256753D+00,
c$$$     -    0.755662801563743066816869966502D+00,
c$$$     -   -0.728893181073822141871191225155D+00,
c$$$     -    0.890352120949259641203425576239D+00,
c$$$     -    0.983768482394454366251318774449D+00,
c$$$     -   -0.881099630592665102691290128583D+00,
c$$$     -   -0.938254611968303487364339023792D+00,
c$$$     -   -0.818072503959121905042384177423D+00,
c$$$     -   -0.699605098071910988400212582964D+00,
c$$$     -    0.380377809890716863114554215598D+00,
c$$$     -    0.701876823794319868632209288516D+00,
c$$$     -    0.922512441980898081256106451930D+00,
c$$$     -    0.312088852070325662571832665018D+00,
c$$$     -    0.409464575912686427596815236802D+00,
c$$$     -    0.935471458253307913795964078839D+00,
c$$$     -    0.180536973692840059604701186057D+00,
c$$$     -   -0.254882234217954227083602228643D+00,
c$$$     -    0.256651183881562284476244668566D+00,
c$$$     -    0.590757774362036772889814959872D+00,
c$$$     -   -0.484127648067765600304856347011D+00,
c$$$     -   -0.315320832661566670016852413906D+00,
c$$$     -   -0.431583010604228764300746544117D-02,
c$$$     -   -0.898182370078303273440198962144D+00,
c$$$     -    0.443683219921730535543288873402D+00,
c$$$     -   -0.934585457473875045694258304872D+00,
c$$$     -   -0.338372338028090124253582118659D+00,
c$$$     -    0.859739870516645506635578825267D+00,
c$$$     -    0.229574735334836015384532216951D+00,
c$$$     -    0.813705754968955761386414286595D+00,
c$$$     -    0.918347935922854845634119508798D+00,
c$$$     -   -0.863372449216151511561984614266D+00,
c$$$     -   -0.972979523109948122687616656014D+00,
c$$$     -   -0.730531871814222705802408358610D+00,
c$$$     -   -0.538233467070080911935183752313D+00,
c$$$     -   -0.433744951160720492513584623634D+00,
c$$$     -    0.491961233650953300785348873058D+00,
c$$$     -   -0.919973053174258516383550200903D+00,
c$$$     -   -0.631832993350726321376846035408D+00,
c$$$     -   -0.614597420166648087605235998476D+00,
c$$$     -   -0.778706398697101333607457081641D+00,
c$$$     -   -0.384711310713382345096775351173D+00,
c$$$     -   -0.301992354047231478827733035938D+00,
c$$$     -    0.933076746705474555811151779519D+00,
c$$$     -    0.516742215995684131300885888414D+00,
c$$$     -   -0.984743330226400918784720287115D+00,
c$$$     -    0.712775567048347239291976978788D+00,
c$$$     -    0.802181206568358148679391864004D+00,
c$$$     -    0.218967157653944152722332329318D+00,
c$$$     -    0.974524379113034143837389288894D+00,
c$$$     -    0.889459742434036328682773954604D+00,
c$$$     -   -0.805114446946467495054184414376D+00/
c$$$        data whts /
c$$$     -    0.460398595861064633418696704358D-02,
c$$$     -    0.761981217725659315597072679948D-02,
c$$$     -    0.155924709094069706701650855913D-01,
c$$$     -    0.716586388891790303799899498080D-02,
c$$$     -    0.526266753746648903021446758392D-02,
c$$$     -    0.107865338911094719843097112108D-01,
c$$$     -    0.122175263684921885035498962934D-01,
c$$$     -    0.130775951976592540191123899761D-01,
c$$$     -    0.121024260530215339155875150335D-01,
c$$$     -    0.576442971891347930299879265007D-01,
c$$$     -    0.149468753071995635499063061831D-01,
c$$$     -    0.742182180618808783789147239361D-01,
c$$$     -    0.690386282414469928791148134137D-01,
c$$$     -    0.150992007043352441065476034512D-01,
c$$$     -    0.569828043914472636259376100818D-01,
c$$$     -    0.220598733776956574239045675450D-01,
c$$$     -    0.305875204898204716488291758425D-01,
c$$$     -    0.659761613211303482784720798000D-02,
c$$$     -    0.794631304025332771485512002248D-01,
c$$$     -    0.482381609615147683981473373810D-01,
c$$$     -    0.560628969869226906439018217012D-01,
c$$$     -    0.370805709261253717277026020856D-01,
c$$$     -    0.832720055534993977223415210205D-01,
c$$$     -    0.111990643751910176714581963227D+00,
c$$$     -    0.968715291537771994559898466139D-01,
c$$$     -    0.928463057447403776714401231772D-01,
c$$$     -    0.220827906907636425080015706828D-01,
c$$$     -    0.571110347426015148553188724989D-01,
c$$$     -    0.785591434368593265064095861479D-01,
c$$$     -    0.665066120437950197703467782347D-01,
c$$$     -    0.349738739205334483340792814121D-01,
c$$$     -    0.195628847391265330117503407276D-01,
c$$$     -    0.609709091573250221404088827577D-01,
c$$$     -    0.263941764820063195742716651007D-01,
c$$$     -    0.277412948144514851853156718941D-01,
c$$$     -    0.930653899270807974590881345928D-01,
c$$$     -    0.505802816084033647114931410026D-01,
c$$$     -    0.902714685338996637312218789196D-01,
c$$$     -    0.734959758039693907967725417393D-02,
c$$$     -    0.902156944357686813537200815775D-01,
c$$$     -    0.752546875047434689306540421734D-01,
c$$$     -    0.645211804593323326523477559890D-02,
c$$$     -    0.224993251302016941087488270874D-01,
c$$$     -    0.122249023275837379907653918609D+00,
c$$$     -    0.123525549055220192558784351722D+00,
c$$$     -    0.841115976235702164349964331024D-01,
c$$$     -    0.120275578509792191701324699799D+00,
c$$$     -    0.729878387319412348545042247352D-01,
c$$$     -    0.125790054971517970891912350695D+00,
c$$$     -    0.347501266599985644230936856815D-01,
c$$$     -    0.398335233888660415541625779497D-01,
c$$$     -    0.780248179369299494120229227259D-02,
c$$$     -    0.817854670335271590780833911506D-01,
c$$$     -    0.612492738964380146350539946757D-01,
c$$$     -    0.100271790169694746696502317604D+00,
c$$$     -    0.285958403173166589883542971699D-01,
c$$$     -    0.451817449534254396623857920046D-01,
c$$$     -    0.596293058565145837147183649744D-01,
c$$$     -    0.272786591846369709736068578805D-01,
c$$$     -    0.887608148317828538233538132608D-01,
c$$$     -    0.970732820662645982129183822257D-01,
c$$$     -    0.402725161009794698084794892423D-01,
c$$$     -    0.116473193198808009930519127783D+00,
c$$$     -    0.438077373993211489642494288781D-01,
c$$$     -    0.567522420499845428874791803417D-01,
c$$$     -    0.852224079854122891845218600231D-02,
c$$$     -    0.286348317795871802983826609974D-01,
c$$$     -    0.309914362637898733638857569559D-01,
c$$$     -    0.884403895910984790510191368257D-01,
c$$$     -    0.283776953197384880624365263482D-01,
c$$$     -    0.132074163939985573458568693335D-01,
c$$$     -    0.183765460047622103909850631511D-01,
c$$$     -    0.669495918528905386322349324354D-01,
c$$$     -    0.595245260233196770255340689194D-01,
c$$$     -    0.663586975360361693276491488040D-01,
c$$$     -    0.259381564635439561146616954376D-01,
c$$$     -    0.561535580929431549219569156586D-01,
c$$$     -    0.590468306586911377731231779118D-01/
c$$$c
c$$$
c$$$!
c$$$        nquad0  = 78
c$$$        do i=1,nquad0
c$$$        xs0(i)   = xs(i)
c$$$        ys0(i)   = ys(i)
c$$$        whts0(i) = whts(i)
c$$$        end do
c$$$        end
c$$$

c$$$
c$$$
c$$$
c$$$        subroutine discquad_info(nquad0)
c$$$        implicit double precision (a-h,o-z)
c$$$        nquad0 = 110
c$$$        end subroutine
c$$$
c$$$
c$$$        subroutine discquad(nquad0,xs0,ys0,whts0)
c$$$        implicit double precision (a-h,o-z)
c$$$        dimension xs(110),ys(110),whts(110)
c$$$        dimension xs0(1),ys0(1),whts0(1)
c$$$c
c$$$c       Return to the user a 110-point quadrature which integrates
c$$$c       polynomials of order 24 on the rectangle [-1,1] x [-1,1].
c$$$c
c$$$c       The quadrature is accurate to around 30 digits.  
c$$$c
c$$$c       It is intended for use in discretizing integral operators by
c$$$!       representing their solutions as 10th order polynomials.
c$$$c
c$$$c                      Output Parameters:
c$$$c
c$$$c   nquad0 - this will be set to the length of the quadrature rule
c$$$c   (xs0,ys0) - these user-supplied arrays will contain the coordinates
c$$$c       of the quadrature nodes
c$$$c   whts0 - this user-supplied array will contain the (positive)
c$$$c       weights
c$$$c   
c$$$        data xs /
c$$$     -    0.578338131358755468487869321594D+00,
c$$$     -    0.989917932324099050193935280546D+00,
c$$$     -    0.935689289820520311538094755839D+00,
c$$$     -    0.983134715216754805520926423259D+00,
c$$$     -    0.655161144093537672312943803856D+00,
c$$$     -   -0.110655216873635288392847688042D+00,
c$$$     -    0.246417459828942791372291267539D+00,
c$$$     -    0.964261169671806015073488453988D+00,
c$$$     -   -0.995285733385343391129818036007D+00,
c$$$     -   -0.996959304328293547935604733653D+00,
c$$$     -    0.291964090180750006566476813415D+00,
c$$$     -    0.436474150480087436291296847570D+00,
c$$$     -   -0.988451227611893992447351666609D+00,
c$$$     -    0.854401157727382250299705972319D+00,
c$$$     -    0.997484321550589929006927719894D+00,
c$$$     -    0.426214296126862241761072890499D+00,
c$$$     -    0.335636159365554643416700936507D+00,
c$$$     -    0.581753765112829079263542995293D+00,
c$$$     -    0.261371395194623688224803542602D+00,
c$$$     -    0.988112658439639188183769353580D+00,
c$$$     -    0.799058303313613530936235971794D+00,
c$$$     -   -0.423919799542618553207034516322D+00,
c$$$     -   -0.449541585027222162092327733859D+00,
c$$$     -    0.979806685879784382380267436753D+00,
c$$$     -   -0.263405985666378587070449056169D+00,
c$$$     -   -0.928450805373773067572289494746D+00,
c$$$     -   -0.920579276148290826392617922165D+00,
c$$$     -   -0.284983631424843087968839785081D+00,
c$$$     -    0.900998414604663572100711353926D+00,
c$$$     -   -0.419509142882353243170099450846D+00,
c$$$     -    0.675822837792625085091399220424D+00,
c$$$     -    0.735681191563475753084312948407D+00,
c$$$     -   -0.453478061289327256192346142702D+00,
c$$$     -    0.555493551187429532058350313632D-01,
c$$$     -   -0.925029461405186415061359901533D+00,
c$$$     -   -0.973940346898449815257749678651D+00,
c$$$     -   -0.992661640104854562230262001359D+00,
c$$$     -   -0.280598552997141223535701410448D+00,
c$$$     -   -0.208322083873715094186667195882D+00,
c$$$     -    0.866237547645044546197848788279D+00,
c$$$     -    0.942959436926699251298770526917D+00,
c$$$     -    0.992544332517046611637362757333D+00,
c$$$     -    0.598372776764861072367919875398D+00,
c$$$     -   -0.397596220040245249340498511314D+00,
c$$$     -    0.888614730314901988671395179838D+00,
c$$$     -   -0.978273253758586041013486821925D+00,
c$$$     -    0.927389020158637257500213299041D+00,
c$$$     -    0.232103426640469657140462847738D+00,
c$$$     -    0.863069122428051282316981950237D+00,
c$$$     -   -0.951821132544472406145221015163D+00,
c$$$     -   -0.591586279399231990293722508995D+00,
c$$$     -   -0.239721511965783405462322952786D-01,
c$$$     -    0.805890687952985921077826935259D+00,
c$$$     -   -0.573316634590597478698706401730D+00,
c$$$     -    0.187974897099413937911881686013D+00,
c$$$     -   -0.989509698301903320846304125640D+00,
c$$$     -    0.830328340431802634804746259359D-01,
c$$$     -    0.520923995230169896390785027991D+00,
c$$$     -    0.722878150557529495361318262679D+00,
c$$$     -    0.765816601077174972405451813030D+00,
c$$$     -    0.951229615832098554558395257173D+00,
c$$$     -   -0.224474035366275805741653041318D+00,
c$$$     -    0.987650080825455349799756428923D+00,
c$$$     -    0.105423826663379597357442149442D+00,
c$$$     -    0.647780389081604045542491539788D+00,
c$$$     -    0.649290070268204179940653614305D-01,
c$$$     -   -0.736276656474779423883100577380D+00,
c$$$     -   -0.850045294127175556695683809176D+00,
c$$$     -    0.481556051687227869629777718461D+00,
c$$$     -   -0.849200801801683920856031068658D-01,
c$$$     -    0.378304037321419630815181545324D+00,
c$$$     -   -0.614255425944048938145133069935D+00,
c$$$     -   -0.951186534829595337625972056628D+00,
c$$$     -   -0.879187890938707422006105340295D+00,
c$$$     -   -0.379455141861730646118569600399D-01,
c$$$     -   -0.255106202264953620380597368537D+00,
c$$$     -   -0.414399313887484609715703345882D+00,
c$$$     -    0.963543994552770888753039928399D+00,
c$$$     -   -0.727602091741514801576195878309D+00,
c$$$     -    0.787486392425711862331676808552D+00,
c$$$     -   -0.972759112566351730906621834530D+00,
c$$$     -   -0.876305298533984338149618694348D+00,
c$$$     -   -0.745504268814309511060346500870D+00,
c$$$     -    0.548211538672122460312125794769D+00,
c$$$     -   -0.737213485817654628361592252519D+00,
c$$$     -   -0.125326434120443540587448473954D+00,
c$$$     -    0.612993333873251227134602949480D+00,
c$$$     -    0.432456244550277552477654644761D+00,
c$$$     -   -0.928017979157815620525813590876D+00,
c$$$     -   -0.554685072129007856874500072879D+00,
c$$$     -   -0.831831291873353540702696777272D+00,
c$$$     -   -0.703512895863964609238058219688D+00,
c$$$     -    0.180972124958704441231501687514D+00,
c$$$     -    0.916070519975212002115586978105D+00,
c$$$     -    0.769047961466390775133842493843D+00,
c$$$     -   -0.842383659619120590325760815849D+00,
c$$$     -    0.398281741509930552640480353081D+00,
c$$$     -   -0.214497141392241964162208667720D+00,
c$$$     -   -0.349074145811955135731180658943D-01,
c$$$     -    0.451057560309855462105171713851D+00,
c$$$     -    0.140000024909246638512141042587D+00,
c$$$     -   -0.573216496218163359969158441824D+00,
c$$$     -   -0.992171640813981647205734356980D+00,
c$$$     -   -0.714281298687964104998705917887D+00,
c$$$     -    0.830580684053784210745149293812D+00,
c$$$     -   -0.578781755885851655359537709569D+00,
c$$$     -   -0.850629641237496353856080928874D+00,
c$$$     -    0.253882448139543283775802815075D+00,
c$$$     -   -0.601035980825033639140114887954D+00,
c$$$     -   -0.836443659001792179101334934813D+00/
c$$$        data ys /                   
c$$$     -    0.987738783592016334452219720187D+00,
c$$$     -   -0.123388216730251988275521550371D+00,
c$$$     -    0.422110729518207013442453559873D+00,
c$$$     -    0.237769962442095514118137254945D+00,
c$$$     -    0.929102264061020734929011246267D+00,
c$$$     -   -0.107260917199191562060753458568D+00,
c$$$     -   -0.959158397915358604981830271782D-01,
c$$$     -    0.982797664268167295738524211060D+00,
c$$$     -    0.313774998985064534494468961096D+00,
c$$$     -   -0.317238564257656165488980039181D+00,
c$$$     -   -0.312543254103086799390898910916D-01,
c$$$     -    0.985691129020807628808962976142D+00,
c$$$     -    0.980919078160093524898693099669D+00,
c$$$     -    0.387175233375450949741751303962D+00,
c$$$     -    0.906038894426891699191143627260D+00,
c$$$     -    0.868306423239176119490301793295D+00,
c$$$     -   -0.928106794565659315583192306837D+00,
c$$$     -    0.744163198498799072850844534611D+00,
c$$$     -    0.937810609614595459929852120836D+00,
c$$$     -    0.570414502339396910399206502924D+00,
c$$$     -   -0.994014922128567818494895923783D+00,
c$$$     -   -0.188477226610190176227740057209D+00,
c$$$     -    0.280914416382924182468088748753D+00,
c$$$     -   -0.988330872888656696583153047318D+00,
c$$$     -    0.475806107325555093866630219529D+00,
c$$$     -    0.223536716121078921610220466969D+00,
c$$$     -   -0.687441241709610641724325809186D+00,
c$$$     -   -0.369994685298573693491972230576D+00,
c$$$     -    0.915530771990217897095448553916D+00,
c$$$     -    0.925913647821213555527251032256D+00,
c$$$     -   -0.944714101932633175341835599528D+00,
c$$$     -    0.581352672207178451671059911405D+00,
c$$$     -    0.646612354548894872457418550110D+00,
c$$$     -   -0.710552837606063190940753149275D+00,
c$$$     -   -0.262480472706355853534623023149D+00,
c$$$     -    0.522500159270853780844318809988D+00,
c$$$     -    0.789453832632186040491295778106D+00,
c$$$     -    0.931287775083829933758388432955D-01,
c$$$     -   -0.985857956607461113249518780559D+00,
c$$$     -   -0.126012238168086025688124284833D+00,
c$$$     -   -0.748732237378739585158648426002D+00,
c$$$     -   -0.584036919585183168808324459813D+00,
c$$$     -    0.396528381830265291520066731885D+00,
c$$$     -   -0.672546115406545413584088341617D+00,
c$$$     -    0.659049056991202310536239244943D+00,
c$$$     -   -0.153363432819385804992589077532D-01,
c$$$     -    0.739454148445587794032387252192D-01,
c$$$     -   -0.822668384654059063177408000164D+00,
c$$$     -   -0.562772420187030489515014425002D+00,
c$$$     -   -0.910960497780475203940899285397D+00,
c$$$     -   -0.993694528707898245844497071740D+00,
c$$$     -   -0.915323382062898979143099818013D+00,
c$$$     -    0.984717518085326622276102045589D+00,
c$$$     -    0.992683889847213119449333755093D+00,
c$$$     -    0.766909998724911116996054531668D+00,
c$$$     -   -0.981075237253124068566503869941D+00,
c$$$     -    0.130552564745526061148162909878D+00,
c$$$     -    0.881994001143697114894829847713D-01,
c$$$     -   -0.716280973472452864819435717504D+00,
c$$$     -   -0.323637496704597571909849704558D+00,
c$$$     -   -0.363113777528609807944340090781D+00,
c$$$     -   -0.819681310895534634596201890011D+00,
c$$$     -   -0.881690905352544661568014188064D+00,
c$$$     -    0.983961474824500893543521508550D+00,
c$$$     -   -0.794648699349394462242807374977D-01,
c$$$     -   -0.318187652351975523326585783255D+00,
c$$$     -    0.942334899824859615071497663191D+00,
c$$$     -    0.449656115591339183077378783485D+00,
c$$$     -   -0.986768549100099562461546323643D+00,
c$$$     -    0.304686756319859013889272141386D+00,
c$$$     -    0.246166961314600731204413582202D+00,
c$$$     -    0.463445760417664210362944425566D+00,
c$$$     -    0.914972558895649833325113600804D+00,
c$$$     -   -0.985037427249410909256584851165D+00,
c$$$     -    0.893963429553482941339319652301D+00,
c$$$     -    0.786960832273606417031100874541D+00,
c$$$     -   -0.938314390943547927378444399693D+00,
c$$$     -    0.783060359432202353852236109814D+00,
c$$$     -    0.219518496620202337653009002605D+00,
c$$$     -    0.818839442076717018050462031341D+00,
c$$$     -   -0.529214869519173086462562343218D+00,
c$$$     -    0.986629502965302453144760327043D+00,
c$$$     -    0.662476113605986935558934010442D+00,
c$$$     -   -0.843254622868324991600838470364D+00,
c$$$     -   -0.940745696197973394055954412527D+00,
c$$$     -   -0.553201285837351425096359394364D+00,
c$$$     -   -0.500766814827076628881161742015D+00,
c$$$     -   -0.677017164449588001756865085874D+00,
c$$$     -    0.682299172236535595676168772503D+00,
c$$$     -   -0.476734657107389201442736668104D+00,
c$$$     -   -0.475258535426877974982331253288D+00,
c$$$     -   -0.252822879633909571313139000577D+00,
c$$$     -    0.436768997253338672918098799180D+00,
c$$$     -   -0.951935465791726957984414602942D+00,
c$$$     -    0.179947740115926907550884899948D+00,
c$$$     -   -0.830364751493477094394829329507D+00,
c$$$     -    0.594572176801124747963056499065D+00,
c$$$     -    0.980897417867983074415155635543D+00,
c$$$     -    0.626912899441285698397178671150D+00,
c$$$     -   -0.284241096717721054709091202614D+00,
c$$$     -   -0.983870030100289638089940861533D+00,
c$$$     -    0.921209208775992472487614855103D-02,
c$$$     -   -0.783562420691980944334249904800D+00,
c$$$     -   -0.672079748346351159104929516004D+00,
c$$$     -   -0.860751470982527741630754268259D+00,
c$$$     -   -0.833059322343442213441093509931D+00,
c$$$     -    0.832692243496268321693103747431D+00,
c$$$     -   -0.507925940827256844289635307401D+00,
c$$$     -    0.822112380956134179088142546192D+00,
c$$$     -   -0.244274942964827857592177959801D-01/
c$$$        data whts /  
c$$$     -    0.467064135341454428114992336163D-02,
c$$$     -    0.116033742684542782857296887627D-01,
c$$$     -    0.198216147550200036448301754826D-01,
c$$$     -    0.135867050988507866717230642088D-01,
c$$$     -    0.249305006440888477783762850573D-01,
c$$$     -    0.821426444360997299417694488346D-01,
c$$$     -    0.554763855636506224675033447356D-01,
c$$$     -    0.409372985712911146433019363735D-02,
c$$$     -    0.629350838190045692503106585095D-02,
c$$$     -    0.548258038243713026134376594405D-02,
c$$$     -    0.349233006344836920918676035553D-01,
c$$$     -    0.941274386432602829720820455335D-02,
c$$$     -    0.203695855368206030548429973710D-02,
c$$$     -    0.306252062145517055801437053442D-01,
c$$$     -    0.295328071434484323265889859264D-02,
c$$$     -    0.355835917978602201358187656056D-01,
c$$$     -    0.265502120189976867147423338148D-01,
c$$$     -    0.470647704580543884195766661599D-01,
c$$$     -    0.219348528296405657314519223871D-01,
c$$$     -    0.101043329772962240497865048031D-01,
c$$$     -    0.542666662903039390922099894915D-02,
c$$$     -    0.663339999614095769219762349826D-01,
c$$$     -    0.672722770258950579806762411118D-01,
c$$$     -    0.204781008460178827406656888109D-02,
c$$$     -    0.686118137615150552549725334458D-01,
c$$$     -    0.328981480943429552981932909667D-01,
c$$$     -    0.234758587114798987520083550203D-01,
c$$$     -    0.686035553921814732099073365547D-01,
c$$$     -    0.157434156966473188537526803315D-01,
c$$$     -    0.330409286480616541804928258746D-01,
c$$$     -    0.213533082528012502898891013147D-01,
c$$$     -    0.464358576541559404183592600424D-01,
c$$$     -    0.596180070958399435009402330678D-01,
c$$$     -    0.601892972898566276125517502087D-01,
c$$$     -    0.338418106912387720239097792629D-01,
c$$$     -    0.149500943473197142998254450259D-01,
c$$$     -    0.643439086210449237461365939570D-02,
c$$$     -    0.740342309167785207112549660994D-01,
c$$$     -    0.133114228978254653024648802795D-01,
c$$$     -    0.396693027656059846502177863921D-01,
c$$$     -    0.214958504437383654343686938041D-01,
c$$$     -    0.872495936303020036064939915323D-02,
c$$$     -    0.686272147418608885705098633484D-01,
c$$$     -    0.628769750420342676446143189957D-01,
c$$$     -    0.274538604906026994689812923336D-01,
c$$$     -    0.178090286596743012428300890543D-01,
c$$$     -    0.291060587890693650171093458640D-01,
c$$$     -    0.407502570690435379135817736924D-01,
c$$$     -    0.405800828630885320790415526068D-01,
c$$$     -    0.119756447561628585208656381230D-01,
c$$$     -    0.796835588542728205231635689346D-02,
c$$$     -    0.356280471389118798011247602815D-01,
c$$$     -    0.866200317996570442940361541139D-02,
c$$$     -    0.878735070049510411972510521057D-02,
c$$$     -    0.613760295240676026817464376080D-01,
c$$$     -    0.198769943863428936273187945770D-02,
c$$$     -    0.759323550922467135131122574168D-01,
c$$$     -    0.514877133792163189719906739681D-01,
c$$$     -    0.439864971677347225148882215996D-01,
c$$$     -    0.521437364681681607056059579294D-01,
c$$$     -    0.284590681893665090793730678549D-01,
c$$$     -    0.509740151218442210212467706468D-01,
c$$$     -    0.615131023107237455146164656561D-02,
c$$$     -    0.129194992330263366681147636423D-01,
c$$$     -    0.643424947247945821580030129662D-01,
c$$$     -    0.832946071491659625154585822257D-01,
c$$$     -    0.220634331299273532277367008174D-01,
c$$$     -    0.439148201848827357892262855227D-01,
c$$$     -    0.111545121521094958148106287014D-01,
c$$$     -    0.674240047388361630786359958184D-01,
c$$$     -    0.698262871162392733371002219968D-01,
c$$$     -    0.580671559991742422335682629943D-01,
c$$$     -    0.114860251852875600151770070292D-01,
c$$$     -    0.714816851758059566773453742668D-02,
c$$$     -    0.432755758481038590681722685371D-01,
c$$$     -    0.565427444805641617432340133316D-01,
c$$$     -    0.299007335265097036587517342118D-01,
c$$$     -    0.145910665251771298786442327725D-01,
c$$$     -    0.599859452983377762465014645347D-01,
c$$$     -    0.314775103167599832401738655891D-01,
c$$$     -    0.159173706993044540244111566836D-01,
c$$$     -    0.673863708634243132734137467944D-02,
c$$$     -    0.458205358007826376209101005778D-01,
c$$$     -    0.381039054455882096016916352383D-01,
c$$$     -    0.217899618504204236156719251201D-01,
c$$$     -    0.720149057681587922155602720458D-01,
c$$$     -    0.622165576935753937790288584924D-01,
c$$$     -    0.577916778755859939304631109829D-01,
c$$$     -    0.225333574669907997810078033788D-01,
c$$$     -    0.690556307583374982005133954438D-01,
c$$$     -    0.462677665976784598537531730984D-01,
c$$$     -    0.655784931456414774724552535121D-01,
c$$$     -    0.789272378651567801053371367378D-01,
c$$$     -    0.107827069022509663796105349386D-01,
c$$$     -    0.597206977409688782545967178656D-01,
c$$$     -    0.268391976152952502987699564456D-01,
c$$$     -    0.699891669161996062731725802168D-01,
c$$$     -    0.164647737830570274586245788109D-01,
c$$$     -    0.728731125791790867808579624736D-01,
c$$$     -    0.823024659193207175505593915908D-01,
c$$$     -    0.140163641492719397904179232859D-01,
c$$$     -    0.690167746540170669695674101967D-01,
c$$$     -    0.657570759620342462685127002476D-02,
c$$$     -    0.490836388040391691340693212726D-01,
c$$$     -    0.257358776454886387730812017386D-01,
c$$$     -    0.430969084144211346794772634065D-01,
c$$$     -    0.262490515662711306014348465039D-01,
c$$$     -    0.763384932026696798361021857891D-01,
c$$$     -    0.435411632062844995369099463939D-01,
c$$$     -    0.516760318065908080296096294818D-01/
c$$$
c$$$!
c$$$        nquad0  = 110
c$$$        do i=1,nquad0
c$$$        xs0(i)   = xs(i)
c$$$        ys0(i)   = ys(i)
c$$$        whts0(i) = whts(i)
c$$$        end do
c$$$        end

