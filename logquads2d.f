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
        ndegree = 16           ! degree of the discretization polynomials
        npoly   = 16           ! degree of polynomials to integrate
        nmax    = 999         ! maximum possible singular quadrature size
!
        nlege   = 20

        allocate(xslege(nlege),whtslege(nlege))
        call legequad(nlege,xslege,whtslege)
c     
c       Fetch the discretization quadrature rule for the interior
c
        call discquad8_info(ndisc)
        allocate(xsdisc(ndisc),ysdisc(ndisc),whtsdisc(ndisc))
        call discquad8(ndisc,xsdisc,ysdisc,whtsdisc)

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
c$$$        write(iw,0100) "subroutine bdyquad(nquad,xs,ys,whts)"
c$$$        write(iw,0100) "implicit double precision (a-h,o-z)"
c$$$        
c$$$        write(iw,0160) "xs",nquadbdy
c$$$        write(iw,0160) "ys",nquadbdy
c$$$        write(iw,0160) "whts",nquadbdy
c$$$        write(iw,0100) "integer          :: nquad"
c$$$
c$$$        write(iw,0300) "nquad     = ",nquadbdy
c$$$
c$$$        do i=1,nquadbdy
c$$$        write(iw,0400) i,xsbdy(i)
c$$$        end do
c$$$c
c$$$        do i=1,nquadbdy
c$$$        write(iw,0500) i,ysbdy(i)
c$$$        end do
c$$$c
c$$$        do i=1,nquadbdy
c$$$        write(iw,0600) i,whtsbdy(i)
c$$$        end do
c$$$c
c$$$        write(iw,0100) "end subroutine"
c$$$        write(iw,*)    ""
c$$$        write(iw,*)    ""
c$$$
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
        call prin2("x1 =*",x1,1)
        call prin2("x2 =*",x2,1)
        call prini("idisc = ",i,1)
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


        subroutine discquad4_info(nquad0)
        implicit double precision (a-h,o-z)
        nquad0 = 16
        end subroutine

        subroutine discquad4(nquad0,xs0,ys0,whts0)
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
        dimension xs(16),ys(16),whts(16)        
        data nquad / 16 /
        data xs /
     -   -0.956027600560075462347694788933D+00,
     -    0.929815286862440803718789192511D+00,
     -    0.892223564678120870833326759751D+00,
     -   -0.875078288406164754745392284601D+00,
     -    0.430042282185442338809611504792D+00,
     -    0.147072117948658060438560112936D+00,
     -    0.517224543236797267186205329340D+00,
     -   -0.617608569146113646697262587649D+00,
     -   -0.679600478458640091783669101011D+00,
     -   -0.119848579525984166507307267422D+00,
     -   -0.905027048066559570244726602685D+00,
     -    0.554885905046915094212302440168D+00,
     -    0.908172066176727446603431456657D+00,
     -   -0.556825827935659979745377926760D-02,
     -   -0.459749208333378670822243225290D+00,
     -    0.886336404883876926806507315815D+00/
        data ys /
     -    0.623891974991350207500175043566D+00,
     -    0.223729040537217545714312033919D+00,
     -   -0.884663190716232820476420983659D+00,
     -   -0.941907019735325735960997323269D+00,
     -    0.962522678718794920356188618345D+00,
     -   -0.951522915712244464114161392368D+00,
     -    0.332419023600223736047390895894D+00,
     -    0.206342440905574014849746201688D+00,
     -    0.933562062541595164333187322135D+00,
     -    0.697892270233722259924964170759D+00,
     -   -0.362979876806572156893862688771D+00,
     -   -0.607267117670451458606579317439D+00,
     -   -0.189661764287239675548585806913D+00,
     -   -0.204821403573513844279392987006D+00,
     -   -0.721687716994303957755070475279D+00,
     -    0.773401355150155816790218102816D+00/
        data whts /
     -    0.121523997711913912329114105370D+00,
     -    0.497351279826182199447311910489D-01,
     -    0.967443574225257284860304103158D-01,
     -    0.730919275495923500605154385337D-01,
     -    0.129854251766320287841926579909D+00,
     -    0.159316412555998718986560041523D+00,
     -    0.480190072341471356298417630500D+00,
     -    0.460269140601132612029640273605D+00,
     -    0.137321249345292210141211045516D+00,
     -    0.432009378449635573581236028211D+00,
     -    0.213626355941694638438908001357D+00,
     -    0.358415508730049096657215384805D+00,
     -    0.192421260606505018380049624193D+00,
     -    0.572852137418635216187785466471D+00,
     -    0.364779124817935795968909362170D+00,
     -    0.157849696758679264667749416472D+00/
!
        nquad0 = nquad
        do i=1,nquad0
        xs0(i)   = xs(i)
        ys0(i)   = ys(i)
        whts0(i) = whts(i)
        end do
        end



        subroutine discquad6_info(nquad0)
        implicit double precision (a-h,o-z)
        nquad0 = 31
        end subroutine



        subroutine discquad6(nquad0,xs0,ys0,whts0)
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
        dimension xs(31),ys(31),whts(31)
        data nquad / 31 /
        data xs /
     -   -0.944625824810328400097220980156D+00,
     -    0.510021108126636252644829927155D+00,
     -   -0.296805388721957857490554864519D+00,
     -   -0.688473828341166336439745504181D+00,
     -    0.977406232488079499172481806285D+00,
     -    0.943398543317638046544969297972D+00,
     -    0.117659617469254773832433254575D+00,
     -    0.945961392527173393143167619773D+00,
     -   -0.703505602747056303862163001455D+00,
     -   -0.940884654269382691616009018980D+00,
     -    0.382612607581477582355242843558D+00,
     -    0.962095182032909545431242281657D+00,
     -   -0.973568804110227815131243759904D+00,
     -    0.805437863558127740617073551153D+00,
     -   -0.339896744906177195227657563378D+00,
     -    0.822159652326719415378815163325D+00,
     -    0.715454737292619740144853902338D+00,
     -    0.746029289208297089824772983478D+00,
     -   -0.899886753137915393595870242441D+00,
     -   -0.750888522730222915731364417309D+00,
     -   -0.775748022508379487394664149306D-01,
     -    0.527734278251722104248379510099D+00,
     -   -0.263699901928653576820552369358D+00,
     -    0.455649276604126756130931793738D-01,
     -    0.112829528041904322700447559160D+00,
     -   -0.507182974573510641685746603058D+00,
     -    0.472165299088584685599061457537D+00,
     -   -0.926696713486177462501339373244D+00,
     -   -0.699367662186085608425821361993D+00,
     -   -0.404030291246551145746063385097D+00,
     -   -0.996636333512829437164945899109D+00/
        data ys /
     -    0.967908142963568102343526659405D+00,
     -    0.964224952716395283727160418983D+00,
     -    0.974984388023193962232338585641D+00,
     -   -0.942828817896867578502635409259D+00,
     -    0.520762979852821162913401281193D+00,
     -   -0.856181415781172113464587198679D+00,
     -    0.819342185786235067085191657510D+00,
     -    0.954055342088526123247372857380D+00,
     -    0.855509650010683431311897012549D+00,
     -    0.620393640485109796051350393915D+00,
     -   -0.831965801581379725335817648934D+00,
     -   -0.298936729484762389677487356149D+00,
     -   -0.952013855657131336815282748659D+00,
     -    0.794719668399697684915162472364D+00,
     -    0.596715394562419413633952012326D+00,
     -    0.147977023679258497320367299712D+00,
     -   -0.986094240408864232403313614683D+00,
     -   -0.589755076067487173957861232201D+00,
     -   -0.730369572500830869437066583266D+00,
     -   -0.409741111077372079247582281358D+00,
     -   -0.955934491414041653786876536896D+00,
     -    0.532477731177373595919621637582D+00,
     -   -0.300108318454116105457807645800D-01,
     -   -0.495322228351262577071274758136D+00,
     -    0.261319407879339506592835297289D+00,
     -   -0.219521416033494450118633573984D+00,
     -   -0.173853560709117142583013543701D+00,
     -    0.360759667502738839830150485878D-01,
     -    0.303517187891841770088856393092D+00,
     -   -0.735116779963307464167268986373D+00,
     -   -0.359228829776019284296084038399D+00/
        data whts /
     -    0.197260226068972666678777683350D-01,
     -    0.626737380671689874555569259489D-01,
     -    0.613518806386189031559385961436D-01,
     -    0.664918858022887350010520325445D-01,
     -    0.515410238020115042249408715363D-01,
     -    0.514506983066765135935574030942D-01,
     -    0.183797853981396252354906341066D+00,
     -    0.212551849820616492525965368049D-01,
     -    0.119737356389486665378887851436D+00,
     -    0.797752310953998215508134079573D-01,
     -    0.168577455075676856293766964045D+00,
     -    0.744515333084728894316050115398D-01,
     -    0.136401611261762472553770038991D-01,
     -    0.103444038954085266223903447835D+00,
     -    0.241035994692925060038052008511D+00,
     -    0.184947178031112158153082434469D+00,
     -    0.372056442638186432885960240786D-01,
     -    0.175572914385626907067874757311D+00,
     -    0.792812849912252355029520690102D-01,
     -    0.155740524061018505746012085882D+00,
     -    0.832653778963040044824730432843D-01,
     -    0.222646146748456416810489593387D+00,
     -    0.233565947177789956381056420356D+00,
     -    0.280940781883578971575834821845D+00,
     -    0.286637715572602681476843429880D+00,
     -    0.139002096048884042371071562367D+00,
     -    0.284196834187627151368808053861D+00,
     -    0.959770424925273440941064509955D-01,
     -    0.203437820629562715978453296124D+00,
     -    0.193270354343497481843318488988D+00,
     -    0.253622784570251659801952974621D-01/

        nquad0 = nquad
        do i=1,nquad0
        xs0(i)   = xs(i)
        ys0(i)   = ys(i)
        whts0(i) = whts(i)
        end do
        end



        subroutine discquad8_info(nquad0)
        implicit double precision (a-h,o-z)
        nquad0 = 52
        end subroutine


        subroutine discquad8(nquad0,xs0,ys0,whts0)
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

        nquad  =  52
        norder =  16
        data xs /
     -    0.957554509813249197942398042693D+00,
     -    0.979330992017623798587011036063D+00,
     -   -0.987085483520060945881713631930D+00,
     -   -0.704279961107218585236737288408D+00,
     -    0.971544394461095061006586686075D+00,
     -    0.985011680722398532724323001832D+00,
     -   -0.219367757975985364794477604470D+00,
     -    0.498944112170079744334245664623D+00,
     -   -0.906520778956612745344348409254D+00,
     -   -0.102528374944668490578399536476D+00,
     -   -0.983729546473877887571822907815D+00,
     -    0.788904256999374427172456961353D+00,
     -    0.973600397370681310282576272887D+00,
     -   -0.689037782434947495784182229481D+00,
     -    0.496182974309953525620254159812D+00,
     -   -0.542087794615958778229182632878D+00,
     -    0.584187533191460054922483278393D+00,
     -    0.881601427067899800407838228161D+00,
     -    0.178304072724625885333130995799D+00,
     -   -0.414374416860569611459932364074D+00,
     -   -0.897193390455823786783714815923D+00,
     -   -0.965570902849928995079983635027D+00,
     -   -0.954114420773726162457426103023D+00,
     -   -0.339698025168779285643452843753D-01,
     -   -0.847569595854621340338256030321D+00,
     -   -0.867696973724043817968767890867D+00,
     -    0.721164245616815290129093670321D+00,
     -   -0.824466970013716910641043608525D+00,
     -    0.448847300931005590594281708273D+00,
     -    0.811552593251842743512852260108D+00,
     -    0.724702309014034732430700909742D+00,
     -   -0.261302838754444370771267240916D+00,
     -    0.819697394267450482658423424565D+00,
     -   -0.432051366493588674992034637738D+00,
     -    0.204390537802567409585801599974D+00,
     -    0.894417630244007432561677504958D+00,
     -   -0.697456051702624454909688171643D+00,
     -   -0.730647516219118592270338903520D+00,
     -   -0.979165343771992634212313588379D+00,
     -   -0.155598525674148427840509337215D+00,
     -   -0.145340340869547273706882121400D+00,
     -    0.152891431394838675855735357685D+00,
     -    0.346422607889851303646861014372D+00,
     -    0.921910418260776252853512082660D+00,
     -    0.588018244422309181961810913250D+00,
     -    0.868518649299822496255326789881D-01,
     -   -0.746463498791393077563927145218D+00,
     -    0.395771100100659619115317212764D+00,
     -   -0.492692596171395354359148206460D+00,
     -   -0.649762628196244183912733491804D+00,
     -   -0.416973868652741496915007666769D+00,
     -    0.634932688888604544723966628519D+00/
        data ys /
     -   -0.787200674099907618380642970506D+00,
     -    0.825764114890971349130844374281D+00,
     -    0.692489976437942170970551924825D+00,
     -    0.514973137146967950162756241337D+00,
     -   -0.977016489760589265254115288390D+00,
     -    0.266283208285061082771419135087D+00,
     -   -0.971309329608316112764831886935D+00,
     -   -0.994488236887355508188687738794D+00,
     -   -0.302447737699151232840224977871D+00,
     -    0.984426629465580940445244945114D+00,
     -   -0.604588873578206040049079947726D+00,
     -    0.925037310677987583005164504050D+00,
     -   -0.374434365986859130565104148056D+00,
     -    0.976663501622027194044850081756D+00,
     -    0.969994386031517553577402154882D+00,
     -   -0.882838481515814798800565931336D+00,
     -   -0.810019463948698877473784023591D+00,
     -   -0.480449025399324351155880948499D-01,
     -    0.869299474738741981350219753831D+00,
     -    0.883758118729167793068891134166D+00,
     -    0.430172239360402875330040276951D+00,
     -   -0.945132730863158878207663832187D+00,
     -    0.955885981643134283302521560305D+00,
     -   -0.755813751930538656880701182721D+00,
     -    0.832226956598935831401762039375D+00,
     -   -0.801235621053988062130185415550D+00,
     -    0.266076326678950577664069564463D+00,
     -   -0.198923237768674753052396081658D+00,
     -    0.669607106869525218327051240712D+00,
     -   -0.936150597872202984698219363879D+00,
     -    0.810249583690233315873922619388D+00,
     -   -0.503866957441679328892427483122D+00,
     -   -0.589653082677998106327590822379D+00,
     -    0.371818871007451841155850924523D+00,
     -   -0.917810503449965337325289922573D+00,
     -    0.577468131106776500935131555491D+00,
     -   -0.562819049983100920605806124737D+00,
     -    0.110041112679988868333171153458D+00,
     -    0.733429828417554744604674411966D-01,
     -    0.733171290929375873456319752226D-01,
     -    0.669017087865854463482874679613D+00,
     -    0.394204222014114932253795604836D+00,
     -   -0.578765025667195974813592014490D+00,
     -    0.977266218867860935786363559807D+00,
     -    0.456132785336496515279850944319D+00,
     -   -0.280669785752676896311970438851D+00,
     -   -0.978925596044971335131754599028D+00,
     -    0.347604910379124927527778479400D-01,
     -   -0.200935882866960343022796883242D+00,
     -    0.664537890403039300313869459391D+00,
     -   -0.763448874769108177738986534326D+00,
     -   -0.294309043922069762683673068508D+00/
        data whts /
     -    0.318099735292673412265282888126D-01,
     -    0.197338355623312648644083020564D-01,
     -    0.202523150183529289328450269193D-01,
     -    0.190617657211528746377447119322D-01,
     -    0.640706942327689571522234001541D-02,
     -    0.304747058661833854020360066254D-01,
     -    0.421787891660591448746810884771D-01,
     -    0.188677697532113423464199051393D-01,
     -    0.701738900041109729852543714606D-01,
     -    0.316165777059210704713477340250D-01,
     -    0.252203770235624471421810479358D-01,
     -    0.229342633719494452970378014252D-01,
     -    0.372303391656891438842608196094D-01,
     -    0.268274883026631729941286890110D-01,
     -    0.381109765360411403082768474440D-01,
     -    0.620021459219487165455934700230D-01,
     -    0.921325348266060918950155285210D-01,
     -    0.942086964151624259734397089700D-01,
     -    0.980360421795991635138662918722D-01,
     -    0.865927962223520739427619862962D-01,
     -    0.794791605654557619208414134131D-01,
     -    0.139832834513040214332132023316D-01,
     -    0.143194829209547315081214148924D-01,
     -    0.111328341194518959986946308875D+00,
     -    0.523865260136279153224882113513D-01,
     -    0.569917889532386935605081046973D-01,
     -    0.115512980797399048573112778658D+00,
     -    0.319506951886203866458590645283D-01,
     -    0.116565480706866905108999647071D+00,
     -    0.380238103246922433912988748651D-01,
     -    0.645711424294891294467356433341D-01,
     -    0.138705532990469960755864658709D+00,
     -    0.903601441996669756511467259062D-01,
     -    0.161438939654590019615712551220D+00,
     -    0.715864653522617508735840557970D-01,
     -    0.755402807932625056500193549114D-01,
     -    0.119107551193515722513064724674D+00,
     -    0.131494023658942462492750676707D+00,
     -    0.359989036433419581190192228070D-01,
     -    0.189804290122648857795632241849D+00,
     -    0.150585001248252059599298870096D+00,
     -    0.181402067366777085348946937591D+00,
     -    0.149702282013638381297968735337D+00,
     -    0.126064783741625001047299421299D-01,
     -    0.652462854414066817479314609180D-01,
     -    0.181308316306324354907027704447D+00,
     -    0.229928635861077217696419637032D-01,
     -    0.182619570392471250045052114741D+00,
     -    0.163709796789011897365691503613D+00,
     -    0.985846965801445441339754567271D-01,
     -    0.607980442358753460681147419771D-01,
     -    0.147423421795519124293651725552D+00/
!
        nquad0 = nquad
        do i=1,nquad0
        xs0(i)   = xs(i)
        ys0(i)   = ys(i)
        whts0(i) = whts(i)
        end do
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine discquad10_info(nquad0)
        implicit double precision (a-h,o-z)
        nquad0 = 78
        end subroutine


        subroutine discquad10(nquad0,xs0,ys0,whts0)
        implicit double precision (a-h,o-z)
        dimension xs(78),ys(78),whts(78)
        dimension xs0(1),ys0(1),whts0(1)
c
c       Return to the user a 78-point quadrature which integrates
c       polynomials of order 20 on the rectangle [-1,1] x [-1,1].
c
c       The quadrature is accurate to around 30 digits.  
c
c       It is intended for use in discretizing integral operators by
!       representing their solutions as 10th order polynomials.
c
c                      Output Parameters:
c
c   nquad0 - this will be set to the length of the quadrature rule
c   (xs0,ys0) - these user-supplied arrays will contain the coordinates
c       of the quadrature nodes
c   whts0 - this user-supplied array will contain the (positive)
c       weights
c   
        data xs /
     -    0.926555975905755513915971020626D+00,
     -   -0.913024526647190478369028753619D+00,
     -   -0.988383918890483819383402716217D+00,
     -   -0.897909422586895589416728775740D+00,
     -    0.918250528282333966739170174725D+00,
     -    0.993143912399823407705012257604D+00,
     -    0.601563907173806964601862595067D+00,
     -    0.591759998035315642428834249077D+00,
     -    0.989273837398036769501509838829D+00,
     -    0.884733440486114138685366797227D+00,
     -    0.987599205146504368430067331523D+00,
     -   -0.540417259542395087683575775109D+00,
     -   -0.778721260506808284854511539355D+00,
     -   -0.889304384214424968764066350354D+00,
     -   -0.839511846147327511251256435860D+00,
     -   -0.612685594580635033592251142141D+00,
     -   -0.941505382885246589850443105423D+00,
     -    0.984048738732380799043733621615D+00,
     -   -0.150010524775347637804222872163D+00,
     -   -0.877201282592437277489172004856D+00,
     -   -0.836917200379954706480569355748D+00,
     -   -0.945990308809321678878553206313D+00,
     -   -0.456193598944486462387801958428D+00,
     -    0.574676524396181193725809306587D+00,
     -   -0.446114543415753682566439186327D+00,
     -   -0.410309846957920086422436798277D+00,
     -   -0.642916721460959665305071493064D+00,
     -    0.807284370961669538482761738836D+00,
     -    0.221915152687594509474002166572D+00,
     -   -0.653296443512545709106091297676D+00,
     -   -0.818353338238848651312936934649D+00,
     -    0.147870626481583843302945488476D+00,
     -    0.123099056887539677578056155333D-01,
     -    0.797550310031424215456924128415D+00,
     -    0.919586820989622809256771707413D+00,
     -   -0.180461794924694458508021614405D+00,
     -   -0.296338789612619901800456842131D+00,
     -   -0.199396171849393871743703895950D+00,
     -   -0.986659079083316701998132121168D+00,
     -   -0.642046340295263110253769590402D+00,
     -    0.643200129865208274234781541213D+00,
     -    0.984615106335576606425482055215D+00,
     -    0.977416813702975436606890120051D+00,
     -    0.344272670575593763050503230385D+00,
     -    0.333447895749259321026740268278D+00,
     -    0.454657277404616748499685536652D+00,
     -    0.754189994881143297810597417724D-01,
     -   -0.253606236270413835788494038084D+00,
     -    0.947241867923324739557767719148D-01,
     -   -0.800249996718763461100192874774D+00,
     -    0.930962764468901076407877765860D+00,
     -   -0.981942622386643710915337719703D+00,
     -   -0.688121808834076768852645297433D+00,
     -   -0.461031209491973810695035990513D+00,
     -   -0.149213663175234621638448027073D+00,
     -    0.924011183305361011036081377753D+00,
     -    0.398529190358020520979892730693D+00,
     -   -0.428407343854623858635750309209D+00,
     -   -0.240172369916866889617323319703D+00,
     -    0.276893339778031360449601997239D+00,
     -    0.524397383009391746127601919110D+00,
     -    0.932996045706880613226821264688D+00,
     -    0.665912752284963279216215882303D-01,
     -    0.414344270105987592196560959318D+00,
     -    0.811236955485436692199089514073D+00,
     -   -0.996366768246143206855615424925D+00,
     -   -0.935245785900152428541161254202D+00,
     -   -0.953404675232094786009673726234D+00,
     -    0.731134516240156292380218685231D+00,
     -    0.791098009763760864638204388590D+00,
     -   -0.991652442640997640154153133121D+00,
     -    0.169860983241556642305791642017D+00,
     -   -0.679450346661891347463572565658D+00,
     -    0.628826776673926873172011341350D+00,
     -    0.795431813190396868786694476380D+00,
     -   -0.279211619091827389292393778096D+00,
     -   -0.359790108177762311921762039076D-01,
     -    0.640016632659688312578046148144D+00/
        data ys /
     -   -0.993648204713861939238672212463D+00,
     -    0.987130947373865906151551613803D+00,
     -   -0.602417236771653182297198661648D-01,
     -   -0.990078266619169016546422409097D+00,
     -    0.992170806540107267359623443082D+00,
     -    0.649923047829817170968616730613D+00,
     -   -0.988992904425822674075772433170D+00,
     -    0.988741360723661080483120087239D+00,
     -   -0.661885609912338212297293041116D+00,
     -   -0.589648677073402235776469461567D-01,
     -   -0.196936485442616018184760091396D+00,
     -   -0.122559070264894807019501360355D+00,
     -    0.662032776392395025277887792183D-01,
     -    0.354671122862718363330850358394D+00,
     -    0.524180112287403829686349341006D+00,
     -   -0.974131307081767833213995019516D+00,
     -    0.749589811011490065755033408537D+00,
     -   -0.932279757497518808467932279191D+00,
     -   -0.196546387804985927687699786736D+00,
     -   -0.158610129929298288536740967721D+00,
     -   -0.570090901257351846150866303548D+00,
     -    0.221950836855751208462952484345D+00,
     -    0.538549475095439271953139876694D+00,
     -    0.496023794868792503828407698271D-02,
     -   -0.525189588437537864756702545208D+00,
     -    0.638128797305795367434545215062D-01,
     -    0.973174302216364122273736045956D+00,
     -    0.636291109246586358535297256753D+00,
     -    0.755662801563743066816869966502D+00,
     -   -0.728893181073822141871191225155D+00,
     -    0.890352120949259641203425576239D+00,
     -    0.983768482394454366251318774449D+00,
     -   -0.881099630592665102691290128583D+00,
     -   -0.938254611968303487364339023792D+00,
     -   -0.818072503959121905042384177423D+00,
     -   -0.699605098071910988400212582964D+00,
     -    0.380377809890716863114554215598D+00,
     -    0.701876823794319868632209288516D+00,
     -    0.922512441980898081256106451930D+00,
     -    0.312088852070325662571832665018D+00,
     -    0.409464575912686427596815236802D+00,
     -    0.935471458253307913795964078839D+00,
     -    0.180536973692840059604701186057D+00,
     -   -0.254882234217954227083602228643D+00,
     -    0.256651183881562284476244668566D+00,
     -    0.590757774362036772889814959872D+00,
     -   -0.484127648067765600304856347011D+00,
     -   -0.315320832661566670016852413906D+00,
     -   -0.431583010604228764300746544117D-02,
     -   -0.898182370078303273440198962144D+00,
     -    0.443683219921730535543288873402D+00,
     -   -0.934585457473875045694258304872D+00,
     -   -0.338372338028090124253582118659D+00,
     -    0.859739870516645506635578825267D+00,
     -    0.229574735334836015384532216951D+00,
     -    0.813705754968955761386414286595D+00,
     -    0.918347935922854845634119508798D+00,
     -   -0.863372449216151511561984614266D+00,
     -   -0.972979523109948122687616656014D+00,
     -   -0.730531871814222705802408358610D+00,
     -   -0.538233467070080911935183752313D+00,
     -   -0.433744951160720492513584623634D+00,
     -    0.491961233650953300785348873058D+00,
     -   -0.919973053174258516383550200903D+00,
     -   -0.631832993350726321376846035408D+00,
     -   -0.614597420166648087605235998476D+00,
     -   -0.778706398697101333607457081641D+00,
     -   -0.384711310713382345096775351173D+00,
     -   -0.301992354047231478827733035938D+00,
     -    0.933076746705474555811151779519D+00,
     -    0.516742215995684131300885888414D+00,
     -   -0.984743330226400918784720287115D+00,
     -    0.712775567048347239291976978788D+00,
     -    0.802181206568358148679391864004D+00,
     -    0.218967157653944152722332329318D+00,
     -    0.974524379113034143837389288894D+00,
     -    0.889459742434036328682773954604D+00,
     -   -0.805114446946467495054184414376D+00/
        data whts /
     -    0.460398595861064633418696704358D-02,
     -    0.761981217725659315597072679948D-02,
     -    0.155924709094069706701650855913D-01,
     -    0.716586388891790303799899498080D-02,
     -    0.526266753746648903021446758392D-02,
     -    0.107865338911094719843097112108D-01,
     -    0.122175263684921885035498962934D-01,
     -    0.130775951976592540191123899761D-01,
     -    0.121024260530215339155875150335D-01,
     -    0.576442971891347930299879265007D-01,
     -    0.149468753071995635499063061831D-01,
     -    0.742182180618808783789147239361D-01,
     -    0.690386282414469928791148134137D-01,
     -    0.150992007043352441065476034512D-01,
     -    0.569828043914472636259376100818D-01,
     -    0.220598733776956574239045675450D-01,
     -    0.305875204898204716488291758425D-01,
     -    0.659761613211303482784720798000D-02,
     -    0.794631304025332771485512002248D-01,
     -    0.482381609615147683981473373810D-01,
     -    0.560628969869226906439018217012D-01,
     -    0.370805709261253717277026020856D-01,
     -    0.832720055534993977223415210205D-01,
     -    0.111990643751910176714581963227D+00,
     -    0.968715291537771994559898466139D-01,
     -    0.928463057447403776714401231772D-01,
     -    0.220827906907636425080015706828D-01,
     -    0.571110347426015148553188724989D-01,
     -    0.785591434368593265064095861479D-01,
     -    0.665066120437950197703467782347D-01,
     -    0.349738739205334483340792814121D-01,
     -    0.195628847391265330117503407276D-01,
     -    0.609709091573250221404088827577D-01,
     -    0.263941764820063195742716651007D-01,
     -    0.277412948144514851853156718941D-01,
     -    0.930653899270807974590881345928D-01,
     -    0.505802816084033647114931410026D-01,
     -    0.902714685338996637312218789196D-01,
     -    0.734959758039693907967725417393D-02,
     -    0.902156944357686813537200815775D-01,
     -    0.752546875047434689306540421734D-01,
     -    0.645211804593323326523477559890D-02,
     -    0.224993251302016941087488270874D-01,
     -    0.122249023275837379907653918609D+00,
     -    0.123525549055220192558784351722D+00,
     -    0.841115976235702164349964331024D-01,
     -    0.120275578509792191701324699799D+00,
     -    0.729878387319412348545042247352D-01,
     -    0.125790054971517970891912350695D+00,
     -    0.347501266599985644230936856815D-01,
     -    0.398335233888660415541625779497D-01,
     -    0.780248179369299494120229227259D-02,
     -    0.817854670335271590780833911506D-01,
     -    0.612492738964380146350539946757D-01,
     -    0.100271790169694746696502317604D+00,
     -    0.285958403173166589883542971699D-01,
     -    0.451817449534254396623857920046D-01,
     -    0.596293058565145837147183649744D-01,
     -    0.272786591846369709736068578805D-01,
     -    0.887608148317828538233538132608D-01,
     -    0.970732820662645982129183822257D-01,
     -    0.402725161009794698084794892423D-01,
     -    0.116473193198808009930519127783D+00,
     -    0.438077373993211489642494288781D-01,
     -    0.567522420499845428874791803417D-01,
     -    0.852224079854122891845218600231D-02,
     -    0.286348317795871802983826609974D-01,
     -    0.309914362637898733638857569559D-01,
     -    0.884403895910984790510191368257D-01,
     -    0.283776953197384880624365263482D-01,
     -    0.132074163939985573458568693335D-01,
     -    0.183765460047622103909850631511D-01,
     -    0.669495918528905386322349324354D-01,
     -    0.595245260233196770255340689194D-01,
     -    0.663586975360361693276491488040D-01,
     -    0.259381564635439561146616954376D-01,
     -    0.561535580929431549219569156586D-01,
     -    0.590468306586911377731231779118D-01/
c

!
        nquad0  = 78
        do i=1,nquad0
        xs0(i)   = xs(i)
        ys0(i)   = ys(i)
        whts0(i) = whts(i)
        end do
        end





        subroutine discquad12_info(nquad0)
        implicit double precision (a-h,o-z)
        nquad0 = 110
        end subroutine


        subroutine discquad12(nquad0,xs0,ys0,whts0)
        implicit double precision (a-h,o-z)
        dimension xs(110),ys(110),whts(110)
        dimension xs0(1),ys0(1),whts0(1)
c
c       Return to the user a 110-point quadrature which integrates
c       polynomials of order 24 on the rectangle [-1,1] x [-1,1].
c
c       The quadrature is accurate to around 30 digits.  
c
c       It is intended for use in discretizing integral operators by
!       representing their solutions as 10th order polynomials.
c
c                      Output Parameters:
c
c   nquad0 - this will be set to the length of the quadrature rule
c   (xs0,ys0) - these user-supplied arrays will contain the coordinates
c       of the quadrature nodes
c   whts0 - this user-supplied array will contain the (positive)
c       weights
c   
        data xs /
     -    0.578338131358755468487869321594D+00,
     -    0.989917932324099050193935280546D+00,
     -    0.935689289820520311538094755839D+00,
     -    0.983134715216754805520926423259D+00,
     -    0.655161144093537672312943803856D+00,
     -   -0.110655216873635288392847688042D+00,
     -    0.246417459828942791372291267539D+00,
     -    0.964261169671806015073488453988D+00,
     -   -0.995285733385343391129818036007D+00,
     -   -0.996959304328293547935604733653D+00,
     -    0.291964090180750006566476813415D+00,
     -    0.436474150480087436291296847570D+00,
     -   -0.988451227611893992447351666609D+00,
     -    0.854401157727382250299705972319D+00,
     -    0.997484321550589929006927719894D+00,
     -    0.426214296126862241761072890499D+00,
     -    0.335636159365554643416700936507D+00,
     -    0.581753765112829079263542995293D+00,
     -    0.261371395194623688224803542602D+00,
     -    0.988112658439639188183769353580D+00,
     -    0.799058303313613530936235971794D+00,
     -   -0.423919799542618553207034516322D+00,
     -   -0.449541585027222162092327733859D+00,
     -    0.979806685879784382380267436753D+00,
     -   -0.263405985666378587070449056169D+00,
     -   -0.928450805373773067572289494746D+00,
     -   -0.920579276148290826392617922165D+00,
     -   -0.284983631424843087968839785081D+00,
     -    0.900998414604663572100711353926D+00,
     -   -0.419509142882353243170099450846D+00,
     -    0.675822837792625085091399220424D+00,
     -    0.735681191563475753084312948407D+00,
     -   -0.453478061289327256192346142702D+00,
     -    0.555493551187429532058350313632D-01,
     -   -0.925029461405186415061359901533D+00,
     -   -0.973940346898449815257749678651D+00,
     -   -0.992661640104854562230262001359D+00,
     -   -0.280598552997141223535701410448D+00,
     -   -0.208322083873715094186667195882D+00,
     -    0.866237547645044546197848788279D+00,
     -    0.942959436926699251298770526917D+00,
     -    0.992544332517046611637362757333D+00,
     -    0.598372776764861072367919875398D+00,
     -   -0.397596220040245249340498511314D+00,
     -    0.888614730314901988671395179838D+00,
     -   -0.978273253758586041013486821925D+00,
     -    0.927389020158637257500213299041D+00,
     -    0.232103426640469657140462847738D+00,
     -    0.863069122428051282316981950237D+00,
     -   -0.951821132544472406145221015163D+00,
     -   -0.591586279399231990293722508995D+00,
     -   -0.239721511965783405462322952786D-01,
     -    0.805890687952985921077826935259D+00,
     -   -0.573316634590597478698706401730D+00,
     -    0.187974897099413937911881686013D+00,
     -   -0.989509698301903320846304125640D+00,
     -    0.830328340431802634804746259359D-01,
     -    0.520923995230169896390785027991D+00,
     -    0.722878150557529495361318262679D+00,
     -    0.765816601077174972405451813030D+00,
     -    0.951229615832098554558395257173D+00,
     -   -0.224474035366275805741653041318D+00,
     -    0.987650080825455349799756428923D+00,
     -    0.105423826663379597357442149442D+00,
     -    0.647780389081604045542491539788D+00,
     -    0.649290070268204179940653614305D-01,
     -   -0.736276656474779423883100577380D+00,
     -   -0.850045294127175556695683809176D+00,
     -    0.481556051687227869629777718461D+00,
     -   -0.849200801801683920856031068658D-01,
     -    0.378304037321419630815181545324D+00,
     -   -0.614255425944048938145133069935D+00,
     -   -0.951186534829595337625972056628D+00,
     -   -0.879187890938707422006105340295D+00,
     -   -0.379455141861730646118569600399D-01,
     -   -0.255106202264953620380597368537D+00,
     -   -0.414399313887484609715703345882D+00,
     -    0.963543994552770888753039928399D+00,
     -   -0.727602091741514801576195878309D+00,
     -    0.787486392425711862331676808552D+00,
     -   -0.972759112566351730906621834530D+00,
     -   -0.876305298533984338149618694348D+00,
     -   -0.745504268814309511060346500870D+00,
     -    0.548211538672122460312125794769D+00,
     -   -0.737213485817654628361592252519D+00,
     -   -0.125326434120443540587448473954D+00,
     -    0.612993333873251227134602949480D+00,
     -    0.432456244550277552477654644761D+00,
     -   -0.928017979157815620525813590876D+00,
     -   -0.554685072129007856874500072879D+00,
     -   -0.831831291873353540702696777272D+00,
     -   -0.703512895863964609238058219688D+00,
     -    0.180972124958704441231501687514D+00,
     -    0.916070519975212002115586978105D+00,
     -    0.769047961466390775133842493843D+00,
     -   -0.842383659619120590325760815849D+00,
     -    0.398281741509930552640480353081D+00,
     -   -0.214497141392241964162208667720D+00,
     -   -0.349074145811955135731180658943D-01,
     -    0.451057560309855462105171713851D+00,
     -    0.140000024909246638512141042587D+00,
     -   -0.573216496218163359969158441824D+00,
     -   -0.992171640813981647205734356980D+00,
     -   -0.714281298687964104998705917887D+00,
     -    0.830580684053784210745149293812D+00,
     -   -0.578781755885851655359537709569D+00,
     -   -0.850629641237496353856080928874D+00,
     -    0.253882448139543283775802815075D+00,
     -   -0.601035980825033639140114887954D+00,
     -   -0.836443659001792179101334934813D+00/
        data ys /                   
     -    0.987738783592016334452219720187D+00,
     -   -0.123388216730251988275521550371D+00,
     -    0.422110729518207013442453559873D+00,
     -    0.237769962442095514118137254945D+00,
     -    0.929102264061020734929011246267D+00,
     -   -0.107260917199191562060753458568D+00,
     -   -0.959158397915358604981830271782D-01,
     -    0.982797664268167295738524211060D+00,
     -    0.313774998985064534494468961096D+00,
     -   -0.317238564257656165488980039181D+00,
     -   -0.312543254103086799390898910916D-01,
     -    0.985691129020807628808962976142D+00,
     -    0.980919078160093524898693099669D+00,
     -    0.387175233375450949741751303962D+00,
     -    0.906038894426891699191143627260D+00,
     -    0.868306423239176119490301793295D+00,
     -   -0.928106794565659315583192306837D+00,
     -    0.744163198498799072850844534611D+00,
     -    0.937810609614595459929852120836D+00,
     -    0.570414502339396910399206502924D+00,
     -   -0.994014922128567818494895923783D+00,
     -   -0.188477226610190176227740057209D+00,
     -    0.280914416382924182468088748753D+00,
     -   -0.988330872888656696583153047318D+00,
     -    0.475806107325555093866630219529D+00,
     -    0.223536716121078921610220466969D+00,
     -   -0.687441241709610641724325809186D+00,
     -   -0.369994685298573693491972230576D+00,
     -    0.915530771990217897095448553916D+00,
     -    0.925913647821213555527251032256D+00,
     -   -0.944714101932633175341835599528D+00,
     -    0.581352672207178451671059911405D+00,
     -    0.646612354548894872457418550110D+00,
     -   -0.710552837606063190940753149275D+00,
     -   -0.262480472706355853534623023149D+00,
     -    0.522500159270853780844318809988D+00,
     -    0.789453832632186040491295778106D+00,
     -    0.931287775083829933758388432955D-01,
     -   -0.985857956607461113249518780559D+00,
     -   -0.126012238168086025688124284833D+00,
     -   -0.748732237378739585158648426002D+00,
     -   -0.584036919585183168808324459813D+00,
     -    0.396528381830265291520066731885D+00,
     -   -0.672546115406545413584088341617D+00,
     -    0.659049056991202310536239244943D+00,
     -   -0.153363432819385804992589077532D-01,
     -    0.739454148445587794032387252192D-01,
     -   -0.822668384654059063177408000164D+00,
     -   -0.562772420187030489515014425002D+00,
     -   -0.910960497780475203940899285397D+00,
     -   -0.993694528707898245844497071740D+00,
     -   -0.915323382062898979143099818013D+00,
     -    0.984717518085326622276102045589D+00,
     -    0.992683889847213119449333755093D+00,
     -    0.766909998724911116996054531668D+00,
     -   -0.981075237253124068566503869941D+00,
     -    0.130552564745526061148162909878D+00,
     -    0.881994001143697114894829847713D-01,
     -   -0.716280973472452864819435717504D+00,
     -   -0.323637496704597571909849704558D+00,
     -   -0.363113777528609807944340090781D+00,
     -   -0.819681310895534634596201890011D+00,
     -   -0.881690905352544661568014188064D+00,
     -    0.983961474824500893543521508550D+00,
     -   -0.794648699349394462242807374977D-01,
     -   -0.318187652351975523326585783255D+00,
     -    0.942334899824859615071497663191D+00,
     -    0.449656115591339183077378783485D+00,
     -   -0.986768549100099562461546323643D+00,
     -    0.304686756319859013889272141386D+00,
     -    0.246166961314600731204413582202D+00,
     -    0.463445760417664210362944425566D+00,
     -    0.914972558895649833325113600804D+00,
     -   -0.985037427249410909256584851165D+00,
     -    0.893963429553482941339319652301D+00,
     -    0.786960832273606417031100874541D+00,
     -   -0.938314390943547927378444399693D+00,
     -    0.783060359432202353852236109814D+00,
     -    0.219518496620202337653009002605D+00,
     -    0.818839442076717018050462031341D+00,
     -   -0.529214869519173086462562343218D+00,
     -    0.986629502965302453144760327043D+00,
     -    0.662476113605986935558934010442D+00,
     -   -0.843254622868324991600838470364D+00,
     -   -0.940745696197973394055954412527D+00,
     -   -0.553201285837351425096359394364D+00,
     -   -0.500766814827076628881161742015D+00,
     -   -0.677017164449588001756865085874D+00,
     -    0.682299172236535595676168772503D+00,
     -   -0.476734657107389201442736668104D+00,
     -   -0.475258535426877974982331253288D+00,
     -   -0.252822879633909571313139000577D+00,
     -    0.436768997253338672918098799180D+00,
     -   -0.951935465791726957984414602942D+00,
     -    0.179947740115926907550884899948D+00,
     -   -0.830364751493477094394829329507D+00,
     -    0.594572176801124747963056499065D+00,
     -    0.980897417867983074415155635543D+00,
     -    0.626912899441285698397178671150D+00,
     -   -0.284241096717721054709091202614D+00,
     -   -0.983870030100289638089940861533D+00,
     -    0.921209208775992472487614855103D-02,
     -   -0.783562420691980944334249904800D+00,
     -   -0.672079748346351159104929516004D+00,
     -   -0.860751470982527741630754268259D+00,
     -   -0.833059322343442213441093509931D+00,
     -    0.832692243496268321693103747431D+00,
     -   -0.507925940827256844289635307401D+00,
     -    0.822112380956134179088142546192D+00,
     -   -0.244274942964827857592177959801D-01/
        data whts /  
     -    0.467064135341454428114992336163D-02,
     -    0.116033742684542782857296887627D-01,
     -    0.198216147550200036448301754826D-01,
     -    0.135867050988507866717230642088D-01,
     -    0.249305006440888477783762850573D-01,
     -    0.821426444360997299417694488346D-01,
     -    0.554763855636506224675033447356D-01,
     -    0.409372985712911146433019363735D-02,
     -    0.629350838190045692503106585095D-02,
     -    0.548258038243713026134376594405D-02,
     -    0.349233006344836920918676035553D-01,
     -    0.941274386432602829720820455335D-02,
     -    0.203695855368206030548429973710D-02,
     -    0.306252062145517055801437053442D-01,
     -    0.295328071434484323265889859264D-02,
     -    0.355835917978602201358187656056D-01,
     -    0.265502120189976867147423338148D-01,
     -    0.470647704580543884195766661599D-01,
     -    0.219348528296405657314519223871D-01,
     -    0.101043329772962240497865048031D-01,
     -    0.542666662903039390922099894915D-02,
     -    0.663339999614095769219762349826D-01,
     -    0.672722770258950579806762411118D-01,
     -    0.204781008460178827406656888109D-02,
     -    0.686118137615150552549725334458D-01,
     -    0.328981480943429552981932909667D-01,
     -    0.234758587114798987520083550203D-01,
     -    0.686035553921814732099073365547D-01,
     -    0.157434156966473188537526803315D-01,
     -    0.330409286480616541804928258746D-01,
     -    0.213533082528012502898891013147D-01,
     -    0.464358576541559404183592600424D-01,
     -    0.596180070958399435009402330678D-01,
     -    0.601892972898566276125517502087D-01,
     -    0.338418106912387720239097792629D-01,
     -    0.149500943473197142998254450259D-01,
     -    0.643439086210449237461365939570D-02,
     -    0.740342309167785207112549660994D-01,
     -    0.133114228978254653024648802795D-01,
     -    0.396693027656059846502177863921D-01,
     -    0.214958504437383654343686938041D-01,
     -    0.872495936303020036064939915323D-02,
     -    0.686272147418608885705098633484D-01,
     -    0.628769750420342676446143189957D-01,
     -    0.274538604906026994689812923336D-01,
     -    0.178090286596743012428300890543D-01,
     -    0.291060587890693650171093458640D-01,
     -    0.407502570690435379135817736924D-01,
     -    0.405800828630885320790415526068D-01,
     -    0.119756447561628585208656381230D-01,
     -    0.796835588542728205231635689346D-02,
     -    0.356280471389118798011247602815D-01,
     -    0.866200317996570442940361541139D-02,
     -    0.878735070049510411972510521057D-02,
     -    0.613760295240676026817464376080D-01,
     -    0.198769943863428936273187945770D-02,
     -    0.759323550922467135131122574168D-01,
     -    0.514877133792163189719906739681D-01,
     -    0.439864971677347225148882215996D-01,
     -    0.521437364681681607056059579294D-01,
     -    0.284590681893665090793730678549D-01,
     -    0.509740151218442210212467706468D-01,
     -    0.615131023107237455146164656561D-02,
     -    0.129194992330263366681147636423D-01,
     -    0.643424947247945821580030129662D-01,
     -    0.832946071491659625154585822257D-01,
     -    0.220634331299273532277367008174D-01,
     -    0.439148201848827357892262855227D-01,
     -    0.111545121521094958148106287014D-01,
     -    0.674240047388361630786359958184D-01,
     -    0.698262871162392733371002219968D-01,
     -    0.580671559991742422335682629943D-01,
     -    0.114860251852875600151770070292D-01,
     -    0.714816851758059566773453742668D-02,
     -    0.432755758481038590681722685371D-01,
     -    0.565427444805641617432340133316D-01,
     -    0.299007335265097036587517342118D-01,
     -    0.145910665251771298786442327725D-01,
     -    0.599859452983377762465014645347D-01,
     -    0.314775103167599832401738655891D-01,
     -    0.159173706993044540244111566836D-01,
     -    0.673863708634243132734137467944D-02,
     -    0.458205358007826376209101005778D-01,
     -    0.381039054455882096016916352383D-01,
     -    0.217899618504204236156719251201D-01,
     -    0.720149057681587922155602720458D-01,
     -    0.622165576935753937790288584924D-01,
     -    0.577916778755859939304631109829D-01,
     -    0.225333574669907997810078033788D-01,
     -    0.690556307583374982005133954438D-01,
     -    0.462677665976784598537531730984D-01,
     -    0.655784931456414774724552535121D-01,
     -    0.789272378651567801053371367378D-01,
     -    0.107827069022509663796105349386D-01,
     -    0.597206977409688782545967178656D-01,
     -    0.268391976152952502987699564456D-01,
     -    0.699891669161996062731725802168D-01,
     -    0.164647737830570274586245788109D-01,
     -    0.728731125791790867808579624736D-01,
     -    0.823024659193207175505593915908D-01,
     -    0.140163641492719397904179232859D-01,
     -    0.690167746540170669695674101967D-01,
     -    0.657570759620342462685127002476D-02,
     -    0.490836388040391691340693212726D-01,
     -    0.257358776454886387730812017386D-01,
     -    0.430969084144211346794772634065D-01,
     -    0.262490515662711306014348465039D-01,
     -    0.763384932026696798361021857891D-01,
     -    0.435411632062844995369099463939D-01,
     -    0.516760318065908080296096294818D-01/

!
        nquad0  = 110
        do i=1,nquad0
        xs0(i)   = xs(i)
        ys0(i)   = ys(i)
        whts0(i) = whts(i)
        end do
        end




        subroutine discquad16_info(nquad0)
        implicit double precision (a-h,o-z)
        nquad0 = 188
        end subroutine


        subroutine discquad16(nquad0,xs0,ys0,whts0)
        implicit double precision (a-h,o-z)
        dimension xs(188),ys(188),whts(188)
        dimension xs0(1),ys0(1),whts0(1)
c
c       Return to the user a 110-point quadrature which integrates
c       polynomials of order 24 on the rectangle [-1,1] x [-1,1].
c
c       The quadrature is accurate to around 30 digits.  
c
c       It is intended for use in discretizing integral operators by
!       representing their solutions as 10th order polynomials.
c
c                      Output Parameters:
c
c   nquad0 - this will be set to the length of the quadrature rule
c   (xs0,ys0) - these user-supplied arrays will contain the coordinates
c       of the quadrature nodes
c   whts0 - this user-supplied array will contain the (positive)
c       weights
c

        nquad  = 188
        norder =  32
        data xs /
     -    0.971628594945475878217412946137D+00,
     -   -0.989957477577230979184700256137D+00,
     -    0.489634458718653300920306180247D+00,
     -    0.168021966344164038010328650093D+00,
     -   -0.797422627051304393546050393859D+00,
     -   -0.887153270040620287117072369148D+00,
     -    0.782579556122000642108923257730D+00,
     -    0.104847376057792407421756462961D+00,
     -   -0.966129747240556198599539322816D+00,
     -   -0.850704767243832917444349359238D+00,
     -   -0.820820145906190269873087266214D+00,
     -    0.532890074127029770849131782606D+00,
     -   -0.822661076255390302520630443130D+00,
     -   -0.827444467748535965551139017024D+00,
     -    0.915696662485388870859520445791D+00,
     -   -0.861625044632796415259599974065D+00,
     -   -0.944678525851313528921309200894D+00,
     -   -0.867497014918292683097686734438D+00,
     -    0.186984330644467867610258141980D+00,
     -    0.986181341622739327840402958934D+00,
     -   -0.997521558492986235257597240527D+00,
     -   -0.957849078475484914602576674380D+00,
     -   -0.703956435329431964378412826757D+00,
     -    0.613941898867860527708150546811D+00,
     -   -0.984994839548056600894791232357D+00,
     -    0.214433883876705379644957320266D+00,
     -    0.577174581833320085391067728895D+00,
     -    0.970025719530933212469215026937D+00,
     -   -0.896400189979763488547877465381D+00,
     -   -0.988643148185537626246520421400D+00,
     -   -0.195611134923602569165293569632D+00,
     -    0.642216854463945321453843044164D+00,
     -   -0.943601726114807499030707363000D+00,
     -   -0.315210680124965302398914276087D+00,
     -   -0.504718613521428541518069771218D-01,
     -   -0.713310122629020968594712598802D+00,
     -    0.823162905785238122891563343147D+00,
     -    0.968114403346128068911140389531D+00,
     -    0.528596119639507515895166087799D+00,
     -    0.578445351490592441145856574113D+00,
     -   -0.965434635674797091715179874284D+00,
     -   -0.593662457059902508767033211148D-01,
     -   -0.687757963818025502837065814672D+00,
     -   -0.461362034456402416577466729830D+00,
     -    0.149779295480257645479320563927D+00,
     -   -0.588780378274470241538904202487D+00,
     -   -0.889132039381478371199882552465D+00,
     -   -0.944267329280955325456788613314D+00,
     -   -0.994456197893206557193015443308D+00,
     -   -0.608415384885246565071615693408D+00,
     -   -0.181968836300819382790322606276D-01,
     -    0.920747128848583213363959936497D+00,
     -    0.669757207752001770432566244505D+00,
     -    0.457056753486206819981725282715D+00,
     -    0.862316088561894314086529652863D+00,
     -    0.662465624049201247339102837505D+00,
     -   -0.819109683427808159925350784177D+00,
     -    0.437964538648999893449563411988D+00,
     -    0.155264650889777232814574343416D+00,
     -   -0.707028316159832771660312191149D+00,
     -    0.535107136612525584792286484415D+00,
     -    0.797767885384697701762578391427D+00,
     -    0.581111937811430956536239695186D+00,
     -    0.720053675341704615543057230156D+00,
     -    0.923233181955500028160166139913D+00,
     -   -0.992198371945881569545555711743D+00,
     -   -0.816094436967410296089848373917D+00,
     -    0.821455371042079624540832083533D+00,
     -   -0.959386340935399320935541140618D+00,
     -    0.261554648570417762971638117538D+00,
     -    0.651419449550449855392687883293D+00,
     -    0.758095608101271644875529884721D+00,
     -    0.224290973806826688748383904751D+00,
     -   -0.993292307307978041364884794811D+00,
     -    0.941425798009569352535286087369D+00,
     -   -0.909201008107934992760918022943D+00,
     -   -0.921334450471203482294218997638D+00,
     -    0.683750122057768565383777918418D+00,
     -    0.118093476677007711023952150133D+00,
     -   -0.277450021110420151976245874231D+00,
     -   -0.593303581110903638932931865237D+00,
     -   -0.199637403994361204158403589472D+00,
     -    0.966397800506158226650123136930D+00,
     -   -0.943433708664524999395692057453D+00,
     -    0.704572353920228357126839955788D+00,
     -   -0.131882705282406958389017245279D-01,
     -    0.801115045442208794245490386549D+00,
     -   -0.972208473175189190233931937170D+00,
     -   -0.417795216363934877142565462894D+00,
     -    0.281171234988311746923731315762D+00,
     -   -0.318435839678858936839383152845D+00,
     -    0.740893791270691654397687041950D+00,
     -   -0.749256690439885625829740973201D+00,
     -    0.657289369541506803579845050815D+00,
     -    0.993518715934057169038732590330D+00,
     -   -0.970958422034872955394053253344D+00,
     -    0.996105638619342266293770688836D+00,
     -   -0.994156926718893020984408421129D+00,
     -   -0.879957105399176962995817681165D+00,
     -   -0.307684414976256442474895181236D+00,
     -    0.798871432899850733592177155369D+00,
     -   -0.864086646480548975123621736943D-01,
     -    0.226130496047335943359969860700D+00,
     -   -0.837393213202779782846753345884D+00,
     -   -0.222940089075422553594829303659D+00,
     -   -0.601357727431579944623771893020D-01,
     -    0.306788534123387280451689818183D+00,
     -    0.780107859499504183277250249027D+00,
     -   -0.769054941273208885532622133056D+00,
     -    0.106727250724502315654667849044D+00,
     -   -0.280072770924797637798508793421D-01,
     -   -0.863775636458261779591631221975D+00,
     -    0.696647737184254464861502762224D+00,
     -   -0.337985224799892109623694215390D+00,
     -    0.899485639342643584450392195432D+00,
     -   -0.245623498318869369552791280414D+00,
     -   -0.437751963115321193429794467979D+00,
     -    0.993854716111344024935787747516D+00,
     -   -0.529757769314950642721317298482D+00,
     -   -0.906686404613953952555580092136D+00,
     -    0.432323575947913547013458166119D+00,
     -   -0.796986486262183835599937884060D+00,
     -    0.686882334654929818288151897570D+00,
     -   -0.370345362016875602423800580837D+00,
     -    0.317728659006174393899737597765D+00,
     -   -0.181736857100818550740410636511D+00,
     -   -0.389695678272760620060948417910D+00,
     -    0.428062706426701146371308326987D+00,
     -   -0.714215388501702311434748974095D+00,
     -    0.926793681205596532000804777859D-01,
     -    0.941169899545314509608477144919D+00,
     -   -0.735613952195008010182295989165D+00,
     -    0.465871631208882568172403204097D+00,
     -    0.958968230897200409749595033876D+00,
     -    0.992278988687072122073179351462D+00,
     -    0.860846062787459140855068498949D+00,
     -   -0.211662607215764527161187272917D+00,
     -    0.299287534790141733501161507607D+00,
     -   -0.557145482221242721489978619870D+00,
     -   -0.848347752210869448695395602845D-01,
     -   -0.479307518063771899609401430521D+00,
     -    0.964701410637196654413482437379D+00,
     -   -0.592893607003086040028036184749D+00,
     -    0.862351094334537097219754557218D+00,
     -    0.560578934666832275603603072592D+00,
     -    0.380463540484907636671523716177D+00,
     -    0.390403016162788072497265585666D-01,
     -    0.991247733324228060590917613876D+00,
     -    0.852224370466867821079605884838D+00,
     -   -0.372216695476573763215306601817D+00,
     -    0.590284888669687236648062851282D+00,
     -    0.143389914929940161316092363906D+00,
     -   -0.684198121468401570211276166163D+00,
     -    0.901096076785672994527076904142D-01,
     -   -0.990818166672458793347984721378D+00,
     -   -0.993345224085789563638089455367D+00,
     -   -0.995335021658593899712348257581D+00,
     -    0.272841284512314757795120519696D+00,
     -    0.343007269691707565131022201340D+00,
     -    0.433039256652871904074617556509D+00,
     -    0.908478753979520399244357002435D+00,
     -   -0.915251749016685705916718381273D+00,
     -   -0.164472311070747477723249095483D+00,
     -    0.916200962748551766583599158973D+00,
     -    0.718274393695268490547241145412D+00,
     -    0.909160493507154488811418353956D+00,
     -    0.813155153966418750291530285391D+00,
     -   -0.789012537736087011890180022155D+00,
     -   -0.472391901433181106500150283779D+00,
     -    0.966534234789935737058112058867D+00,
     -    0.958775206995168294200749088559D+00,
     -   -0.581948349640437978669518266445D+00,
     -   -0.491537420147753848327583640916D+00,
     -   -0.284342557979218527937244071539D-01,
     -   -0.573674229876236924760504701138D+00,
     -    0.886073009953673429876298829881D+00,
     -    0.534227425078654814669135013981D+00,
     -   -0.544704137883904656228873936848D+00,
     -   -0.243482920654032685411832270342D+00,
     -    0.992826346814923716490088149196D+00,
     -    0.994957354795869051189419330148D+00,
     -    0.992592013807231436004652315735D+00,
     -   -0.399655304880205083760818453342D+00,
     -    0.385469346033457881093408739253D+00,
     -    0.776342798009406529804740688577D+00,
     -   -0.633318000149546543582249054534D+00,
     -   -0.661557014572268383907674881177D+00,
     -    0.880495721248826378349493822546D+00/
        data ys /
     -   -0.998369034925239725714792948975D+00,
     -    0.993695929523628105597989226104D+00,
     -    0.997007176764997464115212925469D+00,
     -   -0.179585558485888366906219651852D+00,
     -   -0.690125061379003522030729964475D+00,
     -   -0.787850256809995361533280335517D+00,
     -    0.740666614852151886863060039105D+00,
     -   -0.171500836234048568484682486511D+00,
     -    0.691421959314424377195532490238D+00,
     -    0.796557920455126953350509075732D+00,
     -   -0.837027015708230258478487719984D+00,
     -    0.911614678806193204271245029626D+00,
     -   -0.437137808575668699330585208110D+00,
     -   -0.995318986102753939202608887515D+00,
     -   -0.400766876205617612147926003727D+00,
     -    0.996361522440596968707570451604D+00,
     -   -0.984694141445670442931617983047D+00,
     -   -0.641332477137690840145383734655D+00,
     -    0.786215835447810192866553088130D+00,
     -   -0.229327894921061440882217974631D-02,
     -    0.592665378666608729514526481415D+00,
     -   -0.999311413065990669746291053845D-01,
     -   -0.807580104456997017396616579819D+00,
     -    0.983691607026311375523107126870D+00,
     -   -0.939084228803956092701746349032D+00,
     -    0.993986591359600635982232671580D+00,
     -    0.220966705701514180057769508215D-01,
     -    0.709598009338116530440686323013D+00,
     -    0.762616952231385204369408552579D+00,
     -    0.931150283953223063092637021284D+00,
     -   -0.518486111669525039920057297108D+00,
     -   -0.682418132295616107963592348692D-01,
     -    0.974541917804905567715795033583D+00,
     -   -0.883802643565011673409583077898D+00,
     -    0.986563887450678612624636097526D+00,
     -   -0.280727877597315305116137819826D+00,
     -   -0.727298577975867689174745934996D+00,
     -    0.995660348995830308197082103314D+00,
     -   -0.994980366313441631988588703707D+00,
     -   -0.643927136467449391231517537654D+00,
     -   -0.708392309753493037385166944862D+00,
     -   -0.344255083592265947150951536655D-02,
     -    0.233483466280374049972625418406D-01,
     -   -0.811129010894880398596765334483D+00,
     -   -0.751231871490618122785477169941D+00,
     -   -0.703072706869794789502895452734D+00,
     -    0.566767355706579217170775884202D-01,
     -    0.868164984275312326141952032092D+00,
     -   -0.993136207033875928970849759548D+00,
     -    0.995011094686703251627371918936D+00,
     -   -0.984433624113047335069729139373D+00,
     -    0.568333240572701585411148837603D+00,
     -    0.833352123275563964494552341369D+00,
     -   -0.762102895069957924759935443620D+00,
     -    0.462638514384074441688860287350D-02,
     -   -0.539459446248835082948347548758D+00,
     -    0.654309957287020686667562436499D+00,
     -   -0.139495667690005282036075262493D+00,
     -   -0.499389560306588312037554924090D+00,
     -   -0.569298986849735793746250278605D+00,
     -   -0.905239115344589699874266608650D+00,
     -    0.993809239809809539241584108530D+00,
     -    0.295382805731171870114830799099D+00,
     -   -0.396331346614153281325139200474D+00,
     -    0.984190390129749432782833301211D+00,
     -    0.522497544643483620054469527300D-01,
     -   -0.122774010256915392558463568368D+00,
     -    0.890245681122188838079603708060D+00,
     -    0.189097110261098374173937941540D+00,
     -    0.162730898116184367262302501158D-01,
     -   -0.975072773486698756044851481551D+00,
     -   -0.158119811820270085739759003105D+00,
     -    0.580718923026405645121076169766D+00,
     -    0.314797934782930975746434682535D+00,
     -    0.160488489102592789796091643504D+00,
     -    0.320770619130985899511760772206D+00,
     -   -0.554404273991085091467146506027D+00,
     -    0.611583696159057602503298445472D+00,
     -    0.946059224164432127116139530268D+00,
     -   -0.992427593977159776965091232162D+00,
     -   -0.423627056678482754575264461901D+00,
     -    0.941540761532559179146787543879D+00,
     -   -0.721209917808075310631237713698D+00,
     -   -0.875782266046259408555665267769D+00,
     -   -0.948745057689628367953372049533D+00,
     -   -0.858762100841861769969956978698D+00,
     -    0.494827805519423493114798327748D+00,
     -   -0.414183972115577866273891776498D+00,
     -   -0.306926345273706009058716655720D+00,
     -    0.881613002402191820815615185667D+00,
     -    0.991659139992684611352710882435D+00,
     -   -0.680867163267681900814211661779D+00,
     -    0.971632372195854775557412760561D+00,
     -   -0.811011221217562175705406786310D+00,
     -    0.833072334805514780622724122734D+00,
     -    0.458570111903246609970152599682D+00,
     -    0.572114423332503175019174605077D+00,
     -   -0.237724453553858898766491942139D+00,
     -   -0.948442067442303663962907337127D+00,
     -   -0.426306238262339323730031101785D+00,
     -   -0.993100345100978842316029313877D+00,
     -   -0.316437651397041103602962987067D+00,
     -   -0.993561162763999376618094902163D+00,
     -    0.428088618708999012626565558295D+00,
     -    0.165547583031786912906679250756D+00,
     -    0.645232967900181153515818385901D+00,
     -   -0.861063351228724219412709387013D+00,
     -    0.206519336768185939479344025533D+00,
     -    0.848786142198801892687086868835D+00,
     -    0.746428109873295051023788014878D+00,
     -    0.867703502492413089985457372665D+00,
     -    0.928505335901715867680896433011D+00,
     -    0.954390709264478950448763785820D-01,
     -   -0.673349336601932553920156771130D+00,
     -   -0.819618272700311151578237716435D+00,
     -    0.522444674968332243792166310394D+00,
     -   -0.961645842350324202862921554712D+00,
     -   -0.581604765362862689822101638060D+00,
     -    0.804005875871868047193377867012D+00,
     -   -0.281248793324432517067385376383D+00,
     -    0.800365241608789618544033793542D+00,
     -   -0.895157997228401146102769405210D+00,
     -    0.400168060763316252230986963659D+00,
     -    0.880754334236650498569448394031D+00,
     -   -0.625497738548402177388604860396D+00,
     -   -0.770130776168682332335612456810D+00,
     -    0.676914795804394799758864380860D+00,
     -    0.457974514767841723537277340028D+00,
     -   -0.968460200051585514488819702516D+00,
     -    0.485577353935787326852805498892D+00,
     -   -0.161048017822464609080900943748D+00,
     -    0.528784327864169472688747788326D+00,
     -   -0.470686608273260550326690994911D+00,
     -    0.910076254008774056122238647441D+00,
     -   -0.275279208833346490488142009408D+00,
     -   -0.307010812279389611633194304665D+00,
     -    0.780588682076078347248257805599D+00,
     -   -0.323065439405447418064992151764D+00,
     -   -0.146974975950550515525321463606D+00,
     -    0.356841792894321704615618738696D+00,
     -    0.962096252933476245681189133693D+00,
     -   -0.440190042996473465818060911991D+00,
     -   -0.906235568806002064246546550037D+00,
     -    0.660730376615870848711022929083D+00,
     -    0.711543962794527857201580016534D+00,
     -    0.965771008433301253833219948640D+00,
     -   -0.411191672278298602120647959935D+00,
     -   -0.962661091576205507894856329409D+00,
     -   -0.935062955056663338499866453447D+00,
     -    0.326085497934841300868865197381D+00,
     -   -0.291787002865130029960938889653D+00,
     -   -0.936302787645576088521942596219D+00,
     -    0.728811188081685585312986549480D+00,
     -    0.188505849498732480496399554804D+00,
     -    0.804633708624992320031804782649D+00,
     -   -0.819229139223150644126465912506D+00,
     -   -0.581262203088990913942411805662D+00,
     -    0.333457852052728223640493378460D+00,
     -    0.658722080591985084512696038631D+00,
     -    0.166679881258513117366724887215D+00,
     -   -0.601329480887077494843615838473D+00,
     -    0.566698105362594884737198895500D+00,
     -   -0.941273876574715148475841386776D+00,
     -   -0.976875763592054103557416139682D+00,
     -    0.949509287110179034269363659431D+00,
     -    0.810683922588981878580735640234D+00,
     -   -0.513974849523078541237735362108D+00,
     -    0.207630590401974451871957386718D+00,
     -   -0.558340497533061687019453486859D+00,
     -    0.427010679261117082427500131820D+00,
     -   -0.903731672292979717232206626997D+00,
     -    0.601874335872723193758707844707D+00,
     -    0.470777388467930083685202066940D+00,
     -   -0.634745020502697873282029906496D+00,
     -   -0.994569132370943622719789121285D+00,
     -    0.960265274303890593140257509371D+00,
     -    0.541425662991211728487664745514D+00,
     -    0.183565118838870738714653059057D+00,
     -   -0.158897706582780016695332021897D+00,
     -    0.964732190645674312411756882472D+00,
     -   -0.829348252879587607080404997213D+00,
     -    0.264179820016262327357218431230D+00,
     -    0.128603436215834205731591717017D-01,
     -   -0.963975485827517794503425311884D+00,
     -   -0.876016253590722339109555097228D+00,
     -    0.911283070857786451887744988546D+00,
     -    0.356308569550871739539272296815D+00,
     -    0.337642888445784299911139142814D+00/
        data whts /
     -    0.772164310747589288661764454101D-03,
     -    0.660064450522851867902267474101D-03,
     -    0.236121585627141003940923013671D-02,
     -    0.907197342437337988599817418296D-02,
     -    0.182774237224691795403517140486D-01,
     -    0.138702844675477531364661378761D-01,
     -    0.207301045321937713927340514749D-01,
     -    0.501134313482801605286091402388D-01,
     -    0.929680042268993744966470083592D-02,
     -    0.568915017676560106512058994410D-02,
     -    0.307484114380141223516537831675D-02,
     -    0.199621019377636625398303682257D-01,
     -    0.291814808163075142509262440159D-01,
     -    0.271421409342639952875366856874D-02,
     -    0.787239864700580885983160058716D-02,
     -    0.223245434337196747131035632436D-02,
     -    0.298636596769298680628912023944D-02,
     -    0.915495763477282042655902607265D-02,
     -    0.170382442763766945425754127867D-01,
     -    0.923720412442843846639738650629D-02,
     -    0.279558214259665002720958275310D-02,
     -    0.158119139436708642690825019273D-01,
     -    0.212199049298484397035267720197D-01,
     -    0.629405086613055228383887665945D-02,
     -    0.313775960832353580427079769010D-02,
     -    0.505129805265601944265267596223D-02,
     -    0.309466185915650893733380898530D-01,
     -    0.995359339376932601942855640079D-02,
     -    0.118648324611851254795613404938D-01,
     -    0.281640558965646353661829464055D-02,
     -    0.396467940572824516428137165816D-01,
     -    0.170677948391397183844623031737D-01,
     -    0.418081879629264966001217438112D-02,
     -    0.225379845370345278120799447462D-01,
     -    0.881304396168358488241175308322D-02,
     -    0.364732899188916918609738096213D-01,
     -    0.167852473881966980994788714833D-01,
     -    0.991906230378265769858672922245D-03,
     -    0.399613170831561263948814193166D-02,
     -    0.287356853738861953807316650550D-01,
     -    0.104988817834572888262022570058D-01,
     -    0.580595872789286284100761043298D-01,
     -    0.423564752719225937870641179736D-01,
     -    0.274000203451342428906511522659D-01,
     -    0.392079295203624360055677490263D-01,
     -    0.295490110722369130368978858676D-01,
     -    0.251836158417227082619821077294D-01,
     -    0.876901486252727985487864558992D-02,
     -    0.464570479802088345929875275015D-03,
     -    0.396400454952843847035315944757D-02,
     -    0.947590724660427519415463459159D-02,
     -    0.154583960093190229554973663561D-01,
     -    0.229361562008868050942783425727D-01,
     -    0.323010421940855870765978280750D-01,
     -    0.307445240413363744278804854037D-01,
     -    0.212771001559165029750223398606D-01,
     -    0.199847113148597867316685652335D-01,
     -    0.508212541668541424119443131770D-01,
     -    0.429838091941127150913199866578D-01,
     -    0.298234121523158051430159991525D-01,
     -    0.214287588372754037920955667088D-01,
     -    0.338724260803877808151963182763D-02,
     -    0.378007000092743099558756288962D-01,
     -    0.257385232748294638787124481346D-01,
     -    0.178882115192611764419037812499D-02,
     -    0.562248927759439129030845324250D-02,
     -    0.318558728574239089232315518215D-01,
     -    0.142301011249254299046327072246D-01,
     -    0.127167661332777412560758011132D-01,
     -    0.567206110235633799791637179740D-01,
     -    0.623751907241926010105494559399D-02,
     -    0.362848630684551329605916396804D-01,
     -    0.274924380938066765214453956679D-01,
     -    0.437675067852634472970295761968D-02,
     -    0.194157792541161709328442884109D-01,
     -    0.177794800451406876545241224584D-01,
     -    0.167717905365235090734129984420D-01,
     -    0.269489443206133093889481654482D-01,
     -    0.179580512958202221686426149331D-01,
     -    0.597597124188934637349675040260D-02,
     -    0.361565943625401459200797592720D-01,
     -    0.189992428264509395687200365583D-01,
     -    0.103890557267874952481787519063D-01,
     -    0.850470263319485697915418427992D-02,
     -    0.807469119675315681478621306549D-02,
     -    0.294311075158598968066211946968D-01,
     -    0.256791161557012319504544351505D-01,
     -    0.119592724425706168678544508046D-01,
     -    0.438424573616685574087312270269D-01,
     -    0.243312889100162722975607449983D-01,
     -    0.624534273275028046050520716635D-02,
     -    0.175489870030273215810182038965D-01,
     -    0.922039972728926952421272406402D-02,
     -    0.238814663222319833126582558160D-01,
     -    0.317613242700144539273104512916D-02,
     -    0.103621581257179907676400004453D-01,
     -    0.374114331792550577859698902718D-02,
     -    0.495115638622877489276610285313D-02,
     -    0.774553042732075706474652765961D-02,
     -    0.282124328319861956592751949030D-01,
     -    0.365349700248060555053271819403D-02,
     -    0.480812919728486678276133280079D-01,
     -    0.528667919757964182520823099609D-02,
     -    0.216918730866828843487703145529D-01,
     -    0.540421940616123444582916837087D-01,
     -    0.449372716987458427153500466923D-01,
     -    0.287229162695657045397275427403D-01,
     -    0.314773478171010183940116557820D-01,
     -    0.181202151583241946307157137666D-01,
     -    0.275921788893128499298873093322D-01,
     -    0.295925571109517015336392422853D-01,
     -    0.107767554508033767081094001127D-01,
     -    0.243280024506546197338165556105D-01,
     -    0.352182672141083986598222384499D-01,
     -    0.140429109271858451448642656033D-01,
     -    0.495461341813123564661284734817D-01,
     -    0.143006884390671737467687217322D-01,
     -    0.457498452391486481055596233854D-02,
     -    0.283488049965509395778540465593D-01,
     -    0.232632604235591774485641807986D-01,
     -    0.266005195446347013027049168405D-01,
     -    0.126956486886046105551325740978D-01,
     -    0.282551342334235836058203755033D-01,
     -    0.249560197496227379765552716964D-01,
     -    0.431513879251221154646176456234D-01,
     -    0.342219644134472562497570450491D-01,
     -    0.402090504264477626248135001268D-01,
     -    0.367781276487692740321892826409D-01,
     -    0.102280347428226471541110667459D-01,
     -    0.454049355632879186109627125748D-01,
     -    0.197940661371052570869464378521D-01,
     -    0.250162912367807857906102460424D-01,
     -    0.460875608959172740341260338047D-01,
     -    0.687216458400996940923631206163D-02,
     -    0.611832216917669960453434877783D-02,
     -    0.251501509722114722983315894064D-01,
     -    0.364048163376002486148033051680D-01,
     -    0.530410327942896171489841069415D-01,
     -    0.468624894653097830533743848472D-01,
     -    0.542589706698586184585576109485D-01,
     -    0.137018728216973404810374596694D-01,
     -    0.124150355333534018134758433016D-01,
     -    0.200613245194608587568203451913D-01,
     -    0.140864580165286142032432351290D-01,
     -    0.263190418482717230195807135412D-01,
     -    0.133432165151366038988989193848D-01,
     -    0.261654532385627537505851446813D-01,
     -    0.182002192798788693225638714087D-02,
     -    0.851066613002449563767993753994D-02,
     -    0.468092558349130121122877678510D-01,
     -    0.423549081932679576635810505660D-01,
     -    0.201664795225856256678922853174D-01,
     -    0.256001363375995885703126616173D-01,
     -    0.582958379832489513334030124658D-01,
     -    0.370292051408899565548680477510D-02,
     -    0.320783284581308795812831637165D-02,
     -    0.394651657896258623026873885111D-02,
     -    0.517202094985451276068367185486D-01,
     -    0.344699366754002710257110104703D-01,
     -    0.489245722659052194353294867275D-01,
     -    0.185382739808170225671512010112D-01,
     -    0.161286011490926901192042262438D-01,
     -    0.175938171405805317146924240912D-01,
     -    0.454685828477311453387930641645D-02,
     -    0.111343438348715733050933571722D-01,
     -    0.140581466594114653410712626436D-01,
     -    0.225343643301922410293199701909D-01,
     -    0.346660385089352852382210467986D-01,
     -    0.352728824213369400198618537034D-01,
     -    0.126983975514364290433770860481D-01,
     -    0.706399394990600150134240881347D-02,
     -    0.295657942462505625661990379990D-01,
     -    0.360373560882986231364845861808D-01,
     -    0.447963050109747992180193652814D-01,
     -    0.428244408707459520092456978351D-02,
     -    0.615273146963075859184743337959D-02,
     -    0.267896265682232786671146764045D-01,
     -    0.473661074720743382796651528008D-01,
     -    0.547511093461197519890341294140D-01,
     -    0.155842308740360525044143064010D-02,
     -    0.286915882881202351458551191746D-02,
     -    0.570204131956538147031932417651D-02,
     -    0.528310865304675279102865019509D-01,
     -    0.141311478406941693411606985389D-01,
     -    0.139917114596254999321591576277D-01,
     -    0.185761365466077311888348228768D-01,
     -    0.385495509042297401494724384284D-01,
     -    0.248266274027469610088819811063D-01/
!
        nquad0  = 188
        do i=1,nquad0
        xs0(i)   = xs(i)
        ys0(i)   = ys(i)
        whts0(i) = whts(i)
        end do
        end
