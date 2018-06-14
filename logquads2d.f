c
c
c       Construct tensor product quadrature rules for discretizing integral 
c       operators given on planar domains and which have kernels of the form
c
c             log(|x-y|) \psi(y) + \phi(y)
c
c       with \psi and \phi smooth OR 
c
c            \nabla_x log (|x-y|) \cdot F(y) + \phi(y)
c
c       with \phi a smooth function and F a smooth vector field.
c
        implicit double precision (a-h,o-z)
        double precision, allocatable :: xsdisc(:),ysdisc(:),whtsdisc(:)
        double precision, allocatable :: xsbdy(:),ysbdy(:),whtsbdy(:)

        double precision, allocatable :: xsall(:,:),ysall(:,:)
        double precision, allocatable :: whtsall(:,:)

        double precision, allocatable :: xslege(:),whtslege(:)
        double precision, allocatable :: xs(:),ys(:),whts(:)
c
        integer, allocatable          :: nquadall(:)

        eps   = 1.0d-16
        npoly = 4
c
        nlege = npoly+1
        allocate(xslege(nlege),whtslege(nlege))
        call legequad(nlege,xslege,whtslege)
c
c       Fetch the discretization quadrature for the interior
c
        call discquad_info(ndisc)
        allocate(xsdisc(ndisc),ysdisc(ndisc),whtsdisc(ndisc))
        call discquad(ndisc,xsdisc,ysdisc,whtsdisc)
        call prinf("ndisc    = *",ndisc,1)
        call prin2("xsdisc   = *",xsdisc,ndisc)
        call prin2("ysdisc   = *",ysdisc,ndisc)
        call prin2("whtsdisc = *",whtsdisc,ndisc)
c
c       Put together the boundary quadrature
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
        ysbdy(idx)   = 1
        whtsbdy(idx) = whtslege(i)
        end do

        do i=1,nlege
        idx = idx + 1
        xsbdy(idx)   = 1
        ysbdy(idx)   = xslege(i)
        whtsbdy(idx) = whtslege(i)
        end do

        do i=1,nlege
        idx = idx + 1
        xsbdy(idx)   = xslege(i)
        ysbdy(idx)   = -1
        whtsbdy(idx) = whtslege(i)
        end do
c
        call prinf("nquadbdy = *",nquadbdy,1)
        call prin2("xsbdy    = *",xsbdy,nquadbdy)
        call prin2("ysbdy    = *",ysbdy,nquadbdy)
        call prin2("whtsbdy  = *",whtsbdy,nquadbdy)
c
c       Write out the discretization quadrature rule and the boundary
c       quadrature rule.
c
        iw = 1001
        open(iw,FILE='quads.f90')

 0100 format(A)
 0300 format(A,I3.3)
 0400 format("xs(",I3.3,")   = ",D44.36)
 0500 format("ys(",I3.3,")   = ",D44.36)
 0600 format("whts(",I3.3,") = ",D44.36)
c
        write(iw,0100) "subroutine discquad(nquad,xs,ys,whts)"
        write(iw,0100) "implicit double precision (a-h,o-z)"
        write(iw,0100) "double precision :: xs(:),ys(:),whts(:)"
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
        write(iw,0100) "double precision :: xs(:),ys(:),whts(:)"
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
c       Build the first set of singular quadrature rules
c
        nn = ndisc+nquadbdy

        allocate(xsall(300,nn),ysall(300,nn),whtsall(300,nn))
        allocate(nquadall(nn))
c
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(xs,ys,whts,i,x1,x2,ier,nquad)
        allocate(xs(10000),ys(10000),whts(10000))
!$OMP   DO
        do i=1,ndisc
        x1    = xsdisc(i)
        x2    = ysdisc(i)
        call diagquad(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
        xsall(1:nquad,i)   = xs(1:nquad)
        ysall(1:nquad,i)   = ys(1:nquad)
        whtsall(1:nquad,i) = whts(1:nquad)
        nquadall(i)        = nquad
        end do
!$OMP   END DO
!$OMP   END PARALLEL
c
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(xs,ys,whts,i,x1,x2,ier,nquad)
        allocate(xs(10000),ys(10000),whts(10000))
!$OMP   DO
        do i=1,nquadbdy
        x1    = xsbdy(i)
        x2    = ysbdy(i)
        call diagquad(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
        xsall(1:nquad,i+ndisc)   = xs(1:nquad)
        ysall(1:nquad,i+ndisc)   = ys(1:nquad)
        whtsall(1:nquad,i+ndisc) = whts(1:nquad)
        nquadall(i+ndisc)        = nquad
        end do
!$OMP   END DO
!$OMP   END PARALLEL
c
c       Write the first set of singular rules to the disc
c
        iw = 1001
        open(iw,FILE='quads.f90',STATUS='OLD',ACCESS='append')

        write(iw,0100) "subroutine singquads(nquadsall,xsall," //
     -    "ysall,whtsall)"
        write(iw,0100) "implicit double precision (a-h,o-z)"
        write(iw,0100) "double precision :: xsall(300,:)"
        write(iw,0100) "double precision :: ysall(300,:)"
        write(iw,0100) "double precision :: whtsall(300,:)"
        write(iw,0100) "integer          :: nquadsall(:)"
c
 1100 format("nquadsall(",I3.3,")    = ",I3.3)
 1200 format("xsall(",I3.3,",",I3.3,")   = ",D44.36)
 1300 format("ysall(",I3.3,",",I3.3,")   = ",D44.36)
 1400 format("whtsall(",I3.3,",",I3.3,") = ",D44.36)

        do i=1,ndisc+nquadbdy
        nq = nquadall(i)
        write(iw,1100) i,nq
        do j=1,nq
        write(iw,1200) nq,i,xsall(j,i)
        write(iw,1300) nq,i,ysall(j,i)
        write(iw,1400) nq,i,whtsall(j,i)
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

c
c       Write the second set of singular rules to the disk
c
        nn = ndisc+nquadbdy
c
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(xs,ys,whts,i,x1,x2,ier,nquad)
        allocate(xs(10000),ys(10000),whts(10000))
!$OMP   DO
        do i=1,ndisc
        x1    = xsdisc(i)
        x2    = ysdisc(i)
        call diagquad2(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
        xsall(1:nquad,i)   = xs(1:nquad)
        ysall(1:nquad,i)   = ys(1:nquad)
        whtsall(1:nquad,i) = whts(1:nquad)
        nquadall(i)        = nquad
        end do
!$OMP   END DO
!$OMP   END PARALLEL
c
!$OMP   PARALLEL DEFAULT(SHARED) PRIVATE(xs,ys,whts,i,x1,x2,ier,nquad)
        allocate(xs(10000),ys(10000),whts(10000))
!$OMP   DO
        do i=1,nquadbdy
        x1    = xsbdy(i)
        x2    = ysbdy(i)
        call diagquad2(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
        xsall(1:nquad,i+ndisc)   = xs(1:nquad)
        ysall(1:nquad,i+ndisc)   = ys(1:nquad)
        whtsall(1:nquad,i+ndisc) = whts(1:nquad)
        nquadall(i+ndisc)        = nquad
        end do
!$OMP   END DO
!$OMP   END PARALLEL
c
c       Write the second set of singular rules to the disc
c
        iw = 1001
        open(iw,FILE='quads.f90',STATUS='OLD',ACCESS='append')

        write(iw,0100) "subroutine singquads2(nquadsall,xsall," //
     -    "ysall,whtsall)"
        write(iw,0100) "implicit double precision (a-h,o-z)"
        write(iw,0100) "double precision :: xsall(300,:)"
        write(iw,0100) "double precision :: ysall(300,:)"
        write(iw,0100) "double precision :: whtsall(300,:)"
        write(iw,0100) "integer          :: nquadsall(:)"
c

        do i=1,ndisc+nquadbdy
        nq = nquadall(i)
        write(iw,1100) i,nq
        do j=1,nq
        write(iw,1200) nq,i,xsall(j,i)
        write(iw,1300) nq,i,ysall(j,i)
        write(iw,1400) nq,i,whtsall(j,i)
        end do
        end do
c
        write(iw,0100) "end subroutine"
        write(iw,*)    ""
        write(iw,*)    ""

        close(iw)
        call flush(iw)
c
 2100 format(A," = ",I3.3)
c
        iw = 1001
        open(iw,FILE='quads.f90',STATUS='OLD',ACCESS='append')

        write(iw,0100) "subroutine discquad_info(nquad,nquadbdy," //
     -    "npoly)"
        
        write(iw,0100) "implicit double precision (a-h,o-z)"
        write(iw,2100) "nquad",ndisc
        write(iw,2100) "nquadbdy",nquadbdy
        write(iw,2100) "npoly",npoly
        write(iw,0100) "end subroutine"
        write(iw,*)    ""
        write(iw,*)    ""
        close(iw)
c
        end program



        subroutine diagquad(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension                     :: ab0(2,100)
        double precision, allocatable :: disc(:),coefs(:,:)
        double precision, allocatable :: xs0(:),whts0(:),rints(:)
        double precision, allocatable :: xs1(:),whts1(:),ab(:,:)
        double precision, allocatable :: xs2(:),whts2(:),errs(:)
        double precision              :: xs(1),ys(1),whts(1)
        double precision              :: verts(2,3)
        integer, allocatable          :: idxs(:)
        external funuser1,funuser2,funeval,funquad4,funquad5
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
c
 0100 continue
        a      = -1.0d0
        b      =  1.0d0
        nfuns1 = (npoly+1)*(npoly+2)/2 + npoly+1
c
        nints0    = 1
        ab0(1,1)  = a
        ab0(2,1)  = b
        kdisc     = kdisc+1
c
        call prinf("nfuns1 = *",nfuns1,1)
        call legedisc(ier,nints0,ab0,kdisc,eps,nfuns1,
     -    funuser1,npoly,x1,x2,par4,disc,ldisc,ncoefs,lkeep)
        if (ier .ne. 0) then
        call prinf("after legedisc, ier = *",ier,1)
        stop
        endif
c
c       Gram-Schmidt the input functions.
c
        allocate(coefs(ncoefs,nfuns1))
c
        call legegs(eps,disc,nfuns1,funuser1,npoly,x1,x2,par4,
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
        ngoal = (krank+1)/2+1
        if (nquad1 .gt. ngoal) goto 0100
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
        a      = -1.0d0
        b      =  1.0d0
c
        nints0    = 1
        ab0(1,1)  = a
        ab0(2,1)  = b
        kdisc     = kdisc+1
        nfuns2    = (npoly+1) + (npoly+1)
c
        call legedisc(ier,nints0,ab0,kdisc,eps,nfuns2,
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
c
c
c       Build the Gaussian quadrature.
c
        ifaccept = 0
        ngoal    = 0
c
        call gaussquad(eps,krank,rints,funeval,disc,coefs,krank,
     1    par4,nquad2,xs2,whts2,a,b,ngoal,ifaccept)
c
        deallocate(coefs,xs0,whts0)
        ngoal = (krank+1)/2+1
        if (nquad2 .gt. ngoal) goto 0200
        do i=1,nquad2
        idxs(i) = i
        end do
        call quicksort(nquad2,xs2,idxs)
        whts2(1:nquad2)  = whts2(idxs(1:nquad2))
c        xs2(1:nquad2) = xs2(1:nquad2) + x1

        call prin2("xs2    = *",xs2,nquad2)
        call prin2("whts2  = *",whts2,nquad2)
c
        do i=1,nquad2
        nquad       = nquad+1
        xs(nquad)   = xs2(i)
        ys(nquad)   = xs1(idisc)
        whts(nquad) = whts1(idisc)*whts2(i)
        end do
c
        end do
c
        call prina("*")
        call prina("=========================*")
        call prinf("final nquad = *",nquad,1)
        call prin2("final xs    = *",xs,nquad)
        call prin2("final ys    = *",ys,nquad)
        call prin2("final whts  = *",whts,nquad)
c
c       Test the resulting quadrature rule
c
        
        allocate(errs( (npoly+1)*(npoly+2) ) )
c
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
        call adaptri(ier,eps,verts,funquad4,x1,x2,i1,j1,
     1    val1,nquad01)
c
        verts(1,1) = -1
        verts(2,1) = -1
        verts(1,2) = -1
        verts(2,2) =  1
        verts(1,3) =  1
        verts(2,3) =  1
c
        call adaptri(ier,eps,verts,funquad4,x1,x2,i1,j1,
     1    val2,nquad02)
c
        val0  = val1+val2
c
        val  = 0
        do i=1,nquad
        x   = xs(i)
        y   = ys(i)
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
        call adaptri(ier,eps,verts,funquad5,x1,x2,i1,j1,
     1    val1,nquad01)
c
        verts(1,1) = -1
        verts(2,1) = -1
        verts(1,2) = -1
        verts(2,2) =  1
        verts(1,3) =  1
        verts(2,3) =  1
c
        call adaptri(ier,eps,verts,funquad5,x1,x2,i1,j1,
     1    val2,nquad02)
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
        call prin2("errs = *",errs,idx)
        call prina("=========================*")
        call prina("*")
c
        end



        subroutine diagquad2(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension                     :: ab0(2,100)
        double precision, allocatable :: disc(:),coefs(:,:)
        double precision, allocatable :: xs0(:),whts0(:),rints(:)
        double precision, allocatable :: xs1(:),whts1(:),ab(:,:)
        double precision, allocatable :: xs2(:),whts2(:),errs(:)
        double precision              :: xs(1),ys(1),whts(1)
        double precision              :: verts(2,3)
        integer, allocatable          :: idxs(:)
        external funuser3,funuser4,funeval,funquad6,funquad7,funquad5
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
c       where (x1,x2) is a  user-specified point in [-1,1] x [-1,1], P_k
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
        ier = 0
c
c       First build a quadrature for integrating the functions F_ij, G_ij
c       defined via
c
c                                1
c         F_ij(y2) = P_i(y2) \int   ( (x1-y1) / ( (x1-y1)^2 + (x2-y2)^2 ) P_j(y1) dy1
c                               -1
c
c       and
c
c                                1
c         G_ij(y2) = P_i(y2) \int   ( (x2-y2)/ ( (x1-y1)^2 + (x2-y2)^2 ) P_j(y1)  dy1
c                               -1
c
c       where 0<= i1+j1, j2 <= npoly and P_n denotes the Legendre polynomial of degree n.
c
        call prin2("in diagquad2, x1 = *",x1,1)
        call prin2("in diagquad2, x2 = *",x2,1)
        call prinf("in diagquad2, npoly = *",npoly,1)
c
        ldisc = 100 000 000
        allocate(disc(ldisc))
        allocate(xs1(1000),whts1(1000),rints(1000),idxs(1000))
        allocate(xs2(1000),whts2(1000))
c
        kdisc = 30
        nquad = 0
c
 0100 continue
        a         = -1.0d0
        b         =  1.0d0
        nfuns3    = (npoly+1) + (npoly+1)*(npoly+2)
c
        nints0    = 2
        ab0(1,1)  = a
        ab0(2,1)  = x2-eps
        ab0(1,2)  = x2+eps
        ab0(2,2)  = b

        kdisc     = kdisc+1
c
        call prinf("nfuns3 = *",nfuns3,1)
        call legedisc(ier,nints0,ab0,kdisc,eps,nfuns3,
     -    funuser3,npoly,x1,x2,par4,disc,ldisc,ncoefs,lkeep)
        if (ier .ne. 0) then
        call prinf("after legedisc, ier = *",ier,1)
        stop
        endif
c
c       Gram-Schmidt the input functions.
c
        allocate(coefs(ncoefs,nfuns3))
c
        call legegs(eps,disc,nfuns3,funuser3,npoly,x1,x2,par4,
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
        ngoal = (krank+1)/2+1
        if (nquad1 .gt. ngoal) goto 0100
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
        a      = -1.0d0
        b      =  1.0d0
c
        nints0    = 2
        ab0(1,1)  = a
        ab0(2,1)  = x1-eps
        ab0(1,2)  = x1+eps
        ab0(2,2)  = b
        kdisc     = kdisc+1
        nfuns4    = 3*(npoly+1)
c
        call legedisc(ier,nints0,ab0,kdisc,eps,nfuns4,
     -    funuser4,npoly,x1,x2,y2,disc,ldisc,ncoefs,lkeep)
c
c       Gram-Schmidt the input functions.
c
        allocate(coefs(ncoefs,nfuns4))
c
        call prinf("before legegs, nfuns4 = *",nfuns4,1)
        call legegs(eps,disc,nfuns4,funuser4,npoly,x1,x2,y2,
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
     1    nquad0,xs0,whts0,nquad2,xs2,whts2,rints)
c
c
c       Build the Gaussian quadrature.
c
        ifaccept = 0
        ngoal    = 0
c
        call gaussquad(eps,krank,rints,funeval,disc,coefs,krank,
     1    par4,nquad2,xs2,whts2,a,b,ngoal,ifaccept)
c
        deallocate(coefs,xs0,whts0)
        ngoal = (krank+1)/2+1
        if (nquad2 .gt. ngoal) goto 0200
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
        xs(nquad)   = xs2(i) 
        ys(nquad)   = xs1(idisc) 
        whts(nquad) = whts1(idisc)*whts2(i)
        end do
        end do
c
c
        call prina("*")
        call prina("=========================*")
        call prinf("final nquad = *",nquad,1)
        call prin2("final xs    = *",xs,nquad)
        call prin2("final ys    = *",ys,nquad)
        call prin2("final whts  = *",whts,nquad)
c
c       Test the resulting quadrature rule
c        
        allocate(errs( (npoly+1)*(npoly+2)/2*3 ) )
c
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
        call adaptri(ier,eps,verts,funquad6,x1,x2,i1,j1,
     1    val1,nquad01)
c
        verts(1,1) = -1
        verts(2,1) = -1
        verts(1,2) = -1
        verts(2,2) =  1
        verts(1,3) =  1
        verts(2,3) =  1
c
        call adaptri(ier,eps,verts,funquad6,x1,x2,i1,j1,
     1    val2,nquad02)
c
        val0  = val1+val2
c
        val  = 0
        do i=1,nquad
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
        call funquad6(x,y,x1,x2,i1,j1,val00)
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
        call adaptri(ier,eps,verts,funquad7,x1,x2,i1,j1,
     1    val1,nquad01)
c
        verts(1,1) = -1
        verts(2,1) = -1
        verts(1,2) = -1
        verts(2,2) =  1
        verts(1,3) =  1
        verts(2,3) =  1
c
        call adaptri(ier,eps,verts,funquad7,x1,x2,i1,j1,
     1    val2,nquad02)
c
        val0  = val1+val2
c
        val  = 0
        do i=1,nquad
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
        call funquad7(x,y,x1,x2,i1,j1,val00)
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
        call adaptri(ier,eps,verts,funquad5,x1,x2,i1,j1,
     1    val1,nquad01)
c
        verts(1,1) = -1
        verts(2,1) = -1
        verts(1,2) = -1
        verts(2,2) =  1
        verts(1,3) =  1
        verts(2,3) =  1
c
        call adaptri(ier,eps,verts,funquad5,x1,x2,i1,j1,
     1    val2,nquad02)
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

        call prin2("errs = *",errs,idx)
        call prina("=========================*")
        call prina("*")
c
        stop
c

        end





        subroutine funeval(t,disc,coefs,nfuns,par4,vals,ders,ifders)
        implicit double precision (a-h,o-z)
        dimension vals(1),ders(1)
        call leged_eval(disc,nfuns,coefs,t,vals,ders,ifders)
        end


        subroutine funquad1(y1,val,x1,x2,y2,ideg,par5,par6,par7,
     -    par8)
        implicit double precision (a-h,o-z)
        call lege0(ideg,y1,pol,der)
        dd  = log( (x1-y1)**2 + (x2-y2)**2)
        val = dd*pol
        end subroutine

        subroutine funquad2(y1,val,x1,x2,y2,ideg,par5,par6,par7,
     -    par8)
        implicit double precision (a-h,o-z)
        call lege0(ideg,y1+x1,pol,der)
        dd  = (x1-y1)/((x1-y1)**2 + (x2-y2)**2)
        val = dd*pol
        end subroutine

        subroutine funquad3(y1,val,x1,x2,y2,ideg,par5,par6,par7,
     -    par8)
        implicit double precision (a-h,o-z)
        call lege0(ideg,y1+x1,pol,der)
        dd  = (x2-y2)/((x1-y1)**2 + (x2-y2)**2)
        val = dd*pol
        end subroutine

        subroutine funquad4(y1,y2,x1,x2,i1,j1,val)
        implicit double precision (a-h,o-z)
        dimension pols1(0:100)
        dimension pols2(0:100)
        call lege(i1,y1,pols1)
        call lege(j1,y2,pols2)
        val = log( (x1-y1)**2 + (x2-y2)**2) * pols1(i1)*pols2(j1)
        end subroutine

        subroutine funquad5(y1,y2,x1,x2,i1,j1,val)
        implicit double precision (a-h,o-z)
        dimension pols1(0:100)
        dimension pols2(0:100)
        call lege(i1,y1,pols1)
        call lege(j1,y2,pols2)
        val = pols1(i1)*pols2(j1)
        end subroutine


        subroutine funquad6(y1,y2,x1,x2,i1,j1,val)
        implicit double precision (a-h,o-z)
        dimension pols1(0:100)
        dimension pols2(0:100)
c
        call lege(i1,y1,pols1)
        call lege(j1,y2,pols2)
        val = (x1-y1)/ ((x1-y1)**2 + (x2-y2)**2) * pols1(i1)*pols2(j1)

        end subroutine

        subroutine funquad7(y1,y2,x1,x2,i1,j1,val)
        implicit double precision (a-h,o-z)
        dimension pols1(0:100)
        dimension pols2(0:100)
c
        call lege(i1,y1,pols1)
        call lege(j1,y2,pols2)
c
        val = (x2-y2)/ ((x1-y1)**2 + (x2-y2)**2) * pols1(i1)*pols2(j1)
c
        end subroutine



        subroutine funuser1(y2,vals,npoly,x1,x2,par4)
        implicit double precision (a-h,o-z)
        double precision vals(1),pols(0:100)
        external funquad1,funquad2,funquad3

        a    = -1.0d0
        b    =  1.0d0
        m    = 30
        eps  = epsilon(0.0d0)*100
        ideg = 0
c
        call lege(npoly+1,y2,pols)
c
        idx  = 0 
        do ndeg=0,npoly
        do ideg1=0,ndeg
        ideg2 = ndeg-ideg1
        call adapgauss(ier,a,b,eps,m,val,funquad1,x1,x2,y2,
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
        call lege(npoly+1,y1,pols)
c
        dd = log((x1-y1)**2 + (x2-y2)**2)
        idx   = 0

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





        subroutine funuser3(y2,vals,npoly,x1,x2,par4)
        implicit double precision (a-h,o-z)
        double precision vals(1),pols(0:100)
        external funquad2,funquad3
c
        call lege(npoly+1,y2,pols)
c
        eps  = epsilon(0.0d0)*100

        a1    = -1.0d0
        b1    = x1-eps

        a2    = x1+eps
        b2    = 1.0d0
c
        m    = 30
        idx  = 0 
c
        do ndeg=0,npoly
        do ideg1=0,ndeg
        ideg2 = ndeg-ideg1
        call adapgauss(ier,a1,b1,eps,m,val1,funquad2,x1,x2,y2,
     1    ideg1,par5,par6,par7,par8)
        call adapgauss(ier,a2,b2,eps,m,val2,funquad2,x1,x2,y2,
     1    ideg1,par5,par6,par7,par8)
        val = val1+val2

        idx       = idx + 1
        vals(idx) = val*pols(ideg2)
        end do
        end do
c
        do ndeg=0,npoly
        do ideg1=0,ndeg
        ideg2 = ndeg-ideg1
        call adapgauss(ier,a1,b1,eps,m,val1,funquad3,x1,x2,y2,
     1    ideg1,par5,par6,par7,par8)
        call adapgauss(ier,a2,b2,eps,m,val2,funquad3,x1,x2,y2,
     1    ideg1,par5,par6,par7,par8)
        val = val1+val2

        idx       = idx + 1
        vals(idx) = val*pols(ideg2)
        end do
        end do
c
        do ideg=0,npoly
        idx       = idx + 1
        vals(idx) = pols(ideg)
        end do
c
        return
        end subroutine



        subroutine funuser4(y1,vals,npoly,x1,x2,y2)
        implicit double precision (a-h,o-z)
        dimension ns(1),ys2(1),pols(0:100),vals(1)
        external funquad1
c
c
        call lege(npoly+1,y1,pols)
c
        dd1 = (x1-y1)/((x1-y1)**2 + (x2-y2)**2)
        dd2 = (x2-y2)/((x1-y1)**2 + (x2-y2)**2)
c
        idx   = 0

        do i=0,npoly
        idx       = idx + 1
        vals(idx) = pols(i)*dd1
        end do
c
        do i=0,npoly
        idx       = idx + 1
        vals(idx) = pols(i)*dd2
        end do
c
        do i=0,npoly
        idx       = idx + 1
        vals(idx) = pols(i)
        end do
c
        end subroutine



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






c$$$        subroutine diagquad2(ier,eps,npoly,x1,x2,nquad,xs,ys,whts)
c$$$        implicit double precision (a-h,o-z)
c$$$        dimension                     :: xsgrid(100)
c$$$        double precision, allocatable :: disc(:),coefs(:,:)
c$$$        double precision, allocatable :: xs0(:),whts0(:),rints(:)
c$$$        double precision, allocatable :: xs1(:),whts1(:),ab(:,:)
c$$$        double precision, allocatable :: xs2(:),whts2(:),errs(:)
c$$$        double precision              :: xs(1),ys(1),whts(1)
c$$$        double precision              :: verts(2,3)
c$$$        integer, allocatable          :: idxs(:)
c$$$        external funuser3,funuser4,funeval,funquad6,funquad7
c$$$c
c$$$c       Build a tensor product quadrature on [-1,1] x [-1,1] for integrals of
c$$$c       the form
c$$$c
c$$$c           1     1
c$$$c       \int  \int  (  (x1-y1)/( (x1-y1)^2 + (x2-y2)^2 )  P_i1(y1) P_j1(y2)  +
c$$$c           -1    1
c$$$c
c$$$c                       (x2-y2)/( (x1-y1)^2 + (x2-y2)^2 )  P_i1(y1) P_j1(y2)  +
c$$$c
c$$$c                 P_i2(y1) P_j2(y2) ) dy1 dy2
c$$$c
c$$$c       where (x1,x2) is a  user-specified point in [-1,1] x [-1,1], P_k
c$$$c       denotes the Legendre polynomial of degree k and
c$$$c
c$$$c           0 <= i1+j1 < = npoly,   0 <= i2+j2 <= npoly.
c$$$c
c$$$c       Input parameters:
c$$$c         npoly - the degree of the polynomials 
c$$$c         (x1,x2) - the target node
c$$$c
c$$$c       Output parameters:
c$$$c         ier - an error return code
c$$$c             ier   =  0 indicates successful execution
c$$$c             ier   =  
c$$$c 
c$$$        ier = 0
c$$$
c$$$c
c$$$c       First build a quadrature for integrating the functions F_ij, G_ij and
c$$$c       H_ij defined via
c$$$c
c$$$c                                1
c$$$c         F_ij(y2) = P_iy2) \int   ( log ( (x1-y1)^2 + (x2-y2)^2 ) P_j(y1)
c$$$c                               -1
c$$$c
c$$$c                   + P_i(y1) P_j(y2) ) dy1dy2
c$$$c
c$$$c       where 0<= i1+j1, j2 <= npoly and P_n denotes the Legendre polynomial of degree n.
c$$$c
c$$$        call prin2("in diagquad2, x1 = *",x1,1)
c$$$        call prin2("in diagquad2, x2 = *",x2,1)
c$$$        call prinf("in diagquad2, npoly = *",npoly,1)
c$$$c
c$$$        ldisc = 100 000 000
c$$$        allocate(disc(ldisc))
c$$$        allocate(xs1(1000),whts1(1000),rints(1000),idxs(1000))
c$$$        allocate(xs2(1000),whts2(1000))
c$$$c
c$$$        kdisc = 30
c$$$c
c$$$        nquad=0
c$$$ 0100 continue
c$$$        a      = -1.0d0
c$$$        b      =  1.0d0
c$$$        nfuns1 =  (npoly+1)*(npoly+2) + npoly+1
c$$$c
c$$$        ngrid     = 2
c$$$        xsgrid(1) = a
c$$$        xsgrid(2) = b
c$$$        kdisc     = kdisc+1
c$$$c
c$$$        call prinf("nfuns1 = *",nfuns1,1)
c$$$        call legedisc(ier,ngrid,xsgrid,kdisc,eps,nfuns1,
c$$$     -    funuser3,npoly,x1,x2,par4,disc,ldisc,ncoefs,lkeep)
c$$$        if (ier .ne. 0) then
c$$$        call prinf("after legedisc, ier = *",ier,1)
c$$$        stop
c$$$        endif
c$$$c
c$$$c       Gram-Schmidt the input functions.
c$$$c
c$$$        allocate(coefs(ncoefs,nfuns1))
c$$$c
c$$$        call legegs(eps,disc,nfuns1,funuser3,npoly,x1,x2,par4,
c$$$     1    krank,coefs)
c$$$        call prinf("after legegs, krank=*",krank,1)
c$$$c
c$$$c       Fetch the initial, oversampled, quadrature rule.
c$$$c
c$$$        allocate(xs0(ncoefs),whts0(ncoefs))
c$$$        call legedisc_quad(disc,nquad0,xs0,whts0)
c$$$c
c$$$c       Compute the Chebyshev quadrature.
c$$$c
c$$$        call chebquad(krank,funeval,disc,coefs,krank,par4,
c$$$     1    nquad0,xs0,whts0,nquad1,xs1,whts1,rints)
c$$$c
c$$$c       Build the Gaussian quadrature.
c$$$c
c$$$        ifaccept = 0
c$$$        ngoal    = 0
c$$$c
c$$$        call gaussquad(eps,krank,rints,funeval,disc,coefs,krank,
c$$$     1    par4,nquad1,xs1,whts1,a,b,ngoal,ifaccept)
c$$$c
c$$$        deallocate(coefs,xs0,whts0)
c$$$c
c$$$        ngoal = (krank+1)/2+1
c$$$        if (nquad1 .gt. ngoal) goto 0100
c$$$c
c$$$        do i=1,nquad1
c$$$        idxs(i) = i
c$$$        end do
c$$$        call quicksort(nquad1,xs1,idxs)
c$$$        whts1(1:nquad1)  = whts1(idxs(1:nquad1))
c$$$        call prin2("xs1    = *",xs1,nquad1)
c$$$        call prin2("whts1  = *",whts1,nquad1)
c$$$c
c$$$c
c$$$c       For each point y2 in the preceeding quadrature rule, build
c$$$c       a quadrature for the functions
c$$$c
c$$$c             1
c$$$c        \int    ( log ( (x1-y1)^2 + (x2-y2)^2 ) P_i(y1)  ) dy1 
c$$$c            -1
c$$$c
c$$$c       where 0<= i <= npoly.
c$$$c
c$$$        do idisc=1,nquad1
c$$$        y2    = xs1(idisc)
c$$$        kdisc = 30
c$$$ 0200 continue
c$$$        a      = -1.0d0
c$$$        b      =  1.0d0
c$$$c
c$$$        ngrid     = 2
c$$$        xsgrid(1) = a
c$$$        xsgrid(2) = b
c$$$        kdisc     = kdisc+1
c$$$        nfuns2    = 2*(npoly+1) + (npoly+1)
c$$$c
c$$$        call legedisc(ier,ngrid,xsgrid,kdisc,eps,nfuns2,
c$$$     -    funuser4,npoly,x1,x2,y2,disc,ldisc,ncoefs,lkeep)
c$$$c
c$$$c       Gram-Schmidt the input functions.
c$$$c
c$$$        allocate(coefs(ncoefs,nfuns2))
c$$$c
c$$$        call prinf("before legegs, nfuns2 = *",nfuns2,1)
c$$$        call legegs(eps,disc,nfuns2,funuser4,npoly,x1,x2,y2,
c$$$     1    krank,coefs)
c$$$        call prinf("after legegs, krank=*",krank,1)
c$$$c
c$$$c       Fetch the initial, oversampled, quadrature rule.
c$$$c
c$$$        allocate(xs0(ncoefs),whts0(ncoefs))
c$$$        call legedisc_quad(disc,nquad0,xs0,whts0)
c$$$c
c$$$c       Compute the Chebyshev quadrature.
c$$$c
c$$$        call chebquad(krank,funeval,disc,coefs,krank,par4,
c$$$     1    nquad0,xs0,whts0,nquad2,xs2,whts2,rints)
c$$$c
c$$$c       Build the Gaussian quadrature.
c$$$c
c$$$        ifaccept = 0
c$$$        ngoal    = 0
c$$$c
c$$$        call gaussquad(eps,krank,rints,funeval,disc,coefs,krank,
c$$$     1    par4,nquad2,xs2,whts2,a,b,ngoal,ifaccept)
c$$$c
c$$$        deallocate(coefs,xs0,whts0)
c$$$        ngoal = (krank+1)/2+1
c$$$        if (nquad2 .gt. ngoal) goto 0200
c$$$        do i=1,nquad2
c$$$        idxs(i) = i
c$$$        end do
c$$$        call quicksort(nquad2,xs2,idxs)
c$$$        whts2(1:nquad2)  = whts2(idxs(1:nquad2))
c$$$
c$$$        call prin2("xs2    = *",xs2,nquad2)
c$$$        call prin2("whts2  = *",whts2,nquad2)
c$$$c
c$$$        do i=1,nquad2
c$$$        nquad       = nquad+1
c$$$        xs(nquad)   = xs2(i)
c$$$        ys(nquad)   = xs1(idisc)
c$$$        whts(nquad) = whts1(idisc)*whts2(i)
c$$$        end do
c$$$c
c$$$        end do
c$$$c
c$$$        call prina("*")
c$$$        call prina("=========================*")
c$$$        call prinf("final nquad = *",nquad,1)
c$$$        call prin2("final xs    = *",xs,nquad)
c$$$        call prin2("final ys    = *",ys,nquad)
c$$$        call prin2("final whts  = *",whts,nquad)
c$$$c
c$$$        stop
c$$$c$$$c
c$$$c$$$c
c$$$c$$$c       Test the resulting quadrature rule
c$$$c$$$c
c$$$c$$$        
c$$$c$$$        allocate(errs( (npoly+1)*(npoly+2) ) )
c$$$c$$$c
c$$$c$$$        idx = 0
c$$$c$$$
c$$$c$$$        do nn1=0,npoly
c$$$c$$$        do i1=0,nn1
c$$$c$$$        j1 = nn1-i1
c$$$c$$$c
c$$$c$$$        verts(1,1) = -1
c$$$c$$$        verts(2,1) = -1
c$$$c$$$        verts(1,2) =  1
c$$$c$$$        verts(2,2) = -1
c$$$c$$$        verts(1,3) =  1
c$$$c$$$        verts(2,3) =  1
c$$$c$$$c
c$$$c$$$        call adaptri(ier,eps,verts,funquad4,x1,x2,i1,j1,
c$$$c$$$     1    val1,nquad01)
c$$$c$$$c
c$$$c$$$        verts(1,1) = -1
c$$$c$$$        verts(2,1) = -1
c$$$c$$$        verts(1,2) = -1
c$$$c$$$        verts(2,2) =  1
c$$$c$$$        verts(1,3) =  1
c$$$c$$$        verts(2,3) =  1
c$$$c$$$c
c$$$c$$$        call adaptri(ier,eps,verts,funquad4,x1,x2,i1,j1,
c$$$c$$$     1    val2,nquad02)
c$$$c$$$c
c$$$c$$$        val0  = val1+val2
c$$$c$$$c
c$$$c$$$        val  = 0
c$$$c$$$        do i=1,nquad
c$$$c$$$        x   = xs(i)
c$$$c$$$        y   = ys(i)
c$$$c$$$        wht = whts(i)
c$$$c$$$        call funquad4(x,y,x1,x2,i1,j1,val00)
c$$$c$$$        val = val + wht*val00
c$$$c$$$        end do
c$$$c$$$c
c$$$c$$$        idx = idx + 1
c$$$c$$$        errs(idx) = abs(val-val0)
c$$$c$$$c
c$$$c$$$        end do
c$$$c$$$        end do
c$$$c$$$c
c$$$c$$$        do nn1=0,npoly
c$$$c$$$        do i1=0,nn1
c$$$c$$$        j1 = nn1-i1
c$$$c$$$c
c$$$c$$$        verts(1,1) = -1
c$$$c$$$        verts(2,1) = -1
c$$$c$$$        verts(1,2) =  1
c$$$c$$$        verts(2,2) = -1
c$$$c$$$        verts(1,3) =  1
c$$$c$$$        verts(2,3) =  1
c$$$c$$$c
c$$$c$$$        call adaptri(ier,eps,verts,funquad5,x1,x2,i1,j1,
c$$$c$$$     1    val1,nquad01)
c$$$c$$$c
c$$$c$$$        verts(1,1) = -1
c$$$c$$$        verts(2,1) = -1
c$$$c$$$        verts(1,2) = -1
c$$$c$$$        verts(2,2) =  1
c$$$c$$$        verts(1,3) =  1
c$$$c$$$        verts(2,3) =  1
c$$$c$$$c
c$$$c$$$        call adaptri(ier,eps,verts,funquad5,x1,x2,i1,j1,
c$$$c$$$     1    val2,nquad02)
c$$$c$$$c
c$$$c$$$        val0  = val1+val2
c$$$c$$$c
c$$$c$$$        val  = 0
c$$$c$$$        do i=1,nquad
c$$$c$$$        x   = xs(i)
c$$$c$$$        y   = ys(i)
c$$$c$$$        wht = whts(i)
c$$$c$$$        call funquad5(x,y,x1,x2,i1,j1,val00)
c$$$c$$$        val = val + wht*val00
c$$$c$$$        end do
c$$$c$$$c
c$$$c$$$        idx = idx + 1
c$$$c$$$        errs(idx) = abs(val-val0)
c$$$c$$$c
c$$$c$$$        end do
c$$$c$$$        end do
c$$$c$$$c
c$$$c$$$        call prin2("errs = *",errs,idx)
c$$$c$$$        call prina("=========================*")
c$$$c$$$        call prina("*")
c$$$c$$$
c$$$c$$$c
c$$$
c$$$        end
c$$$
c$$$
c$$$        subroutine funquad2(y1,val,x1,x2,y2,ideg,par5,par6,par7,
c$$$     -    par8)
c$$$        implicit double precision (a-h,o-z)
c$$$        call lege0(ideg,y1,pol,der)
c$$$        dd  = (x1-y1)/((x1-y1)**2 + (x2-y2)**2)
c$$$        val = dd*pol
c$$$        end subroutine
c$$$
c$$$
c$$$        subroutine funquad3(y1,val,x1,x2,y2,ideg,par5,par6,par7,
c$$$     -    par8)
c$$$        implicit double precision (a-h,o-z)
c$$$        call lege0(ideg,y1,pol,der)
c$$$        dd  = (x2-y2)/((x1-y1)**2 + (x2-y2)**2)
c$$$        val = dd*pol
c$$$        end subroutine
c$$$
c$$$
c$$$        subroutine funuser3(y2,vals,npoly,x1,x2,par4)
c$$$        implicit double precision (a-h,o-z)
c$$$        double precision vals(1),pols(0:100)
c$$$        external funquad2,funquad3
c$$$        a    = -1.0d0
c$$$        b    =  1.0d0
c$$$        m    = 30
c$$$        eps  = epsilon(0.0d0)*100
c$$$        ideg = 0
c$$$c
c$$$        call lege(npoly+1,y2,pols)
c$$$        idx = 0
c$$$c
c$$$        do ndeg=0,npoly
c$$$        do ideg1=0,ndeg
c$$$        ideg2 = ndeg-ideg1
c$$$        call adapgauss(ier,a,b,eps,m,val,funquad2,x1,x2,y2,
c$$$     1    ideg1,par5,par6,par7,par8)
c$$$        idx       = idx + 1
c$$$        vals(idx) = val*pols(ideg2)
c$$$        end do
c$$$        end do
c$$$c
c$$$        do ndeg=0,npoly
c$$$        do ideg1=0,ndeg
c$$$        ideg2 = ndeg-ideg1
c$$$        call adapgauss(ier,a,b,eps,m,val,funquad3,x1,x2,y2,
c$$$     1    ideg1,par5,par6,par7,par8)
c$$$        idx       = idx + 1
c$$$        vals(idx) = val*pols(ideg2)
c$$$        end do
c$$$        end do
c$$$c
c$$$        do ideg=0,npoly
c$$$        idx = idx + 1
c$$$        vals(idx) = pols(ideg)
c$$$        end do
c$$$c
c$$$        return
c$$$        end subroutine
c$$$
c$$$
c$$$
c$$$        subroutine funuser4(y1,vals,npoly,x1,x2,y2)
c$$$        implicit double precision (a-h,o-z)
c$$$        dimension ns(1),ys2(1),pols(0:100),vals(1)
c$$$        external funquad1
c$$$        call lege(npoly+1,y1,pols)
c$$$        dd1 = (x1-y1)/((x1-y1)**2 + (x2-y2)**2)
c$$$        dd2 = (x2-y2)/((x1-y1)**2 + (x2-y2)**2)
c$$$c
c$$$        idx   = 0
c$$$c
c$$$        do i=0,npoly
c$$$        idx       = idx + 1
c$$$        vals(idx) = pols(i)*dd1
c$$$        end do
c$$$c
c$$$        do i=0,npoly
c$$$        idx       = idx + 1
c$$$        vals(idx) = pols(i)*dd2
c$$$        end do
c$$$c
c$$$        do i=0,npoly
c$$$        idx       = idx + 1
c$$$        vals(idx) = pols(i)
c$$$        end do
c$$$c
c$$$        end subroutine
c$$$
c$$$
c$$$
c$$$