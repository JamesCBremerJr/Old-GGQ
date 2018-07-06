        implicit double precision (a-h,o-z)
        dimension xslege(1000),whtslege(100)
        dimension xs(10000),whts(10000)
        integer, allocatable          :: nquadall(:)
        double precision, allocatable :: xsall(:,:),whtsall(:,:)
        double precision, allocatable :: xsnear(:),whtsnear(:)
c
        eps   = 1.0d-22
        npoly = 8
c
        call legequad(npoly,xslege,whtslege)
        call prin2("xslege = *",xslege,npoly)
        call prin2("whtslege = *",whtslege,npoly)
c
c       Construct the near quadrature rule
c
        allocate(xsnear(1000),whtsnear(1000))
        call nearquad(eps,npoly,nquadnear,xsnear,whtsnear)
c
c       Construct the singular quadrature rules
c
        allocate(xsall(300,npoly),whtsall(300,npoly))
        allocate(nquadall(npoly))

        do i=1,npoly
        x0 = xslege(i)
        call diagquad(eps,npoly,x0,nquad,xs,whts)
        nquadall(i) = nquad
        do l=1,nquad
        xsall(l,i)   = xs(l)
        whtsall(l,i) = whts(l)
        end do
        end do        
c
        call prinf("nquadall = *",nquadall,npoly)
c
c       Write the quadrature out to the disk
c
 1000 format (A)
 1100 format ("xs(",I2.2,")   = ",D44.36)
 1200 format ("whts(",I2.2,") = ",D44.36)
 1300 format ("nquad    = ",I3.3)
c
 1400 format ("nquadsall(",I2.2,")    = ",I3.3)
 1500 format ("xsall(",I2.2,",",I2.2,")     = ",D44.36)
 1600 format ("whtsall(",I2.2,",",I2.2,")   = ",D44.36)

        
        iw = 1001
        open(iw,FILE='quads1d.f90')
c
        write(iw,1000) "subroutine nearquad(nquad,xs,whts)"
        write(iw,1000) "implicit double precision (a-h,o-z)"
        write(iw,1000) "dimension xs(1),whts(1)"
        write(iw,1300)  nquadnear
        do i=1,nquadnear
        write(iw,1100)  i,xsnear(i)
        write(iw,1200)  i,whtsnear(i)
        end do
c
        write(iw,1000) "end subroutine"
        write(iw,*)    ""
        write(iw,*)    ""

        write(iw,1000) "subroutine singquad(nquadsall,xsall,whtsall)"
        write(iw,1000) "implicit double precision (a-h,o-z)"
        write(iw,1000) "dimension xsall(200,1),whtsall(200,1)"
        write(iw,1000) "dimension nquadsall(1)"

        do i=1,npoly
        write (iw,1400) i,nquadall(i)

        do j=1,nquadall(i)
        write (iw,1500) j,i,xsall(j,i)
        write (iw,1600) j,i,whtsall(j,i)
        end do

        end do
c
        write(iw,1000) "end subroutine"

        close(iw)

        end program




        subroutine diagquad(eps,npoly,x0,nquad,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),whts(1)
c
c       Construct a quadrature rule for integrals of the form
c
c              -1
c          \int    (  log|x0-x| p(x) + q(x) )  dx
c               1
c
c       where p and q are polynomials of degree less than or equal
c       to npoly and x0 is a specified point in the interval [-1,1].
c
c
        double precision, allocatable :: disc(:),coefs(:,:),ab0(:,:)
        double precision, allocatable :: xs0(:),whts0(:),rints(:)
        double precision, allocatable :: errs(:)
        integer                       :: idxs(1000)
        external funuser,funeval,funquad
c
        kdisc  = 31
        nfuns  = 2*(npoly+1)
c
c       Discretize the functions
c
        ldisc = 100 000 000
        allocate(disc(ldisc),rints(10000),ab0(2,100))
c
        nints0    = 2

        a         = -1-x0
        b         =  1-x0
        ab0(1,1)  = a
        ab0(2,1)  = 0.0d0
        ab0(1,2)  = 0.0d0
        ab0(2,2)  = b
        ifadap    = 1
c
        call legedisc(ier,ifadap,nints0,ab0,kdisc,eps,nfuns,
     -    funuser,npoly,x0,par3,par4,disc,ldisc,ncoefs,lkeep)
        call prinf("after legedisc, ncoefs = *",ncoefs,1)
c
c       Gram-Schmidt the input functions.
c
        allocate(coefs(ncoefs,nfuns))
        call legegs(eps,disc,nfuns,funuser,npoly,x0,par3,par4,
     1    krank,coefs)
        call prinf("after legegs, krank=*",krank,1)
c
c       Fetch the initial, extremely oversampled, quadrature rule.
c
        allocate(xs0(ncoefs),whts0(ncoefs))
        call legedisc_quad(disc,nquad0,xs0,whts0)
c
c       Compute the Chebyshev quadrature.
c
        call chebquad(krank,funeval,disc,coefs,krank,par4,
     1    nquad0,xs0,whts0,nquad,xs,whts,rints)
c
c       Build the Gaussian quadrature.
c
        ifaccept = 0
        ngoal    = 0
c
        call gaussquad(eps,krank,rints,funeval,disc,coefs,krank,
     1    par4,nquad,xs,whts,a,b,ngoal,ifaccept)
c
        do i=1,nquad
        xs(i) = xs(i) +x0
        end do
c
        do i=1,nquad
        idxs(i) = i
        end do
        call quicksort(nquad,xs,idxs)
        whts(1:nquad) = whts(idxs(1:nquad))

        call prin2("xs   = *",xs,nquad)
        call prin2("whts = *",whts,nquad)
c
c       Test the quadrature
c
        a = -1.0d0
        b =  1.0d0
        allocate(errs((npoly+1)*(npoly+1)))

        idx = 0
        do i=0,npoly
        do j=0,npoly
c
        m = 30
        call adapgauss(ier,a,b,eps,m,val0,funquad,i,j,x0,
     1      par4,par5,par6,par7,par8)

        val = 0

        do k=1,nquad
        x   = xs(k)
        wht = whts(k)
        call funquad(x,val00,i,j,x0,par4,par5,par6,par7,par8)
        val = val + val00*wht
        end do
c
        idx = idx+1
        errs(idx) = abs(val-val0)

        end do
        end do
c

        call prin2("errs = *",errs,(npoly+1)*(npoly+1))
c
        end subroutine


        subroutine funuser(x,vals,npoly,x0,par3,par4)
        implicit double precision (a-h,o-z)
        dimension vals(1),pols(0:100)
c
        call lege(npoly+1,x+x0,pols)
        dd     = log(abs(x))
        ifun  = 0
c
        do i=0,npoly
        ifun                = ifun+1
        vals(ifun)          = pols(i)
        end do
c
        do i=0,npoly
        ifun                = ifun+1
        vals(ifun)          = pols(i)*dd
        end do
c
        return
        end

        subroutine funeval(t,disc,coefs,nfuns,par4,vals,ders,ifders)
        implicit double precision (a-h,o-z)
        dimension vals(1),ders(1)
        call leged_eval(disc,nfuns,coefs,t,vals,ders,ifders)
        end



        subroutine funquad(x,val,i,j,x0,par4,par5,par6,par7,
     1    par8)
        implicit double precision (a-h,o-z)
        call lege0(i,x,pol1,der)
        call lege0(j,x,pol2,der)
        dd = log(abs(x-x0))
        val = dd*pol1 + pol2
        end


        subroutine nearquad(eps,npoly,nquad,xs,whts)
        implicit double precision (a-h,o-z)
        double precision xs(1),whts(1),zs(10000)
        double precision xslege(1000),whtslege(1000),ab(2,100)
        double precision, allocatable :: disc(:),coefs(:,:),ab0(:,:)
        double precision, allocatable :: xs0(:),whts0(:),rints(:)
        double precision, allocatable :: errs(:)
        integer                       :: idxs(1000)
        external funuser2,funquad2,funeval
c
c       Set the zs
c
        nints   = 5
        k       = 24
        call legequad(k,xslege,whtslege)
        nzs     = 0

c$$$        ab(1,1) = 1.0d0 + 1.0d-7
c$$$        ab(2,1) = 1.0d0 + 1.0d-5

        ab(1,1) = 1.0d0 + 1.0d-5
        ab(2,1) = 1.0d0 + 1.0d-3
c
        ab(1,2) = 1.0d0 + 1.0d-3
        ab(2,2) = 1.0d0 + 1.0d-1
c
        ab(1,3) = 1.0d0 + 1.0d-1
        ab(2,3) = 1.5d0 

        ab(1,4) = 1.5d0
        ab(2,4) = 2.0d0 

        ab(1,5) = 1.5d0
        ab(2,5) = 2.0d0 

        do int=1,nints
        a = ab(1,int)
        b = ab(2,int)
        do i=1,k
        nzs     = nzs+1
        zs(nzs) = (b-a)/2 * xslege(i) + (b+a)/2
        nzs     = nzs+1
        zs(nzs) = -zs(nzs-1)
        end do
        end do
c
        call prin2("zs = *",zs,nzs)
c
        kdisc  = 31
        nfuns  = nzs*(npoly+1) + npoly+1
c
c       Discretize the functions
c
        ldisc = 100 000 000
        allocate(disc(ldisc),rints(10000),ab0(2,100))
c
        nints0    = 1
        a         = -1
        b         =  1
        ab0(1,1)  = a
        ab0(2,1)  = b
        ifadap    = 1
c
        call legedisc(ier,ifadap,nints0,ab0,kdisc,eps,nfuns,
     -    funuser2,npoly,nzs,zs,par4,disc,ldisc,ncoefs,lkeep)
        call prinf("after legedisc, ncoefs = *",ncoefs,1)
c
        allocate(coefs(ncoefs,nfuns))
        call legegs(eps,disc,nfuns,funuser2,npoly,nzs,zs,par4,
     1    krank,coefs)
        call prinf("after legegs, krank=*",krank,1)
c
c       Fetch the initial, extremely oversampled, quadrature rule.
c
        allocate(xs0(ncoefs),whts0(ncoefs))
        call legedisc_quad(disc,nquad0,xs0,whts0)
c
c       Compute the Chebyshev quadrature.
c
        call chebquad(krank,funeval,disc,coefs,krank,par4,
     1    nquad0,xs0,whts0,nquad,xs,whts,rints)
c
c       Build the Gaussian quadrature.
c
        ifaccept = 0
        ngoal    = 0
c
        call gaussquad(eps,krank,rints,funeval,disc,coefs,krank,
     1    par4,nquad,xs,whts,a,b,ngoal,ifaccept)
c
        do i=1,nquad
        idxs(i) = i
        end do
        call quicksort(nquad,xs,idxs)
        whts(1:nquad) = whts(idxs(1:nquad))

        call prin2("xsnear   = *",xs,nquad)
        call prin2("whtsnear = *",whts,nquad)
c
c       Test the quadrature
c
        a = -1.0d0
        b =  1.0d0

        nn = nzs
        allocate(errs(nn))

        idx = 0
        do l=1,nzs
        x0 = zs(l)
        errmax = 0.0d0
        do i=0,npoly
        do j=0,npoly
c
          
        m = 30
        call adapgauss(ier,a,b,eps,m,val0,funquad2,i,j,x0,
     1      par4,par5,par6,par7,par8)

        val = 0

        do k=1,nquad
        x   = xs(k)
        wht = whts(k)
        call funquad2(x,val00,i,j,x0,par4,par5,par6,par7,par8)
        val = val + val00*wht
        end do
c
        errmax = max(errmax,abs(val-val0))

        end do
        end do

        idx = idx+1
        errs(idx) = abs(val-val0)

        end do
c

        call prin2("errs = *",errs,nn)

        end subroutine


        subroutine funuser2(x,vals,npoly,nzs,zs,par4)
        implicit double precision (a-h,o-z)
        double precision :: zs(nzs),vals(1),pols(0:100)
c
        
        call lege(npoly+1,x,pols)
        idx = 0
        do j=1,nzs
        x0 = zs(j)
        dd = log(abs(x-x0))
        do i=0,npoly
        idx = idx + 1
        vals(idx) = pols(i)*dd
        end do
        end do

        do i=0,npoly
        idx       = idx + 1
        vals(idx) = pols(i)
        end do
c
        end subroutine


        subroutine funquad2(x,val,i,j,x0,par4,par5,par6,par7,
     1    par8)
        implicit double precision (a-h,o-z)
        call lege0(i,x,pol1,der)
        call lege0(j,x,pol2,der)
        dd = log(abs(x-x0))
        val = dd*pol1 + pol2
        end
