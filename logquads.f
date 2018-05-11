!
!       Construct a quadrature for evaluating
!
!             1
!       \int     log|x-x0| * p(x) + q(x)
!            -1
!
!       with p and q polynomials of degree less than or equal to
!       norder.
!
        implicit double precision (a-h,o-z)
        dimension                     :: xsgrid(100)
        double precision, allocatable :: disc(:),coefs(:,:)
        double precision, allocatable :: xs0(:),whts0(:),rints(:)
        double precision, allocatable :: xs(:),whts(:)

        external funuser,funeval
c
        eps    = 1.0d-30
        x0     = 0.1d0
        norder = 16
c
        call funuser_init(x0,norder,nfuns)
        call prinf("after funuser_init, nfuns = *",nfuns,1)
c
c       Discretize the functions
c
        ldisc = 100 000 000
        allocate(disc(ldisc))
        allocate(xs(1000),whts(1000),rints(1000))
c
c
        a = -1.0d0   - x0
        b =  1.0d0   - x0
c
c       specify an initial discretization grid
c
        ngrid     = 2
        xsgrid(1) = a
        xsgrid(2) = b
        kdisc     = 30
c
        call legedisc(ier,ngrid,xsgrid,kdisc,eps,nfuns,
     -    funuser,par1,par2,par3,par4,disc,ldisc,ncoefs,lkeep)
c
        if (ier .ne. 0) then
        call prinf("after legedisc, ier = *",ier,1)
        stop
        endif
c
        call prinf("after legedisc, ncoefs = *",ncoefs,1)
c
c       Gram-Schmidt the input functions.
c
        allocate(coefs(ncoefs,nfuns))
c
        call prinf("before legegs, nfuns = *",nfuns,1)
        call legegs(eps,disc,nfuns,funuser,par1,par2,par3,par4,
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
c        if (ifgauss .eq. 1) then
        ifaccept = 0
        ngoal    = 0
c
        call gaussquad(eps,krank,rints,funeval,disc,coefs,krank,
     1    par4,nquad,xs,whts,a,b,ngoal,ifaccept)
c
c        endif
c
        xs = xs + x0
        call prin2("xs = *",xs,nquad)
        call prin2("whts = *",whts,nquad)

c
        end


        subroutine funeval(t,disc,coefs,nfuns,par4,vals,ders,ifders)
        implicit double precision (a-h,o-z)
        dimension vals(1),ders(1)
        call leged_eval(disc,nfuns,coefs,t,vals,ders,ifders)
        end


        subroutine funuser(u,vals,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension vals(1),pols(0:100)
        save
c
        x = u + xx
        call lege(norder,x,pols)
c
        ifun = 0
        dd   = log(abs(u))
c
        do i=0,norder
        ifun                = ifun+1
        vals(ifun)          = pols(i)
        end do
c
        do i=0,norder
        ifun                = ifun+1
        vals(ifun)          = pols(i)*dd
        end do
c
        return
c
        entry funuser_init(x0,norder0,nfuns0)
        xx      = x0
        norder = norder0
        nfuns0 = 2*(norder+1)
c
        end



