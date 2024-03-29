        implicit double precision (a-h,o-z)
        dimension disc(1 000 000),vals(1000),vals0(1000)
        dimension ders(1000),ders0(1000)
c
        dimension xs0(10000),whts0(10000),rints(10000)
        dimension xs(10000),whts(10000)
        dimension ab0(2,100)
c
        
        double precision, allocatable :: coefs(:,:)
        external funuser,funeval
c
        call prini(6,13)
c
        ldisc = 1 000 000
c
        epsdisc = 1.0d-13
        epsgs   = 1.0d-13
        epsnewt = 1.0d-13
c
        a = -1.0d0
        b =  1.0d0
c
        k     = 30
        nfuns = 120
c
        ifadap = 1
        nints0    = 1
        ab0(1,1)  = a
        ab0(2,1)  = b

        call legedisc(ier,ifadap,nints0,ab0,k,epsdisc,nfuns,funuser,
     -    nfuns,par2,par3,par4,disc,ldisc,ncoefs,lkeep)
c
        call prinf("after legedisc, ier = *",ier,1)
        call prinf("after legedisc, ncoefs = *",ncoefs,1)
c
        allocate( coefs(ncoefs,nfuns) )
c
        if (ier .ne. 0) stop
c
c       Test evaluation.
c
        x=1.0d0+1.0d-3
        call legecoefs(disc,nfuns,funuser,nfuns,par2,par3,par4,
     1    coefs)
        call prin2("coefs = *",coefs,ncoefs)
c     
        ifders = 1
        call legeeval(disc,nfuns,coefs,x,vals,ders,ifders)
        call funuser(x,vals0,ders0,nfuns,par2,par3,par4)
c
        errmax=0
        dermax=0
        do 1000 j=1,nfuns
        d = abs(vals(j)-vals0(j))
        dd = abs(ders(j)-ders0(j))
c
        errmax= max(d,errmax)
        dermax=max(dd,dermax)
 1000 continue
c
        call prin2("evaluation errmax = *",errmax,1)
        call prin2("evaluation dermax = *",dermax,1)
c
        call legegs(epsgs,disc,nfuns,funuser,nfuns,par2,par3,
     1    par4,krank,coefs)

        end
c
c
c
        subroutine funuser(x,vals,nfuns,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension vals(1),ders(1)
        do 1000 j=1,nfuns
        vals(j)=x**(j)
 1000 continue
        end
c
c
c
        subroutine funuser2(x,vals,ders,nfuns,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension vals(1),ders(1)
        do 1000 j=1,nfuns
        vals(j)=x**(j)
        ders(j)=j*x**(j-1)
 1000 continue
        end
c
c
c
        subroutine funeval(x,disc,nfuns,coefs,par4,vals,ders,ifders)
        implicit double precision (a-h,o-z)
        dimension vals(1),ders(1),coefs(1)
        call legeeval(disc,nfuns,coefs,x,vals,ders,ifders)
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and begining of the code
c       proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains code for discretizing collections of 
c       piecewise differentiable user-supplied functions via 
c       piecewise Legendre quadrature formulae.  The word ``discretize'' 
c       is used here to mean an embedding of the functions into a 
c       Eucliden space which *preserves inner products*.  
c
c       The following subroutines are user-callable:
c
c   legedisc - this subroutine constructs a piecewise Legendre 
c       quadrature scheme which discretizes of a collection of user-
c       supplied input functions and their derivatves over an interval 
c       [a,b]
c
c   legequad - this subroutine returns the quadrature stored in 
c       a "disc" structure generated by legedisc

c   legegs - this subroutine performs the pivoted Gram-Schmidt with 
c       reorthogonalization algorithm on a collection of input 
c       functions supplied via an external subroutine and for which 
c       a piecewise Legendre discretization scheme is available.  
c       It returns coefficient expansions for the resulting orthonormal 
c       basis functions *and* their derivatives.
c
c   legecoefs - this subroutine computes coefficient expansions for
c       a collection of user-supplied functions given an existing
c       discretization scheme for the input functions
c
c   legeeval  - this subroutine evaluates a collection of functions
c       represented as coefficient expansions at a user-specified
c       point; it can also (optionally) evaluate the derivatives of 
c       the coefficient expansions
c
c   legedisc_ab - return the list of discretization subintervals
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





c
        subroutine legedisc(ier,ifadap,nints0,ab0,k,eps,nfuns,
     -    funuser,par1,par2,par3,par4,disc,ldisc,ncoefs,lkeep)
        implicit double precision (a-h,o-z)
        dimension disc(1),xs0(1),ab0(2,1)
        external funuser
c
c       Construct a piecewise Legendre scheme discretizing a collection
c       of user-supplied input functions given over an interval.  
c 
c       A structure describing the resulting scheme will be returned in
c       the user-supplied array disc.
c
c                          Input Parameters:
c
c   ifadap - a flag indiciating whether or not to conduct adaptive
c       discretization of a collection of functions in order
c       to refine the initial user-specified discretization scheme
c   (nints0,ab0) - an initial discretization scheme
c   k - the number of nodes used in the piecewise expansions on each
c       interval
c   eps - precision for the discretizations
c   nfuns - the number of input functions
c   funuser - an external subroutine supplying the values of the
c       input functions at a point
c
c       funuser(x,vals,par1,par2,par3,par4)
c
c       Return the values of the user-supplied input functions
c       at the point x in the output array vals.
c
c   ldisc - the length of the user-supplied disc array
c
c                         Output Parameters:
c
c   ier - an error return code;
c       ier = 0    indicates successful execution
c       ier = 4    means that the user-supplied disc array is of
c                  insufficient length
c       ier = 8    the internal stack used during discretization 
c                  overflowed before the desired precision could be
c                  achieved
c       ier = 16   means that the maximum number of discretization
c                  intervals was exceeded before the desired precision
c                  could be achieved
c
c   disc - upon return, this user-supplied array will contain a
c       structure describing the discretization
c
        ier = 0
        ncoefs = 0
c
        maxints = 10 000
c      
c       Allocate memory from the disc array for the quadrature and
c       interval list.
c
        ixs = 1000
        lxs = k
c     
        iwhts = ixs+lxs
        lwhts = k
c
        iu = iwhts+lwhts
        lu = k**2
c
        iv = iu+lu
        lv = k**2
c
        iab = iv+lv
        lab = maxints*2
c
        iw = iab+lab
        lw = ldisc-iw
c
        if (lw .le. 0) then
        ier = 4
        return
        endif
c
        lkeep = iw
c
c       Fetch the k-point Chebyshev quadrature data.
c
        call legendre(k,disc(ixs),disc(iwhts),disc(iu),disc(iv))
c
c       Construct the subintervals.
c
        if (ifadap .eq. 1) then
        call discfuns(ier,eps,nints0,ab0,nfuns,funuser,par1,par2,par3,
     -     par4,maxints,disc(iab),nints,k,disc(ixs),disc(iwhts),
     -     disc(iu))
        if (ier .ne. 0) return

        else
        nints = nints0
        do i=1,nints
        disc(iab+2*(i-1))   = ab0(1,i)
        disc(iab+2*(i-1)+1) = ab0(2,i)
        end do
        endif
c
c
        call prin2("in legedisc, ab = *",disc(iab),2*nints)
c
        ncoefs=nints*k
c
c       Build the structure header.
c
        disc(1) = k
        disc(2) = nints
        disc(3) = a
        disc(4) = b
c
        disc(10) = k
        disc(11) = ixs
        disc(12) = iwhts
        disc(13) = iu
c
        disc(20) = nints
        disc(21) = iab
        disc(22) = 1
c
        end
c
c
c
        subroutine discfuns(ier,eps,nints0,ab0,nfuns,
     -     funuser,par1,par2,par3,par4,maxints,ab,nints,k,xs,whts,u)
        implicit double precision (a-h,o-z)
        dimension ab(2,maxints),xs(1),whts(1),u(k,k),xs0(1)
        dimension ab0(2,nints0)
        double precision, allocatable :: vals(:,:),vals0(:),stack(:,:)
        external funuser
c
c       Discretize a collection of input functions given over the
c
        ier    = 0
        nints = 0
c
        maxstack = 20 000
c
        allocate(vals(k,nfuns),vals0(nfuns+100),stack(2,maxstack))
c
c       Initialize the stack.
c     
        istack=0
        do 0100 i=1,nints0
        a = ab0(1,i)
        b = ab0(2,i)
c
        istack=istack+1
        stack(1,istack)=a
        stack(2,istack)=b
 0100 continue
c
c
 1000 continue
c
c       Pop an entry off the stack.
c
        if (istack .eq. 0) goto 3000
c
        aa = stack(1,istack)
        bb = stack(2,istack)
        istack=istack-1
c
        alpha = (bb-aa)/2
        beta  = (bb+aa)/2
c
c       Evaluate the input functions.
c
        do 1100 i=1,k
        x=alpha*xs(i)+beta
        wht = whts(i)*alpha
c
        call funuser(x,vals0,par1,par2,par3,par4)
        do 1200 j=1,nfuns
        vals(i,j)=vals0(j)*sqrt(wht)
 1200 continue
 1100 continue
c
c       Check the trailing coefficients.
c
        do 1300 j=1,nfuns
c
        i = 1
        dd=0
        do 1350 l=1,k
        dd=dd+u(i,l)*vals(l,j)
 1350 continue
        dd = abs(dd)
        dd = max(1.0d0,dd)
c
c        do 1400 i=k/2+1,k
        do 1400 i=k-8,k
        sum=0
        do 1500 l=1,k
        sum=sum+u(i,l)*vals(l,j)
 1500 continue
c
         if( abs(sum) .gt. eps*dd) then
         goto 2000
         endif
c
 1400 continue
 1300 continue
c
c       Accept the interval.
c
        if (nints .eq. maxints) then
        ier=16
        return
        endif
c
        nints=nints+1
        ab(1,nints)=aa
        ab(2,nints)=bb
c
        goto 1000
c
c       Split the interval.
c
 2000 continue
c
        cc=(aa+bb)/2
c
        if (istack+2 .gt. maxstack) then
        ier=8
        return
        endif
c
        istack=istack+1
        stack(1,istack)=aa
        stack(2,istack)=cc
c
        istack=istack+1
        stack(1,istack)=cc
        stack(2,istack)=bb
c
        goto 1000
c
 3000 continue
c
c       Sort the resulting interval list.
c
        call insort0(2*nints,ab)
        end




c
c
c
        subroutine legegs(eps,disc,nfuns,funuser,par1,par2,par3,
     1    par4,krank,coefs)
        implicit double precision (a-h,o-z)
        dimension disc(1),coefs(1)
c
c       Gram-Schmidt a collection of user-supplied input functions for
c       which a discretization has been produced.  This subroutine
c       returns coefficient expansions representing the resulting
c       orthonormal basis.
c
c       NOTE: the user-supplied coefs array must be sufficiently
c       large to hold the coefficient expansions for the resulting
c       orthonormal basis.
c
c                            Input Parameters:
c
c   eps - precision for the GS procedure
c   disc - a structure describing the discretization of the input 
c       functions, presumably produced by legedisc
c   nfuns - the number of user-supplied input functions
c   funuser - an external subroutine supplying the values of the
c       input functions at a point
c
c       funuser(x,vals,par1,par2,par3,par4)
c
c       Return the values of the user-supplied input functions
c       at the point x in the output array vals.
c
c                           Output Parameters:
c
c   krank - the dimension of the orthonormal basis spanning the input
c       functions
c   coefs - a (ncoefs,krank) each column of which contains the
c       coefficient expansion for one function in the orthonormal
c       basis generated by the GS procedure
c

c
c       Fetch data from the disc structure.
c
        k     = disc(1) 
        nints = disc(2)
        a     = disc(3) 
        b     = disc(4)
c
        ixs   = disc(11) 
        iwhts = disc(12) 
        iu    = disc(13) 
c
        nints = disc(20) 
        iab   = disc(21)        
c
        ncoefs = k*nints
c
c       Call an auxillary subroutine to shape arrays.
c
        call legegs0(eps,nints,disc(iab),k,disc(ixs),disc(iwhts),
     1    disc(iu),ncoefs,nfuns,coefs,funuser,par1,par2,par3,
     2    par4,krank,disc)
c
        end
c
c
c
        subroutine legegs0(eps,nints,ab,k,xs,whts,u,ncoefs,nfuns,
     1    coefs,funuser,par1,par2,par3,par4,krank,disc)
        implicit double precision (a-h,o-z)
        dimension ab(2,nints),xs(k),whts(k),u(k,k),xx(k)
c
        dimension coefs(k,nints,nfuns)
c
        double precision, allocatable :: vals0(:),ders0(:)
        double precision, allocatable :: rnorms(:)
        integer, allocatable :: ipivs(:)
c
c        double precision, allocatable w(:)
c
c        dimension vals0(100 000),ders0(100 000),xx(k),rints(1000)
c        dimension w(50 000 000)
c
        allocate(vals0(nfuns+100),ders0(nfuns+100),rnorms(nfuns+k+100))
        allocate(ipivs(nfuns+100))
c
c       Evaluate the functions at the quadrature nodes.
c
        rints=0
c
        do 1000 j=1,nints
c
        a = ab(1,j)
        b = ab(2,j)
c
        alpha = (b-a)/2.0d0
        beta  = (b+a)/2.0d0
c
        do 1100 i=1,k
        x    = xs(i)*alpha+beta
        wht  = whts(i)*alpha
        call funuser(x,vals0,par1,par2,par3,par4)
        do 1200 l=1,nfuns
        coefs(i,j,l)=vals0(l)*sqrt(wht)
 1200 continue
 1100 continue
 1000 continue
c
c       Perform gs.
c
        call legepiv(coefs,ncoefs,nfuns,eps,rnorms,ipivs,krank)
c        call prinf("ipivs=*",ipivs,krank)
c
        call prinf("in legegs, krank =*",krank,1)
        call prin2("in legegs, rnorms =*",rnorms,krank)
c
c       Scale the functions.
c
        do 2000 j=1,nints
        a = ab(1,j)
        b = ab(2,j)
c
        alpha = (b-a)/2.0d0
        beta  = (b+a)/2.0d0
        do 2100 i=1,k
        ww=1/sqrt(alpha)
        do 2200 l=1,krank
        rn=rnorms(l)
        coefs(i,j,l)=coefs(i,j,l)*ww
 2200 continue
 2100 continue
 2000 continue
c
c       Convert the values of those functions to coefficients.
c
        do 3000 ifun=1,krank
        do 3100 int=1,nints
c
        do 3200 i=1,k
        sum=0
        do 3300 j=1,k
        sum=sum+u(i,j)*coefs(j,int,ifun)
 3300 continue
        xx(i)=sum
 3200 continue
c
        do 3400 j=1,k
        coefs(j,int,ifun)=xx(j)
c*sqrt(rnorms(ifun))
 3400 continue
 3100 continue
 3000 continue
        end
c
c
c
c
c
        subroutine legepiv(b,n,m,eps,rnorms,ipivots,ncols)
        implicit double precision (a-h,o-z)
        dimension rnorms(1),ipivots(1)
        double precision b(n,m),cd
c      
c       Set rank to zero just in case.
c     
        ncols=0
c
c        . . . initialize the array of pivots
c 
        do 1100 i=1,m
        ipivots(i)=i
 1100 continue
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
c
        d=0
        do 1200 j=1,n
        d=d+b(j,i)*(b(j,i))
 1200 continue
        rnorms(i)=sqrt(d)
        dtot=dtot+d
 1400 continue
c      
        thresh=dtot*eps**2
        thresh=sqrt(thresh)
c
        if (dtot .eq. 0) then
           ncols = 0
           return
        endif
c 
c       . . . conduct gram-schmidt iterations
c 
        do 4000 i=1,min(m,n)
c 
c       find the pivot
c 
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c 
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c 
        iijj=ipivots(i)
        ipivots(i)=ipivots(ipivot)
        ipivots(ipivot)=iijj
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
c 
c       orthogonalize the i-th column to all preceding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call legerleascap(b(1,i),b(1,j),n,cd)
c
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call legerleascap(b(1,i),b(1,i),n,cd)
c 
        d=cd
        if(d .lt. thresh**2 ) return
c 
        ncols=i
c 
        d=done/sqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
        if(i .eq. m) goto 3400
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call legerleascap(b(1,i),b(1,j),n,cd)
c
        cd=(cd)
c 
        rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*(b(l,j))
 3000 continue
        rnorms(j)=sqrt(rrn)
 3200 continue
c
 3400 continue
c 
 4000 continue
c 
        return
        end
c
c
c
        subroutine legepiv2(b,n,m,rnorms)
        implicit double precision (a-h,o-z)
        dimension rnorms(1)
        double precision b(n,m),cd
c      
c       Set rank to zero just in case.
c     
        ncols=0
        done=1
c 
c       . . . conduct gram-schmidt iterations
c 
        do 4000 i=1,min(m,n)
c 
c       find the pivot
c 
c       orthogonalize the i-th column to all preceding ones
c 
        do 2780 j=1,i-1
c 
        call legerleascap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
c 
c       normalize the i-th column
c 
        call legerleascap(b(1,i),b(1,i),n,cd)
c
        rnorms(i)=sqrt(cd)
c 
        d=cd
        d=done/sqrt(d)
c
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
        if(i .eq. m) goto 3400
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c
        call legerleascap(b(1,i),b(1,j),n,cd)
c 
        cd=(cd)
c 
        rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*(b(l,j))
 3000 continue
 3200 continue
 3400 continue
c 
 4000 continue
c 
        return
        end

c
c
c
        subroutine legerleascap(x,y,n,prod)
        implicit double precision (a-h,o-z)
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*(y(i))
 1200 continue
        return
        end

c
c
c
        subroutine legecoefs(disc,nfuns,funuser,par1,par2,par3,par4,
     1    coefs)
        implicit double precision (a-h,o-z)
        dimension disc(1),coefs(1)
        external funuser
c
c       This subroutine constructs coefficient expansions for a
c       collection of user-supplied input functions for which a
c       discretization scheme is already available.
c
c                            Input Parameters:
c
c   disc - a structure describing the discretization of the input 
c       functions, presumably produced by legedisc
c   nfuns - the number of user-supplied input functions
c   funuser - an external subroutine supplying the values of the
c       input functions at a point
c
c       funuser(x,vals,par1,par2,par3,par4)
c
c       Return the values of the user-supplied input functions
c       at the point x in the output array vals.
c
c                           Output Parameters:
c
c   coefs - a (ncoefs,nfuns) each column of which contains the
c       coefficient expansion for one input function
c
c
c       Fetch data from the disc structure.
c
        k     = disc(1) 
        nints = disc(2)
        a     = disc(3) 
        b     = disc(4)
c
        ixs   = disc(11) 
        iwhts = disc(12) 
        iu    = disc(13) 
c
        nints = disc(20) 
        iab   = disc(21)        
c
        ncoefs = k*nints
c
c       Call an auxillary subroutine to shape arrays.
c
        call legecoefs0(eps,nints,disc(iab),k,disc(ixs),disc(iwhts),
     1    disc(iu),ncoefs,nfuns,coefs,funuser,par1,par2,par3,
     2    par4,krank)
c

        end
c
c
c
        subroutine legecoefs0(eps,nints,ab,k,xs,whts,u,ncoefs,nfuns,
     1    coefs,funuser,par1,par2,par3,par4,krank)
        implicit double precision (a-h,o-z)
        dimension ab(2,nints),xs(k),whts(k),u(k,k),x1(k),x2(k)
        dimension coefs(k,nints,nfuns),vals0(nfuns),ders0(nfuns)
c
c       Evaluate the functions at the quadrature nodes.
c
        do 1000 j=1,nints
        a = ab(1,j)
        b = ab(2,j)
c
        alpha = (b-a)/2.0d0
        beta  = (b+a)/2.0d0
c
        do 1100 i=1,k
        x = xs(i)*alpha+beta
        wht = whts(i)
        call funuser(x,vals0,par1,par2,par3,par4)
        do 1200 l=1,nfuns
        coefs(i,j,l)=vals0(l)*sqrt(wht)
 1200 continue
 1100 continue
 1000 continue
c
c       Convert the values of those functions to coefficients.
c
        do 2000 ifun=1,nfuns
        do 2100 int=1,nints
c
        do 2200 i=1,k
        sum=0
        sum2=0
        do 2300 j=1,k
        sum=sum+u(i,j)*coefs(j,int,ifun)
 2300 continue
        x1(i)=sum
 2200 continue
c
        do 2400 j=1,k
        coefs(j,int,ifun)=x1(j)
 2400 continue
 2100 continue
 2000 continue
c
        end
c
c
c
        subroutine legevals(disc,nfuns,coefs,vals)
        implicit double precision (a-h,o-z)
        dimension disc(1),coefs(1)
        external funuser
c
c
c       Fetch data from the disc structure.
c
        k     = disc(1) 
        nints = disc(2)
        a     = disc(3) 
        b     = disc(4)
c
        ixs   = disc(11) 
        iwhts = disc(12) 
        iu    = disc(13) 
c
        nints = disc(20) 
        iab   = disc(21)        
c
        ncoefs = k*nints
c
c       Call an auxillary subroutine to shape arrays.
c
        call legevals0(nints,disc(iab),k,disc(ixs),disc(iwhts),
     1    disc(iu),disc(iv),ncoefs,nfuns,coefs,vals)
c

        end
c
c
c
        subroutine legevals0(nints,ab,k,xs,whts,u,v,ncoefs,nfuns,
     1    coefs,vals)
        implicit double precision (a-h,o-z)
        dimension ab(2,nints),xs(k),whts(k),u(k,k),x1(k),x2(k)
        dimension coefs(k,nints,nfuns),vals0(nfuns),ders0(nfuns)
        dimension v(k,k),vals(k,nints,nfuns)
c
c
        do 2000 ifun=1,nfuns
        do 2100 int=1,nints
c
        do 2200 i=1,k
        sum=0
        sum2=0
        do 2300 j=1,k
        sum=sum+v(i,j)*coefs(j,int,ifun)
 2300 continue
        x1(i)=sum
 2200 continue
        do 2400 j=1,k
        vals(j,int,ifun)=x1(j)
 2400 continue
 2100 continue
 2000 continue
c
        end

c
c
c
        subroutine leged_eval(disc,nfuns,coefs,x,vals,ders,ifders)
        implicit double precision (a-h,o-z)
        dimension disc(1),vals(1),ders(1)
c
c       Evaluate a collection of functions represented as piecewise
c       coefficient expansions at a single point.
c
c
c
c
c
c       Fetch data from the disc structure.
c
        k     = disc(1) 
        nints = disc(2)
        a     = disc(3) 
        b     = disc(4)
c
        ixs   = disc(11) 
        iwhts = disc(12) 
        iu    = disc(13) 
c
        nints = disc(20) 
        iab   = disc(21)        
        last  = disc(22)
c
c       Call an auxillary routine to perform the evaluations.
c
        call legeeval0(k,nints,disc(iab),last,nfuns,coefs,x,vals,
     1    ders,ifders)
c
        disc(22) = last
        end
c
c
c
        subroutine legeeval0(k,nints,ab,last,nfuns,coefs,x,vals,ders,
     1    ifders)
        implicit double precision (a-h,o-z)
        dimension coefs(k,nints,nfuns),ab(2,nints),vals(1),ders(1)
        dimension vals0(k),ders0(k+1000)
c
c       Find the interval containing the point.
c
        call findint(x,nints,ab,last,int)
c
        a=ab(1,int)
        b=ab(2,int)
c
c       Map the point x back to [-1,1].
c
        alpha = 2/(b-a)
        beta  = (b+a)/(a-b)
c
        u = alpha*x+beta
c
c$$$        do 2000 j=1,nfuns
c$$$        call legefder(u,vals(j),ders(j),coefs(1,int,j),k-1)
c$$$        ders(j)=ders(j)*alpha
c$$$ 2000 continue
c$$$c
c$$$        return
c
c       Sum the expansions.
c
        if (ifders .eq. 0) then
c
        call lege(k-1,u,vals0)
c
        do 1000 j=1,nfuns
        vals(j)=0
 1000 continue
c
        do 1100 i=1,k
        val0=vals0(i)
        do 1200 j=1,nfuns
        vals(j)=vals(j)+coefs(i,int,j)*val0
 1200 continue
 1100 continue
c
        else
c
        call legeders(k-1,u,vals0,ders0)
c
        do 3000 j=1,nfuns
        vals(j)=0
        ders(j)=0
 3000 continue
c
        do 3100 i=1,k
        val0=vals0(i)
        der0=ders0(i)
        do 3200 j=1,nfuns
        vals(j)=vals(j)+coefs(i,int,j)*val0
        ders(j)=ders(j)+coefs(i,int,j)*der0
 3200 continue
 3100 continue
c     
        do 3300 j=1,nfuns
        ders(j)=ders(j)*alpha
 3300 continue
        endif
c
        end
c
c
c
        subroutine findint(x,nints,ab,last,int)
        implicit double precision (a-h,o-z)
        dimension ab(2,1)
c
        ier = 0
c
c       Check to see if we are in the last accessed interval.
c
        a = ab(1,last)
        b = ab(2,last)
c
        if(x .ge. a .AND. x .le. b) then
        int = last
        return
        endif
c
c       Next perform a partial "binary-like" search in an attempt
c       to locate the subinterval containing x.
c
 1000 continue
c
        il=1
        ir=nints
c
        do 1100 j=1,10
        imiddle = (il+ir)/2
        a = ab(1,imiddle)
        if (a .gt. x) ir=imiddle
        if (a .lt. x) il=imiddle
 1100 continue
c
c       Continue the search for interval.
c
        int=il
c        int=1
c
 2000 continue
c
        if (int .eq. nints) goto 2100
c
        a = ab(1,int)
        b = ab(2,int)
c
        if (x .lt. a) goto 2100
        if (x .lt. b) goto 2100
        int=int+1
        goto 2000
c
 2100 continue
        last = int
        end
c
c
        subroutine legedisc_quad(disc,nquad,xs,whts)
        implicit double precision (a-h,o-z)
        dimension disc(1),xs(1),whts(1)
c
c       Fetch data from the disc structure.
c
        k     = disc(1) 
        nints = disc(2)
        a     = disc(3) 
        b     = disc(4)
c
        ixs   = disc(11) 
        iwhts = disc(12) 
        iu    = disc(13) 
c
        nints = disc(20) 
        iab   = disc(21)        
        last  = disc(22)
c
c       Call an auxillary routine to perform the evaluations.
c
        call legedisc_quad0(k,nints,disc(iab),disc(ixs),disc(iwhts),
     1    nquad,xs,whts)
c
        end
c
c
c
        subroutine legedisc_quad0(k,nints,ab,xslege,whtslege,nquad,xs,
     1     whts)
        implicit double precision (a-h,o-z)
        dimension ab(2,nints),xslege(1),whtslege(1)
        dimension xs(1),whts(1)
c     
        nn=1
        do 1000 j=1,nints
        a = ab(1,j)
        b = ab(2,j)
c
        alpha = (b-a)/2.0d0
        beta  = (b+a)/2.0d0
c
        do 1100 i=1,k
        x = xslege(i)*alpha+beta
        ww = whtslege(i)*alpha
        xs(nn)=x
        whts(nn)=ww
        nn=nn+1
 1100 continue
 1000 continue
c
        nquad=nn-1
        end


        subroutine legemove(k,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
        do 1000 j=1,k
        b(j)=a(j)
 1000 continue
        end
