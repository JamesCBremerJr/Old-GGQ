cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains two user-callable procedures for performing
c       the pivoted Gram-Schmidt (GS) algorithm with reorthogonalization.
c
c       The following subroutines are user-callable:
c
c   gspivr - randomized pivoted gram-schmidt with reorthogonalization
c
c   gspiv_skel - compute a column skeleton for a user-supplied input
c       matrix; this version perserves the input matrix
c
c
c   gspiv - perform pivoted Gram-Schmidt with reorthorgonalization
c       on a rank deficient rectangular matrix
c
c   gspiv2 - perform pivoted Gram-Schmidt with reorthogonalization on
c       an (m,n) matrix assumed to have rank m
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine gspiv_skel(ier,n,m,a,eps,rnorms,ipivots,
     1      krank,w,lw)
        implicit double precision (a-h,o-z)
        dimension rnorms(1),ipivots(1),w(1),a(n,1)
c
c       This subroutine uses the "fast randomized skeletonization"
c       to find a column skeleton for a user-supplied input matrix.  If
c       the input matrix is small or its rank is large compared to its 
c       size, it will "fall back" to using standard gram-schmidt.
c
c       NOTE: this procedure does not destroy the input matrix and its
c       output is a list of columns which constitute the skeleton,
c       not the columns themselves
c
c                          Input parameters:
c
c  n,m - the dimensionalities of the matrix to be skeletonized
c  a - the matrix to skeletonize 
c
c  eps - the accuracy to which the matrix is to be skeletonized
c  lw - the length of the user-supplied work array w in double precision words
c
c                         Output parameters:
c
c  ier - error return code;
c        ier = 0  means that the routine was successful
c        ier = 4  means that the rank of the input matrix is too
c                 close to its size for the randomized algorithm to be
c                 efficient
c        ier = 16 means that the work array was of insufficient length
c
c  rnorms - the vector of normalizations produced by the pivoted
c        Gram-Schmidt process. Returned for the user's edification.
c  ipivots - a column skeleton of the matrix A
c  ncols - the rank of A to precision eps
c  lused - the length of the part of the work array w actually used 
c        by the subroutine (in double precision elements)
c  
c
        ier=0
c
c       Check to see if the matrix is small.
c
        goto 3000
        if (n .le. 60 .OR. n*m .le. 2000) goto 3000
c
c       Initialize the fast ``random enough'' transform.
c
        nsteps=6
        call random_transf_init(nsteps,n,w,lkeep)
c
c       Allocate memory
c
        ix = 1+lkeep
        lx = n
c
        iw2=ix+lx
        lw2=lw-iw2
c
        kmax = min(n,(lw2-10 000)/(2*m))
c
        ia=ix+lx
        la=kmax*m
c
        ib=ia+la
        lb=kmax*m
c
        if (kmax .le. 40) then
           ier = 16
           return
        endif
c        
        if (ib+lb .ge. lw) then
           ier = 4
           return
        endif
c
c       Apply the random transform to the subblock, storing as 
c       many output variables as possible.
c
        iout=0
        do 2000 i=1,m
        call random_transf(a(1,i),w(ia+iout),w)
        iout=iout+kmax
 2000 continue
c
        k=40
c     
        do 1000 iter=1,100000
c
        if (lw2 .lt. k*m) then
           ier = 16
           return
        endif
c
c       Call an auxillary routine to do the dirty work.
c
        call gspiv_rskel0(eps,n,m,kmax,k,w(ia),w(ib),krank,rnorms,
     1     ipivots)
c
        if (krank .lt. k-20) goto 1200
c
        k=k*2
c
        if (k .ge. kmax) then
           goto 3000
           return
        endif
c
 1000 continue
 1200 continue
c
c       Compute a skeleton via brute force.
c
 3000 continue
        if (lw .le. n*m) then
           ier = 4
           return
        endif
c
        call gspiv_move(n*m,a,w)
        call gspiv(w,n,m,eps,rnorms,ipivots,krank)
        end
c
c
c
        subroutine gspiv_rskel0(eps,n,m,kmax,k,a,b,krank,rnorms,
     1     ipivots)
        implicit double precision (a-h,o-z)
        dimension a(kmax,m),b(k,m),rnorms(1),ipivots(1)
c
c       Copy k output variables into the work array.
c
        do 1000 j=1,m
        do 1100 i=1,k
        b(i,j)=a(i,j)
 1100 continue
 1000 continue
c
c       Gram-schmidt it.
c
        call gspiv(b,k,m,eps,rnorms,ipivots,krank)
        end
c
c
c
        subroutine gspivr(ier,a,n,m,eps,rnorms,ipivots,ncols,w,lw)
        implicit double precision (a-h,o-z)
        dimension rnorms(1),ipivots(1),w(1),a(n,m)
c
c
c
        ier = 0
c
c       If the matrix is small, just call gspiv.
c     
        if ( n .le. 40 .OR. n*m .le. 2000) goto 5000

        goto 5100
 5000 continue
        call gspiv(a,n,m,eps,rnorms,ipivots,ncols)
        return
c
c       Find the skeleton.
c 
 5100 continue
        call gspiv_skel(ier,n,m,a,eps,rnorms,ipivots,
     1      ncols,w,lw)
c
        if (ier .eq. 32) goto 5000
        if (ier .ne. 0 ) return
c
c       Copy the skeleton into the work array.
c
        ib = 1
        lb = n*ncols
c
        ipivs = ib+lb
        lpivs = ncols
c
        iw2 = ipivs+lpivs
        lw2 = lw-iw2
c
        if (lw2 .le. 0) then
           ier = 4
           return
        endif
c
        ic = ib
        do 1000 j=1,ncols
        do 1100 i=1,n
        w(ic) = a(i,ipivots(j))
        ic=ic+1
 1100 continue
 1000 continue
c      
c       Gram-schmidt it and copy the resultings into the matrix a.
c
        call gspiv2(w(ib),n,ncols,rnorms,w(ipivs))
c
        ic=ib
        do 2000 j=1,ncols
        do 2100 i=1,n
        a(i,j)=w(ic)
        ic=ic+1
 2100 continue
 2000 continue
c
        end
c
c
c
        subroutine gspiv(b,n,m,eps,rnorms,ipivots,ncols)
        implicit double precision (a-h,o-z)
        save
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
c$$$        if(dtot .le. eps**2) then
c$$$           ncols=0
c$$$           return
c$$$        endif

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
        call gsrleascap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call gsrleascap(b(1,i),b(1,i),n,cd)
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
        call gsrleascap(b(1,i),b(1,j),n,cd)
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
 3400 continue
c 
 4000 continue
c 
        return
        end
c
c
c
        subroutine gspiv2(b,n,m,rnorms,ipivots)
        implicit double precision (a-h,o-z)
        dimension rnorms(1),ipivots(1)
        double precision b(n,m),cd
c 
c       This subroutine extracts from the rectangular matrix b
c       (user-supplied) a square full-rank submatrix. The number
c       of columns m of b is expected to be greater than the number
c       n of its columns. The subroutine also returns the integer
c       array ipivots, containing the sequence numbers of columns
c       selected, and the real array rnorms, containing the
c       denominators by which the selected columns had to be
c       divided during the Gram-Schmidt process (a decent measure
c       of stability of the process).
c 
c                    input paramneters:
c 
c  b - the matrix to be gram-schmidt'ed. it is destroyed by this
c       subroutine
c  n,m - dimensionalities of the matrix b
c 
c                     output parameters:
c 
c  b - the matrix of gram-schmidt vectors of the matrix a. note
c        that on exit from this subroutine only the first
c        n columns of  b  are meaningful, the rest making no
c        sense whatsoever, so that effectively on exit b is dimensioned
c        b(n,n)
c  rnorms - the normalizing factors in the gram-schmidt process.
c        Only the first n of them are meaninful, but the array
c        has to be dimensioned by the user to be at least m+1
c        double precision elements long.
c 
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
        do 1400 i=1,m
c 
        d=0
        do 1200 j=1,n
        d=d+b(j,i)*(b(j,i))
 1200 continue
        rnorms(i)=sqrt(d)
 1400 continue
c 
c       . . . conduct gram-schmidt iterations
c
        mm=min(m,n)
        do 4000 i=1,mm
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
        call gsrleascap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call gsrleascap(b(1,i),b(1,i),n,cd)
c 
        d=cd
cccc        ncols=i
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
        call gsrleascap(b(1,i),b(1,j),n,cd)
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
 3400 continue
c 
 4000 continue
c 
        return
        end
c
c
c
        subroutine gsrleascap(x,y,n,prod)
        implicit double precision (a-h,o-z)
        save
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
        subroutine gspiv_move(n,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
        do 1000 j=1,n
        b(j)=a(j)
 1000 continue
        end       
