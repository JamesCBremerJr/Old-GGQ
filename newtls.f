        implicit double precision (a-h,o-z)
c
        call prini(6,13)
c
c
        n = 31
        m = 1010
c
        call test_ls(n,m)

        end
c
c
c
        subroutine test_ls(n,m)       
        implicit double precision (a-h,o-z)
        dimension a(n,m),u(n,n),t(n,n),v(n,m)
        dimension xs(m),whts(m),ut(n,n),a2(n,m)
        dimension y(n),xtrue(m),x(m),y2(n)
c
c       Form the test matrix a.
c
c        itype=1
c        call legeexps(itype,m,xs,u,v,whts)
c
        call legequad(m,xs,whts)
c
        do 1000 j=1,m
        do 1100 i=1,n
        x   = xs(j)
        wht = whts(j)
        call lege0(i-1,x,val,der)
        a(i,j)=val*sqrt(wht)
 1100 continue
 1000 continue
c
c       Form a sum of the columns in y.
c
        do 2000 i=1,n
        y(i)=0
 2000 continue
c
c        call corrand(m,xtrue)
        do 2050 j=1,m
        xtrue(j)=0
 2050 continue
        xtrue(1)=1
c
        do 2100 j=1,m
        do 2200 i=1,n
        y(i)=y(i)+a(i,j)*xtrue(j)
 2200 continue
 2100 continue
c
c       Factor it.
c
        eps=1.0d-15
        call newtls_factor(eps,n,m,a,u,t,v,krank)
c
c       Test the factorization.
c
        call gmatmul(n,krank,n,ut,u,t)
        call gmatmul2(n,krank,m,a2,ut,v)
c
        errl2=0
        do 3000 i=1,n
        do 3100 j=1,m
        d=(a(i,j)-a2(i,j))**2
        errl2=errl2+d
 3100 continue
 3000 continue
        errl2=sqrt(errl2)
        call prin2("factorization errl2=*",errl2,1)
c
        call newtls_solve(n,m,krank,u,t,v,x,y)
c
        stop
c
c       Test the solution.
c
        errl2=0
        do 4000 i=1,n
        sum=0
        do 4100 j=1,m
        sum=sum+a(i,j)*x(j)
 4100 continue
        y2(i)=sum
        dl2=dl2+sum**2
        errl2=errl2+(sum-y(i))**2
 4000 continue
        errl2=sqrt(errl2)
        call prin2("solve errl2=*",errl2,1)
c
c
        dnorm1=0
        dnorm2=0
        do 4200 i=1,m
        dnorm1=dnorm1+xtrue(i)**2
        dnorm2=dnorm2+x(i)**2
 4200 continue
        dnorm1=sqrt(dnorm1)
        dnorm2=sqrt(dnorm2)
        call prin2("original norm = *",dnorm1,1)
        call prin2("least norm = *",dnorm2,1)
c
        end
c
c
c
        subroutine gmatmul(n,k,m,a,b,c)
        implicit double precision (a-j,o-z)
        dimension a(n,m),b(n,k),c(k,m)
c
c       Compute a(n,m) = b(n,k) * c(k,m)
c
        do 1000 i=1,n
        do 1100 j=1,m
        sum=0
        do 1200 l=1,k
        sum=sum+b(i,l)*c(l,j)
 1200 continue
        a(i,j)=sum
 1100  continue
 1000  continue
      
        end
c
c
c
        subroutine gmatmul2(n,k,m,a,b,c)
        implicit double precision (a-j,o-z)
        dimension a(n,m),b(n,k),c(m,k)
c
c       Compute a(n,m) = b(n,k) * c^*(m,k)
c
        do 1000 i=1,n
        do 1100 j=1,m
        sum=0
        do 1200 l=1,k
        sum=sum+b(i,l)*c(j,l)
 1200 continue
        a(i,j)=sum
 1100  continue
 1000  continue
      
        end

c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the code
c       proper.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains utlity subroutines used by the Gaussian
c       quadrature codes newton1d.f and newton2d.f.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine newtls_factor(eps,n,m,a,u,t,v,krank)
        implicit double precision (a-h,o-z)
        dimension a(n,m),u(n,n),t(n,n),v(m,n)
c
        dimension rnorms(m+n),ipivs(m+n)
c
        double precision, allocatable :: q(:,:),rt(:,:)
c
c       Factor a (n,m) matrix a of rank

c           a(n,m) = u(n,krank) t(krank,krank) v^t (krank,m)
c
c       where t is lower triangular and the columns of u and v 
c       comprise of orthonormal bases for the column and row spaces
c       of the input matrix a (respectively).  The decomposition is 
c       obtained via the pivoted Gram-Schmidt with reorthogonalization
c       algorithm.
c                   
c                             Input Parameters:
c
c                            Output Parameters:
c
c
c
c
        allocate(q(n,m))
c
c       Copy the matrix a into q and perform GS.
c
        call newtls_move(n*m,a,q)
c
c        call gspiv2(q,n,m,rnorms,ipivs)
c
        call gspiv(q,n,m,eps,rnorms,ipivs,krank)
        cond = rnorms(1)/rnorms(krank)
c        call prin2("in newtls, cond = *",cond,1)
c        call prin2("in newtls_factor, rnorms = *",rnorms,krank)
c
c       Copy q into u.
c
        call newtls_move(n*krank,q,u)
c
c       Form r'(m,krank) = a'(m,n) q(n,krank)
c
        allocate(rt(m,krank))
c
        do 1000 i=1,m
        do 1100 j=1,krank
        sum=0
        do 1200 l=1,n
        sum=sum+a(l,i)*q(l,j)
 1200 continue
        rt(i,j)=sum
 1100 continue
 1000 continue
c
c       Perform GS on rt to form v.
c      
        call newtls_move(m*krank,rt,v)
        call gspiv2(v,m,krank,rnorms,ipivs)
c
c        call prin2("in newtls,rnorms=*",rnorms,krank)
c
c
c       Form t = r(krank,m) v(m,krank).
c
        do 2000 i=1,krank
        do 2100 j=1,krank
        sum=0
        do 2200 l=1,m
        sum=sum+v(l,i)*rt(l,j)
 2200 continue
        t(j,i)=sum
 2100 continue
 2000 continue
c
        end
c
c
c
        subroutine newtls_solve(n,m,krank,u,t,v,x,b)
        implicit double precision (a-h,o-z)
        dimension u(n,n),t(n,n),v(m,n)
        dimension y(n),x(m),tt(n,n),b(n)
c
c       Find a minimum l^2 norm solution of the nondegenerate linear 
c       system
c
c            a(n,m) * x(m) = y(n)
c
c       given a decomposition of the form (1) for the matrix a.
c
c
c                        Input Parameters:
c
c                        Output Parameters:
c

c
c       Apply u^* to b to form y.
c
        do 1000 i=1,n
        sum=0
        do 1100 l=1,n
        sum=sum+u(l,i)*b(l)
 1100 continue
        y(i)=sum
 1000 continue
c
c       Solve Tz=y.
c
        call newtls_move(n*n,t,tt)
        call qrsolv(tt,n,y,rcond)
c
c       Form x=Vy.
c
        do 2000 i=1,m
        sum=0
        do 2100 l=1,n
        sum=sum+v(i,l)*y(l)
 2100 continue
        x(i)=sum
 2000 continue
        end
c
c
c
        subroutine newtls_solve2(n,m,a,v,x,y)
        implicit double precision (a-h,o-z)
        dimension a(n,m),v(m,n),b(n,n),x(m)
        dimension y(n,n),z(n)
c
c       Find a minimum norm solution of the underdetermined linear
c       system
c
c          a(n,m) x = y 
c     
c       subject to the constraint that x is in the span of a given
c       orthonormal basis of dimension n.

c
c       Form the matrix b(n,n)=a(n,m)*v(m,n).
c
        do 1000 i=1,n
        do 1100 j=1,n
        sum=0
        do 1200 l=1,m
        sum=sum+a(i,l)*v(l,j)
 1200 continue
        b(i,j)=sum
 1100 continue
 1000 continue
c
c       Solve the system bz=y.
c
        call newtls_move(n,y,z)
        call qrsolv(b,n,z,rcond)
c
c       Now form x=Vz.
c
        do 2000 i=1,m
        sum=0
        do 2100 l=1,n
        sum=sum+v(i,l)*z(l)
 2100 continue
        x(i)=sum
 2000 continue
c
        end
c
c
c
        subroutine newtls_move(k,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
        do 1000 j=1,k
        b(j)=a(j)
 1000 continue
        end
