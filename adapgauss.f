        implicit double precision (a-h,o-z)      
        external funuser,funuser_dumb
c
        call prini(6,13)
c
        pi=acos(-1.0d0)
        a=0d0
        b=2d0
c
        m = 20
        eps=1.0d-14
c
        call adapgauss(ier,a,b,eps,m,val,funuser,par1,par2,par3,
     1      par4,par5,par6,par7,par8)
c
        call prinf("after adapgauss, ier = *",ier,1)
        call prin2("after adapgauss, val = *",val,1)
c
        print *,val
c
c       NOTE: answer checked via mathematica.
c

        end
c
c
c
        subroutine funuser_dumb(x,val,n,par2,par3,par4,par5,par6,
     1     par7,par8)
        implicit double precision (a-h,o-z)
        val=1
        end
c
c
c       
        subroutine funuser(x,val,par1,par2,par3,par4,par5,par6,par7,
     1    par8)
        implicit double precision (a-h,o-z)
        val=log(x)*cos(666*x) + sin(34*x)/sqrt(x)
        end
c
c
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and the beginning of the
c       code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is an adaptive Gaussian quadrature code for integrating
c       a user-supplied function specified via an external SUBROUTINE.
c       The following subroutines are user-callable:
c 
c   adapgauss - this procedure integrates a user-supplied function to
c       the prescribed accuracy using adaptive Legendre quadratures
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       
        subroutine adapgauss(ier,a,b,eps,m,val,funuser,par1,par2,par3,
     1    par4,par5,par6,par7,par8)
        implicit double precision (a-h,o-z)
        double precision, allocatable :: w(:),xs(:),whts(:)
        external funuser
c

c
c       This procedure uses adaptive Legendre integration to evaluate
c       an integral of the form
c
c           \int_a^b f(x) dx
c
c       where f(x) is a function supplied by the user via an external
c       subroutine.
c
c       NOTE: the failure of this procedure is usually caused by 
c       asking for an unreasonable level of precision.
c
c                              Input Parameters:
c
c   (a,b) - the domain of integration
c   eps - precision for the computation 
c   m - order of the quadrature formula 
c   funuser - a user-supplied external subroutine with calling
c       sequence:
c
c          subroutine funuser(x,val,par1,par2,par3,par4,par5,par6,par7,
c            par8)
c
c          Return the value of the function f at the point x.
c
c   par? - user-supplied arbitrarily typed variables
c
c                              Output Parameters:
c
c   ier - error return code;
c       ier = 0    means the procedure was successful
c       ier = 4    means that the internal work array was insufficient
c       ier = 512  means that maximum recursion depth was exceeded                   
c
c
        ier = 0 
        done = 1
        dtwo = 2
        MAXDEPTH = 1000
c
        allocate(xs(1000),whts(1000),w(1000000))
c
c       Initialize the quadrature, if necessary.
c
        call legequad(m,xs,whts)
c
        val = 0
c
c       Initialize a stack at the end of the work array.
c
        istack0 = 1 000 000
        istack  = istack0
c
        call oneint(a,b,m,val0,funuser,par1,par2,par3,par4,par5,
     1    par6,par7,par8,xs,whts)
c
        ilevel=0
c
        w(istack)=a
        w(istack-1)=b
        w(istack-2)=val0
        w(istack-3)=ilevel
c
        istack=istack-4
c
 1000 continue
c
c       Pop an element off the stack, assuming it isn't empty.
c
        if (istack .eq. istack0) return
c
        istack=istack+4
c
        a0=w(istack)
        b0=w(istack-1)
        val0=w(istack-2)
        ilevel=w(istack-3)
c
        if (ilevel .gt. MAXDEPTH) then
           ier = 512
           return
        endif
c
c       Compute the integral over the left and right half of the
c       interval ...
c
        c0 = (a0+b0)/dtwo
c
        call oneint(a0,c0,m,vall,funuser,par1,par2,par3,par4,par5,
     1    par6,par7,par8,xs,whts)
c
        call oneint(c0,b0,m,valr,funuser,par1,par2,par3,par4,par5,
     1    par6,par7,par8,xs,whts)
c
c       Compare the two approximations we have obtained ... and
c       if we have integrated with sufficient accuracy, add the
c       contribution to the output value.
c
        val1 = valr+vall
        errabs = abs(val0-val1)
c
        if (errabs .le. eps) then
           val=val+val0
           goto 1000
        endif
c
c       Otherwise push the two children onto the stack and continue.
c
        if (istack .lt. 9) then
           ier = 4
           return
        endif
c
        w(istack)=a0
        w(istack-1)=c0
        w(istack-2)=vall
        w(istack-3)=ilevel+1
        istack=istack-4
c
        w(istack)=c0
        w(istack-1)=b0
        w(istack-2)=valr
        w(istack-3)=ilevel+1
        istack=istack-4
c
        goto 1000
c
        end
c
c
c
        subroutine oneint(a,b,m,val,funuser,par1,par2,par3,par4,par5,
     1    par6,par7,par8,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),whts(1)
        external funuser
c
        dtwo=2
        alpha=(b-a)/dtwo
        beta=(b+a)/dtwo
c
        val=0
c
        do 1000 j=1,m
        x=alpha*xs(j)+beta
        wht=alpha*whts(j)        
        call funuser(x,f,par1,par2,par3,par4,par5,par6,par7,par8)
        val=val+f*wht
 1000 continue
c
c
        end
c
c
c
