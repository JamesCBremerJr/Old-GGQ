        implicit double precision (a-h,o-z)
        dimension xs0(100 000),ys0(100 000),whts0(100 000)
        dimension xs(10 000),ys(10 000),whts(10 000)
        dimension xssym(10 000),yssym(10 000),whtssym(10 000)
c
        dimension w(1000)
        dimension xs1(10000),ys1(10000),whts1(10000)
        dimension xs2(10000),ys2(10000),whts2(10000)
        dimension xs3(10000),ys3(10000),whts3(10000)
c
        dimension xstri(10 000),ystri(10 000),whtstri(10 000)
c
        external funeval,funeval2,funeval3,funeval4,funtest
c
        call prini(6,13)
        call mach_zero(eps)
c
        epscheb = 1.0d-30
        eps     = 1.0d-30
c
        norder = 16
        par1   = norder
c
        npols = (norder+1)*(norder+2)/2
        call funeval(x,y,par1,par2,w,w,w)
        nfuns=par2
c
        d=nfuns
c
        call prinf("norder = *",norder,1)
        call prinf("npols = *",npols,1)
        call prinf("nfuns = *",nfuns,1)
c
        ngoal = (npols+2)/3
        call prinf("ngoal = *",ngoal,1)
c
c       Build an oversampled quadrature for the functions.
c
!        call tensorsquare2(norder,nquad0,xs0,ys0,whts0)
        nnn = norder
 0011 continue
        nnn = nnn + 1
        call tensortri(nnn,nquad0,xs0,ys0,whts0)
        call prinf("after tensorquad, nquad0 = *",nquad0,1)
c
c       Construct a Chebyshev quadrature.
c
        call prinf("nfuns=*",nfuns,1)
        call chebquad(epscheb,nfuns,funeval,par1,par2,nquad0,xs0,ys0,
     1    whts0,nquad1,xs1,ys1,whts1)
c

cccccccccccccccccccccccc
c       Construct Gaussian quadrature.
ccccccccccccccccccccccccc
c 
c       Take the eight-fold path.
c
        ifcauchy     = 0
        ifonlysignifs= 0
        iregion      = 3
c
        d = ngoal
        d = d/8.0
        nn = ceiling(d)
c
        call gaussquad(eps,nfuns,funeval,par1,par2,nquad1,xs1,
     1    ys1,whts1,iregion,ifcauchy,ifonlysignifs,nn)
c
c       Reflect the formula onto the other parts of the square.
c

        nquad2=0
        do 1000 i=1,nquad1
c
        x   = xs1(i)
        y   = ys1(i)
        wht = whts1(i)
c
        if (x .lt. 0 .OR. x .gt. 1.0d0) goto 011
        if (y .lt. 0 .OR. y .gt. x)     goto 011 
        if (wht .lt. 0   )           goto 0011
c
c       Four way-symmetry.
c
        nquad2=nquad2+1
        xs2(nquad2)=x
        ys2(nquad2)=y
        whts2(nquad2)=wht
c
        nquad2=nquad2+1
        xs2(nquad2)=y
        ys2(nquad2)=x
        whts2(nquad2)=wht
c
 1000 continue
c
c
c       Now try 4-way symmetric iterations.
c
        par1 = norder
        call funeval2(x,y,par1,par2,w,w,w)
        nfuns2=par2
c
        ifcauchy=0
        ifonlysignifs=0
        iregion=2
c
        d = ngoal
        d = d/4.0
        nn = ceiling(d)
c
        call gaussquad(eps,nfuns2,funeval2,par1,par2,nquad2,xs2,
     1    ys2,whts2,iregion,ifcauchy,ifonlysignifs,nn,
     1    npts,ilist,ifreplay)
c
c
c       Reflect it (again).
c
        nquad3 = 0
        do 1100 j=1,nquad2
        x = xs2(j)
        y = ys2(j)
        wht=whts2(j)
c

        if (x .lt. 0 .OR. x .gt. 1.0d0)     goto 011
        if (y .lt. 0 .OR. y .gt. 1.0d0)     goto 011 
        if (wht .lt. 0   )                  goto 0011

c
        nquad3=nquad3+1
        xs3(nquad3)=x
        ys3(nquad3)=y
        whts3(nquad3)=wht
c
        nquad3=nquad3+1
        xs3(nquad3)=-x
        ys3(nquad3)=y
        whts3(nquad3)=wht
 1100 continue
c
c
c       Now try 2-way symmetric iterations.
c
        par1 = norder
        call funeval3(x,y,par1,par2,w,w,w)
        nfuns3=par2
c
        d = ngoal
        d = d/2.0
        nn = ceiling(d)
c
        ifcauchy=0
        ifonlysignifs=0
        iregion=1
        call gaussquad(eps,nfuns3,funeval3,par1,par2,nquad3,xs3,
     1    ys3,whts3,iregion,ifcauchy,ifonlysignifs,nn)
c
c       Reflect it (one last time).
c
        nquad = 0
        do 2000 j=1,nquad3
        x = xs3(j)
        y = ys3(j)
        wht=whts3(j)
        nquad=nquad+1
        xs(nquad)=x
        ys(nquad)=y
        whts(nquad)=wht
c
        nquad=nquad+1
        xs(nquad)=x
        ys(nquad)=-y
        whts(nquad)=wht
 2000 continue
c
c
c
c       Conduct Newton iterations on the whole square to polish
c       off remaining points.
c
        nfuns4 = (norder+1)*(norder+2)/2
        par1 = norder
c
        iregion=0
        ifcauchy=0
        ifonlysignifs=0
c
        call gaussquad(eps,nfuns4,funeval4,par1,par2,nquad,xs,
     1    ys,whts,iregion,ifcauchy,ifonlysignifs,ngoal)
c
c       Use them evil Cauchy iterations to try to eliminate a few more
c       points.
c
c$$$        nfuns4 = (norder+1)*(norder+2)/2
c$$$        par1 = norder
c$$$c
c$$$        iregion     =1
c$$$        ifcauchy    =1
c$$$        ifonlysignifs=0
c$$$c
c$$$        call gaussquad(eps,nfuns4,funeval4,par1,par2,nquad,xs,
c$$$     1    ys,whts,iregion,ifcauchy,ifonlysignifs,ngoal)
c
c       We're done.
c
        call prinf("nquad final =*",nquad,1)
        call prin2("xs final= *",xs,nquad)
        call prin2("whs final= *",whts,nquad)
        call prina("*")
        call prina("*")
        call prina("*")
c
        do j=1,nquad
        x=xs(j)
        y=ys(j)
        wht=whts(j)
        if (x .lt. -1 .OR. x. gt. 1) goto 0011
        if (y .lt. -1 .OR. y. gt. 1) goto 0011
        if(wht .lt. 0) goto 011
        end do
c
c       Test the formula.
c
        call tensorsquare(norder,nquad0,xs0,ys0,whts0)
        errl2 = 0
c

        do 3000 i1=0,norder
        do 3100 i2=0,norder-i1
c
        sum0=0
        do 3200 j=1,nquad0
        x=xs0(j)
        y=ys0(j)
        wht=whts0(j)

        call funtest(x,y,val,i1,i2)
        sum0=sum0+wht*val
 3200 continue
c
        sum1=0
        do 3300 j=1,nquad
        x=xs(j)
        y=ys(j)
        wht=whts(j)
        call funtest(x,y,val,i1,i2)
        sum1=sum1+wht*val
 3300 continue
c
c
        errl2=errl2+abs(sum0-sum1)**2
 3100 continue
 3000 continue
c
        errl2=sqrt(errl2)
        call prin2("errl2=*",errl2,1)
c
c       Write the quadrature to disk.
c
 0001 continue

        call write_square(norder,nquad,xs,ys,whts)
        call write_square90(norder,nquad,xs,ys,whts)
c
        end
c
c
c
        subroutine write_square(norder,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
c
        character*12 text
c
c
 1000 format ("square",I2.2,'.txt')
c
 1100 format (8x,'nquad  = ',I3)
 1200 format (8x,'norder = ',I3)
c
 1300 format (5x,3x,'data ',A,' /')
 1400 format ((5x,'-',2x,D38.30,','))
 1450 format ((5x,'-',2x,D38.30,'/'))

 1500 format (5x,3x,'data ',a4,I2,' /')

c
        write(text,1000) norder
c
        iw=100
        open(iw,FILE=text)
c
        write (iw,1100) nquad
        write (iw,1200) norder
c
        write (iw,1300) 'xs'
        write (iw,1400) (xs(i),i=1,nquad-1)
        write (iw,1450) xs(nquad)

c
        write (iw,1300) 'ys'
        write (iw,1400) (ys(i),i=1,nquad-1)
        write (iw,1450) ys(nquad)

c
        write (iw,1300) 'whts'
        write (iw,1400) (whts(i),i=1,nquad-1)
        write (iw,1450) whts(nquad)

c
        close(iw)
c
        end
c
c
c
        subroutine write_square90(norder,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
c
        character*12 text
c
c
 1000 format ("square",I2.2,'.90')
c
 1100 format ('nquad   = ',I3)
 1200 format ('ndegree = ',I3)
c
 1300 format ('data ',A4,' / ',32X,'&')
 1400 format ((2x,D38.30,',   &'))
 1450 format ((2x,D38.30,'    /'))


c
        write(text,1000) norder
c
        iw=100
        open(iw,FILE=text)
c
        write (iw,1100) nquad
        write (iw,1200) norder
c
        write (iw,1300) 'xs'
        write (iw,1400) (xs(i),i=1,nquad-1)
        write (iw,1450) xs(nquad)

c
        write (iw,1300) 'ys'
        write (iw,1400) (ys(i),i=1,nquad-1)
        write (iw,1450) ys(nquad)

c
        write (iw,1300) 'whts'
        write (iw,1400) (whts(i),i=1,nquad-1)
        write (iw,1450) whts(nquad)

c
        close(iw)
c
        end

c
c
c
        subroutine funtest(x,y,val,i,j)
        implicit double precision (a-h,o-z)
        val=x**i*y**j
        return
        end
c
c
c
        subroutine funeval(x,y,par1,par2,vals,dersx,dersy)
        implicit double precision (a-h,o-z)
        dimension vals(1),dersx(1),dersy(1)
c
        dimension xvals(1000),xders(1000)
        dimension yvals(1000),yders(1000)
c      
c       Evaluate the input functions and their derivatives at the point 
c       (x,y).
c
        
c
c       Fetch the values and derivatives of the Legendre polynomials.
c
        norder = par1
c
        call legeders(norder,x,xvals,xders)
        call legeders(norder,y,yvals,yders)
c
c
c       Exploit four-way symmetry.
c
c$$$        nn=1
c$$$        do 1000 j=0,norder,2
c$$$        do 1100 i=0,norder-j,2
c$$$        vals(nn)=xvals(j+1)*yvals(i+1)
c$$$        dersx(nn)=xders(j+1)*yvals(i+1)
c$$$        dersy(nn)=xvals(j+1)*yders(i+1)
c$$$        nn=nn+1
c$$$ 1100 continue
c$$$ 1000 continue
c$$$c
c$$$          par2=nn-1
c
c
c       The eight-fold way.
c
        nn=1
        do 1000 j=0,norder,2
        do 1100 i=j,norder-j,2
        vals(nn)=(xvals(j+1)*yvals(i+1)+xvals(i+1)*yvals(j+1))/2
        dersx(nn)=(xders(j+1)*yvals(i+1)+xders(i+1)*yvals(j+1))/2
        dersy(nn)=(xvals(j+1)*yders(i+1)+xvals(i+1)*yders(j+1))/2
        nn=nn+1
 1100 continue
 1000 continue
c
        par2=nn-1
c
        end
c
c
c
        subroutine funeval2(x,y,par1,par2,vals,dersx,dersy)
        implicit double precision (a-h,o-z)
        dimension vals(1),dersx(1),dersy(1)
c
        dimension xvals(1000),xders(1000)
        dimension yvals(1000),yders(1000)
c      
c       Evaluate the input functions and their derivatives at the point 
c       (x,y).
c
        
c
c       Fetch the values and derivatives of the Legendre polynomials.
c
        norder = par1
c
        call legeders(norder,x,xvals,xders)
        call legeders(norder,y,yvals,yders)
c
c
c       Exploit four-way symmetry.
c
        nn=1
        do 1000 j=0,norder,2
        do 1100 i=0,norder-j,2
        vals(nn)=xvals(j+1)*yvals(i+1)
        dersx(nn)=xders(j+1)*yvals(i+1)
        dersy(nn)=xvals(j+1)*yders(i+1)
        nn=nn+1
 1100 continue
 1000 continue
c
          par2=nn-1
        end

c
c
c
        subroutine funeval3(x,y,par1,par2,vals,dersx,dersy)
        implicit double precision (a-h,o-z)
        dimension vals(1),dersx(1),dersy(1)
c
        dimension xvals(1000),xders(1000)
        dimension yvals(1000),yders(1000),yvals2(1000),yders2(1000)
c      
c       Evaluate the input functions and their derivatives at the point 
c       (x,y).
c
        
c
c       Fetch the values and derivatives of the Legendre polynomials.
c
        norder = par1
c
        call legeders(norder,x,xvals,xders)
        call legeders(norder,y,yvals,yders)
c
c       Explot two-way symmetry.
c
        nn=1
        do 1000 i=0,norder,2
        do 1100 j=0,norder-i
        vals(nn)=xvals(j+1)*yvals(i+1)
        dersx(nn)=xders(j+1)*yvals(i+1)
        dersy(nn)=xvals(j+1)*yders(i+1)
        nn=nn+1
 1100 continue
 1000 continue
c
        par2=nn-1
        end

c
c
c
        subroutine funeval4(x,y,par1,par2,vals,dersx,dersy)
        implicit double precision (a-h,o-z)
        dimension vals(1),dersx(1),dersy(1)
c
        dimension xvals(1000),xders(1000)
        dimension yvals(1000),yders(1000)
c      
c       Evaluate the input functions and their derivatives at the point 
c       (x,y).
c
c
c       Fetch the values and derivatives of the Legendre polynomials.
c
        norder = par1
c
        call legeders(norder,x,xvals,xders)
        call legeders(norder,y,yvals,yders)
c
        nn=1
        do 1000 j=0,norder
        do 1100 i=0,norder-j
        vals(nn)=xvals(j+1)*yvals(i+1)
        dersx(nn)=xders(j+1)*yvals(i+1)
        dersy(nn)=xvals(j+1)*yders(i+1)
        nn=nn+1
 1100 continue
 1000 continue
c
        par2=nn-1
        end

c
c
c
c$$$        subroutine legeders(n,x,pols,ders)
c$$$        implicit double precision (a-h,o-z)
c$$$        dimension pols(1),ders(1)
c$$$c     
c$$$c       Evaluate the Legendre polynomials of degree 0 
c$$$c       through n and their derivatives at the point x.
c$$$c
c$$$        pols(1) = 1
c$$$        ders(1) = 0
c$$$        if (n .eq. 0) goto 2000
c$$$c
c$$$        pols(2) = x
c$$$        ders(2) = 1
c$$$        if (n .eq. 1) goto 2000 
c$$$c
c$$$        do 1000 j=2,n
c$$$        pols(j+1)=((2*j-1)*x*pols(j)-(j-1)*pols(j-1))/j
c$$$ 1000 continue
c$$$ 2000 continue
c$$$c
c$$$c       Compute the derivatives.
c$$$c
c$$$        d=x**2-1
c$$$        do 3000 j=3,n+1
c$$$        ders(j)=(j-1)*(x*pols(j)-pols(j-1))/d
c$$$ 3000 continue
c$$$        end
c$$$c
c
c
c$$$        subroutine lege(n,x,vals)
c$$$        implicit double precision (a-h,o-z)
c$$$        dimension vals(1)
c$$$c     
c$$$c       Evaluate the Legendre polynomials of degree 0 
c$$$c       through n at the point x.  
c$$$c
c$$$c       The polynomials are normalized so as to make their L^2[-1,1]
c$$$c       norms 1.
c$$$c
c$$$        vals(1) = 1
c$$$        if (n .eq. 0) goto 2000
c$$$c
c$$$        vals(2) = x
c$$$        if (n .eq. 1) goto 2000
c$$$c
c$$$        do 1000 j=2,n
c$$$        vals(j+1)=((2*j-1)*x*vals(j)-(j-1)*vals(j-1))/j
c$$$ 1000 continue
c$$$c
c$$$c       Normalize the polynomials.
c$$$c
c$$$ 2000 continue
c$$$        end
c$$$c
c
c
        subroutine tensorsquare(norder,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        dimension xslege(1000),whtslege(1000)
c
c       Return a tensor product quadrature for integrating polynomials
c       of order norder over the triangle with vertices (0,0), (1,0), 
c       and (1,1).
c
c        itype=1
        nlege=(norder+1)
c        call legeexps(itype,nlege,xslege,u,v,whtslege)
        call legequad(nlege,xslege,whtslege)
c
        nquad=0
        do 1000 j=1,nlege
        x=xslege(j)
        xwht=whtslege(j)
c
        do 1100 i=1,nlege
c
        y = xslege(i)
        ywht = whtslege(i)
c
        nquad=nquad+1
        xs(nquad)=x
        ys(nquad)=y
        whts(nquad)=xwht*ywht
 1100 continue
 1000 continue
c        
        end
c
c
c
        subroutine tensorsquare2(norder,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        dimension xslege(1000),whtslege(1000)
c
c       Return a tensor product quadrature for integrating polynomials
c
c        itype=1
        nlege=(norder+1)
        call legequad(nlege,xslege,whtslege)
c        call legeexps(itype,nlege,xslege,u,v,whtslege)
c
        nquad=0
        do 1000 j=1,nlege

        alpha=0.5d0
        beta=0.5d0
c
        x=xslege(j)*alpha+beta
        xwht=whtslege(j)*alpha        
c
        alpha = 0.5d0
        beta  = 0.5d0
c
        do 1100 i=1,nlege
c
        y = xslege(i)*alpha+beta
        ywht = whtslege(i)*alpha
c
        nquad=nquad+1
        xs(nquad)=x
        ys(nquad)=y
        whts(nquad)=xwht*ywht
 1100 continue
 1000 continue
c        
        end
c
c
c
        subroutine tensortri(norder,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        dimension xslege(1000),whtslege(1000)
c
c       Return a tensor product quadrature for integrating polynomials
c       over the triangle with vertices (0,0), (1,0), (1,1)
c
c        itype=1
        nlege=(norder+1)
c        call legeexps(itype,nlege,xslege,u,v,whtslege)
        call legequad(nlege,xslege,whtslege)
c
        nquad=0
        do 1000 j=1,nlege

        alpha = 0.5d0
        beta  = 0.5d0
c
        x=xslege(j)*alpha+beta
        xwht=whtslege(j)*alpha        
c
        alpha = x/2
        beta  = x/2
c
        do 1100 i=1,nlege
c
        y = xslege(i)*alpha+beta
        ywht = whtslege(i)*alpha
c
        nquad=nquad+1
        xs(nquad)=x
        ys(nquad)=y
        whts(nquad)=xwht*ywht
 1100 continue
 1000 continue
c        
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine gaussquad(eps,nfuns,funeval,par1,par2,nquad,xs,ys,
     1    whts,iregion,ifcauchy,ifonlysignifs,ngoal)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        dimension vals0(10 000),dersx0(10 000),dersy0(10 000)
        dimension signifs(100 000),rints(100 000),idxes(10 000)
        dimension ilist(1)
        external funeval
c
c
c       Compute the integrals.
c
        do 1000 j=1,nfuns
        rints(j)=0
 1000 continue
c
        do 1100 j=1,nquad
        x=xs(j)
        y=ys(j)
        wht=whts(j)
c
        call funeval(x,y,par1,par2,vals0,dersx0,dersy0)
        do 1200 i=1,nfuns
        rints(i)=rints(i)+vals0(i)*wht
 1200 continue
 1100 continue
c
c       Remove point, one-by-one (if possible).
c
 2000 continue
c
        if (nquad .eq. ngoal) goto 3000
c
c       Compute signifigances for the points and reorder them.
c
        call compute_signifs(nfuns,funeval,par1,par2,
     1    nquad,xs,ys,whts,signifs,iregion,rints)
c
        if (ifonlysignifs .eq. 1) then
        if(signifs(1) .gt. 0) goto 3000
        endif
c
        call insort4(nquad,signifs,xs,ys,whts)
c
        call prin2("xs=*",xs,nquad)
        call prin2("ys=*",ys,nquad)
        call prin2("whts=*",whts,nquad)
        call prin2("signifs=*",signifs,nquad)
        call prinf("nquad=*",nquad,1)
c
c       Try to remove each point.
c        
        do 2100 ipt=1,nquad/2
c        do 2100 ipt=1,5
c
c        if (ifonlysignifs .eq. 1 .AND. signifs(ipt) .gt. 0) return
c
        if (ifcauchy .eq. 1) then
        call remove_point_cauchy(iresult,ipt,eps,nfuns,funeval,par1,
     1    par2,nquad,xs,ys,whts,rints)
c
        else
        call remove_point(iresult,ipt,eps,nfuns,funeval,par1,
     1    par2,nquad,xs,ys,whts,rints)
        endif
c
        if (iresult .eq. 0) then
        goto 2000
        endif
 2100 continue
c
c$$$        if (ifcauchy .eq. 1) then
c$$$c
c$$$        do 2200 ipt=1,nquad
c$$$c
c$$$        if (ifonlysignifs .AND. signifs(ipt) .gt. 0) return
c$$$c
c$$$        call remove_point_cauchy(iresult,ipt,eps,nfuns,funeval,par1,
c$$$     1    par2,nquad,xs,ys,whts,rints)
c$$$        if (iresult .eq. 0) goto 2000
c$$$ 2200 continue
c$$$
c$$$        endif
c
 3000 continue
c
c       Otherwise, we are done.
c
        end
c
c
c
        subroutine compute_signifs(nfuns,funeval,par1,par2,
     1    nquad,xs,ys,whts,signifs,iregion,rints)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),signifs(1),rints(1)
        dimension vals0(1000),dersx0(1000),dersy0(1000)
c
        dimension w2(10000)
c
        double precision, allocatable :: df(:,:), a(:,:),w(:)
        double precision, allocatable :: ainv(:,:),ainv2(:,:)
        double precision, allocatable :: b(:),u(:),v(:),df2(:,:)
        double precision, allocatable :: rhs(:)
c
        external funeval

c
        nlen = nquad*3
c
c       Construct the linearized system ...
c
        allocate ( df(nfuns,nlen), df2(nfuns,nlen), a(nfuns,nfuns) )
        allocate ( ainv(nfuns,nfuns), ainv2(nfuns,nfuns) )
        allocate ( w(nfuns**2+nfuns+5000) )
        allocate ( rhs(nfuns) )
c
        do 2200 i=1,nquad
c
        x=xs(i)
        y=ys(i)
        wht=whts(i)
c
        call funeval(x,y,par1,par2,vals0,dersx0,dersy0)
c
        do 2300 j=1,nfuns
        df(j,i)=dersx0(j)*wht
        df(j,i+nquad)=dersy0(j)*wht
        df(j,i+nquad)=vals0(j)
 2300 continue
 2200 continue
c
        
c     
c       ... and the normal equations.
c
        call gauss_matmul3(nfuns,nfuns,nlen,df,df,a)
c
c       Compute the inverse.
c
        ainv=a
        call orthom(ainv,nfuns,w,cond)
c
c       Compute the significance of each point.
c
        do 1000 i=1,nquad
c
        x=xs(i)
        y=ys(i)
        wht=whts(i)
c
c       Check that the point is inside the correct region; mark it
c       for removal if not.
c
        if (iregion .eq. 1) then
           
        if (x .lt. -1 .OR. x .gt. 1.0d0) then
        signifs(i)=-1.0d0
        goto 1000
        endif
c
        if (y .lt. 0 .OR. y .gt. 1) then
        signifs(i)=-1.0d0
        goto 1000
        endif
        endif
c
        if (iregion .eq. 2) then
           
        if (x .lt. 0 .OR. x .gt. 1.0d0) then
        signifs(i)=-1.0d0
        goto 1000
        endif
c
        if (y .lt. 0 .OR. y .gt. 1.0d0) then
        signifs(i)=-1.0d0
        goto 1000
        endif
        endif
c
        if (iregion .eq. 3) then
           
        if (x .lt. 0 .OR. x .gt. 1.0d0) then
        signifs(i)=-1.0d0
        goto 1000
        endif
c
        if (y .lt. 0 .OR. y .gt. x) then
        signifs(i)=-1.0d0
        goto 1000
        endif
        endif       
c
        if (iregion .eq. 0) then
        if (x .lt. -1 .OR. x .gt. 1.0d0) then
        signifs(i)=-1.0d0
        goto 1000
        endif
c
        if (y .lt. -1 .OR. y .gt. 1d0) then
        signifs(i)=-1.0d0
        goto 1000
        endif
        endif
c
        
        if (iregion .eq. 2) then
        if (x .lt. 0 .OR. x .gt. 1.0d0) then
        signifs(i)=-1.0d0
        goto 1000
        endif
c
        if (y .lt. 0 .OR. y .gt. 1d0) then
        signifs(i)=-1.0d0
        goto 1000
        endif
        endif
c
c       Check for negative weights.
c
        if (wht .lt. 0.0d0) then
        signifs(i)=0
        goto 1000
        endif
c
c       Construct the neutered normal equations by deleting 3 columns
c       of a^-1.
c
        i1 = i
        i2 = i+nquad
        i3 = i+2*nquad
c
        ainv2=ainv
        ifsign=1
        call gauss_smw(ainv2,nfuns,df(1,i1),df(1,i1),w,ifsign,ifsing)
        if (ifsing .ne. 0) then
        signifs(i)=5000
        goto 1000
        endif
c
        call gauss_smw(ainv2,nfuns,df(1,i2),df(1,i2),w,ifsign,ifsing)
        if (ifsing .ne. 0) then
        signifs(i)=5000
        goto 1000
        endif
c
        call gauss_smw(ainv2,nfuns,df(1,i3),df(1,i3),w,ifsign,ifsing)
        if (ifsing .ne. 0) then
        signifs(i)=5000
        goto 1000
        endif
c
        df2= df
        do 1100 j=1,nfuns
        df2(j,i1)=0
        df2(j,i2)=0
        df2(j,i3)=0
 1100 continue
c
c       Compute the rhs.
c
        do 1200 j=1,nfuns
        rhs(j)=rints(j)
 1200  continue
c
        do 1300 j=1,nquad
        xx = xs(j)
        yy = ys(j)
        whtt = whts(j)
c
        if ( j .eq. i) goto 1300
        call funeval(xx,yy,par1,par2,vals0,dersx0,dersy0)
        do 1400 jj=1,nfuns
        rhs(jj)=rhs(jj)-vals0(jj)*whtt
 1400 continue
c
 1300 continue
c
c       Apply the matrices to compute the solution.
c
        call gauss_apply(nfuns,nfuns,ainv2,rhs,w)
        call gauss_applyt(nfuns,nlen,df2,w,w2)
c
        sum = 0
        do 1500 j=1,nquad
        sum=sum+w2(j)**2
 1500 continue
c
        signifs(i)=sum
c
cc
c$$$        sum=0
c$$$c
c$$$        call funeval(x,y,par1,par2,vals0,dersx0,dersy0)
c$$$c
c$$$        do 1100 j=1,nfuns
c$$$        sum=sum+vals0(j)**2
c$$$ 1100 continue
c$$$        signifs(i)=sum
c
c
c
 1000 continue
c
c       Add some randomization.
c
c$$$        call corrand_norm(nfuns,w,w2)
c$$$        do 3000 j=1,nfuns
c$$$        if( signifs(j) .gt. 0) then
c$$$        signifs(j)=signifs(j)+abs(w2(j))/100
c$$$        endif
c$$$ 3000 continue
        end
c
c
c
c
        subroutine gauss_smw(ainv,n,u,v,work,ifsign,ifsing)
        implicit double precision (a-h,o-z)
c
c       This subroutine applies Sherman-Morrison-Woodbury formula
c       to obtain the rank-1 update of the inverse of the matrix a.
c
c                                          -1
c       More precisely, given the inverse  a  this subroutine 
c       will return either
c
c                         -1
c                (a + uv')  , ifsign=0, or      (1)
c                                       
c                         -1
c                (a - uv')  , ifsign=1,         (2)
c                                       
c     where u and v are two given n-dimensional vectors.
c
c     Please note that u and v are used in the rank-1 update of the 
c     original matrix a, and that all computations in this subroutine
c     are performed with the INVERSE of the matrix a. 
c
c     This is a memory management routine for the routine
c     quaemu0, which performs the actual work.
c
c            Input paramenters:
c
c  a - the matrix dimensioned a(n,n) that contains the inverse of a
c  u,v - the vectors dimensioned u(n) and v(n) 
c
c  ifsign - the parameter that chooses the sign in 
c           Sherman-Morrison-Woodbury formula:
c
c     ifsign=0
c
c              -1    -1    -1              -1   -1     -1
c     (a + uv')   = a   - a   u  (1 + v' a   u)    v' a       (3)
c
c
c     ifsign=1
c
c              -1    -1    -1              -1   -1     -1
c     (a - uv')   = a   + a   u  (1 - v' a   u)    v' a       (4)
c
c            Output parameters:
c
c  a - the rank-1 update of the inverse of a.
c
c            Work arrays:
c
c  work - must be at least 2*n double precision elements long.
c
        double precision ainv(n,n),u(n),v(n),  work(1)
c
        ix=1
        lx=n
c
        iy=ix+lx
        ly=n
c
        call gauss_smw0(ainv,n,u,v,work(ix),work(iy),ifsign,ifsing)
        return
        end
c
c
c
        subroutine gauss_smw0(ainv,n,u,v,x,y,ifsign,ifsing)
        implicit double precision (a-h,o-z)
c
        double precision ainv(n,n),u(n),v(n),x(n),y(n)
c
        call gauss_apply(n,n,ainv,u,x)
        call gauss_applyt(n,n,ainv,v,y)
c
        alpha=0
        do 1200 i=1,n
        alpha=alpha+v(i)*x(i)
 1200   continue
c
        ifsing=0
c
        if (ifsign .eq. 0 .AND. 1+alpha .eq. 0) then
           ifsing=1
           return
        endif
c
        if (ifsign .eq. 1 .AND. 1-alpha .eq. 0) then
           ifsing=1
           return
        endif
c
        if( ifsign.eq.0 ) alpha=+1/(1+alpha)
        if( ifsign.eq.1 ) alpha=-1/(1-alpha)
c
        do 1600 j=1,n
        do 1400 i=1,n
        ainv(j,i)=ainv(j,i)-alpha*x(j)*y(i)
 1400   continue
 1600   continue
c
        return
        end
c
c
c
        subroutine remove_point_cauchy(iresult,ipt,eps,nfuns,funeval,
     1    par1,par2,nquad,xs,ys,whts,rints)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),rints(1)
        dimension vals0(10 000),dersx0(10 000),dersy0(10 000)
        double precision, allocatable :: z0(:)
c        
        external funeval
c
        iresult = 0
c
c       Set algorithm parameters.
c
        maxiters = 20
        maxsteps = 20
c
        ncauchy   = 50
        maxcsteps = 20
c
        dfact    = 0.25d0
c
        nlen = (nquad)*3-3
c
c       Allocate memory for the Gauss-Newton iterations.
c      
        allocate ( z0(nlen) )
c
c       Remove the specified point from the quadrature.
c
        nn=0
        do 1000 j=1,nquad
        if (ipt .eq. j) goto 1000
        nn=nn+1
        z0(nn)=xs(j)
        z0(nn+nquad-1)=ys(j)
        z0(nn+2*nquad-2)=whts(j)
 1000 continue
c
c       Conduct Newton/Cauchy iterations.
c
 1900 continue
        do 2000 iter=1,maxiters       
c        
c
c       Try one Newton step.
c
        call newton_step(iresult,eps,maxsteps,dfact,nquad,z0,
     1    nfuns,funeval,par1,par2,rints,errl2)
c
c        call prin2("Newton errl2=*",errl2,1)
c
        if (errl2 .le. eps**2) goto 3000
c
c
c       Do Cauchy steps if it failed ...
c
        if (iresult .ne. 0 ) then 
        do 2200 j=1,ncauchy
        call cauchy_step(iresult,eps,maxcsteps,dfact,nquad,z0,
     1    nfuns,funeval,par1,par2,rints,errl2)
c
c        call prin2("Cauchy errl2=*",errl2,1)
        if (iresult .ne. 0) return
        if (errl2 .le. eps**2) goto 3000        
 2200 continue
c
        endif
c
 2000 continue
c
        iresult=-1
        goto 9000
 3000 continue
c
c       Copy out the quadrature.
c
        do 5000 j=1,nquad-1
        xs(j)=z0(j)
        ys(j)=z0(j+nquad-1)
        whts(j)=z0(j+2*nquad-2)
 5000 continue
        nquad=nquad-1
c
 9000 continue
        deallocate(z0)
        end

c
c
c
        subroutine remove_point(iresult,ipt,eps,nfuns,funeval,par1,par2,
     1    nquad,xs,ys,whts,rints)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),rints(1)
        dimension vals0(10 000),dersx0(10 000),dersy0(10 000)
        double precision, allocatable :: z0(:)
c        
        external funeval
c
        iresult = 0
c
c       Set algorithm parameters.
c
        maxiters = 16
        maxsteps = 20
        dfact    = 0.25d0
        nextra   = 0
c
        nlen = (nquad)*3-3
c
c       Allocate memory for the Gauss-Newton iterations.
c      
        allocate ( z0(nlen) )
c
c       Remove the specified point from the quadrature.
c
        nn=0
        do 1000 j=1,nquad
        if (ipt .eq. j) goto 1000
        nn=nn+1
        z0(nn)=xs(j)
        z0(nn+nquad-1)=ys(j)
        z0(nn+2*nquad-2)=whts(j)
 1000 continue
c
c       Conduct Newton/Cauchy iterations.
c
 1900 continue
        do 2000 iter=1,maxiters       
        call newton_step(iresult,eps,maxsteps,dfact,nquad,z0,
     1    nfuns,funeval,par1,par2,rints,errl2)
c
        if (iresult .ne. 0)  return
        if (errl2 .le. eps**2) goto 3000
 2000 continue
 2100 continue
        iresult=-1
        goto 9000
 3000 continue
c
c       Perform a few extra newton steps to see if we can refine
c       a bit more.
c
 3200 continue
        errold=errl2
        call newton_step(jresult,eps,maxsteps,dfact,nquad,z0,
     1    nfuns,funeval,par1,par2,rints,errl2)
        if (jresult .eq. 0) goto 3200
c
 3300 continue
c
c       Copy out the quadrature.
c
        do 5000 j=1,nquad-1
        xs(j)=z0(j)
        ys(j)=z0(j+nquad-1)
        whts(j)=z0(j+2*nquad-2)
 5000 continue
        nquad=nquad-1
c
 9000 continue
        deallocate(z0)
        end
c
c
c
        subroutine newton_step(iresult,eps,maxsteps,dfact,
     1    nquad,z0,nfuns,funeval,par1,par2,rints,errnew)
        implicit double precision (a-h,o-z)
        dimension z0(nquad*3-3),rints(1)
        dimension vals0(1000),dersx0(1000),dersy0(1000),rnorms(10 000)
c
        double precision, allocatable :: df(:,:),rhs(:),z1(:)
        double precision, allocatable :: w(:)
c
        double precision, allocatable :: delta_cauchy(:),delta_newton(:)
        external funeval
c
c       Conduct a single Gauss-Newton iteration.
c
        iresult = 0
c
        nlen = nquad*3-3
c
c       Allocate memory for the operations.
c
        allocate (df(nfuns,nlen), rhs(nfuns), z1(nlen) )
        allocate (w(nfuns*nlen*32+5000))
c
        allocate (delta_cauchy(nlen), delta_newton(nlen))
c
c       Form the linear system.
c
        do 1000 j=1,nfuns
        rhs(j)=rints(j)
 1000 continue
c
        do 1200 i=1,nquad-1
        x=z0(i)
        y=z0(i+nquad-1)
        wht=z0(i+2*nquad-2)
c
        call funeval(x,y,par1,par2,vals0,dersx0,dersy0)
c
        do 1300 j=1,nfuns
        rhs(j)=rhs(j)-vals0(j)*wht
        df(j,i)=dersx0(j)*wht
        df(j,i+nquad-1)=dersy0(j)*wht
        df(j,i+2*nquad-2)=vals0(j)
 1300 continue
 1200 continue
c
        errold=0
        do 1400 j=1,nfuns
        errold=errold+rhs(j)**2
 1400 continue
c
c
c       Solve it in a least squares sense to find the Newton direction.
c
        call mach_zero(eps0)
        eps0=eps0*2
        k=1
        ifcheck=1
        call nrleamatll(df,nfuns,nlen,k,rhs,delta_newton,eps0,krank,
     1    rnorms,w,ifcheck,errl2,errmax)
c
c       Compute the Cauchy direction while we're at it.
c
c$$$        do 1500 i=1,nlen
c$$$        sum=0
c$$$        do 1600 j=1,nfuns
c$$$        sum=sum+2*rhs(j)*df(j,i)
c$$$ 1600 continue
c$$$        delta_cauchy(i)=sum
c$$$ 1500 continue
c
c       Perform Newton step-length control.
c        
        alpha=1.0d0
c
        do 2000 istep=1,maxsteps
        z1=z0+alpha*delta_newton
c
        do 2100 j=1,nfuns
        rhs(j)=rints(j)
 2100 continue
c
        do 2200 i=1,nquad-1
        x=z1(i)
        y=z1(i+nquad-1)
        wht=z1(i+2*nquad-2)
c
        call funeval(x,y,par1,par2,vals0,dersx0,dersy0)
c
        do 2300 j=1,nfuns
        rhs(j)=rhs(j)-vals0(j)*wht
 2300 continue
 2200 continue
c
        errnew = 0
        do 2400 j=1,nfuns
        errnew=errnew+rhs(j)**2
 2400 continue
c
c       We will only proceed with Newton steps if we are
c       really in the Newton regime.
c
        if (errnew .lt. errold) goto 3000
        alpha=alpha*dfact
c
 2000 continue
c     
        iresult = -1
        goto 9000
c
 3000 continue
c
c       Accept the new quadrature and return.
c
        z0=z1
c
 9000 continue
c
        deallocate(w,df,z1,delta_newton,delta_cauchy,rhs)
        end
c
c
c
        subroutine cauchy_step(iresult,eps,maxsteps,dfact,
     1    nquad,z0,nfuns,funeval,par1,par2,rints,errnew)
        implicit double precision (a-h,o-z)
        dimension z0(nquad*3-3),rints(1)
        dimension vals0(1000),dersx0(1000),dersy0(1000),rnorms(10 000)
c
        double precision, allocatable :: df(:,:),rhs(:),z1(:),delta(:)
        double precision, allocatable :: w(:)
        external funeval
c
c       Conduct a single Cauchy step.
c
        iresult = 0
c
        nlen = nquad*3-3
c
c       Allocate memory for the operations.
c
        allocate (df(nfuns,nlen), rhs(nfuns), z1(nlen), delta(nlen))
        allocate (w(nfuns*nlen*32+5000))
c
c       Form the Jacobian matrix and the "rhs"
c
        do 1000 j=1,nfuns
        rhs(j)=rints(j)
 1000 continue
c
        do 1200 i=1,nquad-1
        x=z0(i)
        y=z0(i+nquad-1)
        wht=z0(i+2*nquad-2)
c
        call funeval(x,y,par1,par2,vals0,dersx0,dersy0)
c
        do 1300 j=1,nfuns
        rhs(j)=rhs(j)-vals0(j)*wht
        df(j,i)=dersx0(j)*wht
        df(j,i+nquad-1)=dersy0(j)*wht
        df(j,i+2*nquad-2)=vals0(j)
 1300 continue
 1200 continue
c
        errold=0
        do 1400 j=1,nfuns
        errold=errold+rhs(j)**2
 1400 continue
c
c       The step direction is just - 2*rhs*df.
c
        do 1500 i=1,nlen
        sum=0
        do 1600 j=1,nfuns
        sum=sum+2*rhs(j)*df(j,i)
 1600 continue
        delta(i)=sum
 1500 continue
c
c       Perform step-length control.
c        
        alpha=1
c
        do 2000 istep=1,maxsteps
        z1=z0+alpha*delta
c
        do 2100 j=1,nfuns
        rhs(j)=rints(j)
 2100 continue
c
        do 2200 i=1,nquad-1
        x=z1(i)
        y=z1(i+nquad-1)
        wht=z1(i+2*nquad-2)
c
        call funeval(x,y,par1,par2,vals0,dersx0,dersy0)
c
        do 2300 j=1,nfuns
        rhs(j)=rhs(j)-vals0(j)*wht
 2300 continue
 2200 continue
c
        errnew = 0
        do 2400 j=1,nfuns
        errnew=errnew+rhs(j)**2
 2400 continue
c
        if (errnew .lt. errold/2) goto 3000
        alpha=alpha*dfact
c
 2000 continue
c     
        iresult = -1
        goto 9000
c
 3000 continue
c
c       Accept the new quadrature and return.
c
        z0=z1
c
 9000 continue
c
        deallocate(w,df,z1,delta,rhs)
        end

c
c
c
        subroutine chebquad(eps,nfuns,funeval,par1,par2,nquad0,xs0,ys0,
     1    whts0,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs0(1),ys0(1),whts0(1),xs(1),ys(1),whts(1)
        dimension vals(1000),dersx(1000),dersy(1000)
        dimension rnorms(10000),ipivs(10000),rints(10 000)
        double precision, allocatable :: a(:,:),b(:,:),c(:,:),q(:,:)
        double precision, allocatable :: r(:,:)
        double precision, allocatable :: rhs(:)
        external funeval
c
c       Allocate memory for the procedure.
c
        allocate ( a(nfuns,nquad0), b(nfuns,nquad0) )
c        allocate ( rints(nfuns) )
c
c       Form the matrix of function values at quadrature nodes scaled 
c       by weights and compute the integrals of the functions.
c
        do 0900 j=1,nfuns
        rints(j)=0
 0900 continue
c
c        rints(1)=1
c
        do 1000 j=1,nquad0
        x=xs0(j)
        y=ys0(j)
        wht=whts0(j)
        d = sqrt(wht)     
        call funeval(x,y,par1,par2,vals,dersx,dersy)
        do 1100 i=1,nfuns
        a(i,j)=vals(i)*d
        rints(i)=rints(i)+vals(i)*wht
 1100 continue
 1000 continue
c
        call prin2("in chebquad, rints=*",rints,nfuns)
c       
c       Form the QR decomposition of a.
c        
        call move(nfuns*nquad0,a,b)
        call gspiv(b,nfuns,nquad0,eps,rnorms,ipivs,krank)
        call prinf("krank=*",krank,1)
        call prin2("rnorms=*",rnorms,krank)
c
c       Copy out the column skeleton into q.
c
        allocate( q(nfuns,krank), c(nfuns,krank), r(krank,krank) )
        allocate( rhs(krank) )
c
        do 2000 j=1,krank
        do 2100 i=1,nfuns
        q(i,j)=b(i,j)
        c(i,j)=a(i,ipivs(j))
 2100 continue
 2000 continue
c
c       Compute r = q^*c.
c
        do 2200 i=1,krank
        do 2300 j=1,krank
        sum=0
        do 2400 l=1,nfuns
        sum=sum+q(l,i)*c(l,j)
 2400 continue
        r(i,j)=sum
 2300 continue
 2200 continue
c
c       Compute y = b^*rints
c
        do 2500 i=1,krank
        sum=0
        do 2600 j=1,nfuns
        sum=sum+b(j,i)*rints(j)
 2600 continue
        rhs(i)=sum
 2500 continue
c        
c       Solve Rx=y
c     
        call qrsolv(r,krank,rhs,rcond)
        call prin2("rhs=*",rhs,krank)
c
c       Assemble the quadrature,
c
        sum=0
        nquad=krank
        do 3000 j=1,nquad
        jj=ipivs(j)
        xs(j) = xs0(jj)
        ys(j) = ys0(jj)
        whts(j) = rhs(j)*sqrt(whts0(jj))
        sum=sum+whts(j)
 3000 continue
c
        deallocate (a,b,c,q,r,rhs)
        end

c
c
c
        subroutine matmult(n,m,k,q,a,r)
        implicit double precision (a-h,o-z)
        dimension q(m,n),a(m,k),r(n,k)
        do 1000 i=1,n
        do 1100 j=1,k
        sum=0
        do 1200 l=1,m
        sum=sum+q(l,i)*a(l,j)
 1200 continue
        r(i,j)=sum
 1100 continue
 1000 continue
        end

c
c
c
c$$$        subroutine move(k,a,b)
c$$$        implicit double precision (a-h,o-z)
c$$$        dimension a(1),b(1)
c$$$        do 1000 j=1,k
c$$$        b(j)=a(j)
c$$$ 1000 continue
c$$$        end
c
c
c
        subroutine insort4(k,a,b,c,d)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1),c(1),d(1)
c
        do 1000 i=2,k
        val=a(i)
        val2=b(i)
        val3=c(i)
        val4=d(i)
        j=i-1
        do 1100 while (j .ge. 1 .AND. a(j) .gt. val) 
        a(j+1)=a(j)
        b(j+1)=b(j)
        c(j+1)=c(j)
        d(j+1)=d(j)
        j=j-1
 1100 continue
        a(j+1)=val
        b(j+1)=val2
        c(j+1)=val3
        d(j+1)=val4
 1000 continue
        end
c
c
c
        subroutine gauss_apply(n,m,a,x,y)
        implicit double precision (a-h,o-z)
        dimension a(n,m),x(m),y(n)
c
        do 1000 i=1,n
        sum=0
        do 1100 j=1,m
        sum=sum+a(i,j)*x(j)
 1100 continue
        y(i)=sum
 1000 continue
        end
c
c
c
        subroutine gauss_applyt(n,m,u,x,y)
        implicit double precision (a-h,o-z)
        dimension u(n,m),x(n),y(m)
c
        do 1000 i=1,m
        sum=0
        do 1100 j=1,n
        sum=sum+u(j,i)*x(j)
 1100 continue
        y(i)=sum
 1000 continue
        end
c
c
c
        subroutine gauss_matmul2(n,m,k,a,b,c)
        implicit double precision (a-h,o-z)
        dimension c(n,m),a(k,n),b(k,m)
c
c       c = a^*b
c
        do 1000 i=1,n
        do 1100 j=1,m
        sum=0
        do 1200 l=1,k
        sum=sum+a(l,i)*b(l,j)
 1200 continue
        c(i,j)=sum
 1100 continue
 1000 continue
        end
c
c
c
        subroutine gauss_matmul3(n,m,k,a,b,c)
        implicit double precision (a-h,o-z)
        dimension c(n,m),a(n,k),b(m,k)
c
        do 1000 i=1,n
        do 1100 j=1,m
        sum=0
        do 1200 l=1,k
        sum=sum+a(i,l)*b(j,l)
 1200 continue
        c(i,j)=sum
 1100 continue
 1000 continue
        end
