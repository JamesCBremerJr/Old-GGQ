        implicit double precision (a-h,o-z)
        dimension pols(1000),ders(1000),dersx(1000),dersy(1000)
        dimension xs(10 000),ys(10 000),whts(10 000)
        dimension u(1000 000),v(1000000),vals(1000),coefs(1000)
        dimension r(100 00)
        double precision, allocatable :: amatr(:,:),dips(:,:)
        external funwrap1,funwrap2,funuser
c
        call prini(6,13)
c
c       Test the normalization of the Legendre polynomials using
c       adaptive Gaussian quadrature.
c
        n = 6
c
        errl2=0
        do 1000 i=0,n
        do 1100 j=0,n
c
        m   = 20
        call mach_zero(eps)
        eps=eps*10
c
        a=-1
        b=1
        call adapgauss(ier,a,b,eps,m,rint,funwrap1,i,j,
     1      par3,par4,par5,par6,par7,par8)
        d = rint**2
        if (i .eq. j) d = (1-rint)**2
        errl2=errl2+d
 1100 continue
 1000 continue
c     
        errl2=sqrt(errl2)
        call prin2("legendre ip errl2 = *",errl2,1)
c
c       Test the normalization of the Jacobi polynomials using
c       adaptive Gaussian quadrature.
c
        n = 6
        errl2=0
        do 1200 i=0,n
        do 1300 j=0,n
c
        m   = 20
        call mach_zero(eps)
        eps=eps*10
c
        a=-1
        b=1
        aa=3.0d0
        call adapgauss(ier,a,b,eps,m,rint,funwrap2,i,j,
     1      aa,par4,par5,par6,par7,par8)
        d = rint**2
c
        if (i .eq. j) d = (1-rint)**2
        errl2=errl2+d
 1300 continue
 1200 continue
c     
        errl2=sqrt(errl2)
        call prin2("jacobi ip errl2 = *",errl2,1)
c
        x=.1d0
        y=.1d0
        n=6
c
        call clotatim(t1)
        call koorn(n,x,y,vals)
        call clotatim(t2)
        call prin2("koorn time = *",t2-t1,1)
c$$$c
c$$$c       Test the inner products of the koornwinder polynomials
c$$$c       using a tensor product quadrature.
c$$$c
c$$$        n=6
c$$$        npols=(n+1)*(n+2)/2
c$$$c
c$$$        allocate ( dips(npols,npols) )
c$$$        dips=0
c$$$c
c$$$        norder=n*2
c$$$        call tensorsimp(norder,nquad,xs,ys,whts)
c$$$c
c$$$        do 2000 i=1,nquad
c$$$        x=xs(i)
c$$$        y=ys(i)
c$$$        wht=whts(i)
c$$$        call koorn2(n,x,y,pols,dersx,dersy)
c$$$c
c$$$        do 2100 j1=1,npols
c$$$        do 2200 j2=1,npols
c$$$        dips(j1,j2)=dips(j1,j2)+pols(j1)*pols(j2)*wht
c$$$ 2200 continue
c$$$ 2100 continue
c$$$ 2000 continue
c$$$c
c$$$        errl2 = 0
c$$$c
c$$$        do 2300 j=1,npols
c$$$        do 2400 i=1,npols
c$$$        d=dips(i,j)
c$$$        if (i .eq. j) d=dips(i,j)-1
c$$$        errl2=errl2+d**2
c$$$ 2400 continue
c$$$ 2300 continue
c$$$c
c$$$        errl2=sqrt(errl2)
c$$$        call prin2("koorn inner product errl2 = *",errl2,1)
c
c       Test the koornexps routine.
c
        itype  = 1
        norder = 12
c
        call clotatim(t1)
        call koornexps(itype,norder,npols,nquad,xs,ys,whts,u,v)
        call clotatim(t2)
c
        call prin2("koornexps time = *",t2-t1,1)
c     
        call prinf("norder = *",norder,1)
        call prinf("npols  = *",npols,1)
        call prinf("nquad  = *",nquad,1)
c
        allocate ( dips(npols,npols) )
        dips=0
c
        do 2000 i=1,nquad
        x=xs(i)
        y=ys(i)
        wht=whts(i)
        call koorn2(norder,x,y,pols,dersx,dersy)
c
        do 2100 j1=1,npols
        do 2200 j2=1,npols
        dips(j1,j2)=dips(j1,j2)+pols(j1)*pols(j2)*wht
 2200 continue
 2100 continue
 2000 continue
c
        errl2 = 0
c
        do 2300 j=1,npols
        do 2400 i=1,npols
        d=dips(i,j)
        if (i .eq. j) d=dips(i,j)-1
        errl2=errl2+d**2
 2400 continue
 2300 continue
c
        errl2=sqrt(errl2)
        call prin2("koorn inner product errl2 = *",errl2,1)
c
c       Test the products of u and v.
c
        call testprod(npols,nquad,u,v)
        call testprod2(npols,nquad,u,v)
c
c       Test evaluation.
c
        do 3000 i=1,nquad
        x = xs(i)
        y = ys(i)
        call funuser(x,y,val,derx,dery)
        vals(i)=val
 3000 continue
        call apply(npols,nquad,u,vals,coefs)
c        call prin2("coefs=*",coefs,npols)
c
        x=.12323d0
        y=.1d0
c
        call koorn2(norder,x,y,pols,dersx,dersy)
c
        sum0=0
        sum1=0
        sum2=0
        do 4100 i=1,npols        
        sum0=sum0+coefs(i)*pols(i)
        sum1=sum1+coefs(i)*dersx(i)
        sum2=sum2+coefs(i)*dersy(i)
 4100 continue
c
        call funuser(x,y,val,derx,dery)
c
        call prin2("eval err 1 = *",abs(val-sum0),1)
        call prin2("eval err 2 = *",abs(derx-sum1),1)
        call prin2("eval err 3 = *",abs(dery-sum2),1)
c
        sum1=0
        do 1800 i=1,nquad
        x=xs(i)
        y=ys(i)
        wht=whts(i)
        call funuser(x,y,val,derx,dery)
        sum1=sum1+val**2*wht
 1800 continue
c
        sum2=0
        do 1900 i=1,npols
        sum2=sum2+coefs(i)**2
 1900 continue

        print *,sum1
        print *,sum2
c
        stop

        end
c
c
c     
        subroutine funwrap1(x,val,i,j,par3,par4,par5,par6,par7,
     1    par8)
        implicit double precision (a-h,o-z)
        dimension vals0(1000)
        call lege(max(i,j),x,vals0)
        val=vals0(i+1)*vals0(j+1)
        end
c
c
c
        subroutine funwrap2(x,val,i,j,a,par5,par6,par7,par8)
        implicit double precision (a-h,o-z)
        dimension vals0(10 000)
        call jacobi(max(i,j),a,x,vals0)
        val=vals0(j+1)*vals0(i+1)*(1-x)**a
        end

c
c
c
        subroutine funuser(x,y,val,derx,dery)
        implicit double precision (a-h,o-z)
        val = cos(x*y/2)
        derx = -sin(x*y/2)*y/2
        dery = -sin(x*y/2)*x/2
c
c$$$        a=1.0d0/32.0d0
c$$$c
c$$$        x1 = 0
c$$$        y1 = 1-a
c$$$c
c$$$        x2 = a
c$$$        y2 = 1-a
c$$$c
c$$$        x3 = 0
c$$$        y3 = 1
c$$$c
c$$$        a=1.0d0/256.0d0
c$$$c
c$$$        x1 = 1-a
c$$$        y1 = 0
c$$$c
c$$$        x2 = 1
c$$$        y2 = 0
c$$$c
c$$$        x3 = 1-a
c$$$        y3 = a
c$$$
c$$$        a11 = x2-x1
c$$$        a21 = y2-y1
c$$$c
c$$$        a12 = x3-x1
c$$$        a22 = y3-y1
c$$$c
c$$$        a13 = x1
c$$$        a23 = y1
c$$$c
c$$$        s = a11*x+a12*y+a13
c$$$        t = a21*x+a22*y+a23
c$$$c
c$$$        dd = abs(a11*a22-a12*a21)
c$$$        val = s**20
c        
        end
c
c
c
        subroutine testprod(npols,nquad,u,v)
        implicit double precision (a-h,o-z)
        dimension u(npols,nquad),v(nquad,npols)
c
        errl2=0
c
        do 1000 i=1,npols
        do 1100 j=1,npols
        sum=0
        do 1200 l=1,nquad
        sum=sum+u(i,l)*v(l,j)
 1200 continue
c      
        if (i .eq. j) then
        errl2=errl2+(sum-1)**2
        else
        errl2=errl2+sum**2
        endif       
 1100 continue
 1000 continue
c
        errl2=sqrt(errl2)
        call prin2("reconstruction errl2=*",errl2,1)
c
        end

c
c
c
        subroutine testprod2(npols,nquad,u,v)
        implicit double precision (a-h,o-z)
        dimension u(npols,nquad),v(nquad,npols)
        dimension vu(nquad,nquad)
c
        errl2=0
c
c       Form the matrix vu.
c
        do 1000 i=1,nquad
        do 1100 j=1,nquad
        sum=0
        do 1200 l=1,nquad
        sum=sum+v(i,l)*u(l,j)
 1200 continue
        vu(i,j)=sum
 1100 continue
 1000 continue
c
c       Rsest vu*v=v
c
        do 2000 i=1,nquad
        do 2100 j=1,npols
        sum=0
        do 2200 l=1,nquad
        sum=sum+vu(i,l)*v(l,j)
 2200 continue
        errl2=errl2+(sum-v(i,j))**2
 2100 continue
 2000 continue
c
        errl2=sqrt(errl2)
        call prin2("projection errl2=*",errl2,1)
c
        end
c
c
c
        subroutine apply(n,m,a,x,y)
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the code
c       proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains routines for evaluating and manipulating
c       a family of orthonormal polynomials on the two-dimensional 
c       simplex.
c
c       The following subroutines are user-callable:
c
c   koornexps - return a quadrature integrating polynomials of degree
c       less than or equal to 2*norder on the simplex and, optionally, 
c       interpolation matrices for the orthonormal polynomials of 
c       orders 0 through norder.
c
c   koorn - use the Koornwinder representation formula to evaluate the 
c       polynomials, up to a given degree, at a user-specified point 
c       in the simplex.
c
c   koorn2 - use the Koornwinder representation formula to evaluate 
c       the polynomials and their derivatives, up to a given degree, 
c       at a user-specified point in the simplex.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine koornexps(itype,norder,npols,nquad,xs,ys,whts,u,v)
        implicit double precision (a-h,o-z)
c
        parameter ( n16 = 51 )
        dimension xs16(n16),ys16(n16),whts16(n16)
c
        dimension xs(1),ys(1),whts(1),u(1),v(1)
        double precision, allocatable :: eye(:,:),w(:),pols(:),rnorms(:)
c
c       Return a quadrature integrating polynomials of degree 2*norder
c       or less on the simplex and, optionally, interpolation matrices 
c       for  polynomials of degree norder or less on the simplex.
c
c       For certain orders (see below), relatively efficient quadrature
c       formula have been precomputed.  For all other orders, an 
c       inefficient tensor product quadrature is returned.  All
c       quadrature formulae achieve near extended precision accuracy.
c
c       Note that it is the user's responsibility to ensure that the 
c       output arrays are sufficiently large to store the interpolation
c       matrices and quadrature.
c
c                          Input Parameters:
c
c    itype - integer parameter indicating which data should be
c       returned
c
c       itype = 0   means return only quadrature nodes and weights
c       itype = 1   means also construct the matrices u and v
c
c    norder - maximum degree of the polynomials the interpolation 
c       matrices can handles
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c            EFFICIENT QUADRATURES ARE RETURNED FOR ORDERS 
c
c                         8, 12, 10, 16, 20, 30
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c                        Output Parameters:
c  
c   npols - the number of interpolation polynomials
c   nquad - the number of quadrature nodes
c
c   (xs,ys) the x and y coordinates (respectively) of the quadrature
c       nodes
c   whts - the quadrature weights
c
c   u - the (npols,nquad) matrix taking function values to expansion
c       coefficients
c
c   v - the (nquad,npols) matrix taking expansion coefficients to 
c       function values   
c
c
        data xs16/
     -  0.934912492774098837784593173685D+00,
     -  0.496739306154617956990889151711D-01,
     -  0.495159758378218511955073520689D+00,
     -  0.990920398056691395590712408004D+00,
     -  0.980878513784940398479326014232D+00,
     -  0.677795699867181989354429554820D+00,
     -  0.587797377146303664032572333197D+00,
     -  0.982393451950367923184212164375D+00,
     -  0.243191294947147700322444728576D+00,
     -  0.510728615419251912083527682798D+00,
     -  0.534689422065742559136578358140D-01,
     -  0.256185937817400616603851601497D+00,
     -  0.274812949740239493165854580845D+00,
     -  0.454853327691067887439797102742D-01,
     -  0.976478206830688270366401357796D+00,
     -  0.747328570048780024206505463644D+00,
     -  0.305383821076107964254724988928D+00,
     -  0.434203989387440203048785351081D+00,
     -  0.869490497166474845910124824863D+00,
     -  0.935490362913437232107266493781D+00,
     -  0.318111214204234816629475050404D+00,
     -  0.640044381548908745060025120284D-02,
     -  0.614091228795453000519474952363D-02,
     -  0.994348866971379968871641123963D+00,
     -  0.577785623898703274838818811248D+00,
     -  0.561279390133115743502727487204D+00,
     -  0.621329531813481894364963140167D+00,
     -  0.417925486172245941813919086057D+00,
     -  0.436951906847223595816595051967D+00,
     -  0.229920306404883072302375450542D+00,
     -  0.682464734470679440143672795181D+00,
     -  0.202053174513591873834682049890D+00,
     -  0.905385773744742115109131007031D+00,
     -  0.727109931423065515109661643373D+00,
     -  0.811807799737491246802721327002D+00,
     -  0.818639451318818609156980524662D+00,
     -  0.970265667033239599695399556339D-02,
     -  0.945962272674214251423062945647D-02,
     -  0.467915775044694616411163338742D-01,
     -  0.985836874736124696409675807304D+00,
     -  0.376678873559865994092099887664D+00,
     -  0.143306332956560402593956187491D+00,
     -  0.957562046532079730121995933460D+00,
     -  0.580237191008399499843031895761D-01,
     -  0.911290196282975909952216875015D+00,
     -  0.853864687541842614343341551932D+00,
     -  0.140542048023637709388088543720D+00,
     -  0.797843414398625726770594460191D+00,
     -  0.619653660927032687245699996368D+00,
     -  0.124735570248764527966789963276D+00,
     -  0.121925554690411428841332202848D+00/
        data ys16 /
     -  0.591214801642592453141429616216D-01,
     -  0.163937271227152396342982183400D+00,
     -  0.256436779251061029733267234448D+00,
     -  0.707225637315802878941696599280D-02,
     -  0.169040024079034803047481477803D-01,
     -  0.124448212747829247881588356184D+00,
     -  0.337167739110730024963892867581D-02,
     -  0.172098134020616461848080236157D-01,
     -  0.748807756591959804778531292563D+00,
     -  0.787615636502697092737852410231D-01,
     -  0.121106265264441645139865286738D-01,
     -  0.823254762612540145472067112515D-02,
     -  0.117737038632950831805817643418D+00,
     -  0.780591172161901932010722865196D+00,
     -  0.572982945380662216965281184250D-02,
     -  0.236104909670884683098981018729D+00,
     -  0.615730105199813250428710929293D+00,
     -  0.425622225987124508809701616769D+00,
     -  0.128638638744004195964300863140D+00,
     -  0.425702686485749419673482385005D-01,
     -  0.431980150852475671445765348776D+00,
     -  0.924261523987166348410220410395D+00,
     -  0.662747778464560821787847528639D-01,
     -  0.433260498259267136871024749484D-03,
     -  0.365059789413643267252382306876D+00,
     -  0.106351984472890941278412807045D+00,
     -  0.242309477071933262553290977542D+00,
     -  0.315964417232073222216004737727D-01,
     -  0.540682795623167000391201910820D+00,
     -  0.370093550991750677440026655766D+00,
     -  0.283396489106168583315792757081D-01,
     -  0.635617431316561342507840850982D+00,
     -  0.356247606772449221654040434592D-01,
     -  0.202545430027098156610974459909D+00,
     -  0.102718143532693416630125526217D+00,
     -  0.465224587980311425775296799842D-02,
     -  0.325106625080402033943567233230D+00,
     -  0.640197459132655357224868654060D+00,
     -  0.942592407649413564378417622021D+00,
     -  0.698555766814533450739152311902D-02,
     -  0.208557553100046586026544554980D+00,
     -  0.255846639134123815567377518779D+00,
     -  0.664241018769185781544428515012D-03,
     -  0.449947449420473230814193158554D+00,
     -  0.101301192382138132064531841603D-01,
     -  0.118802177759002510867510643620D+00,
     -  0.600297025744256173661557018521D-01,
     -  0.450975120464052299696180323397D-01,
     -  0.377101917394744336133006633858D+00,
     -  0.817161951131468085198511119032D+00,
     -  0.579741251599824248909518430282D+00/
        data whts16 / 
     -  0.732912073834814484181879581851D-03,
     -  0.156871969777336238957122585074D-01,
     -  0.224277769215099974796384106507D-01,
     -  0.555533600399760248732745391015D-04,
     -  0.489684396956794125199908128056D-04,
     -  0.142108669881620506410802760598D-01,
     -  0.331289966815217106846637793103D-02,
     -  0.522208132913360915287364104742D-04,
     -  0.578921902110157950248238209557D-02,
     -  0.106088097762134790879359658332D-01,
     -  0.415773696822125928782147678438D-02,
     -  0.590866267491745440644880391720D-02,
     -  0.245827085928422959447704583227D-01,
     -  0.155437834868859145913309782320D-01,
     -  0.289702552065438535800645759286D-03,
     -  0.524918040934347163938278226866D-02,
     -  0.163210843109416969595652618723D-01,
     -  0.203594789757592057361680935525D-01,
     -  0.888651428232809426443558257471D-03,
     -  0.149176978019345666164997966092D-02,
     -  0.226492671383077475724776195523D-01,
     -  0.341602278485694236183337265737D-02,
     -  0.341915065216762648540849041635D-02,
     -  0.180412036141501391781159056453D-04,
     -  0.135732686964030298393075955840D-01,
     -  0.119490918019100617021462415184D-01,
     -  0.128821311358524969652982025187D-01,
     -  0.130955606843725303653531788895D-01,
     -  0.102214235418538450926481659761D-01,
     -  0.285142704616270156423857860618D-01,
     -  0.829599800698114768147069862989D-02,
     -  0.193627804651181040595760203670D-01,
     -  0.263805221521955988906285658726D-02,
     -  0.763934663669938292110704273302D-02,
     -  0.716627298266479845173768792681D-02,
     -  0.210461990097886904841523785213D-02,
     -  0.774311596552266160167717350343D-02,
     -  0.813030331547507598852421722002D-02,
     -  0.347489826152848974210065687560D-02,
     -  0.151123616692256598374175421577D-03,
     -  0.273060187656058371381336028548D-01,
     -  0.278775467131095411459398988960D-01,
     -  0.176439696082773964155351079644D-03,
     -  0.216282126107872422883263076578D-01,
     -  0.161033784399686196331490867337D-02,
     -  0.353424094176319012979185253538D-02,
     -  0.152809598916178249915399509875D-01,
     -  0.682872244524248252636520137048D-02,
     -  0.282743200148425258518765782930D-02,
     -  0.142771965079444890296504095781D-01,
     -  0.244889698953820012117107312913D-01/
c
c       Copy out the correct quadrature, and call an auxillary
c       routine to construct the interpolation matrices (if the
c       user so desires).
c
        npols = (norder+1)*(norder+2)/2
c
c$$$        if (norder .eq. 8) then
c$$$        nquad = n16
c$$$        call koornmove(nquad,xs16,xs)
c$$$        call koornmove(nquad,ys16,ys)
c$$$        call koornmove(nquad,whts16,whts)
c$$$        goto 1000
c$$$        endif
c
        call koorntensor(norder,nquad,xs,ys,whts)
c
c       Construct the interpolation matrices if the user so desires.
c
 1000 continue
c
        if (itype .eq. 0) return
c
        call koorn0(norder,nquad,npols,xs,ys,whts,u,v)
c
        end
c
c
c
        subroutine koorntensor(norder,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        dimension xslege(1000),whtslege(1000)
c
c       Return a tensor product quadrature for the simplex integrating
c       polynomials of degree norder.
c
        nlege=(norder+1)
c        call legeexps(itype,nlege,xslege,u,v,whtslege)
        call legequad(nlege,xslege,whtslege)
c
        nquad=0
        do 1000 j=1,nlege
        x=(xslege(j)+1)/2
        xwht=whtslege(j)/2
c
        alpha = (1-x)/2
        beta  = (1-x)/2
c
        do 1100 i=1,nlege
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
        subroutine koorn0(norder,nquad,npols,xs,ys,whts,u,v)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),u(npols,nquad),v(nquad,npols)
        dimension eye(npols,npols)
        dimension rnorms(nquad+npols),pols(npols)
        dimension diag(npols),ipivs(nquad+npols)
c
        double precision, allocatable :: w(:)
c
c       Form the interpolation matrices.  
c
        do 1000 i=1,nquad
        x = xs(i)
        y = ys(i)
        wht = whts(i)
c
        call koorn(norder,x,y,pols)
        do 1100 j=1,npols
        v(i,j)=pols(j)
 1100 continue
 1000 continue
c
c       One can construct the coefficient evaluation matrix in the 
c       following fashion; however, it seems to be less accurate
c       then least squares.
c
c
c       Form the matrix u.
c
        do 1300 i=1,nquad
        do 1200 j=1,npols
        u(j,i)=v(i,j)*whts(i)
 1200 continue
 1300 continue
c
        return
c
c$$$c
c$$$c       Construct the identity.
c$$$c
c$$$        allocate(w(nquad*npols*16+5000))
c$$$c
c$$$        do 1200 i=1,npols
c$$$        do 1300 j=1,npols
c$$$        eye(i,j)=0
c$$$ 1300 continue
c$$$        eye(i,i)=1
c$$$ 1200 continue
c$$$
c$$$c
c$$$c       Construct the matrix taking values at quadrature nodes to
c$$$c       expansion coefficients via least squares.
c$$$c
c$$$        call mach_zero(eps)
c$$$        eps=eps*10
c$$$        ifcheck=1
c$$$c
c$$$        call nrleamatrr(v,npols,nquad,npols,eye,u,eps,krank,
c$$$     1    rnorms,w,ifcheck,errl2,errmax)
c$$$c
c$$$        if (ifcheck .eq. 1) then
c$$$        call prin2("in koornexps, errl2=*",errl2,1)
c$$$        endif
c$$$c
c$$$ 2000 continue
c$$$c
c$$$        deallocate(w)
c$$$c
        end
c
c
c
        subroutine koorn(n,x,y,pols)
        implicit double precision (a-h,o-z)
        dimension pols(1)
        dimension plege(n+1),pjac(n+1)
        data sq32 / 5.65685424949238019520675489683879231d0/
c
c       Evaluate the orthonormal polynomials of degree 0 through n
c       at the point (x,y) using the Koornwinder representation.  There 
c       are (n+1)(n+2)/2 such polynomials.
c
c       Note that this routine normalizes the polynomials so that their
c       L^2 norms over the simplex are 1.
c
c       The n+1 polynomials of degree n are given by 
c
c             (2k+1,0)         (0,0)
c         \phi (2x-1)  *   \phi (2y/(1-x)-1) * (1-x)^k 
c             n-k              k
c
c       k=0,...,n, where 
c
c             (a,b)
c         \phi (x)
c             k
c
c       denotes the Jacobi polynomial of degree k with parameters 
c       \alpha=a and \beta=b.       
c
c       The values of the polynomials are sorted first by degree and
c       then by the value of k.        
c
c       
c                          Input Parameters:
c
c   n - maximum order of the polynomials to evaluate
c   (x,y) - point at which to evaluate the polynomials
c
c                         Output Parameters:
c
c   vals - the (n+1)*(n+2)/2 values of the polynomials
c             
c     
c       First, evaluate the scaled Legendre polynomials.
c
        a=0
        b=0
        z = 2*y/(1-x)-1
        call lege(n,z,plege)
c
        alpha=sq32
        do 1000 j=0,n
        plege(j+1)=plege(j+1)*alpha
        alpha=alpha*(1-x)
 1000 continue     
c
        z = 2*x-1
c
        do 1100 k=0,n
c
c       Evaluate the (other) Jacobi polynomials
c
        a = 2*k+1
        call jacobi(n-k,a,z,pjac)
c
        do 1200 l=0,n-k
        m = l+k
        idx = ((m+1)*(m))/2+1+k
c
        pols(idx)=plege(k+1)*pjac(l+1)
c
 1200 continue
 1100 continue
c
        end
c
c
c
        subroutine koorn2(n,x,y,pols,dersx,dersy)
        implicit double precision (a-h,o-z)
        dimension pols(1),dersx(1),dersy(1)
        dimension plege(n+1),pjac(n+1),dlege(n+1),djac(n+1)
        data sq32 / 5.65685424949238019520675489683879231d0/
c
c       Evaluate the orthonormal polynomials of degree 0 through n
c       and their derivatives at the point (x,y) using the Koornwinder 
c       representation.  There are (n+1)(n+2)/2 such polynomials.
c
c       Note that this routine normalizes the polynomials so that their
c       L^2 norms over the simplex are 1.
c
c       The values of the polynomials are sorted first by degree and
c       then by the value of k.        
c
c                          Input Parameters:
c
c   n - maximum order of the polynomials to evaluate
c   (x,y) - point at which to evaluate the polynomials
c
c                         Output Parameters:
c
c   vals - the (n+1)*(n+2)/2 values of the polynomials at (x,y)
c   dersx - the derivatives of said polynomials w.r.t. x at the (x,y)
c   dersy - the derivatives of said polynomials w.r.t. y at the (x,y) 
c
c     
c       Evaluate the scaled Legendre polynomials.
c
        z1    = 2*y/(1-x)-1
        call legeders(n,z1,plege,dlege)
c
        z2 = 2*x-1
c
        nn=0
c
        do 1100 k=0,n
c
c       Evaluate the (other) Jacobi polynomials
c
        a = 2*k+1
        call jacobiders(n-k,a,z2,pjac,djac)
c
        d = sq32*(1-x)**k
c
        do 1200 l=0,n-k
c
        m = l+k
        idx = (m+1)*(m)/2+1+k
c
        pols(idx)=pjac(l+1)*plege(k+1)*d
c
        dd1=djac(l+1)*2*plege(k+1)*d
        dd2=pjac(l+1)*dlege(k+1)*2*y/(1-x)**2*d
c
        if (k .eq. 0) then
        dd3=0
        else
        dd3=-sq32*k*pjac(l+1)*plege(k+1)*(1-x)**(k-1)
        endif
c
        dersx(idx)=(dd1+dd2+dd3)
        dersy(idx)=pjac(l+1)*d*dlege(k+1)*2/(1-x)
c        
 1200 continue
 1100 continue
c
        end
c
c
c
c$$$        subroutine lege(n,x,vals)
c$$$        implicit double precision (a-h,o-z)
c$$$        dimension vals(1)
c$$$c     
c$$$c       Evaluate the normalized Legendre polynomials of degree 0 
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
c$$$c
c$$$        do 2100 j=0,n
c$$$        vals(j+1)=vals(j+1)*sqrt((2*j+1.0d0)/2.0d0)
c$$$ 2100 continue
c$$$        end
c$$$c
c$$$c
c$$$c
c$$$        subroutine legeders(n,x,pols,ders)
c$$$        implicit double precision (a-h,o-z)
c$$$        dimension pols(1),ders(1)
c$$$c     
c$$$c       Evaluate the normalized Legendre polynomials of degree 0 
c$$$c       through n and their derivatives at the point x.
c$$$c
c$$$c       The polynomials are normalized so as to make their L^2[-1,1]
c$$$c       norms 1.
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
c$$$c
c$$$c       Normalize the polynomials and derivatives.
c$$$c
c$$$        do 2100 j=0,n
c$$$        d=sqrt( (2.0d0*j+1.0d0)/2.0d0 ) 
c$$$        pols(j+1)=pols(j+1)*d
c$$$        ders(j+1)=ders(j+1)*d
c$$$ 2100 continue
c$$$c
c$$$        end
c
c
c
        subroutine jacobi(n,a,x,pols)
        implicit double precision (a-h,o-z)
        dimension pols(1)
c     
c       Evaluate the normalized Jacobi polynomials (with parameters
c       \alpha = a and \beta = 0) of degrees 0 through n at the point
c       x.
c
c       The polynomials are normalized so as to make their L^2 norms
c       with respect to the measure (1-x)^a dx equal to 1.
c
        pols(1) = 1
        if (n .eq. 0) goto 2000
c
        pols(2) = (a+(2+a)*x)/2
        if (n .eq. 1) goto 2000
c
        do 1000 l=2,n
c
        a1=2*(l)*(l+a)*(2*l-2+a)
        a2=(2*l-1+a)*a**2
        a3=(2*l-2+a)*(2*l+a-1)*(2*l+a)
        a4=2*(l+a-1)*(l-1)*(2*l+a)
c
        pols(l+1)=((a2+x*a3)*pols(l)-a4*pols(l-1))/a1
 1000 continue
c
c       Normalize the polynomials.
c
 2000 continue
c
        d = 2**(alpha+2)
        do 2100 j=0,n
        pols(j+1)=pols(j+1)*sqrt(2*j+a+1.0d0)/(d)
 2100 continue
        end
c
c
c
        subroutine jacobiders(n,a,x,pols,ders)
        implicit double precision (a-h,o-z)
        dimension pols(1),ders(1)
c     
c       Evaluate the normalized Jacobi polynomials (with parameters
c       \alpha = a and \beta = 0) of degrees 0 through n at the point
c       x and their derivatives.
c
c       The polynomials are normalized so as to make their L^2 norms
c       with respect to the measure (1-x)^a dx equal to 1.
c
        pols(1) = 1
        ders(1) = 0
        if (n .eq. 0) goto 2000
c
        pols(2) = (a+(2+a)*x)/2
        ders(2) = 1+a/2
        if (n .eq. 1) goto 2000
c
        do 1000 l=2,n
c
        a1=2*(l)*(l+a)*(2*l-2+a)
        a2=(2*l-1+a)*a**2
        a3=(2*l-2+a)*(2*l+a-1)*(2*l+a)
        a4=2*(l+a-1)*(l-1)*(2*l+a)
c
        pols(l+1)=((a2+x*a3)*pols(l)-a4*pols(l-1))/a1
        ders(l+1)=((a2+x*a3)*ders(l)+a3*pols(l)-a4*ders(l-1))/a1
 1000 continue
c
c       Normalize the polynomials.
c
 2000 continue
c
        d = 2**(alpha+2)
        do 2100 j=0,n
        pols(j+1)=pols(j+1)*sqrt(2*j+a+1.0d0)/(d)
        ders(j+1)=ders(j+1)*sqrt(2*j+a+1.0d0)/(d)
 2100 continue
        end
c
c
c
        subroutine koornmove(k,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
        do 1000 j=1,k
        b(j)=a(j)
 1000 continue
        end
c
c
c
