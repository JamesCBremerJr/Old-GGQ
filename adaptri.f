        implicit double precision (a-h,o-z)
        dimension verts(2,3),vals(1000),a(2,3)
        external funtest1
c
        verts(1,1) = .001d0
        verts(2,1) = 0
c
        verts(1,2) = 1.0d0
        verts(2,2) = 0
c
        verts(1,3) = 0
        verts(2,3) = 1.0d0
c
        eps=1.0d-14
c
        call adaptri(ier,eps,verts,funtest1,par1,par2,par3,par4,
     1    val,nquad)
c
        print *,"ier = ",ier
        print *,"val = ",val
        print *,"nquad = ",nquad
c         
        end
c        
c
c
        subroutine funtest1(x,y,par1,par2,par3,par4,val)
        implicit double precision (a-h,o-z)
        t   = atan2(y,x)
        val = cos(27*t)+x*y**2*sin(t**2)
        end
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the
c       adaptive integration code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       The following subroutines are user-callable:
c
c   adaptri - adaptively integrate a function supplied by the user via 
c       an external subroutine over a specified triangle in the plane 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c        
        subroutine adaptri(ier,eps,verts,funuser,par1,par2,par3,par4,
     1    val,nquad)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),stack(8,10 000)
        external funuser
c
c       Adaptively integrate a function supplied by the user via an 
c       external subroutine over a specified triangle.  
c
c       This subroutine uses a 79 point quadrature which integrates
c       20th order polynomials on the simplex.
c
c                           Input Parameters:
c
c   eps - precision for the routine
c   verts - a (2,3) array each column of which supplies the coordinates
c       of one vertex of the triangle which constitutes the integration
c       domain
c
c   funuser - user-supplied external subroutine with calling syntax
c
c       subroutine funuser(x,y,par1,par2,par3,par4,val)
c
c       Return the value of the function to integrate at the point
c       (x,y) in the parameter val.
c
c   par? - arbitrarily-typed user-supplied parameters
c
c                          Output Parameters:
c
c   ier - an error return code;
c       ier = 0    indicates success
c       ier = 4    means that stack was exhausted before convergence 
c                  was ovtained
c       ier = 16   means that the maximum recursion depth was exceeded
c                  before convergence was obtained
c
c   val - value of the integral
c   nquad - total number of quadrature nodes necessary to obtain the
c       specified accuracy
c
        ier   = 0
        val   = 0
        nquad = 0
c
c       Set algorithm parameters.
c
        maxstack = 10 000
        maxdepth = 50
c
c       Initialize the stack.
c
        istack=1
c
        x1 = verts(1,1)
        y1 = verts(2,1)
        x2 = verts(1,2)
        y2 = verts(2,2)
        x3 = verts(1,3)
        y3 = verts(2,3)
c
        stack(1,1)=x1
        stack(2,1)=y1
        stack(3,1)=x2
        stack(4,1)=y2
        stack(5,1)=x3
        stack(6,1)=y3
        stack(7,istack)=0
        call adaptri_oneint(stack(1,istack),funuser,par1,par2,par3,par4,
     1    val0)
        stack(8,istack)=val0
c
 1000 continue
c
c       Pop an entry off the stack.
c
        if (istack .eq. 0) goto 2000
c
        x1 = stack(1,istack)
        y1 = stack(2,istack)
        x2 = stack(3,istack)
        y2 = stack(4,istack)
        x3 = stack(5,istack)
        y3 = stack(6,istack)
c
        idepth0 = stack(7,istack)
        val0    = stack(8,istack)
        istack=istack-1
c
        if (idepth0 .gt. maxdepth) then
        ier=16
        return
        endif
c
c       Subdivide the triangle.
c
        u1 = (x1+x2)/2
        v1 = (y1+y2)/2
        u2 = (x1+x3)/2
        v2 = (y1+y3)/2
        u3 = (x2+x3)/2
        v3 = (y2+y3)/2
c
        if (istack+4 .gt. maxstack) then
        ier=4
        return
        endif
c
        istack=istack+1
        stack(1,istack)=x1
        stack(2,istack)=y1
        stack(3,istack)=u1
        stack(4,istack)=v1
        stack(5,istack)=u2
        stack(6,istack)=v2
        stack(7,istack)=idepth0+1
        call adaptri_oneint(stack(1,istack),funuser,par1,par2,par3,par4,
     1    stack(8,istack))
c
        istack=istack+1
        stack(1,istack)=u1
        stack(2,istack)=v1
        stack(3,istack)=u2
        stack(4,istack)=v2
        stack(5,istack)=u3
        stack(6,istack)=v3
        stack(7,istack)=idepth0+1
        call adaptri_oneint(stack(1,istack),funuser,par1,par2,par3,par4,
     1    stack(8,istack))
c
        istack=istack+1
        stack(1,istack)=u1
        stack(2,istack)=v1
        stack(3,istack)=x2
        stack(4,istack)=y2
        stack(5,istack)=u3
        stack(6,istack)=v3
        stack(7,istack)=idepth0+1
        call adaptri_oneint(stack(1,istack),funuser,par1,par2,par3,par4,
     1    stack(8,istack))
c
        istack=istack+1
        stack(1,istack)=u2
        stack(2,istack)=v2
        stack(3,istack)=x3
        stack(4,istack)=y3
        stack(5,istack)=u3
        stack(6,istack)=v3
        stack(7,istack)=idepth0+1
        call adaptri_oneint(stack(1,istack),funuser,par1,par2,par3,par4,
     1    stack(8,istack))
c
        val1 = stack(8,istack)+stack(8,istack-1)+stack(8,istack-2)
     1       + stack(8,istack-3)
        errmax = abs(val1-val0)
 1100 continue
c
c       Accept the parent triangle.
c
        if (errmax .lt. eps) then
        nquad=nquad+79
        val=val+val0
        istack=istack-4
        endif
c
        goto 1000
 2000 continue
        end
c
c
c
        subroutine adaptri_oneint(verts,funuser,par1,par2,par3,par4,sum)
        implicit double precision (a-h,o-z)
c
c       Use a 20th order quadrature on the simplex to integrate the
c       user-supplied function over an arbitrary triangle.
c
        dimension verts(2,3)
        dimension xs(79),ys(79),whts(79)
        data norder /20/
        data nquad  /79/
        data xs /
     -  0.2284735050323727D+00,  0.3319366893275869D+00,
     -  0.2189580470595975D+00,  0.3587736963800968D+00,
     -  0.2395917921515885D+00,  0.3154625530585929D+00,
     -  0.4304399028617899D+00,  0.1275013234218165D+00,
     -  0.4570781952828950D+00,  0.1306995358221298D+00,
     -  0.2030173650969737D+00,  0.3475513891833493D+00,
     -  0.1201791645340710D+00,  0.1345399267611344D+00,
     -  0.3266496842344425D+00,  0.2181075913972982D+00,
     -  0.5324655336566071D+00,  0.4487161705927870D+00,
     -  0.4869424255738578D+00,  0.2120053377439701D+00,
     -  0.1105168556849599D+00,  0.5636132794738498D+00,
     -  0.3033364930769153D+00,  0.5428812476307295D-01,
     -  0.4439932239684421D+00,  0.5206815897991972D-01,
     -  0.5500626889208072D-01,  0.1102961914665522D+00,
     -  0.6239962168981750D+00,  0.3427789398003033D+00,
     -  0.4792491564860371D+00,  0.1155775490615164D+00,
     -  0.1762340695132897D+00,  0.5615923313557501D-01,
     -  0.4803364557762367D-01,  0.6678738031878680D+00,
     -  0.5864034960976094D+00,  0.2168254815877556D+00,
     -  0.6152226766203740D+00,  0.4420803256703856D-01,
     -  0.1168238319185359D+00,  0.3795649922720747D+00,
     -  0.7172476505400961D+00,  0.7544487684646017D-01,
     -  0.2454730551315658D+00,  0.7411881570348956D+00,
     -  0.3621944569638384D-01,  0.4337468822937547D+00,
     -  0.5194720914008033D+00,  0.7670821152355030D+00,
     -  0.2911415011153069D+00,  0.5787450774439998D+00,
     -  0.4633826792163330D-01,  0.1052995941084088D-01,
     -  0.1023661287510408D-01,  0.1292022645374203D+00,
     -  0.1063124167323166D-01,  0.6532706672993116D+00,
     -  0.9554741128271498D-02,  0.1657540557574789D+00,
     -  0.1076925143722874D-01,  0.7151794939817122D+00,
     -  0.8708830955710343D-02,  0.8469498672451297D+00,
     -  0.7700149742155118D+00,  0.6144261677508253D-01,
     -  0.8446686787714371D+00,  0.1488526589464143D-01,
     -  0.7669259085114150D-02,  0.4267133090472286D-01,
     -  0.8335448303747325D+00,  0.2769893611210166D-02,
     -  0.1142531477080084D-01,  0.9180191438309294D+00,
     -  0.8726340631287909D+00,  0.9246137080048001D+00,
     -  0.9447081078109967D+00,  0.1049724260195281D-02,
     -  0.9775786888585845D+00 /
        data ys /
     -  0.3482440773114033D+00,  0.3601815055542639D+00,
     -  0.4830434612220117D+00,  0.2367691028542670D+00,
     -  0.2255504591091730D+00,  0.4906754487510426D+00,
     -  0.3573649814248187D+00,  0.4650433242232578D+00,
     -  0.2391923779685850D+00,  0.3291737077552919D+00,
     -  0.6134802400569791D+00,  0.1312913959604984D+00,
     -  0.6026128197700508D+00,  0.2082159661673957D+00,
     -  0.5709082200316685D+00,  0.1210596031731591D+00,
     -  0.2317350564581508D+00,  0.4341119377715239D+00,
     -  0.1341842478165952D+00,  0.6925469142640721D+00,
     -  0.7263287535250168D+00,  0.3030308406109172D+00,
     -  0.5298069853141057D-01,  0.4489144666660165D+00,
     -  0.5607317421545649D-01,  0.5895170796781314D+00,
     -  0.3139997659041382D+00,  0.1076681310216741D+00,
     -  0.1289547143339058D+00,  0.6168095844894770D+00,
     -  0.4719099926230024D+00,  0.8052875474677370D+00,
     -  0.4664205673493044D-01,  0.1953400526470680D+00,
     -  0.7235922802885636D+00,  0.1921029800646008D+00,
     -  0.5561939131122089D-01,  0.7485084908739520D+00,
     -  0.3283126607270730D+00,  0.8373438528473550D+00,
     -  0.8574753747584495D+00,  0.1056076372483861D-01,
     -  0.5316187443512963D-01,  0.3836800751487691D-01,
     -  0.9585400117180553D-02,  0.2000798838561538D+00,
     -  0.9805645903372592D-01,  0.5577230078004333D+00,
     -  0.1086852850445867D-01,  0.9955269246844377D-01,
     -  0.7023594162138953D+00,  0.4109163850135503D+00,
     -  0.9088856304739723D+00,  0.4339603746564995D+00,
     -  0.5722379991457208D+00,  0.7799647736722636D-02,
     -  0.3023710722654155D+00,  0.1056852154285495D-01,
     -  0.7073998298846009D+00,  0.8301213515986737D+00,
     -  0.1864461911253900D+00,  0.2735342012927543D+00,
     -  0.8287957675869461D+00,  0.9809212602237560D-01,
     -  0.1037153109571940D-01,  0.9325990991612949D+00,
     -  0.3188399213144015D-01,  0.3267404569561113D-01,
     -  0.9246900555593677D+00,  0.5064720872225761D-02,
     -  0.1553532048707776D+00,  0.9286275675075045D-01,
     -  0.9760736456539358D+00,  0.3098669406905532D-01,
     -  0.1516380956493777D-02,  0.6535567038319699D-01,
     -  0.1465415545837219D-02,  0.3640678784650113D-02,
     -  0.1236596945803616D-01 /
        data whts /
     -  0.1373978854437581D-01,  0.1318565866611611D-01,
     -  0.1318323939820397D-01,  0.1378522345886711D-01,
     -  0.1321528996665760D-01,  0.1233783797585319D-01,
     -  0.1173346172189761D-01,  0.1187133011279733D-01,
     -  0.8068797231635116D-02,  0.1164086989568030D-01,
     -  0.1059680034299219D-01,  0.1255421538693777D-01,
     -  0.1048605474861927D-01,  0.1032112216557707D-01,
     -  0.9454530634317015D-02,  0.1053956820449348D-01,
     -  0.1025687013737927D-01,  0.1040996837848747D-01,
     -  0.1292089488682078D-01,  0.8108082597442317D-02,
     -  0.7930755418664014D-02,  0.1048638799851170D-01,
     -  0.8104677939106186D-02,  0.8372987548349911D-02,
     -  0.8968664955997175D-02,  0.7933385739665351D-02,
     -  0.7836517205856028D-02,  0.7597548372004612D-02,
     -  0.1116669174126586D-01,  0.6411466815585067D-02,
     -  0.7735500151664028D-02,  0.5971248755228145D-02,
     -  0.6253125427509120D-02,  0.6764809073387011D-02,
     -  0.6576914398808417D-02,  0.9739934021442645D-02,
     -  0.8540178143756204D-02,  0.5234431998098076D-02,
     -  0.8501495096569582D-02,  0.4795366452222215D-02,
     -  0.3076347135454521D-02,  0.3759825385627354D-02,
     -  0.6945734115124631D-02,  0.3982591840667456D-02,
     -  0.3130122958599739D-02,  0.7625267557360685D-02,
     -  0.4329311350691444D-02,  0.3212906890781465D-02,
     -  0.3855729339422906D-02,  0.7760107791997257D-02,
     -  0.2363315012102513D-02,  0.3782543782314992D-02,
     -  0.3089541649371715D-02,  0.3672815378207577D-02,
     -  0.3613111743952638D-02,  0.2104349423572366D-02,
     -  0.3401183001531508D-02,  0.3424090784971568D-02,
     -  0.3165363267374626D-02,  0.1407804061241869D-02,
     -  0.2907034095167748D-02,  0.3728102829076812D-02,
     -  0.2451186630419030D-02,  0.5289912824386602D-02,
     -  0.2792093351236570D-02,  0.1207912454628245D-02,
     -  0.4402828906502393D-02,  0.1717038416085675D-02,
     -  0.1571283089345151D-02,  0.9822668220272660D-03,
     -  0.2998488249880729D-02,  0.9573988144602159D-03,
     -  0.9393304982541808D-03,  0.2764585208253715D-02,
     -  0.8218146540763565D-03,  0.1878112492476335D-02,
     -  0.5827588575545490D-03,  0.1545398110669673D-03,
     -  0.8195578158890342D-03 /
        external funuser
c
c       Construct an affine mapping take the standard simplex
c       to the user-specified triangle.
c
        b11 = verts(1,2)-verts(1,1)
        b12 = verts(1,3)-verts(1,1)
        b21 = verts(2,2)-verts(2,1)
        b22 = verts(2,3)-verts(2,1)
        b13 = verts(1,1)
        b23 = verts(2,1)
c
        det = abs(b11*b22-b12*b21)
c
        sum=0
        do 1100 i=1,nquad        
        x=b11*xs(i)+b12*ys(i)+b13
        y=b21*xs(i)+b22*ys(i)+b23
        wht=whts(i)*det
        call funuser(x,y,par1,par2,par3,par4,val)
        sum=sum+val*wht
 1100 continue
c
        end
