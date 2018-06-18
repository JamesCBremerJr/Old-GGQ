        implicit double precision (a-h,o-z)
        dimension verts(2,3),vals(1000),a(2,3)
        external funtest1
c
        verts(1,1) = 0.0d0
        verts(2,1) = 0.0d0
c
        verts(1,2) = 1.0d0
        verts(2,2) = 0.0d0
c
        verts(1,3) = 1.0d0
        verts(2,3) = 1.0d0
c
        eps=1.0d-30
c
        call adaptri(ier,eps,verts,funtest1,par1,par2,par3,par4,
     1    val,nquad)
c
        val0 = ( 3 + 147*cos(10.0d0)+470*sin(10.0d0) ) / 3750.0d0

        print *,"ier    =  ",ier
        print *,"val    =  ",val
        print *,"val0   =  ",val0
        print *,"relerr =  ",abs( (val-val0) / val0)
        print *,"nquad  =  ",nquad
c         
        end
c        
c
c
        subroutine funtest1(x,y,par1,par2,par3,par4,val)
        implicit double precision (a-h,o-z)
        val = cos(10*x)*(x**2+y**2)
        end

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
        maxstack = 100 000
        maxdepth = 150
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
c       Use a 14th order quadrature on the simplex to integrate the
c       user-supplied function over an arbitrary triangle.
c
        dimension verts(2,3)

        dimension xs(41),ys(41),whts(41)
        data norder / 14 /
        data nquad  / 41 /
        data xs / 
     -     0.379241137994507245961110233018D-01,
     -     0.388127271554790300546534959810D-02,
     -     0.795002666100827715690234682450D+00,
     -     0.111735737207780419171237103992D+00,
     -     0.964894970344523186234642491377D+00,
     -     0.100022228748153339361157110444D-01,
     -     0.823888165574639200427500724241D-01,
     -     0.369931660325177192240186720635D+00,
     -     0.644370849362085324225385704233D-01,
     -     0.207720673534545277055740252835D+00,
     -     0.716505205526240153467954449054D+00,
     -     0.310657639315710815975193338625D+00,
     -     0.374178432266265041133781813169D+00,
     -     0.498307841762385952187822324630D+00,
     -     0.551215608275970052799400050438D+00,
     -     0.195115900183470666922150017986D+00,
     -     0.872955254613423925061619362702D+00,
     -     0.573000372773246367454665676990D+00,
     -     0.703151857245094710022067810499D-01,
     -     0.606946637256225831075192151426D+00,
     -     0.587811837550479332646742726454D+00,
     -     0.888696807418322544944218914601D+00,
     -     0.225162855190577331406479062465D-01,
     -     0.197035121183009209342568415679D+00,
     -     0.150148823178825671166997224422D-01,
     -     0.241865180190639067346742189284D+00,
     -     0.767794573356033498149425985088D+00,
     -     0.228901000735021918411609081552D+00,
     -     0.202996026529384975653802288854D-01,
     -     0.244683854293849883605338608480D+00,
     -     0.195718646227423605465657307771D+00,
     -     0.378061584517274642333008280137D+00,
     -     0.576528369834951332129920244250D+00,
     -     0.200497605562973018193621541658D-01,
     -     0.403415136352470400711624041701D+00,
     -     0.102917927823445407922212553963D+00,
     -     0.746540692904821697063654528508D+00,
     -     0.104269795407633170742307487750D+00,
     -     0.204094085868971213575340475112D-01,
     -     0.425991541229021324860143447707D+00,
     -     0.103788203949320998847049036350D+00/
        data ys / 
     -     0.102617650294124713662016062442D+01,
     -     0.217093991213185712428498264208D-01,
     -    -0.240555403090534055745584328148D-01,
     -     0.108913635130374533989341061860D+00,
     -     0.153626093673109338382911802633D-01,
     -     0.263993578390463984129886624536D+00,
     -     0.882852876262845569015742243094D+00,
     -     0.611013142124332033892390640965D+00,
     -     0.219027658884281646018104665093D-01,
     -     0.905147802542381889324373787827D-01,
     -     0.428097318557745071028956685434D-01,
     -     0.100824990180069150793518489717D+00,
     -     0.189306018638674294645488565272D-01,
     -     0.820117457441886092159571136763D-01,
     -     0.412605008912271256793460310421D+00,
     -     0.192401919399113452363919165433D-01,
     -     0.256301122353903314296640820007D-01,
     -     0.139217636133773786239694931503D-01,
     -     0.254140056506010585817976102126D+00,
     -     0.284919026740046278279219910273D+00,
     -     0.165333715619788895326254218629D+00,
     -     0.911749071165622071961039549113D-01,
     -     0.112098766554988108940583332864D+00,
     -     0.233104746895981975247702971286D+00,
     -     0.948862022678710137026627826102D+00,
     -     0.385053817074802190418809706183D+00,
     -     0.126886172879332545912243850316D+00,
     -     0.689078355716664662311458812656D+00,
     -     0.646117241341952167189949577688D+00,
     -     0.543024571218556820755321945278D+00,
     -     0.790688324723617330809616759858D+00,
     -     0.213937739969857641358570711119D+00,
     -     0.418756459846227575711455471068D+00,
     -     0.449878051880395178844598100688D+00,
     -     0.492932801160880422227320989106D+00,
     -     0.420559405922335051774577019446D+00,
     -     0.232173919987473530623577544160D+00,
     -     0.597721650658393455491467518692D+00,
     -     0.820796315256715649174256478889D+00,
     -     0.342429246147297121952502462952D+00,
     -     0.753350153662899119449811995806D+00/
        data whts / 
     -     0.461824026594731903180626917301D-04,
     -     0.117253243510077444534257356703D-02,
     -     0.561417558322681032933587027448D-03,
     -     0.132075076763035609842916105746D-01,
     -     0.208424696935672605479324394957D-02,
     -     0.510338425007995062295981702555D-02,
     -     0.771723968315331339543794425125D-02,
     -     0.922943692628572921982954576351D-02,
     -     0.548564775945327899238966732374D-02,
     -     0.739702112541010352991348072097D-02,
     -     0.145704233330180755254365558938D-01,
     -     0.168976356748093057085608507965D-01,
     -     0.930902223068562129787999238049D-02,
     -     0.182460824060165738334834838715D-01,
     -     0.952854797673418443881295118534D-02,
     -     0.780579948823604644324336541272D-02,
     -     0.719822854006889037458126884049D-02,
     -     0.709327423228739340508943094379D-02,
     -     0.154196822451720541401999885077D-01,
     -     0.200292987175634642670001624798D-01,
     -     0.242478875527765661745763037713D-01,
     -     0.570109296782062361976066881508D-02,
     -     0.703192481205137485771045862322D-02,
     -     0.238692546013266078909324132716D-01,
     -     0.336987116182436457972791459148D-02,
     -     0.265921425320554249561174298660D-01,
     -     0.155544041449437858139333846344D-01,
     -     0.158114026039704849449625920888D-01,
     -     0.987126168681624685533893386767D-02,
     -     0.244238487559911536021900495279D-01,
     -     0.550860820292274213900459595299D-02,
     -     0.266147081739593784977931521393D-01,
     -     0.286477470649063734964005438126D-02,
     -     0.100730361992160060104767663198D-01,
     -     0.203841665174592781679241175964D-01,
     -     0.202835197318126782148290709925D-01,
     -     0.903911609493778059625819563275D-02,
     -     0.195896324527121488563941067657D-01,
     -     0.799566440919297315403247728584D-02,
     -     0.278600472699355548709080123195D-01,
     -     0.152110237910669879449917183492D-01/
        external funuser
c

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
