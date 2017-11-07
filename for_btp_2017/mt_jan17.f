c	motion program
	implicit real*8 (a-h,o-z)
	parameter(mxs=50,mxp=200,mxf=100,mxh=20,mxv=20)
c	mxs = max. no. of statins
c	mxp = max. no. of points to define a section
c	mxf = max. no. of frequencies
c	mxh = max. no. of headings
c	mxv = max. no. of velocities

	complex*8 a(mxp,mxp),b(mxp,mxp),f(mxp,4),phi(mxp,4),bb(mxp)
	complex*8 sum,summ,sum_mi,fi(6),fd(6),fex(6)
	complex*8 aa(6,6),zamp(6)
	complex*8 aaa(6,6)
	complex*8 phi1,phi2,phi3,phi4
	complex*8 rca3,rca5
	complex*8 za33,za55,za33n,za55n
	complex*8 q11,q12,q21,q22
	complex*8 rollcm

	common/wksp/pi,grav,rho
	common/space0/volm,xcb,zcb,wparea,xcf,bmt,bml,dels
	common/space1/nstn,npt(mxs),npt1(mxs),npt2(mxs),nsegg(mxs)
	common/space2/al,stn(mxs),xstn(mxs),xjj(mxs,mxp),yjj(mxs,mxp),
     &		      cg1(mxs,mxp),cg2(mxs,mxp),dell(mxs,mxp),
     &		      bn1(mxs,mxp),bn2(mxs,mxp),bnx(mxs,mxp),xt(mxs),
     &		      sa(mxs),by(mxs),sl(mxs),dz(mxs)
	dimension xjt(mxs,mxp),yjt(mxs,mxp)
	dimension admas(4,4),damp(4,4)
	dimension indx(mxp)
	dimension wvlnd(mxf),wvlen(mxf),wegv(mxf),head(mxh),hdrad(mxh)
	dimension vel(mxv),vknot(mxv),frno(mxv)
	dimension amm(18),bmm(18),as33(mxs),bs33(mxs),asl33(mxs)
	dimension av(6,6),bv(6,6)
	dimension fexamp(6),fexph(6),za(6),zph(6)
	dimension a11(mxs,mxf),a12(mxs,mxf),a13(mxs,mxf),a14(mxs,mxf),
     &	 	  a21(mxs,mxf),a22(mxs,mxf),a23(mxs,mxf),a24(mxs,mxf),
     &	 	  a31(mxs,mxf),a32(mxs,mxf),a33(mxs,mxf),a34(mxs,mxf),
     &	 	  a41(mxs,mxf),a42(mxs,mxf),a43(mxs,mxf),a44(mxs,mxf),
     &	 	  b11(mxs,mxf),b12(mxs,mxf),b13(mxs,mxf),b14(mxs,mxf),
     &	 	  b21(mxs,mxf),b22(mxs,mxf),b23(mxs,mxf),b24(mxs,mxf),
     &	 	  b31(mxs,mxf),b32(mxs,mxf),b33(mxs,mxf),b34(mxs,mxf),
     &	 	  b41(mxs,mxf),b42(mxs,mxf),b43(mxs,mxf),b44(mxs,mxf)
	dimension xbk(11),sabk(11),bybk(11),dzbk(11)
	common/space3/weint(mxf),nfr
	common/space4/am0(6,6,mxf),bm0(6,6,mxf)
	common/space5/phi1(mxs,mxf,mxp),phi2(mxs,mxf,mxp),
     &		      phi3(mxs,mxf,mxp),phi4(mxs,mxf,mxp)
	common/space6/aa33(mxs,mxf),bb33(mxs,mxf)

	dimension tmas(6,6),rest(6,6)

	open(1,file='input_folder/new_input_geometry',status="unknown")
	open(2,file='out2')
	open(3,file='out3')
c	open(4,file='out4')
c	open(7,file='out7')
	open(8,file='out8')
	open(11,file='out11')
	open(12,file='out12')
	open(13,file='out13')
	open(14,file='out14')
	open(15,file='out15')
	open(16,file='out16')
c	open(21,file='out21')
c	open(22,file='out22')
	open(23,file='out23')

c---------------------------------------------------------------------
c	File  `input_geometry` is the file sectional geometry data
c	Its structure is as follows:

c	nstn	length	n_ap	n_fp
c	stn. no.	no. of points to define the half station
c	j	y(j)	z(j)
c
c	The above to be repeated for all nstn no. of stations

c	here, nstn = no. of stations
c	length : LBP
c	n_ap : stn. no. at AP
c	stn. no. at FP

c	There can be maximum 50 stations.
c	Stations need not be equispaced
c	station nos are such that stn.no.*ds represent distance from AP
c	here ds = length/(n_fp-n_ap)
c	thus, if n_ap=0,n_fp=20, then ds=L/20.
c	you can have data entered at stations
c	-0.25 0 0.25 .5 1 2 ... 19 19.5 20 20.2  etc.
c	note the first and last stns are aft of ap and ford. of fp.
c

c	the data should be entered for one half,
c	starting from bottommost pt. on centreline to topmost point
c	the top point is usually on waterline, but need not be so for a
c	submerged bulb for example (in which case it will be the topmost
c	point on centreline)

c	max. no of point for half section is 100 (full section is 200)
c

c	the other data are to be ebntered as asked in screen.

c	note the followiong limits (whichj can be easily changed)
c	max. no. of frequency = 100
c	max. no. of headings = 20
c	max. no. of speed = 20
c	files out11,out12,out13,out14 gives zero-speed added mass & damping

c	file out15 gives exciting forces/moments
c	file out16 gives motion raos

c	all are non-dimensionalized.
c	the non-dimensiolanization factors are:
c	frequency : sqrt(L/g)
c	added masses : m, mL and mL^2
c	damp. = addee mass term * sqrt(g/L)
c	forces by m*g*a/L, moments by m*g*a
c	a = wave ampl. is taken as 1
c	RAO for linear are eta/a, for angles are eta/(ka), k=wave-number
c
c	presently no interpolation are done for irregular frequencies
c	that may occur in the sourc method at high frequencies

c	also, no roll damping is yet added

c-------------------------------------------------------------

c	initialization of data

	do i = 1,mxs
		npt(i) = 0.
		npt1(i) = 0.
		npt2(i) = 0.
		nsegg(i) = 0.
		stn(i) = 0.
		xstn(i) = 0.
		xt(i) = 0.
		sa(i)= 0.
		by(i) = 0.
		sl(i) = 0.
		dz(i) = 0.
		do j = 1,mxp
			xjj(i,j) = 0.
			yjj(i,j) = 0.
			cg1(i,j) = 0.
			cg2(i,j) = 0.
			dell(i,j) = 0.
			bn1(i,j) = 0.
			bn2(i,j) = 0.
			bnx(i,j) = 0.
			xjt(i,j) = 0.
			yjt(i,j) = 0.
		enddo
	enddo
	do i = 1,mxf
		wvlnd(i) = 0.
		wvlen(i) = 0.
		wegv(i) = 0.
		weint(i) = 0.
	enddo
	do i = 1,mxh
		head(i) = 0.
		hdrad(i) = 0.
	enddo
	do i = 1,mxv
		vel(i) = 0.
		vknot(i) = 0.
		frno(i) = 0.
	enddo

	do i = 1,mxs
		do j = 1,mxf
			do k = 1,mxp
				phi1(i,j,k) = (0.0,0.0)
				phi2(i,j,k) = (0.0,0.0)
				phi3(i,j,k) = (0.0,0.0)
				phi4(i,j,k) = (0.0,0.0)
			enddo
		enddo
	enddo

	do i = 1,mxs
		do j = 1,mxf
			a11(i,j) = 0.
			a12(i,j) = 0.
			a13(i,j) = 0.
			a14(i,j) = 0.
			a21(i,j) = 0.
			a22(i,j) = 0.
			a23(i,j) = 0.
			a24(i,j) = 0.
			a31(i,j) = 0.
			a32(i,j) = 0.
			a33(i,j) = 0.
			a34(i,j) = 0.
			a41(i,j) = 0.
			a42(i,j) = 0.
			a43(i,j) = 0.
			a44(i,j) = 0.
			b11(i,j) = 0.
			b12(i,j) = 0.
			b13(i,j) = 0.
			b14(i,j) = 0.
			b21(i,j) = 0.
			b22(i,j) = 0.
			b23(i,j) = 0.
			b24(i,j) = 0.
			b31(i,j) = 0.
			b32(i,j) = 0.
			b33(i,j) = 0.
			b34(i,j) = 0.
			b41(i,j) = 0.
			b42(i,j) = 0.
			b43(i,j) = 0.
			b44(i,j) = 0.
		enddo
	enddo

	pi = 3.14159265358979323846
	grav = 9.806
	pi2 = 2.*pi
	rho = 1.025
	tmp = 20.
	g = grav
	radius = 10.
	degrad = pi/180.
	anue = 0.0178/(1.+0.0336*tmp + 0.000221*tmp**2)
	anue = anue/1.0e04
	eproll = pi/720.

	call geom_data

c	write(*,*)' x coordinate of the origin from AP?'
	xorgn = xcb
	write(*,*)' z coor. of cg from origin at waterline'
	read(*,*) zcg
	xg2 = zcg
	og = -zcg
c	determin the stn locations from origin of coordibnate system
	do i = 1,nstn
		xt(i) = xstn(i) - xorgn
	enddo
c	determine maximum breadth
	br = 0.
	do i = 1,nstn
		if (br .lt. by(i) ) br = by(i)
	enddo
	br = 2.*br
	write(*,*)' br =',br

c	determine maximum draft
	dr = 0.
	do i = 1,nstn
		if (dr .lt. dz(i) ) dr = dz(i)
	enddo
	write(*,*)' dr =',dr

c	determine maximum section area (taken as midship)
	sam = 0.
	do i = 1,nstn
		if (sam .lt. sa(i) )sam = sa(i)
	enddo
	write(*,*)' sam =',sam
c	pause

c	reading hydrostatic property values
	vol = volm
	awp = wparea
	gml = zcb + bml - zcg
	cb = vol/(al*br*dr)
	cm = 2.*sam/(br*dr)
	write(*,*)' cb=',cb,' cm=',cm

  123	continue
	gmt = zcb + bmt - zcg
	write(*,*)' values taken are'
	write(*,*)' water plane area =',awp
	write(*,*)' volume =',vol
	write(*,*)' gmt =',gmt
	write(*,*)' gml =',gml

	if (gmt .lt. 0.) then
		write(*,*)' GM_T is negative; re-enter Zcg'
		read(*,*) zcg
		xg2 = zcg
		og = -zcg
		go to 123
	endif

c	set rigid-body mass matrix + restoring coefficients

	write(*,*)' start and end stn. no. of bilge keel'
	read(*,*) stbk1,stbk2
	xbk1 = stbk1*dels
	xbk2 = stbk2*dels
	delbk = (xbk2 - xbk1)/10.
	write(*,*) stbk1,stbk2,xbk1,xbk2,delbk
	write(*,*)' width of bilge keel'
	read(*,*) bk
c	interpolation for bilge keel stns
	do j = 1,11
		xgiven = xbk1 + (j-1)*delbk
		xbk(j) = xgiven - xorgn
c		xbk(j) = xgiven
		call interp(nstn,xstn,sa,by,dz,xgiven,sag,byg,dzg)
		sabk(j) = sag
		bybk(j) = byg
		dzbk(j) = dzg
	enddo
c	do j = 1,11
c	write(4,45) j,xbk(j),bybk(j),dzbk(j),sabk(j)
c	enddo
c   45	format(i6,4f10.4)

  455	continue
	write(*,*)' index for consideration of viscous roll damping'
	write(*,*)' 0 if viscous roll damping is to be ignored'
	write(*,*)' 1 if viscous roll damping is to be considered'
	write(*,*)' enter the value, 0 or 1 only'
	read(*,*) krollv
	if (krollv .eq. 0 .or. krollv .eq. 1) go to 456
	go to 455
  456	continue
c	writing geometric data
	do i = 1,nstn
	write(2,*)' '
	write(2,*)'          station no=',i,' x =',xt(i),xstn(i)
	write(2,*)'    j       y          z        ds         ny
     &nz         nx'
		do j = 1,nsegg(i)
		write(2,33) j,cg1(i,j),cg2(i,j),dell(i,j),
     &			     bn1(i,j),bn2(i,j),bnx(i,j)
		enddo
	enddo
   33	format(i6,3f10.3,3f12.6)

c----------------------------------------------------------------
c	setting rigid body mass matrices

c	setting rigid-body mass matric and restoring matrix
      	do i=1,6
      		do  j=1,6
      			tmas(i,j)=0.0
      			rest(i,j)=0.0
		enddo
	enddo

c	read hydro. + mass properties
	write(*,*)' rx/b,ry/l'
	read(*,*) rx1,ry1
	rxx = rx1 * br
	ryy = ry1 * al
	rzz = dsqrt(rxx**2+ryy**2)
	write(*,*)' rxx=',rxx,' ryy=',ryy,' rzz=',rzz

c	set rigid-body mass matrix + restoring coefficients
	smas = rho*vol
	ti44 = smas * rxx**2
	ti55 = smas * ryy**2
	ti66 = smas * rzz**2
	ti45 = 0.
	ti54 = 0.
	ti46 = 0.
	ti64 = 0.
	ti56 = 0.
	ti65 = 0.
	c33 = rho * grav * awp
	c44 = rho * grav * vol * gmt
	c55 = rho * grav * vol * gml
	c34 = 0.
	c43 = 0.
	c35 = 0.
	c53 = 0.
	c36 = 0.
	c63 = 0.
	c45 = 0.
	c54 = 0.
	c46 = 0.
	c64 = 0.
	c56 = 0.
	c65 = 0.

      	zmcc=smas*zcg
      	tmas(1,5)=zmcc
      	tmas(2,4)=-zmcc
      	tmas(4,2)=-zmcc
      	tmas(5,1)=zmcc
      	tmas(1,1)=smas
      	tmas(2,2)=smas
      	tmas(3,3)=smas
      	tmas(4,4)=ti44
      	tmas(4,6)=-ti46
      	tmas(5,5)=ti55
      	tmas(6,4)=-ti64
      	tmas(6,6)=ti66
      	tmas(4,5)=-ti45
      	tmas(5,4)=-ti45
      	tmas(5,6)=-ti56
      	tmas(6,5)=-ti56
      	rest(3,3)=c33
      	rest(3,5)=c35
      	rest(5,3)=c35
      	rest(4,4)=c44
      	rest(5,5)=c55
      	rest(4,5)=c45
      	rest(5,4)=c45
c-----------------------------------------------------------------
	write(*,*)'do you want to give wavelength/bodylength or period'
	write(*,*)' if former, enter 1, else enter 2'
	read(*,*) kopt
	if (kopt. eq. 1) then
c	here first find the min. and max. freq. limits.
	write(*,*)' no. of wave lengths'
	read(*,*) nfreq
	write(*,*)' wave-length/body length'
	read(*,*) (wvlnd(i),i=1,nfreq)
	else if (kopt. eq. 2) then
	write(*,*)' write min and max. periods, and interval'
	read(*,*) tmin,tmax,delt
	nfreq = (tmax-tmin)/delt+1
	if (nfreq .gt. 100) nfreq = 100
	do i = 1,nfreq
	period = tmin + (i-1)*delt
	wvlnd(i) = (grav*period**2)/(pi2*al)
	enddo
	else
	write(*,*)' wrong option'
	go to 9999
	endif
	do i = 1,nfreq
	write(*,*)'i=',i,' lambda/L=',wvlnd(i)
	enddo

	write(*,*)' no. of headings'
	read(*,*) nhead
	write(*,*)' enter headings. in degree'
	read(*,*) (head(i),i=1,nhead)
	write(*,*)' no. of speeds'
	read(*,*) nspeed
	write(*,*)' speed or Froude no?'
	write(*,*)' for speed, enter 1, else for Fn, enter 0'
	read(*,*) kspeed
	if (kspeed .eq. 1) then
		write(*,*)' enter sepeds in knots'
		read(*,*) (vknot(i),i=1,nspeed)
		do i = 1,nspeed
			vel(i) = vknot(i) * 0.5144
			frno(i) = vel(i)/dsqrt(grav*al)
		enddo
	else
		write(*,*)' enter Froude nos'
		read(*,*) (frno(i),i=1,nspeed)
		do i = 1,nspeed
			vel(i) = frno(i)*dsqrt(grav*al)
			vknot(i) = vel(i)/0.5144
		enddo
	endif
	do i = 1,nspeed
	write(*,*)' vknot=',vknot(i),'vel=',vel(i),' fr=',frno(i)
	enddo
c	pause

c	change heading angle to radian
	do i = 1,nhead
		hdrad(i) = head(i) * degrad
	enddo
c	change speed to m/s
c	do i = 1,nspeed
c		vel(i)=vknot(i) * 0.5144
c	enddo

c	here we find the min and max encouter freq.
	do i = 1,nfreq
		wvlen(i) = wvlnd(i) * al
	enddo
	otmin = 99999.
	otmax = 0.
	do i = 1,nhead
		do j = 1,nspeed
			do k = 1,nfreq
				wnum = pi2/wvlen(k)
		otemp = abs(dsqrt(grav*wnum)-wnum*vel(j)*cos(hdrad(i)))
				if (otemp .lt. otmin) otmin = otemp
				if (otemp .gt. otmax) otmax = otemp
			enddo
		enddo
	enddo


	ttmax = pi2/otmin+1
	wemin = pi2/ttmax
	ttmin = pi2/otmax-1
	if (ttmin .lt. 1) ttmin = 1
	wemax = pi2/ttmin

	write(*,*)' wemin=',wemin,' tmax=',pi2/wemin
	write(*,*)' wemax=',wemax,' tmin=',pi2/wemax

c	stn loop
	do ks = 1,nstn

	write(*,*)' stn =',ks

c	freq. loop
	nfr = 100
	nfr = 50
c	nfr = 10

	delfr = (wemax-wemin)/(nfr-1)
	do kkk = 1,nfr
	weint(kkk) = wemin + (kkk-1)*delfr
	w = weint(kkk)
	t = 2.*pi/w
	freq = w
	anu = w**2/g
	nseg = nsegg(ks)
c	initalize
	do i = 1,mxp
		bb(i) = (0.0,0.0)
		do j = 1,mxp
			a(i,j) = (0.0,0.0)
			b(i,j) = (0.0,0.0)
		enddo
		do j = 1,4
			f(i,j) = (0.0,0.0)
			phi(i,j) = (0.0,0.0)
		enddo
	enddo

	do i = 1,nseg

		x1 = xjj(ks,i)
		y1 = yjj(ks,i)
		x2 = xjj(ks,i+1)
		y2 = yjj(ks,i+1)
		do j = 1,nsegg(ks)
			a1 = xjj(ks,j)
			b1 = yjj(ks,j)
			a2 = xjj(ks,j+1)
			b2 = yjj(ks,j+1)
			ij = 0
			if (i .eq.j) ij = 1

			call assembly(ij,anu,x1,y1,x2,y2,a1,b1,a2,b2,
     &  	term1,term2,term3,term4,bterm1,bterm2,bterm3,bterm4)

	a(i,j) =  dcmplx(( term1- term2+2.* term3)/pi2,-term4)
	b(i,j) =  dcmplx((bterm1-bterm2+2.*bterm3)/pi2,bterm4)

		enddo
	enddo

c	inversion of matrix a
	nnn = nseg
	call ludcmp(a,nnn,mxp,indx)

c	formation of {b} vector

	k = 1
	dumm1 = 0.0
	do i = 1,nseg
		dumm2 = -w*bnx(ks,i)
		bb(i) = dcmplx(dumm1,dumm2)
	enddo
	call lubksb(a,nnn,mxp,indx,bb)
	do i = 1,nseg
		f(i,k) = bb(i)
	enddo

	k = 2
	dumm1 = 0.0
	do i = 1,nseg
		dumm2 = -w*bn1(ks,i)
		bb(i) = dcmplx(dumm1,dumm2)
	enddo
	call lubksb(a,nnn,mxp,indx,bb)
	do i = 1,nseg
		f(i,k) = bb(i)
	enddo

	k = 3
	dumm1 = 0.0
	do i = 1,nseg
		dumm2 = -w*bn2(ks,i)
		bb(i) = dcmplx(dumm1,dumm2)
	enddo
	call lubksb(a,nnn,mxp,indx,bb)
	do i = 1,nseg
		f(i,k) = bb(i)
	enddo

	k = 4
	dumm1 = 0.0
	do i = 1,nseg
		dumm2=-w*(cg1(ks,i)*bn2(ks,i)-(cg2(ks,i)-xg2)*bn1(ks,i))
		bb(i) = dcmplx(dumm1,dumm2)
	enddo
	call lubksb(a,nnn,mxp,indx,bb)
	do i = 1,nseg
		f(i,k) = bb(i)
	enddo

	do k = 1,4
		do i = 1,nseg
			summ = (0.0,0.0)
			do j = 1,nseg
				summ = summ + b(i,j)*f(j,k)
			enddo
			phi(i,k) = summ
		enddo
	enddo

c	storing the 2d sectional potentials

	do i = 1,nseg
		phi1(ks,kkk,i) = phi(i,1)
		phi2(ks,kkk,i) = phi(i,2)
		phi3(ks,kkk,i) = phi(i,3)
		phi4(ks,kkk,i) = phi(i,4)
	enddo
c
c	calculation of sectional added mass and damping
c	initialize
	do k = 1,4
		do j = 1,4
			admas(k,j) = 0.0
			damp(k,j) = 0.0
		enddo
	enddo

	do j = 1,4
		do k = 1,4
			summ = (0.0,0.0)
			do i = 1,nseg
				if (k .eq. 1) then
					dn = bnx(ks,i)
				else if (k .eq. 2) then
					dn = bn1(ks,i)
				else if (k .eq. 3) then
					dn = bn2(ks,i)
				else if (k .eq. 4) then
			dn=cg1(ks,i)*bn2(ks,i)-(cg2(ks,i)-xg2)*bn1(ks,i)
				endif
				summ = summ + phi(i,j)*dell(ks,i)*dn
			enddo

			ax1 = aimag(summ)
			ax2 = real(summ)

			admas(k,j) = rho*ax1/freq
			damp(k,j)  = rho*ax2

		enddo
	enddo

c	here we store the sectional aded masses for each freq. and station

	a11(ks,kkk) = admas(1,1)
	a12(ks,kkk) = admas(1,2)
	a13(ks,kkk) = admas(1,3)
	a14(ks,kkk) = admas(1,4)
	a21(ks,kkk) = admas(2,1)
	a22(ks,kkk) = admas(2,2)
	a23(ks,kkk) = admas(2,3)
	a24(ks,kkk) = admas(2,4)
	a31(ks,kkk) = admas(3,1)
	a32(ks,kkk) = admas(3,2)
	a33(ks,kkk) = admas(3,3)
	a34(ks,kkk) = admas(3,4)
	a41(ks,kkk) = admas(4,1)
	a42(ks,kkk) = admas(4,2)
	a43(ks,kkk) = admas(4,3)
	a44(ks,kkk) = admas(4,4)

	b11(ks,kkk) = damp(1,1)
	b12(ks,kkk) = damp(1,2)
	b13(ks,kkk) = damp(1,3)
	b14(ks,kkk) = damp(1,4)
	b21(ks,kkk) = damp(2,1)
	b22(ks,kkk) = damp(2,2)
	b23(ks,kkk) = damp(2,3)
	b24(ks,kkk) = damp(2,4)
	b31(ks,kkk) = damp(3,1)
	b32(ks,kkk) = damp(3,2)
	b33(ks,kkk) = damp(3,3)
	b34(ks,kkk) = damp(3,4)
	b41(ks,kkk) = damp(4,1)
	b42(ks,kkk) = damp(4,2)
	b43(ks,kkk) = damp(4,3)
	b44(ks,kkk) = damp(4,4)

	aa33(ks,kkk) = a33(ks,kkk)
	bb33(ks,kkk) = b33(ks,kkk)

 	enddo

 	enddo

c	calculation for zero speed added masses
	do kkk = 1,nfr
		do ii = 1,6
			do jj = 1,6
				am0(ii,jj,kkk) = 0.0
				bm0(ii,jj,kkk) = 0.0
			enddo
		enddo
	enddo

	do kkk = 1,nfr
c	write(*,*)' kkk =',kkk
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
		sum1 = sum1 + 0.5*( a11(ks-1,kkk) + a11(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		sum2 = sum2 + 0.5*( b11(ks-1,kkk) + b11(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(1,1,kkk) = sum1
		bm0(1,1,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
		sum1 = sum1 + 0.5*( a13(ks-1,kkk) + a13(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		sum2 = sum2 + 0.5*( b13(ks-1,kkk) + b13(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(1,3,kkk) = sum1
		bm0(1,3,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
		sum1 = sum1 + 0.5*( a31(ks-1,kkk) + a31(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		sum2 = sum2 + 0.5*( b31(ks-1,kkk) + b31(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(3,1,kkk) = sum1
		bm0(3,1,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
	sum1 = sum1-0.5*(xt(ks-1)*a13(ks-1,kkk) + xt(ks)*a13(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
	sum2 = sum2-0.5*(xt(ks-1)*b13(ks-1,kkk) + xt(ks)*b13(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(1,5,kkk) = sum1
		bm0(1,5,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
	sum1 = sum1-0.5*(xt(ks-1)*a31(ks-1,kkk) + xt(ks)*a31(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
	sum2 = sum2-0.5*(xt(ks-1)*b31(ks-1,kkk) + xt(ks)*b31(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(5,1,kkk) = sum1
		bm0(5,1,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
		sum1 = sum1 + 0.5*( a33(ks-1,kkk) + a33(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		sum2 = sum2 + 0.5*( b33(ks-1,kkk) + b33(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(3,3,kkk) = sum1
		bm0(3,3,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
	sum1 = sum1-0.5*(xt(ks-1)*a33(ks-1,kkk) + xt(ks)*a33(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
	sum2 = sum2-0.5*(xt(ks-1)*b33(ks-1,kkk) + xt(ks)*b33(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(3,5,kkk) = sum1
		bm0(3,5,kkk) = sum2

		am0(5,3,kkk) = am0(3,5,kkk)
		bm0(5,3,kkk) = bm0(3,5,kkk)
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
	sum1=sum1+0.5*(xt(ks-1)**2*a33(ks-1,kkk)+xt(ks)**2*a33(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
	sum2=sum2+0.5*(xt(ks-1)**2*b33(ks-1,kkk)+xt(ks)**2*b33(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(5,5,kkk) = sum1
		bm0(5,5,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
		sum1 = sum1 + 0.5*( a22(ks-1,kkk) + a22(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		sum2 = sum2 + 0.5*( b22(ks-1,kkk) + b22(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(2,2,kkk) = sum1
		bm0(2,2,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
		sum1 = sum1 + 0.5*( a24(ks-1,kkk) + a24(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		sum2 = sum2 + 0.5*( b24(ks-1,kkk) + b24(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(2,4,kkk) = sum1
		bm0(2,4,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
		sum1 = sum1 + 0.5*( a42(ks-1,kkk) + a42(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		sum2 = sum2 + 0.5*( b42(ks-1,kkk) + b42(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(4,2,kkk) = sum1
		bm0(4,2,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
		sum1 = sum1 + 0.5*( a44(ks-1,kkk) + a44(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		sum2 = sum2 + 0.5*( b44(ks-1,kkk) + b44(ks,kkk) )
     &			       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(4,4,kkk) = sum1
		bm0(4,4,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
	sum1 = sum1+0.5*(xt(ks-1)*a24(ks-1,kkk) + xt(ks)*a24(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
	sum2 = sum2+0.5*(xt(ks-1)*b24(ks-1,kkk) + xt(ks)*b24(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(4,6,kkk) = sum1
		bm0(4,6,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
	sum1 = sum1+0.5*(xt(ks-1)*a42(ks-1,kkk) + xt(ks)*a42(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
	sum2 = sum2+0.5*(xt(ks-1)*b42(ks-1,kkk) + xt(ks)*b42(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(6,4,kkk) = sum1
		bm0(6,4,kkk) = sum2
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
	sum1 = sum1+0.5*(xt(ks-1)*a22(ks-1,kkk) + xt(ks)*a22(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
	sum2 = sum2+0.5*(xt(ks-1)*b22(ks-1,kkk) + xt(ks)*b22(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(2,6,kkk) = sum1
		bm0(2,6,kkk) = sum2

		am0(6,2,kkk) = am0(2,6,kkk)
		bm0(6,2,kkk) = bm0(2,6,kkk)
	enddo

	do kkk = 1,nfr
		sum1 = 0.
		sum2 = 0.
		do ks = 2,nstn
	sum1=sum1+0.5*(xt(ks-1)**2*a22(ks-1,kkk)+xt(ks)**2*a22(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
	sum2=sum2+0.5*(xt(ks-1)**2*b22(ks-1,kkk)+xt(ks)**2*b22(ks,kkk))
     &	       * ( xt(ks)-xt(ks-1) )
		enddo
		am0(6,6,kkk) = sum1
		bm0(6,6,kkk) = sum2
	enddo

c	write zero speed added mass and damping
	write(11,*)'      we(nd)      a(1,1)       a(2,2)        a(3,3)
     & 	   a(4,4)        a(5,5)         a(6,6)'
	write(13,*)'      we(nd)      b(1,1)       b(2,2)        b(3,3)
     & 	   b(4,4)        b(5,5)         b(6,6)'
	write(12,*)'      we(nd)      a(3,5)       a(5,3)        a(2,6)
     & 	   a(6,2)        a(2,4)         a(4,2)       a(4,6)        a(6,4
     &)'
	write(14,*)'      we(nd)      b(3,5)       b(5,3)        b(2,6)
     & 	   b(6,2)        b(2,4)         b(4,2)       b(4,6)        b(6,4
     &)'
	af0 = dsqrt(al/grav)
	af1 = smas
	af2 = smas*al
	af3 = smas*al**2
	bf1 = af1/af0
	bf2 = af2/af0
	bf3 = af3/af0
	do kkk = 1,nfr
		freq = weint(kkk)
		fnd = freq*af0
		write(11,81) fnd,am0(1,1,kkk)/af1,am0(2,2,kkk)/af1,
     &				 am0(3,3,kkk)/af1,am0(4,4,kkk)/af3,
     &				 am0(5,5,kkk)/af3,am0(6,6,kkk)/af3
		write(12,82) fnd,am0(3,5,kkk)/af2,am0(5,3,kkk)/af2,
     &			         am0(2,6,kkk)/af2,am0(6,2,kkk)/af2,
     &			   	 am0(2,4,kkk)/af2,am0(4,2,kkk)/af2,
     &				 am0(4,6,kkk)/af3,am0(6,4,kkk)/af3
		write(13,81) fnd,bm0(1,1,kkk)/bf1,bm0(2,2,kkk)/bf1,
     &				 bm0(3,3,kkk)/bf1,bm0(4,4,kkk)/bf3,
     &				 bm0(5,5,kkk)/bf3,bm0(6,6,kkk)/bf3
		write(14,82) fnd,bm0(3,5,kkk)/bf2,bm0(5,3,kkk)/bf2,
     &			         bm0(2,6,kkk)/bf2,bm0(6,2,kkk)/bf2,
     &			   	 bm0(2,4,kkk)/bf2,bm0(4,2,kkk)/bf2,
     &				 bm0(4,6,kkk)/bf3,bm0(6,4,kkk)/bf3

	enddo
   81	format(f12.4,6e14.4)
   82	format(f12.4,8e14.4)

c	the zero-speed 3d added mass damping is now over

c	calculation for motions

c	non-dimensionalization factors
	wamp = 1.0
	cf0 = af0
	cf1 = smas*grav*wamp/al
	cf2 = smas*grav*wamp
	arf = rho*grav*br**2*wamp**2/al
	factdm = sqrt(br/(2.*grav)) /(rho*vol*br**2)

	write(15,*)'         non-dimensional exciting forces'
	write(15,*)'    (non-dimensional amplitudes and phases in deg.)'

	write(16,*)'         non-dimensional motion amplitudes'
	write(16,*)'    (non-dimensional RAOs and phases in deg.)'

	write(8,*)'       roll convergence results'

c	speed loop
	do i = 1,nspeed
		u0 = vel(i)

c		heading loop
		do j = 1,nhead
			hd = hdrad(j)
		write(15,*)' '
		write(15,*)' '
		write(15,*)'speed =',u0/0.5144,'knots',
     &			' Fn=',frno(i),
     &		' heading =',hd/degrad,'deg.'
		write(15,*)' '
		write(15,*)' we    w0    period  l/L        surge
     &sway           heave              roll               pitch
     &    yaw'
		write(15,*)' '

		write(16,*)' '
		write(16,*)' '
		write(16,*)'speed =',u0/0.5144,'knots',
     &			' Fn=',u0/dsqrt(al*grav),
     &		' heading =',hd/degrad,'deg.'
		write(16,*)' '
		write(16,*)' we    w0    period  l/L        surge
     &sway           heave              roll               pitch
     &    yaw'
		write(16,*)' '


	if (abs(hd-pi) .gt. eps) go to 9987
		write(23,*)' '
		write(23,*)' '
		write(23,*)'speed =',u0/0.5144,'knots',
     &			' Fn=',u0/dsqrt(al*grav),
     &		' heading =',hd/degrad,'deg.'
	write(23,*)' 1 : Havelock, 2: Joosen,  3: Gerritsma & Bukelman,
     &  4: Salvesen, 5 : Boese'
	write(23,*)' '
		write(23,*)' we        l/L  sqrt(L/l)     1         2
     &3        4         5'
 9987	continue

		write(8,*)' '
		write(8,*)' '
		write(8,*)'speed =',u0/0.5144,'knots',
     &			' Fn=',frno(i),
     &		' heading =',hd/degrad,'deg.'
		write(8,*)' '

c			frequency loop
			do k = 1,nfreq
				wnum = pi2/wvlen(k)
				w0 = dsqrt(grav*wnum)
				we = w0 - wnum*u0*cos(hd)
				we = abs(we)
				do ii = 1,6
					fi(ii) = (0.,0.)
					fd(ii) = (0.,0.)
					fex(ii) = (0.,0.0)
					zamp(ii) = (0.0,0.0)
					fexamp(ii) = 0.
					fexph(ii) = 0.
					za(ii) = 0.
					zph(ii) = 0.
					do jj = 1,6
						av(ii,jj) = 0.0
						bv(ii,jj) = 0.0
					enddo
				enddo

c	here we call subroutine to interpolate for added masses
				call inter_amass(we,amm,bmm)

				av(1,1) = amm(1)
				bv(1,1) = bmm(1)

				av(1,3) = amm(2)
				bv(1,3) = bmm(2)

				av(3,1) = amm(3)
				bv(3,1) = bmm(3)

				av(1,5) = amm(4) - (u0/we**2) * bv(1,3)
				bv(1,5) = bmm(4) +  u0        * av(1,3)

				av(5,1) = amm(5) + (u0/we**2) * bv(3,1)
				bv(5,1) = bmm(5) -  u0        * av(3,1)

				av(3,3) = amm(6)
				bv(3,3) = bmm(6)

				av(3,5) = amm(7) - (u0/we**2) * bv(3,3)
				bv(3,5) = bmm(7) +  u0        * av(3,3)

				av(5,3) = amm(8) + (u0/we**2) * bv(3,3)
				bv(5,3) = bmm(8) -  u0        * av(3,3)

				av(5,5) = amm(9) + (u0**2/we**2)*av(3,3)
				bv(5,5) = bmm(9) + (u0**2/we**2)*bv(3,3)

				av(2,2) = amm(10)
				bv(2,2) = bmm(10)

				av(2,4) = amm(11)
				bv(2,4) = bmm(11)

				av(4,2) = amm(12)
				bv(4,2) = bmm(12)

				av(4,4) = amm(13)
				bv(4,4) = bmm(13)

				av(4,6) = amm(14) + (u0/we**2) * bv(2,4)
				bv(4,6) = bmm(14) -  u0        * av(2,4)

				av(6,4) = amm(15) - (u0/we**2) * bv(4,2)
				bv(6,4) = bmm(15) +  u0        * av(4,2)

				av(2,6) = amm(16) + (u0/we**2) * bv(2,2)
				bv(2,6) = bmm(16) -  u0        * av(2,2)

				av(6,2) = amm(17) - (u0/we**2) * bv(2,2)
				bv(6,2) = bmm(17) +  u0        * av(2,2)

				av(6,6) = amm(18) + (u0**2/we**2)*av(2,2)
				bv(6,6) = bmm(18) + (u0**2/we**2)*bv(2,2)

c	write(21,7777) we*af0,av(3,3)/af1,av(5,5)/af3,av(3,5)/af2,
c     &			      av(5,3)/af2,
c     &			      bv(3,3)/bf1,bv(5,5)/bf3,bv(3,5)/bf2,
c     &			      bv(5,3)/bf2

c-----------------------------------------------------
c	first, we find FK forces

				call fkforce(xg2,hd,wnum,fi)

c	now we find diffraction forces

				call fdforce(xg2,hd,w0,we,wnum,u0,fd)

c				form external force matrix

				do ii = 1,6
					fex(ii) = fi(ii) + fd(ii)
					fexamp(ii) = cabs(fex(ii))
					qq1 = aimag(fex(ii))
					qq2 = real(fex(ii))
					fexph(ii) = datan2(qq1,qq2)/degrad
				enddo

	write(15,55)we*cf0,w0*cf0,pi2/w0,wvlen(k)/al,
     &	fexamp(1)/cf1,fexph(1),fexamp(2)/cf1,fexph(2),
     &	fexamp(3)/cf1,fexph(3),fexamp(4)/cf2,fexph(4),
     &	fexamp(5)/cf2,fexph(5),fexamp(6)/cf2,fexph(6)

c--------------------------------------------------------------
c	here we calculate the viscous roll damping coefficients
	omega = we
	call roll_damp_fric(al,br,dr,cb,og,omega,u0,anue,bf)
	if (k .eq. 1) call roll_damp_lift(al,br,dr,cm,og,u0,bl)
	bfnd = bf * factdm
	if (k .eq. 1) blnd = bl * factdm

	write(8,*)'   we    iter    theta1     theta2'
	write(8,*)' (rad/s)         (deg)      (deg)'
	iter = 1
	theta = 0.
    6	continue
	if (bk .ne. 0.) then
	call roll_damp_bk(xbk,bybk,dzbk,sabk,br,dr,og,omega,
     &	bk,delbk,theta,bbkn,bbkh)
	else
	bbkn = 0.
	bbkh = 0.
	endif
	bbknd = bbkn * factdm
	bbkhd = bbkh * factdm
	bknd = bbknd + bbkhd
	call roll_damp_ed(xbk,bybk,dzbk,sabk,al,br,dr,og,omega,
     &	u0,delbk,theta,be)
	bend = be * factdm
	brnd = bv(4,4)*factdm
	b44tot = bfnd + blnd + bknd + bend
	b44v = b44tot/factdm
	fn = u0/dsqrt(grav*al)
c	write(7,445)iter,hd,fn,we,brnd,bfnd,blnd,bknd,bend,b44tot
c  445	format(i4,3f8.4,6e14.4)
c	find sdof roll
	t1_44 = - (we**2)*( tmas(4,4) + av(4,4) )
	t2_44 = we*( bv(4,4) + b44v )
	t3_44 = rest(4,4)
	t4_44 = t1_44 + t3_44
	rollcm = fex(4)/dcmplx(t4_44,t2_44)
	rollrl = cabs(rollcm)
	rollcv = abs(rollrl - theta) / eproll
	write(8,446) we,iter,theta*180./pi,
     &		    rollrl*180./pi
  446	format(f8.4,i6,2f10.2)
	iter = iter+1
	theta = rollrl
	if (rollcv .lt. 1. .or. iter. gt. 10) go to 7
	go to 6
    7	continue
	brtot = bv(4,4) + b44v
	termvv = (2.*sqrt( ( tmas(4,4)+av(4,4) ) * rest(4,4) ) )
	dmprad = bv(4,4)/termvv
	dmpvis = b44v/termvv
	dmptot = brtot/termvv

	write(8,*)' radiation damping  =',dmprad*100.00,
     &	'% of critical damping'
	write(8,*)' viscous damping    =',dmpvis*100.00,
     &	'% of critical damping'
	write(8,*)' total roll damping =',dmptot*100.00,
     &	'% of critical damping'
	write(8,*)' '

c--------------------------------------------------------

c				form the mass matrix

				do ii = 1,6
				do jj = 1,6
				t1 = - (we**2)*(tmas(ii,jj) + av(ii,jj))
				t2 = we*bv(ii,jj)
				if (ii .eq. 4 .and. jj .eq. 4)
     &					t2 = we*(bv(ii,jj)+krollv*b44v)
				t3 = rest(ii,jj)
				t4 = t1+t3
				aa(ii,jj) = dcmplx(t4,t2)
				aaa(ii,jj) = aa(ii,jj)
				enddo
				enddo

c				inversion of matrix a

				call invert1(aa)

c	here we solve for motions to get complex motion amplitudes
				do ii = 1,6
					sum = (0.0,0.0)
					do jj = 1,6
					sum = sum + aa(ii,jj)*fex(jj)
					enddo
					zamp(ii) = sum
					za(ii) = cabs(zamp(ii))
					qq1 = aimag(zamp(ii))
					qq2 = real(zamp(ii))
					zph(ii) = datan2(qq1,qq2)/degrad
				enddo

				rca3 = zamp(3)
				rca5 = zamp(5)

	write(16,55)we*cf0,w0*cf0,pi2/w0,wvlen(k)/al,
     &	za(1),zph(1),za(2),zph(2),za(3),zph(3),za(4)/wnum,zph(4),
     &	za(5)/wnum,zph(5),za(6)/wnum,zph(6)

c	write(22,7778) we*cf0,wvlen(k)/al,dsqrt(al/wvlen(k)),
c     &			za(3),zph(3),za(5)*al/pi2,zph(5)

c----------------------------------------------------------

c	write(*,*)' added res. starts'

c	added resistance computation done here
c	these shoulld be done only for 180 deg. heading

	if (abs(hd-pi) .gt. eps) go to 9988

c	option 1 : Havelock`s method as derived by Joosen

	b033 = bmm(6)
	b055 = bmm(9)
	b035 = bmm(7)

	bv33 = bv(3,3)
	bv55 = bv(5,5)

	ra3 = za(3)
	ra5 = za(5)

c	write(*,*)' ra3=',ra3,' ra5=',ra5

	e3 = zph(3)*degrad
	e5 = zph(5)*degrad

	raw1 = (we**2/(2.*grav)) * (b033*ra3**2 + b055*ra5**2)

c	option 2 :

	raw2 = (we**2/(2.*grav)) * (b033*ra3**2 + b055*ra5**2
     &				  - 2.*b035*ra3*ra5*cos(e3-e5))


c	option 3 Gerritsma and Beulelman

c	first we need to interpolate and get sectional added masses

	call inter_sec_mass(nstn,we,as33,bs33)

c	find slope of as33
	call slope(nstn,xt,as33,asl33)

c	write(*,*)' ra3=',ra3,' ra5=',ra5
	call rawgb(we,wnum,u0,bs33,asl33,ra3,ra5,e3,e5,wamp,raw)
	raw3 = raw

	call rawsv(w0,we,wnum,u0,as33,bs33,rca3,rca5,wamp,raw)
	raw4 = raw

	call rawbs(we,wnum,volm,ra3,ra5,e3,e5,wamp,raw)
	raw5 = raw


	write(23,7779) we*cf0,wvlen(k)/al,dsqrt(al/wvlen(k)),
     &		       raw1/arf,raw2/arf,raw3/arf,raw4/arf,raw5/arf

 9988	continue

c			freq. loop ends
			enddo

c		heading loop ends
		enddo

c	speed loop ends
	enddo
   55	format(2f6.3,2f6.2,6(f11.4,f7.1))
   56	format(f6.3,3f6.2,6(f11.4,f7.1))
 7777	format(f8.4,8e14.4)
 7778	format(3f8.4,2(f11.4,f7.1))
 7779	format(f8.4,2f8.2,5f10.3)


 9999	continue
	stop
	end

c--------------------------------------------

	subroutine geom_data
	implicit real*8 (a-h,o-z)

	parameter(mxs=50,mxp=200)
	common/space1/nstn,npt(mxs),npt1(mxs),npt2(mxs),nsegg(mxs)
	common/space2/al,stn(mxs),xstn(mxs),xjj(mxs,mxp),yjj(mxs,mxp),
     &		      cg1(mxs,mxp),cg2(mxs,mxp),dell(mxs,mxp),
     &		      bn1(mxs,mxp),bn2(mxs,mxp),bnx(mxs,mxp),xt(mxs),
     &		      sa(mxs),by(mxs),sl(mxs),dz(mxs)
	common/space0/volm,xcb,zcb,wparea,xcf,bmt,bml,dels
	dimension xjt(mxs,mxp),yjt(mxs,mxp),vm(mxs)

	eps = 1.0e-08

	read(1,*) nstn,al,stap,stfp
	dels = al/(stfp-stap)
	do i = 1,nstn
		read(1,*) stn(i),npt(i)
		xstn(i) = stn(i)*dels
		do j = 1,npt(i)
			read(1,*) jt,xjt(i,j),yjt(i,j)
		enddo
		sum1 = 0.
		sum2 = 0.
		nptt = npt(i)
		do j = 2,nptt
			y1 = xjt(i,j-1)
			y2 = xjt(i,j)
			x1 = yjt(i,j-1)
			x2 = yjt(i,j)
			sum1 = sum1 + 0.5*(y1+y2)*(x2-x1)
			sum2 = sum2 + 0.5*(y1*x1+y2*x2)*(x2-x1)
		enddo
		by(i) = xjt(i,nptt)
		dz(i) = -yjt(i,1)
		sa(i) = sum1
		vm(i) = sum2
	enddo
c	do i = 1,nstn
c		write(4,44) i,xstn(i),by(i),dz(i),sa(i)
c	enddo
   44	format(i6,4f10.4)
c	find slope of waterline
	do i = 1,nstn
		if (i .eq. 1) then
			ig = 2
			xg = xstn(1)
		else if (i .eq. nstn) then
			ig = nstn-1
			xg = xstn(nstn)
		else
			ig = i
			xg = xstn(i)
		endif
		a2 = ((by(ig)-by(ig-1))/(xstn(ig)-xstn(ig-1)) -
     & (by(ig+1)-by(ig))/(xstn(ig+1)-xstn(ig)))/(xstn(ig-1)-xstn(ig+1))
		a1 = (by(ig)-by(ig-1))/(xstn(ig)-xstn(ig-1))
     &		   - a2*(xstn(ig-1)+xstn(ig))
		sl(i) = a1 + 2.*a2*xg

c	write(23,*)'i=',i,' ig=',ig,' xg=',xg
c	write(23,*)' x=',xstn(ig-1),xstn(ig),xstn(ig+1)
c	write(23,*)' y=',by(ig-1),by(ig),by(ig+1)
c	write(23,*)' a2=',a2,' a1=',a1,' sl=',sl(i)

	enddo

c	find volume etc.
	sumvl = 0.
	sumlm = 0.
	sumvm = 0.

	sumwp = 0.
	sumcf = 0.
	sumix = 0.
	sumiy = 0.
	do i = 2,nstn
		x1 = xstn(i-1)
		x2 = xstn(i)

		y1 = sa(i-1)
		y2 = sa(i)
		sumvl = sumvl + 0.5*(y1+y2)*(x2-x1)

		y1 = vm(i-1)
		y2 = vm(i)
		sumvm = sumvm + 0.5*(y1+y2)*(x2-x1)

		y1 = sa(i-1)*xstn(i-1)
		y2 = sa(i)  *xstn(i)
		sumlm = sumlm + 0.5*(y1+y2)*(x2-x1)

		y1 = by(i-1)
		y2 = by(i)
		sumwp = sumwp + 0.5*(y1+y2)*(x2-x1)

		y1 = by(i-1)*xstn(i-1)
		y2 = by(i)  *xstn(i)
		sumcf = sumcf + 0.5*(y1+y2)*(x2-x1)

		y1 = by(i-1)**3
		y2 = by(i)**3
		sumix = sumix + 0.5*(y1+y2)*(x2-x1)

		y1 = by(i-1)*xstn(i-1)**2
		y2 = by(i)  *xstn(i)  **2
		sumiy = sumiy + 0.5*(y1+y2)*(x2-x1)
	enddo
	volm = 2.*sumvl
	zcb = sumvm/sumvl
	xcb = sumlm/sumvl

	wparea = 2.*sumwp
	xcf = sumcf/sumwp
	aixx = (2./3.)*sumix
	aiyy = 2.*sumiy - wparea*xcf**2
	bmt = aixx/volm
	bml = aiyy/volm

	write(*,*)' vol=',volm,' kb=',zcb,' lcb=',xcb
	write(*,*)' awp=',wparea,' lcf =',xcf
	write(*,*)' bmt=',bmt,' bml=',bml

c	pause


	do i = 1,nstn

		do j = 1,npt(i)
			xjj(i,j) = -xjt(i,npt(i)-j+1)
			yjj(i,j) =  yjt(i,npt(i)-j+1)
		enddo

		npt1(i) = npt(i) + 1
		npt2(i) = 2.*npt(i) - 1
		do j = npt1(i),npt2(i)
			xjj(i,j) = xjt(i,j-npt(i)+1)
			yjj(i,j) = yjt(i,j-npt(i)+1)
		enddo

		npt(i) = npt2(i)
		nsegg(i) = npt(i) - 1

	enddo

c	write(2,*)' nstn=',nstn
c	do i = 1,nstn
c	write(2,*)' stn no =',stn(i),'xstn=',xstn(i),' njt=',npt(i),
c    &	' nseg=',nsegg(i)
c		do j = 1,npt(i)
c			write(2,*) j,xjj(i,j),yjj(i,j)
c		enddo
c	enddo

	do i = 1,nstn
		do j = 1,nsegg(i)
			cg1(i,j) = 0.5*(xjj(i,j)+xjj(i,j+1))
			cg2(i,j) = 0.5*(yjj(i,j)+yjj(i,j+1))
			ylen = yjj(i,j+1)-yjj(i,j)
			xlen = xjj(i,j+1)-xjj(i,j)
			dell(i,j) = sqrt(xlen**2+ylen**2)
			bn1(i,j) = -ylen/dell(i,j)
			bn2(i,j) =  xlen/dell(i,j)
		enddo
	enddo

c	calculation for N_x

	do i = 1,nstn
		if (i .eq. 1) then
			i1 = 1
			i2 = 2
			i3 = 3
			xgv = xstn(1)
			kstype = 1
		else if (i .eq. nstn) then
			i1 = nstn-2
			i2 = nstn-1
			i3 = nstn
			xgv = xstn(nstn)
			ktype = 2
		else
			i1 = i-1
			i2 = i
			i3 = i+1
			xgv = xstn(i2)
			ktype = 0
		endif

		do j = 1,nsegg(i)/2
c			check if segment horizontal
			if (abs(yjj(i,j+1)-yjj(i,j)) .lt. eps) go to 1
c			panel not horizontal
			ygv = cg1(i,j)
			zgv = cg2(i,j)

			call inter1(ktype,i1,i2,i3,xgv,ygv,zgv,anx)
			bnx(i,j) = anx
			go to 2

    1		continue
c			panel is horizontal
			ygv = cg1(i,j)
			zgv = cg2(i,j)
			call inter2(ktype,i1,i2,i3,xgv,ygv,zgv,anx)
			bnx(i,j) = anx
    2			continue
		enddo

		do j = nsegg(i)/2,nsegg(i)
			bnx(i,j) = bnx(i,nsegg(i)+1-j)
		enddo
	enddo

	return
	end

c-----------------------------------------------------------------------

	subroutine inter1(ktype,i1,i2,i3,xgvn,ygvn,zgvn,anxx)
	implicit real*8 (a-h,o-z)

	parameter(mxs=50,mxp=200)
	common/space1/nstn,npt(mxs),npt1(mxs),npt2(mxs),nsegg(mxs)
	common/space2/al,stn(mxs),xstn(mxs),xjj(mxs,mxp),yjj(mxs,mxp),
     &		      cg1(mxs,mxp),cg2(mxs,mxp),dell(mxs,mxp),
     &		      bn1(mxs,mxp),bn2(mxs,mxp),bnx(mxs,mxp),xt(mxs),
     &		      sa(mxs),by(mxs),sl(mxs),dz(mxs)

	x1 = xstn(i1)
	x2 = xstn(i2)
	x3 = xstn(i3)

c	find ig1,ig2,ig3

	if (ktype .eq. 1) then
		ig1 = i2
		ig2 = i3
	else if (ktype .eq. 2) then
		ig1 = i1
		ig2 = i2
	else
		ig1 = i1
		ig2 = i3
	endif

c	for stns.nos ig1 and ug2, find the ys for zgvn

c	first, we do for stn. ig1

	j = 1
    1	continue
	if ( (zgvn .ge. yjj(ig1,j) .and. zgvn .le. yjj(ig1,j+1) ) .or.
     &	(zgvn .le. yjj(ig1,j) .and. zgvn .ge. yjj(ig1,j+1) ) )  go to 2
	j = j+1
	if (j .gt. npt(ig1)) go to 3
	go to 1
    2	continue
	y1 = xjj(ig1,j)
	y2 = xjj(ig1,j+1)
	z1 = yjj(ig1,j)
	z2 = yjj(ig1,j+1)
	yg1 = (y2-y1)/(z2-z1)*(zgvn-z1) + y1
	go to 4
   3	continue
	yg1 = 0
   4	continue

c	we now find for stn ig2
	j = 1
   11	continue
	if ( (zgvn .ge. yjj(ig2,j) .and. zgvn .le. yjj(ig2,j+1) ) .or.
     &	(zgvn .le. yjj(ig2,j) .and. zgvn .ge. yjj(ig2,j+1) ) ) go to 12
	j = j+1
	if (j .gt. npt(ig2)) go to 13
	go to 11
   12  	continue
	y1 = xjj(ig2,j)
	y2 = xjj(ig2,j+1)
	z1 = yjj(ig2,j)
	z2 = yjj(ig2,j+1)
	yg2 = (y2-y1)/(z2-z1)*(zgvn-z1) + y1
	go to 14
   13   continue
	yg2 = 0
   14	continue

c	here we find the slope
	if (ktype .eq. 1) then
		yy1 = ygvn
		yy2 = yg1
		yy3 = yg2
	else if (ktype .eq. 2) then
		yy1 = yg1
		yy2 = yg2
		yy3 = ygvn
	else
		yy1 = yg1
		yy2 = ygvn
		yy3 = yg2
	endif

	a2 = ((yy2-yy1)/(x2-x1) - (yy3-yy2)/(x3-x2))/(x1-x3)
	a1 = (yy2-yy1)/(x2-x1) - a2*(x1+x2)

	slope = a1 + 2.*a2*xgvn
	anxx = slope/sqrt(1+slope**2)

	return
	end

c---------------------------------------------------------------------

	subroutine inter2(ktype,i1,i2,i3,xgvn,ygvn,zgvn,anxx)
	implicit real*8 (a-h,o-z)

	parameter(mxs=50,mxp=200)
	common/space1/nstn,npt(mxs),npt1(mxs),npt2(mxs),nsegg(mxs)
	common/space2/al,stn(mxs),xstn(mxs),xjj(mxs,mxp),yjj(mxs,mxp),
     &		      cg1(mxs,mxp),cg2(mxs,mxp),dell(mxs,mxp),
     &		      bn1(mxs,mxp),bn2(mxs,mxp),bnx(mxs,mxp),xt(mxs),
     &		      sa(mxs),by(mxs),sl(mxs),dz(mxs)

	x1 = xstn(i1)
	x2 = xstn(i2)
	x3 = xstn(i3)

c	find ig1,ig2,ig3

	if (ktype .eq. 1) then
		ig1 = i2
		ig2 = i3
	else if (ktype .eq. 2) then
		ig1 = i1
		ig2 = i2
	else
		ig1 = i1
		ig2 = i3
	endif

c	for stns.nos ig1 and ug2, find the ys for zgvn

c	first, we do for stn. ig1

	j = 1
    1	continue
	if ( (ygvn .ge. xjj(ig1,j) .and. ygvn .le. xjj(ig1,j+1) ) .or.
     &	(ygvn .le. xjj(ig1,j) .and. ygvn .ge. xjj(ig1,j+1) ) )  go to 2
	j = j+1
	if (j .gt. npt(ig1)) go to 3
	go to 1
    2	continue
	y1 = yjj(ig1,j)
	y2 = yjj(ig1,j+1)
	z1 = xjj(ig1,j)
	z2 = xjj(ig1,j+1)
	yg1 = (y2-y1)/(z2-z1)*(ygvn-z1) + y1
	go to 4
   3	continue
	yg1 = 0
   4	continue

c	we now find for stn ig2
	j = 1
   11	continue
	if ( (ygvn .ge. xjj(ig2,j) .and. ygvn .le. xjj(ig2,j+1) ) .or.
     &	(ygvn .le. xjj(ig2,j) .and. ygvn .ge. xjj(ig2,j+1) ) ) go to 12
	j = j+1
	if (j .gt. npt(ig2)) go to 13
	go to 11
   12  	continue
	y1 = yjj(ig2,j)
	y2 = yjj(ig2,j+1)
	z1 = xjj(ig2,j)
	z2 = xjj(ig2,j+1)
	yg2 = (y2-y1)/(z2-z1)*(ygvn-z1) + y1
	go to 14
   13   continue
	yg2 = 0
   14	continue

c	here we find the slope
	if (ktype .eq. 1) then
		yy1 = zgvn
		yy2 = yg1
		yy3 = yg2
	else if (ktype .eq. 2) then
		yy1 = yg1
		yy2 = yg2
		yy3 = zgvn
	else
		yy1 = yg1
		yy2 = zgvn
		yy3 = yg2
	endif

	a2 = ((yy2-yy1)/(x2-x1) - (yy3-yy2)/(x3-x2))/(x1-x3)
	a1 = (yy2-yy1)/(x2-x1) - a2*(x1+x2)

	slope = a1 + 2.*a2*xgvn
	anxx = slope/sqrt(1+slope**2)

	return
	end

c-----------------------------------------------

	subroutine assembly(ij,anu,x1,y1,x2,y2,a1,b1,a2,b2,
     &  term1,term2,term3,term4,bterm1,bterm2,bterm3,bterm4)

	implicit real*8 (a-h,o-z)
	common/wksp/pi,grav,rho

	gamma = 0.577215664901532860606512
	pi2 = 2.*pi


c	term 1

	ylen = y2-y1
	xlen = x2-x1
	alfai = datan2(ylen,xlen)
	if (alfai .lt. 0.) alfai = alfai + pi2
	blen = b2-b1
	alen = a2-a1
	alfaj = datan2(blen,alen)
	deli = dsqrt(xlen**2 + ylen**2)
	delj = dsqrt(alen**2 + blen**2)
	anxi =  ylen/deli
	anyi = -xlen/deli
	anxj =  blen/delj
	anyj = -alen/delj

	xi = 0.5*(x1+x2)
	yi = 0.5*(y1+y2)

	r1 = dsqrt( (xi-a1)**2 + (yi-b1)**2 )
	r2 = dsqrt( (xi-a2)**2 + (yi-b2)**2 )

	si = sin(alfai)
	sj = sin(alfaj)
	ci = cos(alfai)
	cj = cos(alfaj)
	s11 = si*cj - ci*sj
	c11 = ci*cj + si*sj
c	s11 = sin(alfai - alfaj)
c	c11 = cos(alfai - alfaj)
	alr1 = dlog(r1)
	alr2 = dlog(r2)
	term11 = alr1 - alr2
	ang1 = datan2( (yi-b1),(xi-a1) )
	ang2 = datan2( (yi-b2),(xi-a2) )
	ang = ang1 - ang2
	if (ang .le. -pi) then
		ang = ang + pi2
	else if (ang .gt. pi) then
		ang = ang - pi2
	endif
	term12 = ang
c	write(*,*)' ang1=',ang1*180./pi,' ang2=',ang2*180./pi
c	write(*,*)' ang =', term12*180./pi

	term1 = s11*term11 + c11*term12
	if (ij .eq. 1) term1 = pi

c	cal. for bterm1

	xx1 = a1
	yy1 = b1
	xx2 = a2
	yy2 = b2
	alen = dsqrt( (xx2-xx1)**2 + (yy2-yy1)**2 )
	d1 = ( (xi-xx1)*(xx2-xx1) + (yi-yy1)*(yy2-yy1) )/alen
	d2 = alen - d1
	b = ( -(xi-xx1)*(yy2-yy1) + (yi-yy1)*(xx2-xx1) )/alen
	r1 = dsqrt(d1**2 + b**2)
	r2 = dsqrt(d2**2 + b**2)
	theta = datan(d1/b) + datan(d2/b)
	bterm1 = d1*dlog(r1) + d2*dlog(r2) - alen + b*theta

c	term 2 : image segment

	r1 = dsqrt( (xi-a1)**2 + (yi+b1)**2 )
	r2 = dsqrt( (xi-a2)**2 + (yi+b2)**2 )

	s12 = si*cj + sj*ci
	c12 = ci*cj - si*sj
	alr1 = dlog(r1)
	alr2 = dlog(r2)
	term11 = alr1 - alr2
	ang1 = datan2( (yi+b1),(xi-a1) )
	ang2 = datan2( (yi+b2),(xi-a2) )
	ang = ang1 - ang2
	if (ang .le. -pi) then
		ang = ang + pi2
	else if (ang .gt. pi) then
		ang = ang - pi2
	endif
	term12 = ang

	ta12 = term12
	alr12 = term11

	term2 = s12*term11 + c12*term12

c	bterm2

	xx1 = a1
	yy1 = -b1
	xx2 = a2
	yy2 = -b2
	alen = dsqrt( (xx2-xx1)**2 + (yy2-yy1)**2 )
	d1 = ( (xi-xx1)*(xx2-xx1) + (yi-yy1)*(yy2-yy1) )/alen
	d2 = alen - d1
	b = ( -(xi-xx1)*(yy2-yy1) + (yi-yy1)*(xx2-xx1) )/alen
	r1 = dsqrt(d1**2 + b**2)
	r2 = dsqrt(d2**2 + b**2)
	theta = datan(d1/b) + datan(d2/b)
	bterm2 = d1*dlog(r1) + d2*dlog(r2) - alen + b*theta

c	term4

	c1 = anu*(yi + b1)
	c2 = anu*(yi + b2)
	d1 = anu*(xi - a1)
	d2 = anu*(xi - a2)
	e1 = exp(c1)
	e2 = exp(c2)
	term4 = -s12*(e1*cos(d1) - e2*cos(d2))
     &		+c12*(e1*sin(d1) - e2*sin(d2))

	bterm4 = ( e1*sin(d1-alfaj) - e2*sin(d2-alfaj) )/anu

c	term 3
c	This is the Principal Value (PV) integration

	call pv(xi,yi,a1,b1,anu,pvc1,pvs1)

	call pv(xi,yi,a2,b2,anu,pvc2,pvs2)

	term3 = s12*(pvc1 - pvc2)
     &        - c12*(pvs1 - pvs2)

	bterm3 = ( sin(alfaj)*( alr12 + pvc2 - pvc1 )
     &	       +   cos(alfaj)*( ta12 + pvs1 - pvs2 ) )/anu

	return
	end

c-------------------------------------

	subroutine pv(x,y,a,b,anu,pv1,pv2)
	implicit real*8 (a-h,o-z)

	common/wksp/pi,grav,rho

	complex*8 z,xi,rr

	gamma = 0.577215664901532860606512
	eps = 1.0d-20

	z = dcmplx(x,y)
	xi = dcmplx(a,b)

	rr1 = anu*(y+b)
	rr2 = -anu*(x-a)
	rr = dcmplx(rr1,rr2)

	r = cabs(rr)

	qq3 = a-x
	qq4 = y+b

	if (abs (qq3) .lt. eps) then
			theta = 0.
	else if (qq3 .lt. 0.) then
			theta = datan2(qq3,qq4) + pi
	else if (qq3 .gt. 0.) then
			theta = datan2(qq3,qq4) - pi
	endif

c	calculating c(r,theta) & s(r,theta)

	n = 0
	term = 1
	sum1 = 0.
	sum2 = 0.
    1	n = n+1
	thn = n*theta
	term = term*(r/n)
	sum1 = sum1 + term * cos(thn)/n
	sum2 = sum2 + term * sin(thn)/n
	if (term .gt. eps) go to 1

	crt = gamma + dlog(r) + sum1
	srt = theta          + sum2

	aa = exp(anu*(y+b))
	cc = cos(anu*(x-a))
	ss = sin(anu*(x-a))

	pv1 = aa*(crt*cc + srt*ss)
	pv2 = aa*(crt*ss - srt*cc)

	return
	end

c---------------------------------------------------------------------
      subroutine ludcmp(a,n,np,indx)
c---------------------------------------------------------------------

c	program for matrix soln.
c-------------------------------------------
	implicit real*8 (a-h,o-z)

        parameter (mxp=200)
        complex*8 a(np,np),dum,sum
        dimension indx(mxp),vv(mxp)

	tiny = 1.0d-20

      	d = 1.
      	do i = 1,n
        	aamax = 0.
        	do j = 1,n
        	if (cabs(a(i,j)) .gt. aamax) aamax = cabs(a(i,j))
		enddo
      		if (aamax .eq. 0.) write(*,*) 'singular matrix'
      		vv(i) = 1./aamax
	enddo

      	do j = 1,n
		do i = 1,j-1
			sum = a(i,j)
			do k = 1,i-1
            			sum = sum - a(i,k) * a(k,j)
			enddo
        		a(i,j) = sum
		enddo
        	aamax = 0.0
        	do i = j,n
			sum = a(i,j)
			do k = 1,j-1
				sum = sum - a(i,k) * a(k,j)
			enddo
			a(i,j) = sum
			dum = dcmplx(vv(i) * cabs(sum),0.0)
			if ( cabs(dum) .ge. aamax ) then
				imax = i
				aamax = cabs(dum)
			endif
		enddo
		if ( j .ne. imax ) then
			do k = 1,n
				dum = a(imax,k)
				a(imax,k) = a(j,k)
				a(j,k) = dum
			enddo
			d = -d
			vv(imax) = vv(j)
		endif
		indx(j) = imax
		if (cabs(a(j,j)) .eq. 0.) a(j,j) = tiny

			if ( j .ne. n ) then
				dum = 1./a(j,j)
				do i = j+1,n
					a(i,j) = a(i,j) * dum
c				write(*,*)' i=',i,' j=',j,' a=',a(i,j)
				enddo
			endif
	enddo

	return
	end


c------------------ e n d -----------------------------------------------

c------------------------------------------------------------------------
      subroutine lubksb(a,n,np,indx,b)
c------------------------------------------------------------------------
c	program for matrix soln.
c     solve the system ax = b, where a = lu decomposition of orig. a
c------------------------------------------------------------
	implicit real*8 (a-h,o-z)

        parameter (mxp=200)
	complex*8 a(mxp,mxp),b(mxp),sum
	dimension indx(mxp)

	ii = 0
	do i = 1,n
		ll = indx(i)
		sum = b(ll)
		b(ll) = b(i)
		if ( ii .ne. 0 ) then
			do j = ii,i-1
				sum = sum - a(i,j) * b(j)
			enddo
		else if ( sum .ne. 0. ) then
			ii = i
		endif
		b(i) = sum
	enddo

	do i = n,1,-1
		sum = b(i)
		do j = i+1,n
			sum = sum - a(i,j) * b(j)
		enddo
		b(i) = sum/a(i,i)
	enddo

	return
	end

c-------------------------- e n d -------------------------------------------


	subroutine inter_amass(we,amm,bmm)

	implicit real*8 (a-h,o-z)
	parameter(mxf=100)
	common/space3/weint(mxf),nfr
	common/space4/am0(6,6,mxf),bm0(6,6,mxf)
	dimension amm(18),bmm(18)

c	initialization
	do i = 1,18
		amm(i) = 0.
		bmm(i) = 0.
	enddo

c	interpolation here
c	first find within which frequencies the two are

	k = 2
    2	continue
    	if (we .ge. weint(k-1) .and. we .le. weint(k)) go to 1
	k = k+1
	if (k .lt. nfr) go to 2
    1	continue
	k1 = k-1
	k2 = k
	dx = weint(k2)-weint(k1)
	dy = we - weint(k1)
	dyx = dy/dx

	amm(1) = (am0(1,1,k2)-am0(1,1,k1))*dyx + am0(1,1,k1)
	bmm(1) = (bm0(1,1,k2)-bm0(1,1,k1))*dyx + bm0(1,1,k1)

	amm(2) = (am0(1,3,k2)-am0(1,3,k1))*dyx + am0(1,3,k1)
	bmm(2) = (bm0(1,3,k2)-bm0(1,3,k1))*dyx + bm0(1,3,k1)

	amm(3) = (am0(3,1,k2)-am0(3,1,k1))*dyx + am0(3,1,k1)
	bmm(3) = (bm0(3,1,k2)-bm0(3,1,k1))*dyx + bm0(3,1,k1)

	amm(4) = (am0(1,5,k2)-am0(1,5,k1))*dyx + am0(1,5,k1)
	bmm(4) = (bm0(1,5,k2)-bm0(1,5,k1))*dyx + bm0(1,5,k1)

	amm(5) = (am0(5,1,k2)-am0(5,1,k1))*dyx + am0(5,1,k1)
	bmm(5) = (bm0(5,1,k2)-bm0(5,1,k1))*dyx + bm0(5,1,k1)

	amm(6) = (am0(3,3,k2)-am0(3,3,k1))*dyx + am0(3,3,k1)
	bmm(6) = (bm0(3,3,k2)-bm0(3,3,k1))*dyx + bm0(3,3,k1)

	amm(7) = (am0(3,5,k2)-am0(3,5,k1))*dyx + am0(3,5,k1)
	bmm(7) = (bm0(3,5,k2)-bm0(3,5,k1))*dyx + bm0(3,5,k1)

	amm(8) = (am0(5,3,k2)-am0(5,3,k1))*dyx + am0(5,3,k1)
	bmm(8) = (bm0(5,3,k2)-bm0(5,3,k1))*dyx + bm0(5,3,k1)

	amm(9) = (am0(5,5,k2)-am0(5,5,k1))*dyx + am0(5,5,k1)
	bmm(9) = (bm0(5,5,k2)-bm0(5,5,k1))*dyx + bm0(5,5,k1)

	amm(10) = (am0(2,2,k2)-am0(2,2,k1))*dyx + am0(2,2,k1)
	bmm(10) = (bm0(2,2,k2)-bm0(2,2,k1))*dyx + bm0(2,2,k1)

	amm(11) = (am0(2,4,k2)-am0(2,4,k1))*dyx + am0(2,4,k1)
	bmm(11) = (bm0(2,4,k2)-bm0(2,4,k1))*dyx + bm0(2,4,k1)

	amm(12) = (am0(4,2,k2)-am0(4,2,k1))*dyx + am0(4,2,k1)
	bmm(12) = (bm0(4,2,k2)-bm0(4,2,k1))*dyx + bm0(4,2,k1)

	amm(13) = (am0(4,4,k2)-am0(4,4,k1))*dyx + am0(4,4,k1)
	bmm(13) = (bm0(4,4,k2)-bm0(4,4,k1))*dyx + bm0(4,4,k1)

	amm(14) = (am0(4,6,k2)-am0(4,6,k1))*dyx + am0(4,6,k1)
	bmm(14) = (bm0(4,6,k2)-bm0(4,6,k1))*dyx + bm0(4,6,k1)

	amm(15) = (am0(6,4,k2)-am0(6,4,k1))*dyx + am0(6,4,k1)
	bmm(15) = (bm0(6,4,k2)-bm0(6,4,k1))*dyx + bm0(6,4,k1)

	amm(16) = (am0(2,6,k2)-am0(2,6,k1))*dyx + am0(2,6,k1)
	bmm(16) = (bm0(2,6,k2)-bm0(2,6,k1))*dyx + bm0(2,6,k1)

	amm(17) = (am0(6,2,k2)-am0(6,2,k1))*dyx + am0(6,2,k1)
	bmm(17) = (bm0(6,2,k2)-bm0(6,2,k1))*dyx + bm0(6,2,k1)

	amm(18) = (am0(6,6,k2)-am0(6,6,k1))*dyx + am0(6,6,k1)
	bmm(18) = (bm0(6,6,k2)-bm0(6,6,k1))*dyx + bm0(6,6,k1)

	return
	end

c-------------------------------------------------

	subroutine fkforce(xg2,hd,wnum,fi)

	implicit real*8 (a-h,o-z)
	parameter(mxs=50,mxp=200)

	common/wksp/pi,grav,rho
	common/space1/nstn,npt(mxs),npt1(mxs),npt2(mxs),nsegg(mxs)
	common/space2/al,stn(mxs),xstn(mxs),xjj(mxs,mxp),yjj(mxs,mxp),
     &		      cg1(mxs,mxp),cg2(mxs,mxp),dell(mxs,mxp),
     &		      bn1(mxs,mxp),bn2(mxs,mxp),bnx(mxs,mxp),xt(mxs),
     &		      sa(mxs),by(mxs),sl(mxs),dz(mxs)
	complex*8 term,summ,sum1,sum2,sum3,sum4,y1,y2,tt,t(mxs)
	complex*8 fis1(mxs),fis2(mxs),fis3(mxs),fis4(mxs)
	complex*8 fi(6)

	sinmu = sin(hd)
	cosmu = cos(hd)
c	write(2,*)' sinmu=',sinmu,' cosmu=',cosmu
	rhog = rho*grav
c	write(2,*)' rhog=',rhog
c	initialize
	do i = 1,mxs
		fis1(i) = (0.0,0.0)
		fis2(i) = (0.0,0.0)
		fis3(i) = (0.0,0.0)
		fis4(i) = (0.0,0.0)
	enddo

	do ks = 1,nstn
		nseg = nsegg(ks)
c		write(2,*)' nseg=',nseg
		sum1 = (0.0,0.0)
		sum2 = (0.0,0.0)
		sum3 = (0.0,0.0)
		sum4 = (0.0,0.0)
		do j = 1,nseg
c			write(2,*)' j=',j
c****		check here for axis system***
			y = -cg1(ks,j)
			z =  cg2(ks,j)
			an1 =  bnx(ks,j)
			an2 = -bn1(ks,j)
			an3 =  bn2(ks,j)
c*** check here
			an4 = y*an3 - (z-xg2)*an2
			q1 = wnum*z
			q2 = -wnum*y*sinmu
			tt = dcmplx(q1,q2)
			term = exp(tt)
			sum1 = sum1 + an1*term*dell(ks,j)
			sum2 = sum2 + an2*term*dell(ks,j)
			sum3 = sum3 + an3*term*dell(ks,j)
			sum4 = sum4 + an4*term*dell(ks,j)
		enddo
		fis1(ks) = rhog*sum1
		fis2(ks) = rhog*sum2
		fis3(ks) = rhog*sum3
		fis4(ks) = rhog*sum4
c		write(2,*)' fis1=',fis1(ks),' fis2=',fis2(ks)
c		write(2,*)' fis3=',fis3(ks),' fis4=',fis4(ks)
c		write(2,*)' k=',ks,' fis=',fis3(ks)
	enddo

c	integrate to get total for ship
	do ks = 1,nstn
		xx = xt(ks)
		q1 = 0.0
		q2 = -wnum*xx*cosmu
		tt = dcmplx(q1,q2)
		t(ks) = exp(tt)
	enddo

	do i = 1,6
		summ = (0.0,0.0)
		do ks = 2,nstn
			x1 = xt(ks-1)
			x2 = xt(ks)
			if (i .eq. 1) then
				y1 = t(ks-1)*fis1(ks-1)
				y2 = t(ks)*fis1(ks)
			else if (i .eq. 2) then
				y1 = t(ks-1)*fis2(ks-1)
				y2 = t(ks)*fis2(ks)
			else if (i .eq. 3) then
				y1 = t(ks-1)*fis3(ks-1)
				y2 = t(ks)*fis3(ks)
			else if (i .eq. 4) then
				y1 = t(ks-1)*fis4(ks-1)
				y2 = t(ks)*fis4(ks)
			else if (i .eq. 5) then
				y1 = -t(ks-1)*x1*fis3(ks-1)
				y2 = -t(ks)*x2*fis3(ks)
			else if (i .eq. 6) then
				y1 =  t(ks-1)*x1*fis2(ks-1)
				y2 =  t(ks)*x2*fis2(ks)
			endif
			summ = summ + 0.5*(x2-x1)*(y1+y2)
		enddo
		fi(i) = summ
c		write(2,*)' i=',i,' fi=',fi(i)
	enddo

	return
	end

c--------------------------------------------------------

	subroutine fdforce(xg2,hd,w0,we,wnum,u0,fd)

	implicit real*8 (a-h,o-z)
	parameter(mxs=50,mxp=200,mxf=100)

	common/wksp/pi,grav,rho
	common/space1/nstn,npt(mxs),npt1(mxs),npt2(mxs),nsegg(mxs)
	common/space2/al,stn(mxs),xstn(mxs),xjj(mxs,mxp),yjj(mxs,mxp),
     &		      cg1(mxs,mxp),cg2(mxs,mxp),dell(mxs,mxp),
     &		      bn1(mxs,mxp),bn2(mxs,mxp),bnx(mxs,mxp),xt(mxs),
     &		      sa(mxs),by(mxs),sl(mxs),dz(mxs)
	common/space3/weint(mxf),nfr
	common/space5/phi1(mxs,mxf,mxp),phi2(mxs,mxf,mxp),
     &		      phi3(mxs,mxf,mxp),phi4(mxs,mxf,mxp)

	complex*8 term1,term2,term3
	complex*8 summ,sum1,sum2,sum3,sum4,y1,y2,tt,t(mxs)
	complex*8 fds1(mxs),fds2(mxs),fds3(mxs),fds4(mxs)
	complex*8 fd(6),ph1(mxp),ph2(mxp),ph3(mxp),ph4(mxp)
	complex*8 phi1,phi2,phi3,phi4

	sinmu = sin(hd)
	cosmu = cos(hd)
	rhow = rho*w0

	do ks = 1,nstn
		nseg = nsegg(ks)
		sum1 = (0.0,0.0)
		sum2 = (0.0,0.0)
		sum3 = (0.0,0.0)
		sum4 = (0.0,0.0)
c	here we have to interpolate in freq. band the sectional potentials
c	for the given frequency
c		write(2,*)' ks=',ks,' we =',we
		call inter_pot(ks,nseg,we,ph1,ph2,ph3,ph4)
c		do j = 1,nseg
c			write(2,*) j,real(ph2(j)),aimag(ph2(j))
c		enddo
		do j = 1,nseg
c****		check here for axis system***
			y = -cg1(ks,j)
			z =  cg2(ks,j)
			an1 =  bnx(ks,j)
			an2 = -bn1(ks,j)
			an3 =  bn2(ks,j)
c*** check here
			an4 = y*an3 - (z-xg2)*an2
			q1 = wnum*z
			q2 = -wnum*y*sinmu
			tt = dcmplx(q1,q2)
			term1 = exp(tt)
			term2 = dcmplx((an1*cosmu+an2*sinmu),an3)
			term3 = term1*term2
			sum1 = sum1 + term3*ph1(j)*dell(ks,j)
			sum2 = sum2 + term3*ph2(j)*dell(ks,j)
			sum3 = sum3 + term3*ph3(j)*dell(ks,j)
			sum4 = sum4 + term3*ph4(j)*dell(ks,j)
		enddo
		fds1(ks) = rhow*sum1
		fds2(ks) = rhow*sum2
		fds3(ks) = rhow*sum3
		fds4(ks) = rhow*sum4
c		write(2,*)' fds1=',fds1(ks),' fds2=',fds2(ks)
c		write(2,*)' fds3=',fds3(ks),' fds4=',fds4(ks)
	enddo

c	integrate to get total for ship
c	integrate to get total for ship
	do ks = 1,nstn
		xx = xt(ks)
		q1 = 0.0
		q2 = -wnum*xx*cosmu
		tt = dcmplx(q1,q2)
		t(ks) = exp(tt)
	enddo

	do i = 1,6
		summ = (0.0,0.0)
		do ks = 2,nstn
			x1 = xt(ks-1)
			x2 = xt(ks)
			if (i .eq. 1) then
				y1 = t(ks-1)*fds1(ks-1)
				y2 = t(ks)*fds1(ks)
			else if (i .eq. 2) then
				y1 = t(ks-1)*fds2(ks-1)
				y2 = t(ks)*fds2(ks)
			else if (i .eq. 3) then
				y1 = t(ks-1)*fds3(ks-1)
				y2 = t(ks)*fds3(ks)
			else if (i .eq. 4) then
				y1 = t(ks-1)*fds4(ks-1)
				y2 = t(ks)*fds4(ks)
			else if (i .eq. 5) then
				y1=-t(ks-1)*dcmplx(x1,-u0/we)*fds3(ks-1)
				y2=-t(ks)*dcmplx(x2,-u0/we)*fds3(ks)
			else if (i .eq. 6) then
				y1= t(ks-1)*dcmplx(x1,-u0/we)*fds2(ks-1)
				y2= t(ks)*dcmplx(x2,-u0/we)*fds2(ks)
			endif
			summ = summ + 0.5*(x2-x1)*(y1+y2)
		enddo
		fd(i) = summ
		if (i .eq. 2) fd(i) = -fd(i)
		if (i .eq. 4) fd(i) = -fd(i)
		if (i .eq. 4) fd(i) = -fd(i)
	enddo

	return
	end

c-------------------------------------------------------

	subroutine inter_pot(ks,nseg,we,ph1,ph2,ph3,ph4)
	implicit real*8 (a-h,o-z)
	parameter(mxs=50,mxp=200,mxf=100)

	complex*8 ph1(mxp),ph2(mxp),ph3(mxp),ph4(mxp)
	complex*8 phi1,phi2,phi3,phi4

	common/space3/weint(mxf),nfr
	common/space5/phi1(mxs,mxf,mxp),phi2(mxs,mxf,mxp),
     &		      phi3(mxs,mxf,mxp),phi4(mxs,mxf,mxp)


	do i = 1,nseg

	k = 2
    2	continue
	if (we .ge. weint(k-1) .and. we  .le. weint(k)) go to 1
	k = k+1
	if (k .lt. nfr) go to 2
    1	continue
	k1 = k-1
	k2 = k
	dx = weint(k2)-weint(k1)
	dy = we - weint(k1)
	dyx = dy/dx

	ph1(i) = (phi1(ks,k2,i)-phi1(ks,k1,i))*dyx + phi1(ks,k1,i)
	ph2(i) = (phi2(ks,k2,i)-phi2(ks,k1,i))*dyx + phi2(ks,k1,i)
	ph3(i) = (phi3(ks,k2,i)-phi3(ks,k1,i))*dyx + phi3(ks,k1,i)
	ph4(i) = (phi4(ks,k2,i)-phi4(ks,k1,i))*dyx + phi4(ks,k1,i)

	enddo

	return
	end

c------------------------------------------------------

      	subroutine invert1(aa)
      	implicit real*8 (a-h,o-z)
      	dimension m(6)
      	complex*8 aa(6,6),cc(6),d,temp,de

	nn = 6
	n = 6
      	de=(1.0e0,0.0e0)
  600 	do 610 i=1,nn
      	m(i)=-i
  610 	continue
      	do 620 i=1,nn
      	x=0.0e0
      	do 630 l=1,nn
      	if(m(l) .gt. 0) go to 630
      	do 640 k=1,nn
      	if (m(k) .gt. 0) go to 640
      	d=aa(l,k)
	y = cabs(d)
      	if (x .gt. y) go to 640
      	ld=l
      	kd=k
      	x=y
  640 	continue
  630 	continue
      	d=aa(ld,kd)
      	de=d
      	l=-m(ld)
      	m(ld)=m(kd)
      	m(kd)=l
      	do 660 j=1,nn
      	cc(j)=aa(ld,j)
      	aa(ld,j)=aa(kd,j)
  660 	aa(kd,j)=cc(j)
      	do 670 k=1,nn
      	aa(k,kd)=aa(k,kd)/d
  670 	continue
      	do 700 j=1,n
      	if (j .eq. kd) go to 700
      	do 710 k=1,nn
      	aa(k,j)=aa(k,j)-cc(j)*aa(k,kd)
  710 	continue
  700 	continue
      	cc(kd)=(-1.0e0,0.0e0)
      	do 780 k=1,nn
      	aa(kd,k)=-cc(k)/d
  780 	continue
  620 	continue
      	do 840 i=1,nn
      	l=0
  820 	l=l+1
      	if (m(l) .ne. i) go to 820
      	m(l)=m(i)
      	m(i)=i
      	do 840 k=1,nn
      	temp=aa(k,l)
      	aa(k,l)=aa(k,i)
  840 	aa(k,i)=temp
      	det=cabs(de)
  900 	continue

	return
	end
c---------------------------------------------------------

	subroutine inter_sec_mass(nstn,we,as33,bs33)

	implicit real*8 (a-h,o-z)
	parameter(mxs=50,mxf=100)
	common/space3/weint(mxf),nfr
	common/space6/aa33(mxs,mxf),bb33(mxs,mxf)
	dimension as33(mxs),bs33(mxs)

c	write(*,*)' we=',we,' nstn=',nstn
c	initialization
	do i = 1,nstn
		as33(i) = 0.
		bs33(i) = 0.
	enddo
c	write(*,*)' we are here'

c	interpolation here
c	first find within which frequencies the two are

	k = 2
    2	continue
    	if (we .ge. weint(k-1) .and. we .le. weint(k)) go to 1
	k = k+1
	if (k .lt. nfr) go to 2
    1	continue
	k1 = k-1
	k2 = k
	dx = weint(k2)-weint(k1)
	dy = we - weint(k1)
	dyx = dy/dx

	do i = 1,nstn
		as33(i) = (aa33(i,k2) - aa33(i,k1))*dyx + aa33(i,k1)
		bs33(i) = (bb33(i,k2) - bb33(i,k1))*dyx + bb33(i,k1)
	enddo

	return
	end

c-------------------------------------------------

	subroutine slope(nstn,xt,as33,asl33)
	implicit real*8 (a-h,o-z)
	parameter(mxs=50)
	dimension xt(mxs),as33(mxs),asl33(mxs)

	do i = 1,nstn
		if (i .eq. 1) then
			x1 = xt(1)
			x2 = xt(2)
			x3 = xt(3)
			y1 = as33(1)
			y2 = as33(2)
			y3 = as33(3)
			xg = x1
		else if (i .eq. nstn) then
			x1 = xt(nstn-2)
			x2 = xt(nstn-1)
			x3 = xt(nstn)
			y1 = as33(nstn-2)
			y2 = as33(nstn-1)
			y3 = as33(nstn)
			xg = x3
		else
			x1 = xt(i-1)
			x2 = xt(i)
			x3 = xt(i+1)
			y1 = as33(i-1)
			y2 = as33(i)
			y3 = as33(i+1)
			xg = x2
		endif

		a2 = ((y2-y1)/(x2-x1) - (y3-y2)/(x3-x2))/(x1-x3)
		a1 = (y2-y1)/(x2-x1) - a2*(x1+x2)
		slpe = a1 + 2.*a2*xg
		asl33(i) = slpe
c		write(23,*) i,xt(i),asl33(i)
	enddo

	return
	end



c---------------------------------------------------------------
	subroutine rawgb(we,wnum,u0,bs33,asl33,r3,r5,e3,e5,wamp,raw)

	implicit real*8 (a-h,o-z)
	parameter(mxs=50,mxp=200)
	common/space1/nstn,npt(mxs),npt1(mxs),npt2(mxs),nsegg(mxs)
	common/space2/al,stn(mxs),xstn(mxs),xjj(mxs,mxp),yjj(mxs,mxp),
     &		      cg1(mxs,mxp),cg2(mxs,mxp),dell(mxs,mxp),
     &		      bn1(mxs,mxp),bn2(mxs,mxp),bnx(mxs,mxp),xt(mxs),
     &		      sa(mxs),by(mxs),sl(mxs),dz(mxs)
	dimension bs33(mxs),asl33(mxs),ord(mxs)

c	write(*,*)' r3=',r3,' r5=',r5,' e3=',e3,' e5=',e5

	eps = 1.0e-08
	vel = u0

	do i = 1,nstn

c	find relative velocity

c	first we find the attenuation factor for inc. wave
		sum = 0.
		do j = npt1(i)-1,npt2(i)-1
			y1 = xjj(i,j)*exp(wnum*yjj(i,j))
			y2 = xjj(i,j+1)*exp(wnum*yjj(i,j+1))
			x1 = yjj(i,j)
			x2 = yjj(i,j+1)
			sum = sum + 0.5*(y1+y2)*(x2-x1)
		enddo
		if (by(i) .lt. eps) then
			fact = 1.
		else
			fact = 1. - (wnum/by(i)) * sum
		endif

		xb = xt(i)

		term1 = - we*r3*cos(e3) + we*xb*r5*cos(e5)
     &		- vel*r5*sin(e5) + fact*wamp*we*cos(wnum*xb*cos(wang))
		term2 = -we*r3*sin(e3) + we*xb*r5*sin(e5)
     &		+ vel*r5*cos(e5) - fact*wamp*we*sin(wnum*xb*cos(wang))

		vzamp = dsqrt(term1**2 + term2**2)

		ord(i) = (bs33(i) - u0*asl33(i)) * (vzamp**2)
c		ord(i) = (bs33(i) + u0*asl33(i)) * (vzamp**2)

	enddo

   22	format(i6,7f11.4)

c	summing for added resistance

	sum = 0.
	do i = 2,nstn
		sum = sum + 0.5*(ord(i)+ord(i-1))*(xt(i)-xt(i-1))
	enddo

	raw = wnum/(2.*we) * sum

	return
	end

c-------------------------


	subroutine rawsv(w0,we,wnum,u0,as33,bs33,rca3,rca5,wamp,raw)
	implicit real*8 (a-h,o-z)
	parameter(mxs=50,mxp=200)
	common/wksp/pi,grav,rho
	common/space1/nstn,npt(mxs),npt1(mxs),npt2(mxs),nsegg(mxs)
	common/space2/al,stn(mxs),xstn(mxs),xjj(mxs,mxp),yjj(mxs,mxp),
     &		      cg1(mxs,mxp),cg2(mxs,mxp),dell(mxs,mxp),
     &		      bn1(mxs,mxp),bn2(mxs,mxp),bnx(mxs,mxp),xt(mxs),
     &		      sa(mxs),by(mxs),sl(mxs),dz(mxs)
	dimension as33(mxs),bs33(mxs),ord7(mxs)
	complex*8 rca3,rca5,t3,t5,tt,et,sum1,sum2,f3,f5,rw
	complex*8 ord3(mxs),ord5(mxs)

	eps = 1.0e-08

	do i = 1,nstn
		c33 = 2.*rho*grav*by(i)
		t31 = c33 - w0*we*as33(i)
		t32 = w0*bs33(i)
		t3 = dcmplx(t31,t32)
		t51 = xt(i)*c33 - xt(i)*w0*we*as33(i)-u0*bs33(i)*w0/we
		t52 = xt(i)*w0*bs33(i) - u0*w0*as33(i)
		t5 = dcmplx(t51,t52)

		if (by(i) .lt. eps) then
			tt1 = 0.
			tt7 = 0.
		else
			tt1 = wnum*sa(i)/by(i)
			tt7 = -2.*wnum*sa(i)/by(i)
			tt7 = -2.*tt1
		endif
		tt2 = wnum*xt(i)
		tt = dcmplx(tt1,tt2)
		et = exp(-tt)
		ord3(i) = et * t3
		ord5(i) = et * t5
		ord7(i) = exp(tt7) * bs33(i)

	enddo
  234	format(i5,9f10.4)

c	summing
	sum1 = (0.0,0.0)
	sum2 = (0.0,0.0)
	sum7 = 0.0

	do i = 2,nstn
		dx = 0.5*(xt(i) - xt(i-1))
		sum1 = sum1 + (ord3(i)+ord3(i-1))*dx
		sum2 = sum2 + (ord5(i)+ord5(i-1))*dx
		sum7 = sum7 + (ord7(i)+ord7(i-1))*dx
	enddo
	f3 =  wamp*sum1
	f5 = -wamp*sum2
	r7 = 0.5*(wamp**2)*wnum*((w0**2)/we) * sum7

	rw = 0.5*dcmplx(0.0,1.0)*wnum*(rca3*f3 + rca5*f5) + r7
	raw = cabs(rw)

	return
	end

c--------------------------------------------------------------------------

	subroutine rawbs(we,wnum,volm,ra3,ra5,e3,e5,wamp,raw)
	implicit real*8 (a-h,o-z)
	parameter(mxs=50,mxp=200)

	common/wksp/pi,grav,rho
	common/space1/nstn,npt(mxs),npt1(mxs),npt2(mxs),nsegg(mxs)
	common/space2/al,stn(mxs),xstn(mxs),xjj(mxs,mxp),yjj(mxs,mxp),
     &		      cg1(mxs,mxp),cg2(mxs,mxp),dell(mxs,mxp),
     &		      bn1(mxs,mxp),bn2(mxs,mxp),bnx(mxs,mxp),xt(mxs),
     &		      sa(mxs),by(mxs),sl(mxs),dz(mxs)
	dimension rmot(mxs)

	do i = 1,nstn
		r1 = ra3*cos(e3)-xt(i)*ra5*cos(e5)-wamp*cos(wnum*xt(i))
		r2 = ra3*sin(e3)-xt(i)*ra5*sin(e5)+wamp*sin(wnum*xt(i))
		rmot(i) = dsqrt(r1**2+r2**2)
c		write(*,*)' i=',i,' rmot=',rmot(i),' sl=',sl(i)
	enddo

	sum = 0.
	do i = 2,nstn
		y1 = (rmot(i-1)**2) * sl(i-1)
		y2 = (rmot(i)  **2) * sl(i)
		x1 = xt(i-1)
		x2 = xt(i)
		sum = sum + 0.5*(y1+y2)*(x2-x1)
	enddo

c	term1 = 0.5*rho*grav*sum
	term1 = - 0.5*rho*grav*sum
	term2 =   0.5*rho*volm*we**2*ra3*ra5*cos(e3-e5)

c	write(*,*)' ra3=',ra3,' ra5=',ra5*180./pi,
c     &	' e3=',e3*180/pi,' e5=',e5*180/pi
c	write(*,*)' we=',we,' cos=',cos(e3-e5)
c	write(*,*)' term1=',term1,' term2=',term2

	raw = term1 + term2

	return
	end

c------------------------------------------------

	subroutine interp(nstn,xstn,sa,by,dz,xgiven,sag,byg,dzg)
	implicit real*8 (a-h,o-z)
	parameter(mxs=50)
	dimension xstn(mxs),sa(mxs),by(mxs),dz(mxs)

	if (xgiven .le. xstn(1) ) go to 1
	j = 2
    3	continue
	if (xgiven .ge. xstn(j-1) .and. xgiven .le. xstn(j)) go to 2
	j = j+1
	if (j .lt. nstn) go to 3
    2	continue
	j1 = j-1
	j2 = j
	go to 5
    1	j1 = 1
	j2 = 2
    5	continue

	dx = xstn(j2) - xstn(j1)
	dy = xgiven - xstn(j1)
	dyx = dy/dx
	sag = (sa(j2)-sa(j1))*dyx + sa(j1)
	byg = (by(j2)-by(j1))*dyx + by(j1)
	dzg = (dz(j2)-dz(j1))*dyx + dz(j1)

	return
	end

c------------------------------------------------------------------

	subroutine roll_damp_fric(al,br,dr,cb,og,omega,u0,anue,bf)
	implicit real*8 (a-h,o-z)
	common/wksp/pi,grav,rho

	pi2 = 2.*pi
	u = u0

c	friction damping
	s = al*(1.7*dr + cb*br)
	r = ( (0.887+0.145*cb)*(s/al) - 2.*og )/pi

	bf0 = 0.787*rho*s*(r**2)*sqrt(omega*anue)
	vfact = 1.+4.1*u/(omega*al)
	bf = bf0*vfact

c	write(*,*)' bf0=',bf0,' bf=',bf

	return
	end

c------------------------------------------------------

	subroutine roll_damp_lift(al,br,dr,cm,og,u0,bl)
	implicit real*8 (a-h,o-z)
	common/wksp/pi,grav,rho

	pi2 = 2.*pi
	u = u0

c	lift damping
	if (cm .le. 0.92) akapa = 0.
	if (cm .le. 0.97 .and. cm .gt. 0.92) akapa = 0.1
	if (cm .gt. 0.97) akapa = 0.3
	akn = pi2 * (dr/al) + akapa * (4.1*br/al - 0.045)
	al0 = 0.3*dr
	alr = 0.5*dr
c	write(*,*)'cm=',cm,' kapa =',akapa,' u=',u
c	write(*,*)'kn=',akn,' l0=',al0,' lr=',alr,' og=',og

	bl = 0.5*rho*u*al*dr*akn*al0*alr
     &	   * ( 1. - 1.4*og/alr + 0.7*og**2/(al0*alr) )
c	write(*,*)' bl=',bl

	return
	end

c------------------------------------------------------

	subroutine roll_damp_bk(xbk,bybk,dzbk,sabk,br,dr,og,omega,
     &	bk,delbk,theta,bbkn,bbkh)
	implicit real*8 (a-h,o-z)

	dimension xbk(11),bybk(11),dzbk(11),sabk(11),b1(11),b2(11)
	common/wksp/pi,grav,rho

	do i = 1,11
		by = bybk(i)
		dz = dzbk(i)
		sa = sabk(i)
		ogd = og/dr
		ca = sa/(by*dz)
		h0 = by/dz
c		bilge circle radius
		rr = 2.*dz*sqrt(h0*(ca-1)/(pi-4.))
		rd = rr/dz
		if (h0 .le. 1. .and. rd .ge. h0) rr = by
		if (h0 .gt. 1. .and. rd .ge. 1.) rr = dz
		rd = rr/dz

		f = 1. + 0.3*exp(-160.*(1.-ca))
		r = dz * sqrt( h0-0.292893*rd + (1.-ogd-0.292893*rd)**2)

		am1 = rd
		am2 = ogd
		am3 = 1.0 - am1 - am2
		am4 = h0 - am1
		am5 = (0.414*h0 + 0.0651*am1**2-(0.382*h0+0.0106)*am1)/
     &	      	      ((h0-0.215*am1)*(1.0-0.215*am1))
		am6 = (0.414*h0 + 0.0651*am1**2-(0.382+0.0106*h0)*am1)/
     &	      		((h0-0.215*am1)*(1.0-0.215*am1))
		s0 = 0.3*(pi*f*r*theta) + 1.95*bk
		am7 = s0/dz - 0.25*pi*am1
		r1 = 0.25*pi*rr
		if (s0 .lt. r1) am7 = 0.
		am8 = am7 + 0.414*am1
		if (s0 .lt. r1) am8 = am7 + 1.414214*(1.-cos(s0/rr))*am1

		a = (am3+am4)*am8 - am7**2
		b = am4**3/(3.*h0-0.215*am1) + (1-am1**2)*(2.*am3-am2)
     &	    	    /(6.*(1.-0.215*am1)) + am1*(am3*am5 + am4*am6)
		cpplus = 1.2
		cpmins = 0.
		if (theta .ne. 0.) cpmins = -22.5*bk/(pi*r*f*theta) - 1.2
		cd = cpplus - cpmins

		ai = -a*cpmins + b*cpplus
		ratio = 2.*r*bk*cd/(dz**2 * ai)

		b1(i) = (8./(3.*pi))*rho*(r**2)*(bk**2)*omega*(f**2)*
     &		       (22.5/(pi*f)+2.40*r*theta/bk)

		b2(i) = (4./(3.*pi))*rho*(r**2)*(dz**2)
     &		        *omega*theta*ai*(f**2)

	enddo

c	summing up over length

	bbkn = (delbk/3.)*(b1(1) + 4.*b1(2) + 2.*b1(3) + 4.*b1(4) +
     &		        2.*b1(5) + 4.*b1(6) + 2.*b1(7) + 4.*b1(8) +
     &		        2.*b1(9) + 4.*b1(10)+    b1(11))

	bbkh = (delbk/3.)*(b2(1) + 4.*b2(2) + 2.*b2(3) + 4.*b2(4) +
     &		        2.*b2(5) + 4.*b2(6) + 2.*b2(7) + 4.*b2(8) +
     &		        2.*b2(9) + 4.*b2(10)+    b2(11))

c	write(*,*)' bbkn =',bbkn,' bbkh =',bbkh

	return
	end

c--------------------------------------------------------------

	subroutine roll_damp_ed(xbk,bybk,dzbk,sabk,al,br,dr,og,omega,
     &	u0,delbk,theta,be)
	implicit real*8 (a-h,o-z)

	dimension xbk(11),bybk(11),dzbk(11),sabk(11),cr(11)
	common/wksp/pi,grav,rho

	do i = 1,11

	by = bybk(i)
	dz = dzbk(i)
	sa = sabk(i)
	ogd = og/dr
	ca = sa/(by*dz)
	h0 = by/dz

	sig = ca
	ah0 = h0/(1.-ogd)
	sigma = (sig-ogd)/(1.-ogd)
	e = (ah0-1)/(ah0+1)
	e2 = e**2
	a = 4.*sigma*(1-e2)/pi + e2
	o = -a/(a+3)
	o2 = sqrt(o**2 - (a-1.)/(a+3.))
	a3 = o + o2
	a1 = e*(1+a3)
	am = 0.5*by/(1.+a1+a3)
	am = by/(1.+a1+a3)
	aa1 = 0.25*a1*(1.+a3)/a3
	if (aa1 .gt. 1.0) aa1 = 1.0
	if (aa1 .lt. -1.0) aa1 = -1.0

	psi1 = 0.0
	psi2 = 0.5*acos(aa1)

	ss1 = sin(psi1)
	cs1 = cos(psi1)
	s3s1 = sin(3.*psi1)
	c3s1 = cos(3.*psi1)
	ss2 = sin(psi2)
	cs2 = cos(psi2)
	s3s2 = sin(3.*psi2)
	c3s2 = cos(3.*psi2)

	rmax1 = am*sqrt(((1+a1)*ss1-a3*s3s1)**2+((1-a1)*cs1+a3*c3s1)**2)
	rmax2 = am*sqrt(((1+a1)*ss2-a3*s3s2)**2+((1-a1)*cs2+a3*c3s2)**2)

	if (rmax1 .ge. rmax2) then
		psi = psi1
	else if (rmax1 .lt. rmax2) then
		psi = psi2
	endif
	ss = sin(psi)
	cs = cos(psi)
	s3s = sin(3.*psi)
	c3s = cos(3.*psi)
	rmax = am*sqrt(((1+a1)*ss-a3*s3s)**2+((1-a1)*cs+a3*c3s)**2)

	c2s = cos(2.*psi)
	c4s = cos(4.*psi)
	c5s = cos(5.*psi)
	s2s = sin(2.*psi)
	s4s = sin(4.*psi)
	s5s = sin(5.*psi)

	hh = 1. + a1**2 + 9.*a3**2 + 2.*a1*(1-3.*a3)*c2s - 6.*a3*c4s
	aa = -2.*a3*c5s + a1*(1.-a3)*c3s +
     &	      ( (6.-3.*a1)*a3**2 + (a1**2 - 3.*a1)*a3 + a1**2 ) * cs
	bb = -2.*a3*s5s + a1*(1.-a3)*s3s +
     &	      ( (6.+3.*a1)*a3**2 + (3.*a1 + a1**2)*a3 + a1**2 ) * ss

	f3 = 1. + 4.*exp( -165000.*(1.-sigma)**2 )

	rmean = 2.*dz*(1.-ogd)*sqrt(ah0*sigma/pi)
c	write(*,*)' rmean =',rmean
	vmax = 2.*am*sqrt(aa**2 + bb**2) / hh

	gama = f3*(rmax+vmax)/rmean

	cp = 0.5*(0.87*exp(-gama) - 4.*exp(-0.187*gama) + 3)
	f1 = 0.5*( 1.+tanh( 20.*(sig -0.7) ) )
	f2 = 0.5*( 1.-cos(pi*sig) ) - 1.5*( 1.-exp(-5.*(1-sig)) )
     &		*(sin(pi*sig))**2

c	bilge circle radius
	rr = 2.*dz*sqrt(h0*(ca-1)/(pi-4.))
	rd = rr/dz
	if (h0 .le. 1. .and. rd .ge. h0) rr = by
	if (h0 .gt. 1. .and. rd .ge. 1.) rr = dz
	rd = rr/dr

	cr(i) = (rmax/dz)**2*cp*((1.-f1*rd)
     &	      * (1.-ogd-f1*rd)+f2*(h0-f1*rd)**2)
	cr(i) = (rmax/dz)**2*cp*((1.-f1*rd)
     &	      *(1.-ogd)+f2*(h0-f1*rd)**2)

	enddo

c	integrate for cr(i)

c	summing up over length

	crr = (delbk/3.)*( cr(1) + 4.*cr(2) + 2.*cr(3) + 4.*cr(4) +
     &		        2.*cr(5) + 4.*cr(6) + 2.*cr(7) + 4.*cr(8) +
     &		        2.*cr(9) + 4.*cr(10)+    cr(11) )

	be0 = 4.*rho*dr**4*omega*theta/(3.*pi) * crr

	term = 0.04 * omega**2 * al**2
	factsp = term/(u0**2+term)
	be = be0 * factsp
c	write(*,*)' be0=',be0,' factsp=',factsp,' be=',be

	return
	end

c----------------------------

