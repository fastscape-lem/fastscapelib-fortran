program flexure_test

	! test example to demonstrate how to use the flexure routine
	! see relevant section below for explanations

	implicit none

	integer :: nx,ny,istep,nstep,nn,nfreq,ibc,nreflector
	double precision, dimension(:), allocatable :: u,h,b,hp,h0,rhos,p,x,y,kf,kd
	real :: time_in,time_out
	double precision :: kfsed,m,n,kdsed,g,expp
	double precision xl,yl,dx,dy,dt,pi,vex
	double precision rhoa, eet, tflux, eflux, bflux

	integer i,j,ij

	nx = 201
	ny = 401
	nn = nx*ny

	pi=atan(1.d0)*4.d0

	vex = 50.d0

	call FastScape_Init ()
	call FastScape_Set_NX_NY (nx,ny)
	call FastScape_Setup ()

	xl=200.d3
	yl=400.d3
	dx=xl/(nx-1)
	dy=yl/(ny-1)
	call FastScape_Set_XL_YL (xl,yl)

	allocate (x(nx*ny),y(nx*ny))
	x = (/((xl*float(i-1)/(nx-1), i=1,nx),j=1,ny)/)
	y = (/((yl*float(j-1)/(ny-1), i=1,nx),j=1,ny)/)

	dt=20.d3
	call FastScape_Set_DT (dt)

	allocate (kf(nn),kd(nn))

	kf=1.d-5
	kfsed=1.d-5
	m=0.4d0
	n=1.d0
	kd=1.d-2
	kdsed=1.d-2
	kd=0.d0
	kdsed=0.d0
	g=1.d0
	expp=-1.d0
	!expp=10.d0
	call FastScape_Set_Erosional_Parameters (kf,kfsed,m,n,kd,kdsed,g,g,expp)

	ibc=1000
	call FastScape_Set_BC (ibc)

	allocate (h(nn),b(nn),u(nn),hp(nn),h0(nn),rhos(nn),p(nn))
	call random_number (h)
	!where (y.lt.3.d0*yl/4.d0) h = h + 10.d0
	where (y.gt.2.d0*yl/4.d0) h = h + 100.d0
	h=h+1.*cos(x/xl*2.d0*pi)
	!h=h+10.*y/yl

	call FastScape_Init_H (h)

	nstep = 1000
	nfreq = 10
	call FastScape_View()
	istep = nstep

	rhoa = 3250.d0
	eet = 10.d3
	rhos = 2400.d0

	u = 3.d-4
	where (y.lt.2.d0*yl/4.d0) u =0.d0
	call FastScape_Set_U (u)

	p = 1.d0
	call FastScape_Set_Precip (p)

	nreflector = 5
	call FastScape_Strati (nstep, nreflector, nfreq, vex)

	call cpu_time (time_in)

	open (88,file='Fluxes.txt',status='unknown')

	do while (istep.le.nstep)

		! stores topography at time t
		call FastScape_Copy_H (hp)

		! execute an erosion/deposition step
		call FastScape_Execute_Step ()
		call FastScape_Get_Step (istep)

		! get solution at time t+Dt and stores it in h0
		call FastScape_Copy_H (h)
		h0 = h

		! apply flexure
		call flexure (h,hp,nx,ny,xl,yl,rhos,rhoa,eet,ibc)
		h0 = h-h0
		call FastScape_Set_All_Layers (h0)

		! when needed saves solution to VTK file
		if (mod(istep,nfreq)==0) then
			call FastScape_Copy_H (h)
			call FastScape_Copy_Erosion_Rate (b)
			print*,istep
			print*,'topo',minval(h),sum(h)/nn,maxval(h)
			print*,'erosion rate',minval(b),sum(b)/nn,maxval(b)
			call FastScape_VTK (b,-vex)
		endif

		call FastScape_Get_Fluxes (tflux, eflux, bflux)

		write (88,*) tflux, eflux, bflux

	enddo

	close (88)

	call cpu_time (time_out)
	print*,'Total run time',time_out-time_in

	call FastScape_Destroy ()

	deallocate (u,h,b,hp,h0,rhos,x,y,kf,kd)

end program flexure_test
