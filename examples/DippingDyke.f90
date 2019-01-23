program DippingDyke

! problem to test the variability in erodibility (kf)
! we assume that a dyke dipping at 30 degrees is buried beneath the landscape
! and is progressively exhumed by erosion; for this we use the total erosion
! to define the erodibility array kf

implicit none

integer :: nx, ny, istep, nstep, i, j
double precision :: xl, yl, dt, m, n, g, xDyke, dxDyke, angle, cotana
double precision, dimension(:), allocatable :: h, u, chi, kf, kd, x, y, e

! initialize FastScape
call FastScape_Init ()

! set grid size
nx = 201
ny = 201
call FastScape_Set_NX_NY (nx,ny)

! allocate memory
call FastScape_Setup ()

! set model dimensions
xl = 100.d3
yl = xl
call FastScape_Set_XL_YL (xl,yl)

! computes x and y coordinate arrays
allocate (x(nx*ny),y(nx*ny))
x = (/((float(i-1)/(nx-1)*xl,i=1,nx),j=1,ny)/)
y = (/((float(j-1)/(ny-1)*yl,i=1,nx),j=1,ny)/)

! set time step
dt = 1.d5
call FastScape_Set_DT (dt)

! set random initial topography
allocate (h(nx*ny))
call random_number (h)
call FastScape_Init_H (h)

! set erosional parameters
allocate (kf(nx*ny),kd(nx*ny))
kf = 2.d-5
m = 0.4d0
n = 1.d0
kd = 1.d-1
g = 0.d0
call FastScape_Set_Erosional_Parameters (kf, -1.d0, m, n, kd, -1.d0, g, -1.d0, 1.d0)

! set uplift rate
allocate (u(nx*ny))
u = 1.d-3
u(1:nx)=0.d0
u(nx:nx*ny:nx)=0.d0
u(1:nx*ny:nx)=0.d0
u(nx*(ny-1)+1:nx*ny)=0.d0
call FastScape_Set_U (u)

! set boundary conditions
call FastScape_Set_BC (1111)

! set number of time steps and initialize counter istep
nstep = 500
call FastScape_Get_Step (istep)

allocate (chi(nx*ny))

xDyke = xl/10.d0
dxDyke = xl/50.d0
angle = 30.d0
cotana = 1.d0/tan(angle*3.141592654d0/180.d0)

allocate (e(nx*ny))

! loop on time stepping
	do while (istep < nstep)
	call FastScape_Copy_Total_Erosion (e)
	kf = 2.d-5
	where ((x-xDyke-e*cotana-dxDyke)*(x-xDyke-e*cotana+dxDyke).le.0.d0) kf = 1.d-5
	call FastScape_Set_Erosional_Parameters (kf, -1.d0, m, n, kd, -1.d0, g, -1.d0, 1.d0)
	! execute step
	call FastScape_Execute_Step()
	! get value of time step counter
	call FastScape_Get_Step (istep)
	! extract solution
	call FastScape_Copy_Chi (chi)
	! create VTK file
	call FastScape_VTK (chi, 2.d0)
	! outputs h values
	call FastScape_Copy_h (h)
	print*,'h range:',minval(h),sum(h)/(nx*ny),maxval(h)
	enddo

! output timing
call FastScape_Debug()

! end FastScape run
call FastScape_Destroy ()

deallocate (h,u,kf,kd,chi,x,y,e)

end program DippingDyke
