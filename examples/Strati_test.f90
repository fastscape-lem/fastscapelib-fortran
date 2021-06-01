program Strati_test

  ! Example of the use of the FastScapeInterface
  ! where a square domain (100x100km) is subjected to constant and uniform uplift
  ! of 1 mm/yr while adjacent area (100x100 km) is kept 1000 m below sea level
  ! all boundaries are at base level
  ! initial random topography
  ! nonlinear erosion law (n=2, m=0.8)
  ! transport coefficient g = 1
  ! marine tranport is activated too

  ! note that this example introduces 2 arrays x and y containing the position
  ! of the nodes on the grid; these are used to define the uplift function and
  ! the initial topography

  use FastScapeAPI

  implicit none

  integer :: nx, ny, istep, nstep, nfreq, i, j, nreflector
  double precision :: xl, yl, dt, kfsed, m, n, kdsed, g, sealevel, poro, zporo, ratio, L, kds, vex, pi
  double precision, dimension(:), allocatable :: h, u, x, y, kf, kd

  pi=atan(1.d0)*4.d0

  ! initialize FastScape
  call FastScape_Init ()

  ! set grid size
  nx = 101
  ny = 101
  call FastScape_Set_NX_NY (nx,ny)

  ! allocate memory
  call FastScape_Setup ()

  ! set model dimensions
  xl = 400.d3
  yl = 400.d3
  call FastScape_Set_XL_YL (xl,yl)

  ! construct nodal coordinate arrays x and y
  allocate (x(nx*ny),y(nx*ny))
  x = (/((xl*float(i-1)/(nx-1), i=1,nx),j=1,ny)/)
  y = (/((yl*float(j-1)/(ny-1), i=1,nx),j=1,ny)/)

  ! set time step
  dt = 1.d3
  call FastScape_Set_DT (dt)

  ! set random initial topography
  allocate (h(nx*ny),kf(nx*ny),kd(nx*ny))
  call random_number (h)
  where (y<2.d0*yl/3.d0) h = h - 1000.d0*(2.d0*yl/3.d0-y)/(2.d0*yl/3.d0)
  h=h+10.*cos(x/xl*2.d0*pi)
  call FastScape_Init_H (h)

  ! set erosional parameters
  kf = 1.d-5
  kfsed = 1.d-5
  m = 0.4d0
  n = 1.d0
  kd = 1.d-2
  kdsed = 1.d-2
  g = 1.d0
  call FastScape_Set_Erosional_Parameters (kf, kfsed, m, n, kd, kdsed, g, g, 1.d0)

  ! set marine transport parameters
  sealevel = 0.d0
  poro = 0.d0
  zporo = 1.d3
  ratio = 0.5d0
  L = 1.d2
  kds = 5.d2
  call FastScape_Set_Marine_Parameters &
       (sealevel, poro, poro, zporo, zporo, ratio, L, kds, kds/2.d0)

  ! set uplift rate
  allocate (u(nx*ny))
  u = 3.d-4
  where (y<2.d0*yl/3.d0) u = 0.d0
  call FastScape_Set_U (u)

  ! set boundary conditions
  call FastScape_Set_BC (1000)

  ! set number of time steps and initialize counter istep
  nstep = 5000
  call FastScape_Get_Step (istep)

  nreflector = 5
  nfreq = 500
  vex = 100.d0
  call FastScape_Strati (nstep, nreflector, nfreq, vex)

  ! loop on time stepping
  do while (istep<nstep)
   ! execute step
    call FastScape_Execute_Step()
    ! get value of time step counter
    call FastScape_Get_Step (istep)
    ! message to the screen
    if (mod(istep,nfreq).eq.0) then
      call FastScape_Copy_H (h)
      print*,'Calling Strati','step',istep,'h-range:',minval(h),sum(h)/(nx*ny),maxval(h)
    endif
  enddo

  call FastScape_Debug()
  ! end FastScape run
  call FastScape_Destroy ()

end program Strati_test
