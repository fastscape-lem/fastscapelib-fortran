program Margin

  ! Example of the use of the FastScapeInterface
  ! where a square domain (100x100km) is subjected to constant and uniform uplift
  ! of 1 mm/yr while adjacent area (100x100 km) is kept 1000 m below sea level
  ! bottom boundary is at base level, top is nu flux and left and rigth are cyclic
  ! initial random topography
  ! nonlinear erosion law (n=2, m=0.8)
  ! transport coefficient g = 1
  ! marine tranport is activated too

  ! note that this example introduces 2 arrays x and y containing the position
  ! of the nodes on the grid; these are used to define the uplift function and
  ! the initial topography

  use FastScapeAPI

  implicit none

  integer :: nx, ny, istep, nstep, nfreq, i, j
  double precision :: xl, yl, dt, kfsed, m, n, kdsed, g, sealevel, poro, zporo, ratio, L, kds
  double precision, dimension(:), allocatable :: h, u, x, y, kf, kd, fd

  ! initialize FastScape
  call FastScape_Init ()

  ! set grid size
  nx = 101
  ny = 151
  call FastScape_Set_NX_NY (nx,ny)

  ! allocate memory
  call FastScape_Setup ()

  ! set model dimensions
  xl = 100.d3
  yl = 150.d3
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
  where (y<yl/2.d0) h = h - 1000.d0
  call FastScape_Init_H (h)

  ! set erosional parameters
  kf = 1.d-5
  kfsed = 1.d-5
  m = 0.8d0
  n = 2.d0
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
  kds = 3.d2
  call FastScape_Set_Marine_Parameters &
       (sealevel, poro, poro, zporo, zporo, ratio, L, kds, kds/2.d0)

  ! set uplift rate
  allocate (u(nx*ny))
  u = 1.d-3
  where (y<yl/2.d0) u = 0.d0
  call FastScape_Set_U (u)

  ! set boundary conditions
  call FastScape_Set_BC (1000)

  ! set number of time steps and initialize counter istep
  nstep = 1000
  nfreq = 10
  call FastScape_Get_Step (istep)

  allocate (fd(nx*ny))

  ! loop on time stepping
  do while (istep<nstep)
    ! execute step
    call FastScape_Execute_Step()
    ! get value of time step counter
    call FastScape_Get_Step (istep)
    if ((istep/nfreq)*nfreq.eq.istep) then
      call FastScape_Debug()
      call FastScape_Copy_F (fd)
      call FastScape_VTK (fd, -10.d0)
    endif
  enddo

  ! end FastScape run
  call FastScape_Destroy ()

end program Margin
