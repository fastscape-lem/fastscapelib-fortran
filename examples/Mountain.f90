program Mountain

  ! simple example of the use of the FastScapeLib
  ! where a square domain (100x100km) is subjected to constant and uniform uplift
  ! of 1 mm/yr
  ! all boundaries are at base level
  ! initial random topography
  ! nonlinear erosion law (n=1.5, m=0.6)
  ! transport coefficient g = 1

  implicit none

  integer :: nx, ny, istep, nstep
  double precision :: xl, yl, dt, kfsed, m, n, kdsed, g
  double precision, dimension(:), allocatable :: h, u, chi, kf, kd, b

  ! initialize FastScape
  call FastScape_Init ()

  ! set grid size
  nx = 401
  ny = 401
  call FastScape_Set_NX_NY (nx,ny)

  ! allocate memory
  call FastScape_Setup ()

  ! set model dimensions
  xl = 100.d3
  yl = xl
  call FastScape_Set_XL_YL (xl,yl)

  ! set time step
  dt = 1.d5
  call FastScape_Set_DT (dt)

  ! set random initial topography
  allocate (h(nx*ny))
  call random_number (h)
  call FastScape_Init_H (h)

  ! set erosional parameters
  allocate (kf(nx*ny),kd(nx*ny))
  kf = 2.d-6
  kfsed = -1.d0
  m = 0.6d0
  n = 1.5d0
  kd = 1.d-1
  kd = 0.d0
  kdsed = -1.d0
  g = 0.d0
  call FastScape_Set_Erosional_Parameters (kf, kfsed, m, n, kd, kdsed, g, g, -2.d0)

  ! set uplift rate (uniform while keeping bounaries at base level)
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
  nstep = 100
  call FastScape_Get_Step (istep)

  !allocate memory to extract chi
  allocate (chi(nx*ny),b(nx*ny))

  ! loop on time stepping
  do while (istep<nstep)
    ! execute step
    call FastScape_Execute_Step()
    ! get value of time step counter
    call FastScape_Get_Step (istep)
    ! extract solution
    call FastScape_Copy_Chi (chi)
    ! create VTK file
    call FastScape_VTK (chi, -2.d0)
    ! outputs h values
    call FastScape_Copy_h (h)
    print*,'step',istep
    print*,'h range:',minval(h),sum(h)/(nx*ny),maxval(h)
    call FastScape_Copy_Basement (b)
    print*,'b range:',minval(b),sum(b)/(nx*ny),maxval(b)
  enddo

  ! output timing
  call FastScape_Debug()

  ! end FastScape run
  call FastScape_Destroy ()

  deallocate (h,u,kf,kd,chi,b)

end program Mountain
