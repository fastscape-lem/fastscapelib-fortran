program Fan

  ! Example of a small fan

  ! This fan is created by relaxation/erosion of a pre-existing plateau

  ! Note how sediments are progressively deposited in the valleys until
  ! the fan is coimpletely full (it has reached steady-state) at which stage
  ! the valleys are progressively emptied

  ! The evolution of the sedimentary flux out of the system and out of the
  ! plateau only is stored in the Fluxes.txt file

  use FastScapeAPI

  implicit none

  integer :: nx,ny,istep,nstep,nn,ibc,ierr
  double precision, dimension(:), allocatable :: h,hp,x,y,kf,kd,b
  real :: time_in,time_out
  double precision :: kfsed,m,n,kdsed,g1,g2,expp
  double precision xl,yl,dt,pi,vex

  integer i,j

  ! set model resolution
  nx = 101
  ny = 201
  nn = nx*ny

  pi=atan(1.d0)*4.d0

  ! initialize FastScape
  call FastScape_Init ()
  call FastScape_Set_NX_NY (nx,ny,ierr)
  call FastScape_Setup ()

  ! set model dimensions
  xl=10.d3
  yl=20.d3
  call FastScape_Set_XL_YL (xl,yl,ierr)

  ! construct nodal coordinate arrays x and y
  allocate (x(nx*ny),y(nx*ny))
  x = (/((xl*float(i-1)/(nx-1), i=1,nx),j=1,ny)/)
  y = (/((yl*float(j-1)/(ny-1), i=1,nx),j=1,ny)/)

  ! set time step
  dt=2.d3
  call FastScape_Set_DT (dt,ierr)

  ! we make the sediment slightly more easily erodible
  allocate (kf(nn),kd(nn))
  kf=1.d-4
  kfsed=1.5d-4
  m=0.4d0
  n=1.d0
  kd=1.d-2
  kdsed=1.5d-2
  g1=1.d0
  g2=1.d0
  expp=1.d0
  call FastScape_Set_Erosional_Parameters (kf,kfsed,m,n,kd,kdsed,g1,g2,expp,ierr)

  ! bottom side is fixed only
  ibc=1000
  call FastScape_Set_BC (ibc,ierr)

  ! initial topography is a 1000 m high plateau
  allocate (h(nn),b(nn),hp(nn))
  call random_number (h)
  where (y.gt.yl/2.d0) h=h+1000.d0
  call FastScape_Init_H (h, ierr)

  ! set number of time steps
  nstep = 200

  ! echo model setup
  call FastScape_View (ierr)

  ! initializes time step
  call FastScape_Get_Step (istep,ierr)

   ! set vertical exaggeration
  vex = 3.d0

  ! start of time loop
  call cpu_time (time_in)
  do while (istep.lt.nstep)

    ! execute FastScape step
    call FastScape_Execute_Step (ierr)
    call FastScape_Get_Step (istep,ierr)

    ! output vtk with sediment thickness information
    call FastScape_Copy_H (h,ierr)
    call FastScape_Copy_Basement (b,ierr)
    call FastScape_VTK (h-b, vex, ierr)

  enddo

  ! display timing information
  call FastScape_Debug(ierr)
  call cpu_time (time_out)
  print*,'Total run time',time_out-time_in

  ! exits FastScape
  call FastScape_Destroy (ierr)

  ! deallocate memory
  deallocate (h,x,y,kf,kd,b,hp)

end program Fan
