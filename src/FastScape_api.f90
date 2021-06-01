! FastScape API

! -----------------------------------------------------------------------------------------

! contains a series of subroutines that can be accessed from outside FastScape
! to setup, run and close FastScape when used as a subroutine from another program

! subroutines and their functions

! FastScape_Init ()
! Must be called before any other routine to initialize nx, ny and step

! FastScape_SetUp ()
! Must be called to allocate memory for all internal arrays
! can only be called once FastScapeSetNXNY has been used to set up nx and ny

! FastScape_Execute_Step ()
! Executes a single step solving the SPL and diffusion equations

! FastScape_Destroy ()
! Must be called to deallocate memory

! FastScape_View ()
! prints to standard output the value of nx,ny,nn,step,xl,yl,dt,K,m,n,kd,ibc
! as well as min, mean and max values of h and u

! FastScape_Set_NX_NY (nx,ny)
! sets the value of nx and ny, the rectangular grid dimensions
! nx and ny are integer

! FastScape_Set_XL_YL (xl,yl)
! sets the value of xl,yl, the rectangular grid extent (in m)
! xl and yl are double precision

! FastScape_Set_Erosional_Parameters (k1,k2,m,n,kd1,kd2,g1,g2)
! sets the value of the erosional parameters
! k1,k2 are rate coefficient in the stream power law (in m^(1-2*m)/yr) for bedrock and sediment respectively
! m is the area exponent in the stream power law
! kd1, kd2 are the hillslope transport coefficient or diffusivity (in m^2/yr) for bedrock and sediment respectively
! g1, g2 are the sediment fluvial transport/deposition coefficients (dimensionless) for bedrock and sediment respectively
! all parameters are double precision

! FastScape_Set_Marine_Parameters (sealevel, poro1, poro2, zporo1, zporo2, ratio, length, kds1, kds2)
! sets the value of the marine transport parameters
! sl is sea level (in m)
! poro1 is surface porosity for shale (dimensionless)
! poro2 is surface porosity for sand (dimensionless)
! zporo1 is e-folding porosity depth for shale (in m)
! zporo2 is e-folding porosity depth for sand (in m)
! ratio is the ratio of sand in the incoming flux from the continent (dimensionless)
! length is the thickness of the "mixed" surface layer (in m) at the bottom of the ocean
! kds1 and kds2 are the marine transport coefficients (diffusivities) for shale and sand respectively (in m^2/yr)

! FastScape_Set_DT (dt)
! sets the time step length (in yr)
! dt is double precision

! FastScape_Set_BC (ibc)
! sets the boundary conditions
! two types are allowed (0 is reflective bc and 1 is fixed base level)
! ibc should be an integer made of 0 and 1 corresponding to the four boundaries in the
! following order: bottom, right, top and left
! ibc is integer

! FastScape_Set_U (u)
! sets the uplift velocity/rate (in m/yr)
! an array of dimension nn(=nx*ny) should be passed
! u is double precision of size nn

! FastScape_Set_V (vx,vy)
! sets the x- and y-direction advection velocities/rates (in m/yr)
! two array of dimension nn(=nx*ny) should be passed
! vx and vy are double precision of size nn

! FastScape_Init_H (h)
! sets the initial topography (in m) as well as the basement heigh to h
! an array of dimension nn(=nx*ny) should be passed
! h is double precision of size nn

! FastScape_Init_F (F)
! sets the initial sand to shale ratio to F
! an array of dimension nn(=nx*ny) should be passed
! F is double precision of size nn

! FastScape_Copy_H (h)
! returns the current topographic height (in m)
! as an array of dimension nn(=nx*ny)
! h is double precision of size nn

! FastScape_Copy_Basement (b)
! returns the current basement height (in m)
! as an array of dimension nn(=nx*ny)
! b is double precision of size nn

! FastScape_Copy_F (F)
! returns the current surface sand-to-shale ratio (in m)
! as an array of dimension nn(=nx*ny)
! F is double precision of size nn

! FastScape_Copy_Etot (etot)
! returns the current cumulative erosion (in m)
! as an array of dimension nn(=nx*ny)
! etot is double precision of size nn

! FastScape_Reset_Cumulative_Erosion ()
! resets current cumulative erosion

! FastScape_Copy_Area (area)
! returns the drainage area at each point (in m^2)
! as an array of dimension nn(=nx*ny)
! area is double precision of size nn

! FastScape_Copy_Erate (erate)
! returns the current erosion rate (in m/yr)
! as an array of dimension nn(=nx*ny)
! erate is double precision of size nn

! FastScape_Get_Sizes (nx,ny)
! returns the value of the grid size
! nx and ny are integer

! FastScape_Get_Step (step)
! returns the value of the current time step
! step is integer

! FastScape_Set_H (h)
! resets the surface topography (in m)
! an array of dimension nn(=nx*ny) should be passed
! h is double precision of size nn

! FastScape_Set_Basement (b)
! resets the basement topography (in m)
! an array of dimension nn(=nx*ny) should be passed
! b is double precision of size nn

! FastScape_Set_Precip (p)
! resets the precipitation rate (in m/yr)
! an array of dimension nn(=nx*ny) should be passed
! p is double precision of size nn

! FastScape_Debug()
! writes debugging information to the default output

! -----------------------------------------------------------------------------------------
#include "Error.fpp"

module FastScapeAPI

  use iso_c_binding

  contains

  subroutine FastScape_Init(ierr) bind(C, name='fastscape_init')

    use FastScapeContext
    implicit none

    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Init()')

    call Init()

    !return

  end subroutine FastScape_Init

  !--------------------------------------------------------------------------

  subroutine FastScape_Setup(ierr) bind(C, name='fastscape_setup')

    use FastScapeContext
    implicit none

    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Setup()')

    if (nx.eq.0) then
      FSCAPE_RAISE_MESSAGE('FastScape_Setup(): nx cannot be zero',ERR_ParameterInvalid,ierr_)
    end if
    if (nx.le.0) then
       FSCAPE_RAISE_MESSAGE('FastScape_Setup(): nx cannot be negative',ERR_ParameterOutOfRange,ierr_)
    end if
    if (ny.eq.0) then
       FSCAPE_RAISE_MESSAGE('FastScape_Setup(): ny cannot be zero',ERR_ParameterInvalid,ierr_)
    end if
    if (ny.le.0) then
      FSCAPE_RAISE_MESSAGE('FastScape_Setup(): ny cannot be negative',ERR_ParameterOutOfRange,ierr_)
    end if
    FSCAPE_CHKERR_OPT(ierr, ierr_) ! Call FSCAPE_CHKERR() so that all possible exceptions above will be displayed

    call SetUp()

    return

  end subroutine FastScape_Setup

  !--------------------------------------------------------------------------

  subroutine FastScape_Destroy(ierr) bind(C, name='fastscape_destroy')

    use FastScapeContext

    implicit none

    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Destroy()')

    call Destroy()

    return

  end subroutine FastScape_Destroy

  !--------------------------------------------------------------------------

  subroutine FastScape_View(ierr) bind(C, name='fastscape_view')

    use FastScapeContext

    implicit none

    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_View()')

    call View()

    return

  end subroutine FastScape_View

  !--------------------------------------------------------------------------
  subroutine FastScape_Execute_Step(ierr) bind(C, name='fastscape_execute_step')

    use FastScapeContext

    implicit none

    real :: time_in, time_out
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Execute_Step()')

    if (runAdvect) then
      call cpu_time (time_in)
      call Advect ()
      call cpu_time (time_out)
      timeAdvect = timeAdvect + time_out-time_in
    endif

    if (runUplift) then
      call cpu_time (time_in)
      call Uplift()
      call cpu_time (time_out)
      timeUplift = timeUplift + time_out-time_in
    endif

    if (runSPL) then
      call cpu_time (time_in)
      if (SingleFlowDirection) then
        call FlowRoutingSingleFlowDirection (ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
        call FlowAccumulationSingleFlowDirection ()
        call StreamPowerLawSingleFlowDirection ()
      else
        call FlowRouting (ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
        call FlowAccumulation ()
        call StreamPowerLaw ()
      endif
      call cpu_time (time_out)
      timeSPL = timeSPL + time_out-time_in
    endif

    if (runDiffusion) then
      call cpu_time (time_in)
      call Diffusion (ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
      call cpu_time (time_out)
      timeDiffusion = timeDiffusion + time_out-time_in
    endif

    if (runMarine) then
       call cpu_time (time_in)
       call Marine (ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
       call cpu_time (time_out)
       timeMarine = timeMarine + time_out-time_in
    endif

    if (runStrati) then
       call cpu_time (time_in)
       call Run_Strati ()
       call cpu_time (time_out)
       timeStrati = timeStrati + time_out-time_in
    endif

    step=step+1

    return

  end subroutine FastScape_Execute_Step

  !--------------------------------------------------------------------------

  subroutine FastScape_Init_H(hp,ierr) bind(C, name='fastscape_init_h')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: hp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Init_H()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call InitH(hp)

    return

  end subroutine FastScape_Init_H

  !--------------------------------------------------------------------------

  subroutine FastScape_Init_F(Fmixp,ierr) bind(C, name='fastscape_init_f')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: Fmixp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Init_F()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call InitF (Fmixp)

    return

  end subroutine FastScape_Init_F

  !--------------------------------------------------------------------------

  subroutine FastScape_Copy_H(hp,ierr) bind(C, name='fastscape_copy_h')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: hp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Copy_H()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call CopyH(hp)

    return

  end subroutine FastScape_Copy_H

  !--------------------------------------------------------------------------

  subroutine FastScape_Copy_Basement(bp,ierr) bind(C, name='fastscape_copy_basement')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: bp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Copy_Basement()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call CopyBasement(bp)

    return

  end subroutine FastScape_Copy_Basement

  !--------------------------------------------------------------------------

  subroutine FastScape_Copy_Total_Erosion (etotp,ierr)  bind(C, name='fastscape_copy_total_erosion')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: etotp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Copy_Total_Erosion()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call CopyEtot(etotp)

    return

  end subroutine FastScape_Copy_Total_Erosion

  !--------------------------------------------------------------------------

  subroutine FastScape_Copy_Drainage_Area (ap,ierr) bind(C, name='fastscape_copy_drainage_area')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: ap
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Copy_Drainage_Area()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call CopyArea(ap)

    return

  end subroutine FastScape_Copy_Drainage_Area

  !--------------------------------------------------------------------------

  subroutine FastScape_Copy_Erosion_Rate (eratep,ierr) bind(C, name='fastscape_copy_erosion_rate')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: eratep
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Copy_Erosion_Rate()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call CopyERate(eratep)

    return

  end subroutine FastScape_Copy_Erosion_Rate

  !--------------------------------------------------------------------------

  subroutine FastScape_Copy_Chi (chip,ierr) bind(C, name='fastscape_copy_chi')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: chip
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Copy_Chi()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call CopyChi(chip)

    return

  end subroutine FastScape_Copy_Chi

  !--------------------------------------------------------------------------

  subroutine FastScape_Copy_Slope (slopep,ierr) bind(C, name='fastscape_copy_slope')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: slopep
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Copy_Slope()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call CopySlope(slopep)

    return

  end subroutine FastScape_Copy_Slope

  !--------------------------------------------------------------------------

  subroutine FastScape_Copy_Curvature (curvaturep,ierr) bind(C, name='fastscape_copy_curvature')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: curvaturep
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Copy_Curvature()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call CopyCurvature(curvaturep)

    return

  end subroutine FastScape_Copy_Curvature

  !--------------------------------------------------------------------------

  subroutine FastScape_Copy_Catchment (catchp,ierr) bind(C, name='fastscape_copy_catchment')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: catchp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Copy_Catchment()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call CopyCatchment (catchp)

    return

  end subroutine FastScape_Copy_Catchment

  !--------------------------------------------------------------------------

  subroutine FastScape_Copy_F(Fmixp,ierr) bind(C, name='fastscape_copy_f')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: Fmixp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Copy_F()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call CopyF(Fmixp)

    return

  end subroutine FastScape_Copy_F

  !--------------------------------------------------------------------------

  subroutine FastScape_Copy_Lake_Depth(Lp,ierr) bind(C, name='fastscape_copy_lake_depth')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: Lp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Copy_Lake_Depth()')

    if (.not.setup_has_been_run) then
      FSCAPE_RAISE(ERR_SetupNotRun,ierr_);FSCAPE_CHKERR_OPT(ierr, ierr_)
    end if

    call CopyLakeDepth(Lp)

    return

  end subroutine FastScape_Copy_Lake_Depth

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_NX_NY (nnx,nny,ierr) bind(C, name='fastscape_set_nx_ny')

    use FastScapeContext

    implicit none

    integer(c_int), intent(in) :: nnx,nny
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_NX_NY()')

    call SetNXNY (nnx,nny)

    return

  end subroutine FastScape_Set_NX_NY

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_XL_YL (xxl,yyl,ierr) bind(C, name='fastscape_set_xl_yl')

    use FastScapeContext

    implicit none

    real(c_double), intent(in) :: xxl,yyl
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_XL_YL()')

    call SetXLYL (xxl,yyl)

    return

  end subroutine FastScape_Set_XL_YL

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_DT (dtt,ierr) bind(C, name='fastscape_set_dt')

    use FastScapeContext

    implicit none

    real(c_double), intent(in) :: dtt
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_DT()')

    call SetDT (dtt)

    return

  end subroutine FastScape_Set_DT

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_Erosional_Parameters (kkf,kkfsed,mm,nnn,kkd,kkdsed,gg1,gg2,pp,ierr)  bind(C, name='fastscape_set_erosional_parameters')

    use FastScapeContext

    implicit none

    real(c_double), intent(in), dimension(*) :: kkf,kkd
    real(c_double), intent(in) :: kkfsed,mm,nnn,kkdsed,gg1,gg2,pp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_Erosional_Parameters()')

    call SetErosionalParam (kkf,kkfsed,mm,nnn,kkd,kkdsed,gg1,gg2,pp)

    return

  end subroutine FastScape_Set_Erosional_Parameters

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_Marine_Parameters (sl, p1, p2, z1, z2, r, l, kds1, kds2,ierr) bind(C, name='fastscape_set_marine_parameters')

    use FastScapeContext

    implicit none

    real(c_double), intent(in) :: sl, p1, p2, z1, z2, r, l, kds1, kds2
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_Marine_Parameters()')

    call SetMarineParam (sl, p1, p2, z1, z2, r, l, kds1, kds2)

    return

  end subroutine FastScape_Set_Marine_Parameters

  !--------------------------------------------------------------------------

  subroutine FastScape_Get_Sizes (nnx,nny,ierr) bind(C, name='fastscape_get_sizes')

    use FastScapeContext

    implicit none

    integer(c_int), intent(out) :: nnx,nny
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Get_Sizes()')

    call GetSizes (nnx,nny)

    return

  end subroutine FastScape_Get_Sizes

  !--------------------------------------------------------------------------

  subroutine FastScape_Get_Step (sstep,ierr) bind(C, name='fastscape_get_step')

    use FastScapeContext

    implicit none

    integer(c_int), intent(out) :: sstep
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Get_Step()')

    call GetStep (sstep)

    return

  end subroutine FastScape_Get_Step

  !--------------------------------------------------------------------------

  subroutine FastScape_Debug(ierr) bind(C, name='fastscape_debug')

    use FastScapeContext

    implicit none

    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Debug()')

    call Debug()

    return

  end subroutine FastScape_Debug

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_BC(jbc,ierr) bind(C, name='fastscape_set_bc')

    use FastScapeContext

    implicit none

    integer(c_int), intent(in) :: jbc
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_BC()')

    call SetBC (jbc)

    return

  end subroutine FastScape_Set_BC

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_U (up,ierr) bind(C, name='fastscape_set_u')

    use FastScapeContext

    implicit none

    real(c_double), intent(in), dimension(*) :: up
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_U()')

    call SetU(up)

    return

  end subroutine FastScape_Set_U

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_V (ux,uy,ierr) bind(C, name='fastscape_set_v')

    use FastScapeContext

    implicit none

    real(c_double), intent(in), dimension(*) :: ux,uy
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_V()')

    call SetV(ux,uy)

    return

  end subroutine FastScape_Set_V

  !--------------------------------------------------------------------------

  subroutine FastScape_Reset_Cumulative_Erosion (ierr) bind(C, name='fastscape_reset_cumulative_erosion')

    use FastScapeContext

    implicit none

    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Reset_Cumulative_Erosion()')

    call ResetCumulativeErosion ()

    return

  end subroutine FastScape_Reset_Cumulative_Erosion

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_H(hp,ierr) bind(C, name='fastscape_set_h')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: hp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_H()')

    call SetH(hp)

    return

  end subroutine FastScape_Set_H

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_All_Layers (dhp,ierr) bind(C, name='fastscape_set_all_layers')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: dhp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_All_Layers()')

    call SetAllLayers(dhp)

    return

  end subroutine FastScape_Set_All_Layers

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_Basement(bp,ierr) bind(C, name='fastscape_set_basement')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: bp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_Basement()')

    call SetBasement(bp)

    return

  end subroutine FastScape_Set_Basement

  !--------------------------------------------------------------------------

  subroutine FastScape_Set_Precip (precipp,ierr) bind(C, name='fastscape_set_precip')

    use FastScapeContext

    implicit none

    real(c_double), intent(inout), dimension(*) :: precipp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Set_Precip()')

    call SetPrecip (precipp)

    return

  end subroutine FastScape_Set_Precip

  !--------------------------------------------------------------------------

  subroutine FastScape_VTK (fp, vexp, ierr) bind(C, name='fastscape_vtk')

    use FastScapeContext

    implicit none

    real(c_double), intent(in), dimension(*) :: fp
    real(c_double), intent(in) :: vexp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_VTK()')

    call Make_VTK (fp, vexp)

    return

  end subroutine FastScape_VTK

  !--------------------------------------------------------------------------

  subroutine FastScape_Strati (nstepp, nreflectorp, nfreqp, vexp, ierr) bind(C, name='fastscape_strati')

    use FastScapeContext

    implicit none

    integer(c_int), intent(inout) :: nstepp, nreflectorp, nfreqp
    real(c_double), intent(inout) :: vexp
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Strati()')

    call Activate_Strati (nstepp, nreflectorp, nfreqp, vexp)

    return

  end subroutine FastScape_Strati

  !--------------------------------------------------------------------------

  subroutine FastScape_Get_Fluxes (ttectonic_flux, eerosion_flux, bboundary_flux, ierr) bind(C, name='fastscape_get_fluxes')

    use FastScapeContext

    implicit none

    real(c_double), intent(out) :: ttectonic_flux, eerosion_flux, bboundary_flux
    integer(c_int), optional, intent(out) :: ierr
    integer :: ierr_

    FSCAPE_INITERR(ierr, ierr_, 'FastScape_Get_Fluxes()')

    call compute_fluxes (ttectonic_flux, eerosion_flux, bboundary_flux)

    return

  end subroutine FastScape_Get_Fluxes

end module FastScapeAPI
