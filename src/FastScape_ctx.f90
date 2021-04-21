
#include "Error.fpp"

module FastScapeContext

  ! Context module for FastScape api
  ! should not be accessed or changed
  ! see API for name of routines and externally accessible variables
  use FastScapeErrorCodes

  implicit none

  integer :: nx, ny, nn, nstack
  integer :: bounds_ibc
  integer :: bounds_i1, bounds_i2, bounds_j1, bounds_j2
  logical :: bounds_xcyclic, bounds_ycyclic
  logical, dimension(:), allocatable :: bounds_bc
  integer :: step
  integer :: nGSStreamPowerLaw, nGSMarine
  logical :: setup_has_been_run
  double precision, target, dimension(:), allocatable :: h,u,vx,vy,length,a,erate,etot,catch,catch0,b,precip,kf,kd
  double precision, target, dimension(:), allocatable :: Sedflux, Fmix
  double precision, target, dimension(:), allocatable :: g
  double precision, target, dimension(:), allocatable :: p_mfd_exp
  double precision, dimension(:,:), pointer, contiguous :: h2, vx2, vy2, etot2, b2
  double precision :: xl, yl, dt, kfsed, m, n, kdsed, g1, g2, p
  double precision :: sealevel, poro1, poro2, zporo1, zporo2, ratio, layer, kdsea1, kdsea2
  integer, dimension(:), allocatable :: stack, ndon, rec
  integer, dimension(:,:), allocatable :: don
  logical :: runSPL, runAdvect, runDiffusion, runStrati, runUplift, runMarine
  real :: timeSPL, timeAdvect, timeDiffusion, timeStrati, timeUplift, timeMarine
  double precision, dimension(:,:), allocatable :: reflector
  double precision, dimension(:,:,:), allocatable :: fields
  integer nfield, nfreq, nreflector, nfreqref, ireflector
  double precision :: vexref
  logical :: SingleFlowDirection
  double precision, dimension(:), allocatable :: lake_depth, hwater
  integer, dimension(:), allocatable :: mnrec,mstack
  integer, dimension(:,:), allocatable :: mrec
  double precision, dimension(:,:), allocatable :: mwrec,mlrec



  contains

  subroutine Init()
    implicit none

    nx=0
    ny=0
    step=0
    setup_has_been_run = .false.
    timeSPL = 0.
    timeAdvect = 0.
    timeDiffusion = 0.
    timeStrati = 0.
    timeMarine = 0.
    timeUplift = 0.

  end subroutine Init

  !---------------------------------------------------------------

  subroutine SetUp()
    implicit none

    nn=nx*ny

    call Destroy()
    allocate (h(nn),u(nn),vx(nn),vy(nn),stack(nn),ndon(nn),rec(nn),don(8,nn),catch0(nn),catch(nn),precip(nn))
    allocate (g(nn))
    allocate (bounds_bc(nn))
    allocate (p_mfd_exp(nn))
    allocate (length(nn),a(nn),erate(nn),etot(nn),b(nn),Sedflux(nn),Fmix(nn),kf(nn),kd(nn))
    allocate (lake_depth(nn),hwater(nn),mrec(8,nn),mnrec(nn),mwrec(8,nn),mlrec(8,nn),mstack(nn))

    h2(1:nx,1:ny) => h
    b2(1:nx,1:ny) => b
    vx2(1:nx,1:ny) => vx
    vy2(1:nx,1:ny) => vy
    etot2(1:nx,1:ny) => etot

    call SetBC (1111)
    call random_number (h)
    h(1:nx) = 0.d0
    h(nx*(ny-1)+1:nx*ny) = 0.d0
    h(1:nx*ny:nx) = 0.d0
    h(nx:nx*ny:nx) = 0.d0
    u = 0.d0
    vx = 0.d0
    vy = 0.d0
    etot = 0.d0
    b = h
    precip = 1.d0
    p_mfd_exp(1:nn) = 1.d0
    call random_number (catch0)
    sealevel = 0.d0
    Fmix = 0.5d0
    lake_depth = 0.d0

    runSPL = .false.
    runAdvect = .false.
    runDiffusion = .false.
    runStrati = .false.
    runMarine = .false.
    runUplift = .false.

    nGSStreamPowerLaw = 0
    nGSMarine = 0

    setup_has_been_run = .true.

    return

  end subroutine SetUp

  !---------------------------------------------------------------

  subroutine Destroy()

    if (allocated(h)) deallocate(h)
    if (allocated(u)) deallocate(u)
    if (allocated(vx)) deallocate(vx)
    if (allocated(vy)) deallocate(vy)
    if (allocated(stack)) deallocate(stack)
    if (allocated(ndon)) deallocate(ndon)
    if (allocated(rec)) deallocate(rec)
    if (allocated(don)) deallocate(don)
    if (allocated(catch0)) deallocate(catch0)
    if (allocated(catch)) deallocate(catch)
    if (allocated(length)) deallocate (length)
    if (allocated(a)) deallocate (a)
    if (allocated(b)) deallocate (b)
    if (allocated(sedflux)) deallocate (sedflux)
    if (allocated(Fmix)) deallocate (Fmix)
    if (allocated(erate)) deallocate(erate)
    if (allocated(etot)) deallocate(etot)
    if (allocated(precip)) deallocate(precip)
    if (allocated(kd)) deallocate(kd)
    if (allocated(kf)) deallocate(kf)
    if (allocated(reflector)) deallocate(reflector)
    if (allocated(fields)) deallocate(fields)
    if (allocated(lake_depth)) deallocate(lake_depth)
    if (allocated(hwater)) deallocate(hwater)
    if (allocated(mrec)) deallocate(mrec)
    if (allocated(mnrec)) deallocate(mnrec)
    if (allocated(mwrec)) deallocate(mwrec)
    if (allocated(mlrec)) deallocate(mlrec)
    if (allocated(mstack)) deallocate(mstack)
    if (allocated(g)) deallocate(g)
    if (allocated(p_mfd_exp)) deallocate(p_mfd_exp)
    if (allocated(bounds_bc)) deallocate(bounds_bc)

    return

  end subroutine Destroy

  !---------------------------------------------------------------

  subroutine CopyH (hp)

    double precision, intent(out), dimension(*) :: hp

    hp(1:nn)=h

    return

  end subroutine CopyH

  !---------------------------------------------------------------

  subroutine CopyBasement (bp)

    double precision, intent(out), dimension(*) :: bp

    bp(1:nn)=b

    return

  end subroutine CopyBasement

  !---------------------------------------------------------------

  subroutine CopyEtot (etotp)

    double precision, intent(inout), dimension(*) :: etotp

    etotp(1:nn)=etot

    return

  end subroutine CopyEtot

  !---------------------------------------------------------------

  subroutine CopyArea (ap)

    double precision, intent(inout), dimension(*) :: ap

    ap(1:nn)=a

    return

  end subroutine CopyArea

  !---------------------------------------------------------------

  subroutine CopyErate (eratep)

    double precision, intent(inout), dimension(*) :: eratep

    eratep(1:nn)=erate

    return

  end subroutine CopyErate

  !---------------------------------------------------------------

  subroutine Copychi (chip)

    double precision, intent(inout), dimension(*) :: chip
    double precision, dimension(:), allocatable :: chi
    integer ij,ijk
    double precision dx,dy,a0

    allocate (chi(nn))
    chi=0.d0
    dx=xl/(nx-1)
    dy=yl/(ny-1)
    a0=dx*dy*10.d0
    do ij=1,nn
      ijk=stack(ij)
      if (a(ijk).gt.a0) chi(ijk)=chi(rec(ijk))+(a0/a(ijk))**(m/n)*length(ijk)
    enddo
    chip(1:nn)=chi
    deallocate(chi)

    return

  end subroutine CopyChi

  !---------------------------------------------------------------

  subroutine CopySlope (slopep)

    double precision, intent(inout), dimension(*) :: slopep
    double precision, dimension(:), allocatable :: s
    double precision dx,dy

    allocate (s(nn))
    dx=xl/(nx-1)
    dy=yl/(ny-1)
    call slope (h,s,nx,ny,dx,dy)
    slopep(1:nn)=s
    deallocate(s)

    return

  end subroutine CopySlope

  !---------------------------------------------------------------

  subroutine CopyCurvature (curvaturep)

    double precision, intent(inout), dimension(*) :: curvaturep
    double precision, dimension(:), allocatable :: c
    double precision dx,dy

    allocate (c(nn))
    dx=xl/(nx-1)
    dy=yl/(ny-1)
    call curvature (h,c,nx,ny,dx,dy)
    curvaturep(1:nn)=c
    deallocate(c)

    return

  end subroutine CopyCurvature

  !---------------------------------------------------------------

  subroutine CopyCatchment (catchp)

    double precision, intent(inout), dimension(*) :: catchp

    catchp(1:nn)=catch

    return

  end subroutine CopyCatchment

  !---------------------------------------------------------------

  subroutine CopyF (Fmixp)

    double precision, intent(out), dimension(*) :: Fmixp

    Fmixp(1:nn) = Fmix

    return

  end subroutine CopyF

  !---------------------------------------------------------------

  subroutine CopyLakeDepth (Lp)

    double precision, intent(out), dimension(*) :: Lp

    Lp(1:nn) = lake_depth

    return

  end subroutine CopyLakeDepth

  !---------------------------------------------------------------

  subroutine InitH (hp)

    double precision, intent(in), dimension(*) :: hp

    h = hp(1:nn)
    b = h

    return

  end subroutine InitH

  !---------------------------------------------------------------

  subroutine InitF (Fmixp)

    double precision, intent(in), dimension(*) :: Fmixp

    Fmix = Fmixp(1:nn)

    return

  end subroutine InitF

  !---------------------------------------------------------------

  subroutine ResetCumulativeErosion ()

    etot = 0.d0

    return

  end subroutine ResetCumulativeErosion

  !---------------------------------------------------------------

  subroutine View()

    write (*,*) 'FastScapeContext:'
    write (*,*) 'nx,ny',nx,ny
    write (*,*) 'nn',nn
    write (*,*) 'step',step
    write (*,*) 'xl,yl',xl,yl
    write (*,*) 'dt',dt
    write (*,*) 'Kf,Kfsed,,m,n,Kd,Kdsed,G1,G2',sum(kf)/nn,kfsed,m,n,sum(kd)/nn,kdsed,g1,g2
    write (*,*) 'ibc',bounds_ibc
    write (*,*) 'h',minval(h),sum(h)/nn,maxval(h)
    write (*,*) 'u',minval(u),sum(u)/nn,maxval(u)

    return

  end subroutine View

  !---------------------------------------------------------------

  subroutine SetNXNY (nnx,nny)

    integer, intent(in) :: nnx,nny

    nx = nnx
    ny = nny

    return

  end subroutine SetNXNY

  !---------------------------------------------------------------

  subroutine SetXLYL (xxl,yyl)

    double precision, intent(in) :: xxl,yyl

    xl = xxl
    yl = yyl

    return

  end subroutine SetXLYL

  !---------------------------------------------------------------

  subroutine SetErosionalParam (kkf,kkfsed,mm,nnn,kkd,kkdsed,gg1,gg2,pp)

    double precision, intent(in), dimension(*) :: kkf,kkd
    double precision, intent(in) :: kkfsed,mm,nnn,kkdsed,gg1,gg2,pp

    runSPL = .true.

    kf(1:nn) = kkf(1:nn)
    kfsed = kkfsed
    m = mm
    n = nnn
    kd(1:nn) = kkd(1:nn)
    kdsed = kkdsed
    g1 = gg1
    g2 = gg2
    p = pp
    p_mfd_exp(1:nn) = pp
    SingleFlowDirection = .false.
    if (pp.lt.-1.5d0) then
      SingleFlowDirection = .true.
      p = 1.d0
    endif

    if (maxval(kd).gt.tiny(kd).or.kdsed.gt.tiny(kdsed)) runDiffusion = .true.

    return

  end subroutine SetErosionalParam

  !---------------------------------------------------------------

  subroutine SetMarineParam (sl, p1, p2, z1, z2, r, l, kds1, kds2)

    double precision, intent(in) :: sl, p1, p2, z1, z2, r, l, kds1, kds2

    runMarine = .true.

    sealevel = sl
    poro1 = p1
    poro2 = p2
    zporo1 = z1
    zporo2 = z2
    ratio = r
    layer = l
    kdsea1 = kds1
    kdsea2 = kds2

    return

  end subroutine SetMarineParam

  !---------------------------------------------------------------

  subroutine SetDT (dtt)

    double precision, intent(in) :: dtt

    dt = dtt

    return

  end subroutine SetDT

  !---------------------------------------------------------------

  subroutine GetSizes (nnx,nny)

    integer, intent(out) :: nnx,nny

    nnx = nx
    nny = ny

    return

  end subroutine GetSizes

  !---------------------------------------------------------------

  subroutine GetStep (nstep)

    integer, intent(out) :: nstep

    nstep = step

    return

  end subroutine GetStep

  !---------------------------------------------------------------

  subroutine Debug ()

    implicit none

    integer i,j,ij,counter

    write (*,*) '--------------------------------------------------------'

    write (*,*) 'Time step', step

    write (*,*) 'Debug information'

    write (*,*) 'Total number of nodes (nx*ny)',nn

    write (*,*) 'Stack size',nstack

    write (*,*) 'Number of nil elements in stack',count(stack==0)

    write (*,*) 'Total number of donors',sum(ndon)

    counter=0
    do ij=1,nn
      if (rec(ij)==ij) counter=counter+1
    enddo

    write (*,*) 'Total number of self donors',counter

    counter=0
    do j=bounds_j1,bounds_j2
      do i=bounds_i1,bounds_i2
        ij=(j-1)*nx+i
        if (rec(ij)==ij) counter=counter+1
      enddo
    enddo

    write (*,*) 'Total number of local minima',counter

    write (*,*) 'Number of Gauss-Siedel iterations (SPL)',nGSStreamPowerLaw
    write (*,*) 'Number of Crank-Nicholson iterations (Marine)',nGSMarine

    write (*,*) 'Timing:'
    if (runSPL) write (*,*) 'SPL:',timeSPL
    if (runDiffusion) write (*,*) 'Diffusion:',timeDiffusion
    if (runMarine) write (*,*) 'Marine:',timeMarine
    if (runAdvect) write (*,*) 'Advection:',timeAdvect
    if (runUplift) write (*,*) 'Uplift:',timeUplift
    if (runStrati) write (*,*) 'Strati:',timeStrati

  end subroutine Debug

  !---------------------------------------------------------------

  subroutine SetBC (jbc)

    implicit none

    integer, intent(in) :: jbc
    character*4 :: cbc

    bounds_ibc = jbc

    write (cbc,'(i4)') jbc
    bounds_bc=.FALSE.
    bounds_i1=1
    bounds_i2=nx
    bounds_j1=1
    bounds_j2=ny
    if (cbc(4:4).eq.'1') bounds_i1=2
    if (cbc(2:2).eq.'1') bounds_i2=nx-1
    if (cbc(1:1).eq.'1') bounds_j1=2
    if (cbc(3:3).eq.'1') bounds_j2=ny-1
    if (cbc(4:4).eq.'1') bounds_bc(1:nn:nx)=.TRUE.
    if (cbc(2:2).eq.'1') bounds_bc(nx:nn:nx)=.TRUE.
    if (cbc(1:1).eq.'1') bounds_bc(1:nx)=.TRUE.
    if (cbc(3:3).eq.'1') bounds_bc(nx*(ny-1)+1:nn)=.TRUE.
    bounds_xcyclic=.FALSE.
    bounds_ycyclic=.FALSE.
    if (cbc(4:4).ne.'1'.and.cbc(2:2).ne.'1') bounds_xcyclic=.TRUE.
    if (cbc(1:1).ne.'1'.and.cbc(3:3).ne.'1') bounds_ycyclic=.TRUE.

  end subroutine SetBC

  !---------------------------------------------------------------

  subroutine SetU (up)

    implicit none

    double precision, intent(in) :: up(*)
    integer i

    runUplift = .true.

    do i=1,nn
      u(i) = up(i)
    enddo

    return

  end subroutine SetU

  !---------------------------------------------------------------

  subroutine SetV (ux,uy)

    implicit none

    double precision, intent(in) :: ux(*),uy(*)
    integer i

    runAdvect = .true.

    do i=1,nn
      vx(i) = ux(i)
      vy(i) = uy(i)
    enddo

    return

  end subroutine SetV

  !---------------------------------------------------------------

  subroutine SetH (hp)

    double precision, intent(in), dimension(*) :: hp

    h = hp(1:nn)

    return

  end subroutine SetH

  !---------------------------------------------------------------

  subroutine SetPrecip (precipp)

    double precision, intent(in), dimension(*) :: precipp

    precip = precipp(1:nn)

    return

  end subroutine SetPrecip

  !---------------------------------------------------------------

  subroutine SetAllLayers (dh)

    double precision, intent(in), dimension(*) :: dh

    integer i

    h = h + dh(1:nn)
    b = b + dh(1:nn)
    if (runStrati) then
      do i = 1, nreflector
        reflector(:,i) = reflector(:,i) + dh(1:nn)
      enddo
    endif

    return

  end subroutine SetAllLayers

  !---------------------------------------------------------------

  subroutine SetBasement (bp)

    double precision, intent(in), dimension(*) :: bp

    b = bp(1:nn)

    return

  end subroutine SetBasement

  !---------------------------------------------------------------

  subroutine Make_VTK (f, vex)

    ! subroutine to create a simple VTK file for plotting

    implicit none

    double precision, intent(in) :: vex
    double precision, intent(in), dimension(*) :: f

    integer nheader,nfooter,npart1,npart2
    character header*1024,footer*1024,part1*1024,part2*1024,nxc*6,nyc*6,nnc*12
    integer i,j
    character*7 cstep
    double precision dx,dy

    dx = xl/(nx - 1)
    dy = yl/(ny - 1)

    write (cstep,'(i7)') step
    if (step.lt.10) cstep(1:6)='000000'
    if (step.lt.100) cstep(1:5)='00000'
    if (step.lt.1000) cstep(1:4)='0000'
    if (step.lt.10000) cstep(1:3)='000'
    if (step.lt.100000) cstep(1:2)='00'
    if (step.lt.1000000) cstep(1:1)='0'

#ifdef ON_WINDOWS
    call system ('if not exist "VTK" mkdir VTK')
#else
      call system ("mkdir -p VTK")
#endif

      write (nxc,'(i6)') nx
      write (nyc,'(i6)') ny
      write (nnc,'(i12)') nn

      header(1:1024)=''
      header='# vtk DataFile Version 3.0'//char(10)//'FastScape'//char(10) &
      //'BINARY'//char(10)//'DATASET STRUCTURED_GRID'//char(10) &
      //'DIMENSIONS '//nxc//' '//nyc//' 1'//char(10)//'POINTS' &
      //nnc//' float'//char(10)
      nheader=len_trim(header)
      footer(1:1024)=''
      footer='POINT_DATA'//nnc//char(10)
      nfooter=len_trim(footer)
      part1(1:1024)=''
      part1='SCALARS '
      npart1=len_trim(part1)+1
      part2(1:1024)=''
      part2=' float 1'//char(10)//'LOOKUP_TABLE default'//char(10)
      npart2=len_trim(part2)

      open(unit=77,file='VTK/Topography'//cstep//'.vtk',status='unknown',form='unformatted',access='direct', &
      recl=nheader+3*4*nn+nfooter+(npart1+1+npart2+4*nn) &
      +(npart1+5+npart2+4*nn),convert='big_endian')
      write (77,rec=1) &
      header(1:nheader), &
      ((sngl(dx*(i-1)),sngl(dy*(j-1)),sngl(h(i+(j-1)*nx)*abs(vex)),i=1,nx),j=1,ny), &
      footer(1:nfooter), &
      part1(1:npart1)//'H'//part2(1:npart2),sngl(h(1:nn)), &
      part1(1:npart1)//'HHHHH'//part2(1:npart2),sngl(f(1:nn))
      close(77)

      if (vex.lt.0.d0) then
        open(unit=77,file='VTK/Basement'//cstep//'.vtk',status='unknown',form='unformatted',access='direct', &
        recl=nheader+3*4*nn+nfooter+(npart1+1+npart2+4*nn) &
        +(npart1+5+npart2+4*nn),convert='big_endian')
        write (77,rec=1) &
        header(1:nheader), &
        ((sngl(dx*(i-1)),sngl(dy*(j-1)),sngl(b(i+(j-1)*nx)*abs(vex)),i=1,nx),j=1,ny), &
        footer(1:nfooter), &
        part1(1:npart1)//'B'//part2(1:npart2),sngl(b(1:nn)), &
        part1(1:npart1)//'HHHHH'//part2(1:npart2),sngl(f(1:nn))
        close(77)
        open(unit=77,file='VTK/SeaLevel'//cstep//'.vtk',status='unknown',form='unformatted',access='direct', &
        recl=nheader+3*4*nn+nfooter+(npart1+2+npart2+4*nn),convert='big_endian')
        write (77,rec=1) &
        header(1:nheader), &
        ((sngl(dx*(i-1)),sngl(dy*(j-1)),sngl(sealevel*abs(vex)),i=1,nx),j=1,ny), &
        footer(1:nfooter), &
        part1(1:npart1)//'SL'//part2(1:npart2),(sngl(sealevel),i=1,nn)
        close(77)
      endif

      return
    end subroutine Make_VTK

    !---------------------------------------------------------------

    subroutine Activate_Strati (nstepp, nreflectorp, nfreqp, vexp)

      implicit none

      double precision, intent(in) :: vexp
      integer, intent(in) :: nstepp, nreflectorp, nfreqp

      nfield = 10

      nfreqref = nstepp/nreflectorp
      ireflector = 1
      vexref = vexp
      nreflector = nreflectorp
      nfreq = nfreqp

      allocate (reflector(nn,nreflector),fields(nn,nfield,nreflector))

      fields=0.d0

      !call Strati (b, Fmix, nx, ny, xl, yl, reflector, nreflector, ireflector, 0, &
      !fields, nfield, vexref, dt*nfreqref, stack, rec, length, sealevel)

      runStrati = .true.

    end subroutine Activate_Strati

    !---------------------------------------------------------------

    subroutine run_Strati ()

      implicit none

      integer i

      ! uplift reflectors
      do i = 1, nreflector
        reflector(:,i) = reflector(:,i) + u*dt
      enddo

      ! updates erosion below each reflector
      do i= 1, ireflector
        fields(:,10,i) = fields(:,10,i)+max(0.,reflector(:,i)-h)
      enddo

      do i = 1, ireflector - 1
        reflector(:,i) = min(reflector(:,i),h)
      enddo

      do i = ireflector, nreflector
        reflector(:,i) = h
      enddo

      if (((step+1)/nfreq)*nfreq.eq.(step+1)) then
        if (((step+1)/nfreqref)*nfreqref.eq.(step+1)) ireflector = ireflector + 1
        call Strati (b, Fmix, nx, ny, xl, yl, reflector, nreflector, ireflector, step + 1, &
        fields, nfield, vexref, dt*nfreqref, stack, rec, length, sealevel)
      endif

    end subroutine run_Strati

    !---------------------------------------------------------------

    subroutine compute_fluxes (tectonic_flux, erosion_flux, boundary_flux)

      implicit none

      double precision, intent(out) :: tectonic_flux, erosion_flux, boundary_flux
      double precision :: surf
      double precision, dimension(:), allocatable :: flux
      integer ij,ijk,k

      surf = xl*yl/(nx - 1)/(ny - 1)

      tectonic_flux = sum(u)*surf
      erosion_flux = sum(erate)*surf

      ! computes receiver and stack information for multi-direction flow
      !allocate (mrec(8,nn), mnrec(nn), mwrec(8,nn), mlrec(8,nn), mstack(nn), hwater(nn)
      !call find_mult_rec (h, rec, stack, hwater, mrec, mnrec, mwrec, mlrec, mstack, nx, ny, xl/(nx-1), yl/(ny-1), p, bounds, p_mfd_exp)
      ! computes sediment flux
      allocate (flux(nn))

      flux = erate

      if (SingleFlowDirection) then
        do ij = nn ,1, -1
          ijk = stack(ij)
          flux(ijk)=max(0.d0,flux(ijk))
          flux(rec(ijk)) = flux(rec(ijk)) + flux(ijk)
        enddo
      else
        do ij = 1, nn
          ijk = mstack(ij)
          flux(ijk)=max(0.d0,flux(ijk))
          do k = 1, mnrec(ijk)
            flux(mrec(k,ijk)) = flux(mrec(k,ijk)) + flux(ijk)*mwrec(k,ijk)
          enddo
        enddo
      endif

      ! compute boundary flux
      boundary_flux = sum(flux,bounds_bc)*surf

      deallocate (flux)

    end subroutine compute_fluxes

  end module FastScapeContext
