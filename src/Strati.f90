subroutine Strati (h,b,F,nx,ny,xl,yl,reflector,nreflector,ireflector,istep,fields,nfield,vex,dt, &
  stack,rec,length,sealevel)

  ! this routine tracks information (fields) on a set of reflectors (reflector)
  ! and outputs it to a set of VTKs.

  implicit none

  integer :: nx, ny, nn, nreflector, ireflector, i, nfield, istep
  double precision reflector(nx*ny,nreflector), h(nx*ny), b(nx*ny), F(nx*ny), dx, dy, vex, dt
  double precision fields(nx*ny,nfield,nreflector), xl, yl,sealevel
  double precision, dimension(nx*ny) :: length
  integer, dimension(nx*ny) :: stack,rec
  character*30 names(nfield)

  double precision, dimension(:), allocatable :: s,dist
  character*3 :: ref

  ! 1: current depth
  ! 2: current slope
  ! 3: current thickness between me and next horizon
  ! 4: current thickness between me and basement
  ! 5: paleo bathymetry
  ! 6: paleo slope
  ! 7: paleo distance to shore
  ! 8: paleo F
  ! 9: age of reflector (always with respect to present)
  ! 10: how much has been eroded below

  nn = nx*ny
  allocate (s(nn),dist(nn))

  dx = xl/(nx - 1)
  dy = yl/(ny - 1)

  do i = 1, nreflector
    fields(:,1,i) = reflector(:,i)
    call slope (reflector(:,i),s,nx,ny,dx,dy)
    fields(:,2,i) = s
    if (rec(1).ne.0) call distance_to_shore (reflector(:,i),dist,nx,ny,stack,rec,length)
    if (i.gt.1) then
      fields(:,3,i) = reflector(:,i)-reflector(:,i-1)
    else
      fields(:,3,i) = reflector(:,i)-b
    endif
    fields(:,4,i) = reflector(:,i)-b
    if (i.ge.ireflector) then
      fields(:,5,i) = reflector(:,i) + sealevel
      fields(:,6,i) = s
      fields(:,7,i) = dist
      fields(:,8,i) = F
    endif
    fields(:,9,i) = max(dt*(nreflector-i),0.)
  enddo

  names(1) = '1.CurrentDepth(m)'
  names(2) = '2.CurrentSlope(Deg)'
  names(3) = '3.ThicknessToNextReflector(m)'
  names(4) = '4.ThicknessToBasement(m)'
  names(5) = '5.DepositionalBathymetry(m)'
  names(6) = '6.DepositionalSlope(Deg)'
  names(7) = '7.DistanceToShore(m)'
  names(8) = '8.Sand/ShaleRatio'
  names(9) = '9.ReflectorAge(yr)'
  names(10) = 'A.ThicknessErodedBelow(m)'

  do i = 1, nreflector
    write (ref,'(i3)') i
    if (i.lt.10) ref(1:2)='00'
    if (i.lt.100) ref(1:1)='0'
    call VTK (reflector(:,i),'Horizon'//ref//'-',nfield,fields(:,1:nfield,i),names, nx,ny,dx,dy,istep,vex)
  enddo

  deallocate (s,dist)

  if (ireflector.eq.nreflector) call VTK_CUBE (fields, nx, ny, nfield, nreflector, xl, yl, names)

  return

end subroutine Strati

!-------------------------------------------------------------------

subroutine distance_to_shore (h,d,nx,ny,stack,rec,length)

  implicit none

  integer nx, ny
  double precision, dimension(nx*ny) :: h, d
  integer, dimension(nx*ny) :: stack, rec
  double precision, dimension(nx*ny) :: length

  integer :: i, ij

  d = length
  do i = nx*ny, 1, -1
    ij = stack(i)
    if (h(rec(ij)).gt.0.d0) then
      d(ij) = 1.d0
      elseif (h(ij).gt.0.d0) then
        d(ij) = length(ij)*(-h(rec(ij))/(h(ij)-h(rec(ij))))
      endif
      d(rec(ij)) = d(rec(ij)) + d(ij)
    enddo

    return

  end subroutine distance_to_shore
