subroutine Strati (b,F,nx,ny,xl,yl,reflector,nreflector,ireflector,istep,fields,nfield,vex,dt, &
  stack,rec,length,sealevel)

  ! this routine tracks information (fields) on a set of reflectors (reflector)
  ! and outputs it to a set of VTKs.

  implicit none

  integer :: nx, ny, nn, nreflector, ireflector, i, nfield, istep
  double precision reflector(nx*ny,nreflector), b(nx*ny), F(nx*ny), dx, dy, vex, dt
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

subroutine slope (h,s,nx,ny,dx,dy)

  implicit none

  double precision h(nx*ny),s(nx*ny),dx,dy
  integer nx,ny

  integer i,j,ij,ia,ib,ic,id,ie,if,ig,ih,ii
  double precision dzdx,dzdy,con

  con=45.d0/atan(1.d0)

  s=0.d0
  do j=2,ny-1
    do i=2,nx-1
      ij=i+(j-1)*nx
      ia=ij+nx-1
      ib=ia+1
      ic=ib+1
      id=ij-1
      ie=ij
      if=ij+1
      ig=ij-nx-1
      ih=ig+1
      ii=ih+1
      dzdx=((h(ic)+2.d0*h(if)+h(ii))-(h(ia)+2.d0*h(id)+h(ig)))/8.d0/dx
      dzdy=((h(ig)+2.d0*h(ih)+h(ii))-(h(ia)+2.d0*h(ib)+h(ic)))/8.d0/dy
      s(ij)=dzdx**2+dzdy**2
      if (s(ij).gt.tiny(s(ij))) s(ij)=atan(sqrt(s(ij)))*con
    enddo
  enddo

end subroutine slope

!-------------------------------------------------------------------

subroutine curvature (h,curv,nx,ny,dx,dy)

  implicit none

  double precision h(nx*ny),curv(nx*ny),dx,dy
  integer nx,ny

  integer i,j,ij,i1,i2,i3,i4,i5,i6,i7,i8,i9
  double precision a,b,c,d,e,f

  curv=0.d0
  do j=2,ny-1
    do i=2,nx-1
    ij=i+(j-1)*nx
    i1=ij+nx-1
    i2=i1+1
    i3=i2+1
    i4=ij-1
    i5=ij
    i6=ij+1
    i7=ij-nx-1
    i8=i7+1
    i9=i8+1
    a=(h(i1)+h(i3)+h(i4)+h(i6)+h(i7)+h(i9))/dx/dx/12.d0-(h(i2)+h(i5)+h(i8))/dx/dx/6.d0
    b=(h(i1)+h(i2)+h(i3)+h(i7)+h(i8)+h(i9))/dy/dy/12.d0-(h(i4)+h(i5)+h(i6))/dy/dy/6.d0
    c=(h(i3)+h(i7)-h(i1)-h(i9))/dx/dy/4.d0
    d=(h(i3)+h(i6)+h(i9)-h(i1)-h(i4)-h(i7))/dx/6.d0
    e=(h(i1)+h(i2)+h(i3)-h(i7)-h(i8)-h(i9))/dy/6.d0
    f=(2.d0*(h(i2)+h(i4)+h(i6)+h(i8))-(h(i1)+h(i3)+h(i7)+h(i9))+5.d0*h(i5))/9.d0
    curv(ij)=1.d0+d**2+e**2
    if (curv(ij).gt.tiny(curv(ij))) curv(ij)=(a*(1.d0+e**2)+b*(1.d0+d**2)-c*d*e)/(curv(ij)**(3.d0/2.d0))
    enddo
  enddo

end subroutine curvature

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
