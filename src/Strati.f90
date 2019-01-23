subroutine Strati (h,F,nx,ny,xl,yl,reflector,nreflector,ireflector,istep,fields,nfield,vex,dt, &
			stack,rec,length)

! this routine tracks information (fields) on a set of reflectors (reflector)
! and outputs it to a set of VTKs.

implicit none

integer :: nx, ny, nn, nreflector, ireflector, i, nfield, istep
double precision reflector(nx*ny,0:nreflector), h(nx*ny), F(nx*ny), dx, dy, vex, dt
double precision fields(nx*ny,nfield,0:nreflector), xl, yl
double precision, dimension(nx*ny) :: length
integer, dimension(nx*ny) :: stack,rec
character*30 names(nfield)
integer ibc

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

do i = 0, ireflector - 1
reflector(:,i) = min(reflector(:,i),h)
enddo

do i = ireflector, nreflector
reflector(:,i) = h
enddo

do i = 0, nreflector
fields(:,1,i) = reflector(:,i)
call slope (reflector(:,i),s,nx,ny,dx,dy)
fields(:,2,i) = s
call distance_to_shore (reflector(:,i),dist,nx,ny,stack,rec,length)
	if (i.gt.0) then
	fields(:,3,i) = reflector(:,i)-reflector(:,i-1)
	fields(:,4,i) = reflector(:,i)-reflector(:,0)
	endif
	if (i.ge.ireflector) then
	fields(:,5,i) = reflector(:,i)
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

	do i = 0, nreflector
	write (ref,'(i3)') i
	if (i.lt.10) ref(1:2)='00'
	if (i.lt.100) ref(1:1)='0'
	call VTK (reflector(:,i),'Horizon'//ref//'-',nfield,fields(:,1:nfield,i),names, nx,ny,dx,dy,istep,vex)
	enddo

deallocate (s,dist)

if (ireflector.eq.nreflector) call VTK_CUBE (fields, nx, ny, nfield, nreflector, xl, yl, names)

return

end

!-------------------------------------------------------------------

subroutine slope (h,s,nx,ny,dx,dy)

implicit none

double precision h(nx*ny),s(nx*ny),dx,dy
integer nx,ny,nn

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

subroutine distance_to_shore (h,d,nx,ny,stack,rec,length)

implicit none

integer nx, ny
double precision, dimension(nx*ny) :: h, d
double precision dx, dy
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
