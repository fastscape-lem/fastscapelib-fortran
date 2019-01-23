subroutine VTK (h,name,nf,f,fname,nx,ny,dx,dy,istep,vex)

! subroutine to create a simple VTK file for plotting

implicit none

double precision h(nx,ny),f(nx,ny,nf),dx,dy,vex
integer nx,ny,nf,istep,nheader,nfooter,npart1,npart2,nn,namel
character header*1024,footer*1024,part1*1024,part2*1024,nxc*6,nyc*6,nnc*12
character*(*) name,fname(nf)

integer i,j,k
character*7 cstep

namel = sum(len_trim(fname(:)))

write (cstep,'(i7)') istep
if (istep.lt.10) cstep(1:6)='000000'
if (istep.lt.100) cstep(1:5)='00000'
if (istep.lt.1000) cstep(1:4)='0000'
if (istep.lt.10000) cstep(1:3)='000'
if (istep.lt.100000) cstep(1:2)='00'
if (istep.lt.1000000) cstep(1:1)='0'

!if (nf.gt.10) stop 'too many fields to be displayed by VTK, maximum is 10'
!  do i=1,nf
!  write (CI(i),'(i1)') i
!  enddo
  
call system ("mkdir -p VTK")

nn=nx*ny
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

call system ('rm -f VTK/'//trim(name)//cstep//'.vtk')

open(unit=77,file='VTK/'//trim(name)//cstep//'.vtk',status='unknown',form='unformatted',access='direct', &
      	recl=nheader+3*4*nn+nfooter+(npart1+1+npart2+4*nn) &
		+nf*(npart1+npart2+4*nn)+namel,convert='big_endian')
write (77,rec=1) &
	header(1:nheader), &
	((sngl(dx*(i-1)),sngl(dy*(j-1)),sngl(h(i,j)*vex),i=1,nx),j=1,ny), &
	footer(1:nfooter), &
	part1(1:npart1)//'H'//part2(1:npart2),sngl(h), &
	(part1(1:npart1)//trim(fname(k))//part2(1:npart2),sngl(f(:,:,k)),k=1,nf)
	close(77)

return
end subroutine VTK

!------------------------------------------------

subroutine VTK_CUBE (fields, nx, ny, nf, nreflector, xl, yl, fname)

implicit none

double precision, dimension(nx,ny,nf,0:nreflector) :: fields
integer :: nx, ny, nf, nreflector
double precision :: xl, yl
integer nheader,nfooter,npart1,npart2,nn,namel
character header*1024,footer*1024,part1*1024,part2*1024,nxc*6,nyc*6,nnc*12,nrefc*3
character*(*) fname(nf)
double precision :: dx, dy, dz
integer :: i, j, k

dx = xl/(nx - 1)
dy = yl/(ny - 1)
dz = (dx+dy)/2.d0

namel = sum(len_trim(fname(:)))

nn=nx*ny*(nreflector + 1)
write (nxc,'(i6)') nx
write (nyc,'(i6)') ny
write (nrefc,'(i3)') nreflector + 1
write (nnc,'(i12)') nn

header(1:1024)=''
header='# vtk DataFile Version 3.0'//char(10)//'FastScape'//char(10) &
//'BINARY'//char(10)//'DATASET STRUCTURED_GRID'//char(10) &
//'DIMENSIONS '//nxc//' '//nyc//' '//nrefc//char(10)//'POINTS' &
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

open(unit=77,file='VTK/CUBE.vtk',status='unknown',form='unformatted',access='direct', &
	recl=nheader+3*4*nn+nfooter+nf*(npart1+npart2+4*nn)+namel,convert='big_endian')
write (77,rec=1) &
	header(1:nheader), &
	(((sngl(dx*(i-1)),sngl(dy*(j-1)),sngl(dz*(k-1)),i=1,nx),j=1,ny),k=1,nreflector+1), &
	footer(1:nfooter), &
	(part1(1:npart1)//trim(fname(k))//part2(1:npart2),sngl(fields(:,:,k,:)),k=1,nf)
close(77)

end subroutine VTK_CUBE
