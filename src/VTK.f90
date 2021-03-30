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

#ifdef ON_WINDOWS
    call system ('if not exist "VTK" mkdir VTK')
#else
    call system ("mkdir -p VTK")
#endif

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

#ifdef ON_WINDOWS
    call system ('del VTK/'//trim(name)//cstep//'.vtk')
#else
    call system ('rm -f VTK/'//trim(name)//cstep//'.vtk')
#endif

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

    double precision, dimension(nx,ny,nf,nreflector) :: fields
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
    write (nrefc,'(i3)') nreflector
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
    (((sngl(dx*(i-1)),sngl(dy*(j-1)),sngl(dz*(k-1)),i=1,nx),j=1,ny),k=1,nreflector), &
    footer(1:nfooter), &
    (part1(1:npart1)//trim(fname(k))//part2(1:npart2),sngl(fields(:,:,k,:)),k=1,nf)
    close(77)

end subroutine VTK_CUBE

!------------------------------------------------

subroutine VTK_filled (basement, nreflector, reflector, nfield, fields, names, nx, ny, dx, dy, istep, vex, distb)

  ! This routine produces a 3D volume containing all the information produced by Strati in vtk format
  ! The resulting file(s) called Strati-number (where number is the step number) can be loaded in
  ! Paraview for viewing and data extraction (along synthetic wells)
  ! Note however that the data is only stored at the location of the reflectors and that the information
  ! in between reflectors is the result of interpolation made in Paraview
  ! The larger the number of reflectors, the less the interpolation

  ! Note that contrary to all other VTK producing routines, this one outputs the geometric information
  ! in ASCII format (not binary); this means that the files it generates are larger and take more time
  ! to load into Paraview

  integer :: nx, ny, nreflector, nfield, istep, i, j, k, l, ij 
  double precision basement(nx*ny), reflector(nx*ny,nreflector), dx, dy, vex, distb(nx*ny)
  double precision fields(nx*ny,nfield,nreflector)
  character*30 names(nfield)

  character cstep*7, name*128
  integer iunit, nnode, nelem

  write (cstep,'(i7)') istep
  if (istep.lt.10) cstep(1:6)='000000'
  if (istep.lt.100) cstep(1:5)='00000'
  if (istep.lt.1000) cstep(1:4)='0000'
  if (istep.lt.10000) cstep(1:3)='000'
  if (istep.lt.100000) cstep(1:2)='00'
  if (istep.lt.1000000) cstep(1:1)='0'

  iunit=30

#ifdef ON_WINDOWS
  call system ('if not exist "VTK" mkdir VTK')
#else
  call system ("mkdir -p VTK")
#endif

name='Strati-'

#ifdef ON_WINDOWS
  call system ('del VTK/'//trim(name)//cstep//'.vtk')
#else
  call system ('rm -f VTK/'//trim(name)//cstep//'.vtk')
#endif

  nnode = nx*ny*(nreflector+1)
  nelem = (nx-1)*(ny-1)*nreflector

  open(unit=iunit,file='VTK/'//trim(name)//cstep//'.vtk')
  write(iunit,'(a)')'# vtk DataFile Version 3.0'
  write(iunit,'(a)')'FilledStratigraphy'
  write(iunit,'(a)')'ASCII'
  write(iunit,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(iunit,'(a7,i10,a6)')'POINTS ',nnode,' float'

  do k = 0, nreflector
    do j = 1, ny
      do i = 1, nx
        ij = (j - 1)*nx + i
        if (k.eq.0) then
          write(iunit,'(3f16.4)') dx*(i - 1), dy*(j - 1), basement(ij)*vex
        else
          write(iunit,'(3f16.4)') dx*(i - 1), dy*(j - 1), reflector(ij, k)*vex
        endif
      enddo
    enddo
  enddo

  write(iunit,'(A6, 2I10)') 'CELLS ',nelem,9*nelem
  do k=1,nreflector
    do j=1,ny-1
      do i=1,nx-1
        ij = (k-1)*nx*ny+(j-1)*nx+i-1
        write(iunit,'(9I10)') 8 , ij, ij+1, ij+1+nx, ij+nx, &
                                  ij+nx*ny, ij+1+nx*ny, ij+1+nx+nx*ny, ij+nx+nx*ny
      enddo
    enddo
  enddo

  write(iunit,'(A11, I10)') 'CELL_TYPES ',nelem
  do k=1,nelem
  write(iunit,'(I2)') 12 ! octree  (8 nodes)
  enddo

  write(iunit,'(a11,i10)')'POINT_DATA ',nnode

  write(iunit,'(a)')'SCALARS 0.Reflector float 1'
  write(iunit,'(a)')'LOOKUP_TABLE default'
  do k = 0, nreflector
    do j = 1, ny
      do i = 1, nx
        ij = (j - 1)*nx + i
        write(iunit,'(e10.4)') float(k)
      enddo
    enddo
  enddo

  do l = 1, nfield
    write(iunit,'(a)')'SCALARS '//names(l)//' float 1'
    write(iunit,'(a)')'LOOKUP_TABLE default'
    do k = 0, nreflector
      do j = 1, ny
        do i = 1, nx
          ij = (j - 1)*nx + i
          if (k.eq.0) then
            if (l.eq.1.or.l.eq.2) then
              write(iunit,'(e10.4)') fields(ij,l,k+1)
            elseif (l.eq.7) then
              write(iunit,'(e10.4)') distb(ij)
            elseif (l.eq.9) then
              write(iunit,'(e10.4)') 2*fields(ij,l,k+1)-fields(ij,l,k+2)
            else
              write(iunit,'(e10.4)') 0.
            endif
          else
            write(iunit,'(e10.4)') fields(ij,l,k)
          endif
        enddo
      enddo
    enddo
  enddo

  close (iunit)

end subroutine VTK_filled
