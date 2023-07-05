subroutine Fastscape_Named_VTK (f, vex, istep, foldername, k)

    ! Subroutine modified from VTK.f90 to create a simple
    ! VTK file inside of a specified folder path.
    ! This additionally outputs the topography, basement
    ! erosion rate, total erosion, drainage area, and
    ! catchment on top of the given HHHHH field.
    ! This assumes the output folder is created beforehand.

    use FastScapeContext
    implicit none

    integer, intent(in) :: k, istep
    double precision, intent(in) :: vex
    double precision, intent(in), dimension(*) :: f
    character(len=k), intent(in) :: foldername
    character cstep*7

    integer nheader,nfooter,npart1,npart2, ftime
    character header*1024,footer*1024,part1*1024,part2*1024,nxc*6,nyc*6,nnc*12
    integer i,j
    double precision dx,dy

    dx = xl/(nx - 1)
    dy = yl/(ny - 1)

    write (cstep,'(i7)') istep
    if (istep.lt.10) cstep(1:6)='000000'
    if (istep.lt.100) cstep(1:5)='00000'
    if (istep.lt.1000) cstep(1:4)='0000'
    if (istep.lt.10000) cstep(1:3)='000'
    if (istep.lt.100000) cstep(1:2)='00'
    if (istep.lt.1000000) cstep(1:1)='0'

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

    ! Formatting to output the correct size for each visualization field.
    open(unit=77,file=(foldername//'/Topography'//cstep//'.vtk'),status='unknown', &
    form='unformatted',access='direct', recl=nheader+3*4*nn+nfooter+(npart1+10+npart2+4*nn) &
    +(npart1+5+npart2+4*nn)+(npart1+8+npart2+4*nn)+(npart1+12+npart2+4*nn)+(npart1+13+npart2+4*nn) &
    +(npart1+13+npart2+4*nn)+(npart1+9+npart2+4*nn))
    write (77,rec=1) &
    header(1:nheader), &
    ((sngl(dx*(i-1)),sngl(dy*(j-1)),sngl(h(i+(j-1)*nx)*abs(vex)),i=1,nx),j=1,ny), &
    footer(1:nfooter), &
    part1(1:npart1)//'topography'//part2(1:npart2),sngl(h(1:nn)), &
    part1(1:npart1)//'HHHHH'//part2(1:npart2),sngl(f(1:nn)), &
    part1(1:npart1)//'basement'//part2(1:npart2),sngl(b(1:nn)), &
    part1(1:npart1)//'erosion_rate'//part2(1:npart2),sngl(erate(1:nn)), &
    part1(1:npart1)//'total_erosion'//part2(1:npart2),sngl(etot(1:nn)), &
    part1(1:npart1)//'drainage_area'//part2(1:npart2),sngl(a(1:nn)) , &
    part1(1:npart1)//'catchment'//part2(1:npart2),sngl(catch(1:nn))
    close(77)

    if (vex.lt.0.d0) then
      open(unit=77,file=(foldername//'/Basement'//cstep//'.vtk'),status='unknown',form='unformatted',access='direct', &
      recl=nheader+3*4*nn+nfooter+(npart1+1+npart2+4*nn) &
      +(npart1+5+npart2+4*nn))
      write (77,rec=1) &
      header(1:nheader), &
      ((sngl(dx*(i-1)),sngl(dy*(j-1)),sngl(b(i+(j-1)*nx)*abs(vex)),i=1,nx),j=1,ny), &
      footer(1:nfooter), &
      part1(1:npart1)//'B'//part2(1:npart2),sngl(b(1:nn)), &
      part1(1:npart1)//'HHHHH'//part2(1:npart2),sngl(f(1:nn))
      close(77)
      open(unit=77,file=(foldername//'/SeaLevel'//cstep//'.vtk'),status='unknown',form='unformatted',access='direct', &
      recl=nheader+3*4*nn+nfooter+(npart1+2+npart2+4*nn))
      write (77,rec=1) &
      header(1:nheader), &
      ((sngl(dx*(i-1)),sngl(dy*(j-1)),sngl(sealevel*abs(vex)),i=1,nx),j=1,ny), &
      footer(1:nfooter), &
      part1(1:npart1)//'SL'//part2(1:npart2),(sngl(sealevel),i=1,nn)
      close(77)
    endif


    return
  end subroutine Fastscape_Named_VTK
