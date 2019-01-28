subroutine Advect ()

  use FastScapeContext

  ! routine to advect the topography (h), basement (hb) and total erosion(etot)
  ! by a velocity field vx, vy, vz known at the same locations on a rectangular
  ! grid of nx by ny points separated by dx, dy over a time step dt

  implicit none

  double precision, dimension(:), allocatable :: diag,sup,inf,rhs,res
  double precision dx,dy
  integer i,j,k

  !print*,'Advect'

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  ! x-advection using an implicit, second-order scheme to solve advection eq.

  allocate (diag(nx),sup(nx),inf(nx),rhs(nx),res(nx))

  do j=1,ny

    diag=1.d0
    sup=0.d0
    inf=0.d0

    do i=1,nx
      if (vx2(i,j).gt.0.d0) then
        diag(i)=1.d0+vx2(i,j)*dt/dx
        inf(i)=-vx2(i,j)*dt/dx
      elseif (vx2(i,j).lt.0.d0) then
        diag(i)=1.d0-vx2(i,j)*dt/dx
        sup(i)=vx2(i,j)*dt/dx
      endif
    enddo
    sup(1)=0.d0
    diag(1)=1.d0
    diag(nx)=1.d0
    inf(nx)=0.d0

    rhs=h2(:,j)
    call tridag (inf,diag,sup,rhs,res,nx)
    h2(:,j)=res

    rhs=b2(:,j)
    call tridag (inf,diag,sup,rhs,res,nx)
    b2(:,j)=res

    rhs=etot2(:,j)
    call tridag (inf,diag,sup,rhs,res,nx)
    etot2(:,j)=res

  enddo

  deallocate (diag,sup,inf,rhs,res)

  ! y-advection using an implicit, second-order scheme to solve advection eq.

  allocate (diag(ny),sup(ny),inf(ny),rhs(ny),res(ny))

  do i=1,nx

    diag=1.d0
    sup=0.d0
    inf=0.d0

    do j=1,ny
      if (vy2(i,j).gt.0.d0) then
        diag(j)=1.d0+vy2(i,j)*dt/dy
        inf(j)=-vy2(i,j)*dt/dy
      elseif (vy2(i,j).lt.0.d0) then
        diag(j)=1.d0-vy2(i,j)*dt/dy
        sup(j)=vy2(i,j)*dt/dy
      endif
    enddo
    sup(1)=0.d0
    diag(1)=1.d0
    diag(ny)=1.d0
    inf(ny)=0.d0

    rhs=h2(i,:)
    call tridag (inf,diag,sup,rhs,res,ny)
    h2(i,:)=res

    rhs=b2(i,:)
    call tridag (inf,diag,sup,rhs,res,ny)
    b2(i,:)=res

    rhs=etot2(i,:)
    call tridag (inf,diag,sup,rhs,res,ny)
    etot2(i,:)=res

  enddo

  deallocate (diag,sup,inf,rhs,res)

  b=min(b,h)

  return
  end
