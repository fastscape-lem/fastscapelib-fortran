subroutine Diffusion ()

  ! subroutine to solve the diffusion equation by ADI

  use FastScapeContext

  implicit none

  double precision, dimension(:), allocatable :: f,diag,sup,inf,res
  double precision, dimension(:,:), allocatable :: zint,kdint,zintp
  integer i,j,ij
  double precision factxp,factxm,factyp,factym,dx,dy
  character cbc*4

  !print*,'Diffusion'

  write (cbc,'(i4)') bounds_ibc

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  ! creates 2D internal arrays to store topo and kd

  allocate (zint(nx,ny),kdint(nx,ny),zintp(nx,ny))

  do j=1,ny
    do i=1,nx
      ij=(j-1)*nx+i
      zint(i,j)=h(ij)
      kdint(i,j)=kd(ij)
      if (kdsed.gt.0.d0 .and. (h(ij)-b(ij)).gt.1.d-6) kdint(i,j)=kdsed
    enddo
  enddo

  zintp = zint

  ! first pass along the x-axis

  allocate (f(nx),diag(nx),sup(nx),inf(nx),res(nx))
  f=0.d0
  diag=0.d0
  sup=0.d0
  inf=0.d0
  res=0.d0
  do j=2,ny-1
    do i=2,nx-1
    factxp=(kdint(i+1,j)+kdint(i,j))/2.d0*(dt/2.)/dx**2
    factxm=(kdint(i-1,j)+kdint(i,j))/2.d0*(dt/2.)/dx**2
      factyp=(kdint(i,j+1)+kdint(i,j))/2.d0*(dt/2.)/dy**2
    factym=(kdint(i,j-1)+kdint(i,j))/2.d0*(dt/2.)/dy**2
    diag(i)=1.d0+factxp+factxm
    sup(i)=-factxp
    inf(i)=-factxm
    f(i)=zintp(i,j)+factyp*zintp(i,j+1)-(factyp+factym)*zintp(i,j)+factym*zintp(i,j-1)
    enddo
! left bc
    if (cbc(4:4).eq.'1') then
    diag(1)=1.
    sup(1)=0.
    f(1)=zintp(1,j)
    else
    factxp=(kdint(2,j)+kdint(1,j))/2.d0*(dt/2.)/dx**2
    factyp=(kdint(1,j+1)+kdint(1,j))/2.d0*(dt/2.)/dy**2
    factym=(kdint(1,j-1)+kdint(1,j))/2.d0*(dt/2.)/dy**2
    diag(1)=1.d0+factxp
    sup(1)=-factxp
    f(1)=zintp(1,j)+factyp*zintp(1,j+1)-(factyp+factym)*zintp(1,j)+factym*zintp(1,j-1)
    endif
! right bc
    if (cbc(2:2).eq.'1') then
    diag(nx)=1.
    inf(nx)=0.
    f(nx)=zintp(nx,j)
    else
    factxm=(kdint(nx-1,j)+kdint(nx,j))/2.d0*(dt/2.)/dx**2
    factyp=(kdint(nx,j+1)+kdint(nx,j))/2.d0*(dt/2.)/dy**2
    factym=(kdint(nx,j-1)+kdint(nx,j))/2.d0*(dt/2.)/dy**2
    diag(nx)=1.d0+factxm
    inf(nx)=-factxm
    f(nx)=zintp(nx,j)+factyp*zintp(nx,j+1)-(factyp+factym)*zintp(nx,j)+factym*zintp(nx,j-1)
    endif
  call tridag (inf,diag,sup,f,res,nx)
    do i=1,nx
    zint(i,j)=res(i)
    enddo
  enddo
  deallocate (f,diag,sup,inf,res)

  ! second pass along y-axis

  allocate (f(ny),diag(ny),sup(ny),inf(ny),res(ny))
  f=0.d0
  diag=0.d0
  sup=0.d0
  inf=0.d0
  res=0.d0
  do i=2,nx-1
    do j=2,ny-1
    factxp=(kdint(i+1,j)+kdint(i,j))/2.d0*(dt/2.)/dx**2
    factxm=(kdint(i-1,j)+kdint(i,j))/2.d0*(dt/2.)/dx**2
    factyp=(kdint(i,j+1)+kdint(i,j))/2.d0*(dt/2.)/dy**2
    factym=(kdint(i,j-1)+kdint(i,j))/2.d0*(dt/2.)/dy**2
    diag(j)=1.d0+factyp+factym
    sup(j)=-factyp
    inf(j)=-factym
    f(j)=zint(i,j)+factxp*zint(i+1,j)-(factxp+factxm)*zint(i,j)+factxm*zint(i-1,j)
    enddo
! bottom bc
    if (cbc(1:1).eq.'1') then
    diag(1)=1.
    sup(1)=0.
    f(1)=zint(i,1)
    else
    factxp=(kdint(i+1,1)+kdint(i,j))/2.d0*(dt/2.)/dx**2
    factxm=(kdint(i-1,1)+kdint(i,1))/2.d0*(dt/2.)/dx**2
    factyp=(kdint(i,2)+kdint(i,1))/2.d0*(dt/2.)/dy**2
    diag(1)=1.d0+factyp
    sup(1)=-factyp
    f(1)=zint(i,1)+factxp*zint(i+1,1)-(factxp+factxm)*zint(i,1)+factxm*zint(i-1,1)
    endif
! top bc
    if (cbc(3:3).eq.'1') then
    diag(ny)=1.
    inf(ny)=0.
    f(ny)=zint(i,ny)
    else
    factxp=(kdint(i+1,ny)+kdint(i,ny))/2.d0*(dt/2.)/dx**2
    factxm=(kdint(i-1,ny)+kdint(i,ny))/2.d0*(dt/2.)/dx**2
    factym=(kdint(i,ny-1)+kdint(i,ny))/2.d0*(dt/2.)/dy**2
    diag(ny)=1.d0+factym
    inf(ny)=-factym
    f(ny)=zint(i,ny)+factxp*zint(i+1,ny)-(factxp+factxm)*zint(i,ny)+factxm*zint(i-1,ny)
    endif
  call tridag (inf,diag,sup,f,res,ny)
    do j=1,ny
    zintp(i,j)=res(j)
    enddo
  enddo
  deallocate (f,diag,sup,inf,res)

  ! stores result in 1D array

  do j=1,ny
    do i=1,nx
      ij=(j-1)*nx+i
      etot(ij)=etot(ij)+h(ij)-zintp(i,j)
      erate(ij)=erate(ij)+(h(ij)-zintp(i,j))/dt
      h(ij)=zintp(i,j)
    enddo
  enddo

  b=min(h,b)

  deallocate (zint,kdint,zintp)

  return

end subroutine Diffusion

!----------

! subroutine to solve a tri-diagonal system of equations (from Numerical Recipes)

      SUBROUTINE tridag(a,b,c,r,u,n)

      implicit none

      INTEGER n
      double precision a(n),b(n),c(n),r(n),u(n)
      INTEGER j
      double precision bet
      double precision,dimension(:),allocatable::gam

      allocate (gam(n))

      if(b(1).eq.0.d0) stop 'in tridag'

! first pass

bet=b(1)
u(1)=r(1)/bet
do 11 j=2,n
  gam(j)=c(j-1)/bet
  bet=b(j)-a(j)*gam(j)
  if(bet.eq.0.) then
    print*,'tridag failed'
    stop
  endif
  u(j)=(r(j)-a(j)*u(j-1))/bet
  11    continue

  ! second pass

  do 12 j=n-1,1,-1
    u(j)=u(j)-gam(j+1)*u(j+1)
    12    continue

    deallocate (gam)

    return

    END
