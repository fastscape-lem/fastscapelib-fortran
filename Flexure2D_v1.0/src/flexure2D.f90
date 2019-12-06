subroutine flexure (hh2,hp2,nx,ny,xl,yl,rhos2,rhoa,eet,ibc)

  ! Routine to compute the flexural response of erosion
  ! in input:
  ! hp(nx,ny), the topography (in isostatic equilibrium) before erosion, in m,
  ! h(nx,ny), the topography (out of isostatic equilibrium) after erosion, in m,
  ! rhos(nx,ny), surface rock density, in kg/m^3
  ! rhoa, asthenospheric density, in kg/m^3
  ! eet, effective elastic thickness, in m
  ! nx,ny, resolution of the input topography
  ! xl,yl, horizontal dimensions of the input topography, in m

  ! Here fixed values are assumed for:
  ! Young modulus, 1.d11 Pa
  ! Poisson ratio, 0.25
  ! g, 9.81 m/s^2

  ! the flexural, biharmonic equation is solved by FFT method (see Nunn and Aires, 1988)
  ! on a 1024x1024 mesh

  implicit none

  integer, intent(in) :: nx,ny,ibc
  double precision, intent(inout), dimension(*) :: hh2
  double precision, intent(in), dimension(*) :: hp2,rhos2
  double precision, intent(in) :: xl,yl,rhoa,eet

  integer nxflex,nyflex,i,j,ii,jj,ij
  double precision, dimension(:,:), allocatable :: w,h,hp,rhos
  double precision hx,hy,dflex,d,xk,pihx,pihy,g,fi,fj,tij,dx,dy,r,s,h1,h2,h3,h4
  double precision ddxf,ddyf,xloc,yloc,dw,xflexloc,yflexloc
  integer iflexmin,iflexmax,jflexmin,jflexmax
  character*4 cbc

  double precision, dimension(:,:,:), allocatable :: hw
  integer, dimension(:,:), allocatable :: iiw,jjw

  allocate (h(nx,ny),hp(nx,ny),rhos(nx,ny))

  do j=1,ny
    do i=1,nx
      ij=(j-1)*nx+i
      h(i,j)=hh2(ij)
      hp(i,j)=hp2(ij)
      rhos(i,j)=rhos2(ij)
    enddo
  enddo

  write (cbc,'(i4)') ibc
  if (cbc(1:1).eq.'') cbc(1:1)='0'
  if (cbc(2:2).eq.'') cbc(2:2)='0'
  if (cbc(3:3).eq.'') cbc(3:3)='0'
  if (cbc(4:4).eq.'') cbc(4:4)='0'

  ! allocate memory

  nxflex=1
  do while (nxflex.lt.nx)
    nxflex=nxflex*2
  enddo

  nyflex=1
  do while (nyflex.lt.ny)
    nyflex=nyflex*2
  enddo

  allocate (hw(4,nx,ny), iiw(nx,ny), jjw(nx,ny))
  allocate (w(nxflex,nyflex))

  ! compute relevant geometrical, flexural and spectral parameters

  iflexmin=nxflex/2-nxflex/8
  iflexmax=nxflex/2+nxflex/8
  jflexmin=nyflex/2-nyflex/8
  jflexmax=nyflex/2+nyflex/8

  dx=xl/(nx-1)
  dy=yl/(ny-1)
  ddxf=xl/(iflexmax-iflexmin)
  ddyf=yl/(jflexmax-jflexmin)
  hx=ddxf*(nxflex-1)
  hy=ddyf*(nyflex-1)
  dflex=1.d11/12.d0/(1.d0-0.25d0**2)
  d=dflex*eet**3
  g=9.81d0
  xk=rhoa*g
  pihx=3.141592654d0/hx
  pihy=3.141592654d0/hx

  ! compute weigths corresponding to the increase in topography by interpolation
  ! from the nx,ny grid to the nflex, nflex grid, using a bilinear interpolation scheme

  w=0.d0
  jj=jflexmin
  yflexloc=0.d0
    do j=1,ny
    yloc=(j-1)*dy
    if (yloc.gt.yflexloc+ddyf) jj=jj+1
    yflexloc=(jj-jflexmin)*ddyf
    ii=iflexmin
    xflexloc=0.d0
    do i=1,nx
      xloc=(i-1)*dx
      if (xloc.gt.xflexloc+ddxf) ii=ii+1
      xflexloc=(ii-iflexmin)*ddxf
      r=(xloc-xflexloc)/ddxf*2.d0-1.d0
      s=(yloc-yflexloc)/ddyf*2.d0-1.d0
      h1=(1.d0-r)*(1.d0-s)/4.d0
      h2=(1.d0+r)*(1.d0-s)/4.d0
      h3=(1.d0-r)*(1.d0+s)/4.d0
      h4=(1.d0+r)*(1.d0+s)/4.d0
      iiw(i,j)=ii
      jjw(i,j)=jj
      hw(1,i,j)=h1
      hw(2,i,j)=h2
      hw(3,i,j)=h3
      hw(4,i,j)=h4
      dw=(hp(i,j)-h(i,j))*rhos(i,j)*dx*dy*g
      w(ii,jj) = w(ii,jj) + dw*h1
      w(ii+1,jj) = w(ii+1,jj) + dw*h2
      w(ii,jj+1) = w(ii,jj+1) + dw*h3
      w(ii+1,jj+1) = w(ii+1,jj+1) + dw*h4
    enddo
  enddo

  call addw (w,nxflex,nyflex,iflexmin,iflexmax,jflexmin,jflexmax,cbc)

  ! compute FFT of weights

  do j=1,nyflex
    call sinft (w(:,j),nxflex)
  enddo

  w=transpose(w)

  do i=1,nxflex
    call sinft (w(:,i),nyflex)
  enddo

  ! apply filter to FFT of weights to simulated flexure (see Nunn and Aires, 1988)

  w=w*4./hx/hy

  do j=1,nyflex
    fj=(j*pihx)**2
    do i=1,nxflex
      fi=(i*pihx)**2
      tij=d/xk*(fi**2+2.d0*fi*fj+fj**2)+1.d0
      w(j,i)=w(j,i)/xk/tij
    enddo
  enddo

  ! compute inverse FFT of filtered weights to obtain deflection

  do i=1,nxflex
    call sinft (w(:,i),nyflex)
  enddo

  w=transpose(w)

  do j=1,nyflex
    call sinft (w(:,j),nxflex)
  enddo

  ! add  deflection by interpolation from the nflex,nflex grid to the nx,ny grid
  ! by bilinear interpolation

  do j=1,ny
    do i=1,nx
      ii=iiw(i,j)
      jj=jjw(i,j)
      h1=hw(1,i,j)
      h2=hw(2,i,j)
      h3=hw(3,i,j)
      h4=hw(4,i,j)
      h(i,j)=h(i,j)+w(ii,jj)*h1+w(ii+1,jj)*h2+w(ii,jj+1)*h3+w(ii+1,jj+1)*h4
    enddo
  enddo

      ! deallocate memory

  do j=1,ny
    do i=1,nx
      ij=(j-1)*nx+i
      hh2(ij)=h(i,j)
    enddo
  enddo

  deallocate (w,h,hp,rhos,iiw,jjw,hw)

end subroutine flexure

!-----------

subroutine addw (w,nxflex,nyflex,iflexmin,iflexmax,jflexmin,jflexmax,cbc)

  implicit none

  double precision w(nxflex,nyflex)
  integer :: nxflex,nyflex,i,j,iflexmin,iflexmax,jflexmin,jflexmax
  character*4 cbc

  if (cbc(1:1).eq.'0') w(:,jflexmin)=w(:,jflexmin+1)
  if (cbc(2:2).eq.'0') w(iflexmax,:)=w(iflexmax-1,:)
  if (cbc(3:3).eq.'0') w(:,jflexmax)=w(:,jflexmax-1)
  if (cbc(4:4).eq.'0') w(iflexmin,:)=w(iflexmin+1,:)

  do j=jflexmin,jflexmax
    do i=iflexmin,iflexmax
      if (cbc(1:1).eq.'0') w(i,jflexmin-(j-jflexmin+1))=w(i,jflexmin-(j-jflexmin+1))+w(i,j)
      if (cbc(3:3).eq.'0') w(i,jflexmax+(jflexmax-j+1))=w(i,jflexmax+(jflexmax-j+1))+w(i,j)
      if (cbc(4:4).eq.'0') w(iflexmin-(i-iflexmin+1),j)=w(iflexmin-(i-iflexmin+1),j)+w(i,j)
      if (cbc(2:2).eq.'0') w(iflexmax+(iflexmax-i+1),j)=w(iflexmax+(iflexmax-i+1),j)+w(i,j)
      if (cbc(1:1).eq.'0'.and.cbc(2:2).eq.'0') w(iflexmax+(iflexmax-i+1),jflexmin-(j-jflexmin+1))= &
      w(iflexmax+(iflexmax-i+1),jflexmin-(j-jflexmin+1))+w(i,j)
      if (cbc(2:2).eq.'0'.and.cbc(3:3).eq.'0') w(iflexmax+(iflexmax-i+1),jflexmax+(jflexmax-j+1))= &
      w(iflexmax+(iflexmax-i+1),jflexmax+(jflexmax-j+1))+w(i,j)
      if (cbc(3:3).eq.'0'.and.cbc(4:4).eq.'0') w(iflexmin-(i-iflexmin+1),jflexmax+(jflexmax-j+1))= &
      w(iflexmin-(i-iflexmin+1),jflexmax+(jflexmax-j+1))+w(i,j)
      if (cbc(4:4).eq.'0'.and.cbc(1:1).eq.'0') w(iflexmin-(i-iflexmin+1),jflexmin-(j-jflexmin+1))= &
      w(iflexmin-(i-iflexmin+1),jflexmin-(j-jflexmin+1))+w(i,j)
    enddo
  enddo

  return

end subroutine addw
