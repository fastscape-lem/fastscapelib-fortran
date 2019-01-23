subroutine flexure (h,hp,nx,ny,xl,yl,rhos,rhoa,eet,ibc)

! Routine to compute the flexural response of erosin
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

integer nx,ny,ibc
double precision h(nx,ny),hp(nx,ny),xl,yl,rhos(nx,ny)
integer nflex,i,j,ii,jj
double precision, dimension(:,:), allocatable :: w
double precision hx,dflex,d,rhoa,xk,pihx,g,fi,fj,tij,rat,dx,dy,eet,r,s,h1,h2,h3,h4,pi
double precision ddxf,ddyf,xloc,yloc,hh,hhp,rrhos,dw
integer iflexmin,iflexmax,jflexmin,jflexmax
character*4 cbc

write (cbc,'(i4)') ibc

! allocate memory

nflex=1
do while (nflex.lt.max(nx,ny))
nflex=nflex*2
enddo

allocate (w(nflex,nflex))

! compute relevant geometrical, flexural and spectral parameters

dx=xl/(nx-1)
dy=yl/(ny-1)
hx=max(7.d0*xl,7.d0*yl)
ddxf=hx/(nflex-1)
ddyf=hx/(nflex-1)
dflex=1.d11/12.d0/(1.d0-0.25d0**2)
d=dflex*eet**3
g=9.81d0
xk=rhoa*g
pihx=3.141592654d0/hx

! compute weigths corresponding to the increase in topography by interpolation
! from the nx,ny grid to the nflex, nflex grid, using a bilinear interpolation scheme

iflexmin=1+(0.d0-xl/2.d0+hx/2.d0)/ddxf
iflexmax=1+(xl-xl/2.d0+hx/2.d0)/ddxf
jflexmin=1+(0.d0-yl/2.d0+hx/2.d0)/ddyf
jflexmax=1+(yl-yl/2.d0+hx/2.d0)/ddyf

w=0.d0
  do j=jflexmin,jflexmax
  yloc=(j-1)*hx/(nflex-1)+yl/2.d0-hx/2.d0
  jj=1+(ny-1)*yloc/yl
  jj=min(jj,ny-1)
  jj=max(1,jj)
    do i=iflexmin,iflexmax
    xloc=(i-1)*hx/(nflex-1)+xl/2.d0-hx/2.d0
    ii=1+(nx-1)*xloc/xl
    ii=min(ii,nx-1)
    ii=max(1,ii)
    r=(xloc-(ii-1)*dx)/dx*2.d0-1.d0
    s=(yloc-(jj-1)*dy)/dy*2.d0-1.d0
    h1=(1.d0-r)*(1.d0-s)/4.d0
    h2=(1.d0+r)*(1.d0-s)/4.d0
    h3=(1.d0-r)*(1.d0+s)/4.d0
    h4=(1.d0+r)*(1.d0+s)/4.d0
    hh=h(ii,jj)*h1+h(ii+1,jj)*h2+h(ii,jj+1)*h3+h(ii+1,jj+1)*h4
    hhp=hp(ii,jj)*h1+hp(ii+1,jj)*h2+hp(ii,jj+1)*h3+hp(ii+1,jj+1)*h4
    rrhos=rhos(ii,jj)*h1+rhos(ii+1,jj)*h2+rhos(ii,jj+1)*h3+rhos(ii+1,jj+1)*h4
	dw=(hhp-hh)*rrhos*ddxf*ddyf*g
    w(i,j)=w(i,j)+dw
! modifications brought by Jean to account for boundary conditions (request from Laure) Dec 2018
!goto 1111
	if (cbc(1:1).eq.'0') w(i,jflexmin-(j-jflexmin+1))=w(i,jflexmin-(j-jflexmin+1))+dw
	if (cbc(3:3).eq.'0') w(i,jflexmax+(jflexmax-j+1))=w(i,jflexmax+(jflexmax-j+1))+dw
	if (cbc(4:4).eq.'0') w(iflexmin-(i-iflexmin+1),j)=w(iflexmin-(i-iflexmin+1),j)+dw
	if (cbc(2:2).eq.'0') w(iflexmax+(iflexmax-i+1),j)=w(iflexmax+(iflexmax-i+1),j)+dw
	if (cbc(1:1).eq.'0'.and.cbc(2:2).eq.'0') w(iflexmax+(iflexmax-i+1),jflexmin-(j-jflexmin+1))= &
		w(iflexmax+(iflexmax-i+1),jflexmin-(j-jflexmin+1))+dw
	if (cbc(2:2).eq.'0'.and.cbc(3:3).eq.'0') w(iflexmax+(iflexmax-i+1),jflexmax+(jflexmax-j+1))= &
		w(iflexmax+(iflexmax-i+1),jflexmax+(jflexmax-j+1))+dw
	if (cbc(3:3).eq.'0'.and.cbc(4:4).eq.'0') w(iflexmin-(i-iflexmin+1),jflexmax+(jflexmax-j+1))= &
		w(iflexmin-(i-iflexmin+1),jflexmax+(jflexmax-j+1))+dw
	if (cbc(4:4).eq.'0'.and.cbc(1:1).eq.'0') w(iflexmin-(i-iflexmin+1),jflexmin-(j-jflexmin+1))= &
		w(iflexmin-(i-iflexmin+1),jflexmin-(j-jflexmin+1))+dw
1111 continue
    enddo
  enddo

! compute FFT of weights

  do j=1,nflex
  call sinft (w(:,j),nflex)
  enddo

w=transpose(w)

  do i=1,nflex
  call sinft (w(:,i),nflex)
  enddo

! apply filter to FFT of weights to simulated flexure (see Nunn and Aires, 1988)

w=w*4./hx/hx

  do j=1,nflex
  fj=(j*pihx)**2
    do i=1,nflex
    fi=(i*pihx)**2
    tij=d/xk*(fi**2+2.d0*fi*fj+fj**2)+1.d0
    w(j,i)=w(j,i)/xk/tij
    enddo
  enddo

! compute inverse FFT of filtered weights to obtain deflection

  do i=1,nflex
  call sinft (w(:,i),nflex)
  enddo

w=transpose(w)

  do j=1,nflex
  call sinft (w(:,j),nflex)
  enddo

! add  deflection by interpolation from the nflex,nflex grid to the nx,ny grid
! by bilinear interpolation

  do j=1,ny
  yloc=(j-1)*dy+hx/2.d0-yl/2.d0
  jj=1+(nflex-1)*yloc/hx
  jj=min(jj,nflex-1)
  jj=max(1,jj)
    do i=1,nx
    xloc=(i-1)*dx+hx/2.d0-xl/2.d0
    ii=1+(nflex-1)*xloc/hx
    ii=min(ii,nflex-1)
    ii=max(1,ii)
    r=(xloc-(ii-1)*ddxf)/ddxf*2.d0-1.d0
    s=(yloc-(jj-1)*ddyf)/ddyf*2.d0-1.d0
    h1=(1.d0-r)*(1.d0-s)/4.d0
    h2=(1.d0+r)*(1.d0-s)/4.d0
    h3=(1.d0-r)*(1.d0+s)/4.d0
    h4=(1.d0+r)*(1.d0+s)/4.d0
    h(i,j)=h(i,j)+w(ii,jj)*h1+w(ii+1,jj)*h2+w(ii,jj+1)*h3+w(ii+1,jj+1)*h4
    enddo
  enddo

! deallocate memory

deallocate (w)

end subroutine flexure
