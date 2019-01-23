subroutine StreamPowerLaw ()

! subroutine to solve the stream power law equation following the FastScape method described
! in Braun and Willett, Geomorphology, 2015

use FastScapeContext

implicit none

integer :: i,j,ij,ii,jj,iii,jjj,ijk,ijr,i1,i2,j1,j2,np,n_local_minima,k
double precision :: dx,dy,smax,l,slope,fact,diff,tol,err,w_rcv
double precision :: f,df,errp,h0,hn,omega,tolp,depos
character cbc*4
logical xcyclic,ycyclic
double precision, dimension(:), allocatable :: rhs,ht,g,kfint,bt,dh,hp
double precision, dimension(:), allocatable :: elev,hwater
logical, dimension(:), allocatable :: bc
integer, dimension(:), allocatable :: mnrec,mstack
integer, dimension(:,:), allocatable :: mrec
double precision, dimension(:,:), allocatable :: mwrec,mlrec

allocate (rhs(nn),ht(nn),g(nn),kfint(nn),bc(nn),bt(nn),dh(nn),hp(nn))
allocate (elev(nn))

! sets boundary conditions

write (cbc,'(i4)') ibc
bc=.FALSE.
i1=1
i2=nx
j1=1
j2=ny
if (cbc(4:4).eq.'1') i1=2
if (cbc(2:2).eq.'1') i2=nx-1
if (cbc(1:1).eq.'1') j1=2
if (cbc(3:3).eq.'1') j2=ny-1
if (cbc(4:4).eq.'1') bc(1:nn:nx)=.TRUE.
if (cbc(2:2).eq.'1') bc(nx:nn:nx)=.TRUE.
if (cbc(1:1).eq.'1') bc(1:nx)=.TRUE.
if (cbc(3:3).eq.'1') bc(nx*(ny-1)+1:nn)=.TRUE.
xcyclic=.FALSE.
ycyclic=.FALSE.
if (cbc(4:4).ne.'1'.and.cbc(2:2).ne.'1') xcyclic=.TRUE.
if (cbc(1:1).ne.'1'.and.cbc(3:3).ne.'1') ycyclic=.TRUE.

dx=xl/(nx-1)
dy=yl/(ny-1)

! defines g, dimensionless parameter for sediment transport and deposition
g=g1
kfint=kf
if (g2.gt.0.d0) where ((h-b).gt.1.d-6) g=g2
if (kfsed.gt.0.d0) where ((h-b).gt.1.d-6) kfint=kfsed

! uplift

h=h+u*dt

! resets receiver and distance between node and its receiver

  do ij=1,nn
  rec(ij)=ij
  length(ij)=0.d0
  enddo

! finds receiver using steepest descent/neighbour method

  do j=j1,j2
    do i=i1,i2
    ij=i+(j-1)*nx
      smax=tiny(smax)
      do jj=-1,1
        do ii=-1,1
        iii=i+ii
        if (iii.lt.1.and.xcyclic) iii=iii+nx
        iii=max(iii,1)
        if (iii.gt.nx.and.xcyclic) iii=iii-nx
        iii=min(iii,nx)
        jjj=j+jj
        if (jjj.lt.1.and.ycyclic) jjj=jjj+ny
        jjj=max(jjj,1)
        if (jjj.gt.ny.and.ycyclic) jjj=jjj-ny
        jjj=min(jjj,ny)
        ijk=iii+(jjj-1)*nx
          if (ijk.ne.ij) then
          l=sqrt((dx*ii)**2+(dy*jj)**2)
          slope=(h(ij)-h(ijk))/l
            if (slope.gt.smax) then
            smax=slope
            rec(ij)=ijk
            length(ij)=l
            endif
          endif
        enddo
      enddo
    enddo
  enddo

! inverts receiver array to compute donor arrays
ndon=0
  do ij=1,nn
    if (rec(ij).ne.ij) then
    ijk=rec(ij)
    ndon(ijk)=ndon(ijk)+1
    don(ndon(ijk),ijk)=ij
    endif
  enddo

! computes stack by recursion
  nstack=0
  catch=catch0
    do ij=1,nn
      if (rec(ij).eq.ij) then
      nstack=nstack+1
      stack(nstack)=ij
      call find_stack (ij,don,ndon,nn,stack,nstack,catch)
    endif
  enddo

! removes local minima
call LocalMinima (stack,rec,bc,ndon,don,h,length,nx,ny,dx,dy)

! computes receiver and stack information for mult-direction flow
allocate (mrec(8,nn),mnrec(nn),mwrec(8,nn),mlrec(8,nn),mstack(nn),hwater(nn))
call find_mult_rec (h,rec,stack,hwater,mrec,mnrec,mwrec,mlrec,mstack,nx,ny,dx,dy,p,ibc)

if (count(mstack==0).ne.0) print*,'incomplete stack',count(mstack==0),nn

! computes drainage area
a=dx*dy*precip
  do ij=1,nn
  ijk=mstack(ij)
        do k =1,mnrec(ijk)
    a(mrec(k,ijk))=a(mrec(k,ijk))+a(ijk)*mwrec(k,ijk)
        enddo
  enddo

! calculate the elevation / SPL, including sediment flux
tol=1.d-3
err=2*tol

! store the elevation at t
ht=hwater
bt=b

! Gauss-Seidel iteration
nGSStreamPowerLaw=0

do while (err.gt.tol)
123 continue
   nGSStreamPowerLaw=nGSStreamPowerLaw+1
   ! guess/update the elevation at t+Δt (k)
   hp=hwater
   ! update the base elevation
   b=min(hp,bt+u*dt)
   ! calculate erosion/deposition at each node
   dh=ht-hp

   ! sum the erosion in stack order
   do ij=1,nn
      ijk=mstack(ij)
      do k=1,mnrec(ijk)
         ijr=mrec(k,ijk)
         dh(ijr)=dh(ijr)+dh(ijk)*mwrec(k,ijk)
      enddo
   enddo

   where (bc)
      elev=ht
   elsewhere
      elev=ht+(dh-(ht-hp))*g*dx*dy/a
   endwhere

   ! apply modified stream power law using lake surface (hwater)

   if (abs(n-1.d0).lt.tiny(n)) then

      do ij=nn,1,-1
         ijk=mstack(ij)
         if (mnrec(ijk).gt.0) then
            if (ht(ijk).ge.sealevel) then
               f = elev(ijk)
               df = 1.d0
               do k=1,mnrec(ijk)
                  if (ht(ijk).ge.ht(mrec(k,ijk))) then
                     fact = kfint(ijk)*dt*(a(ijk)*mwrec(k,ijk))**m/mlrec(k,ijk)
                     f = f + fact*hwater(mrec(k,ijk))
                     df = df + fact
                  endif
               enddo
               hwater(ijk)=f/df
            endif
         endif
      enddo

   else

      do ij=nn,1,-1
         ijk=mstack(ij)
         if (mnrec(ijk).gt.0) then
            if (ht(ijk).ge.sealevel) then
               omega=0.875d0/n
               tolp=1.d-3
               errp=2.d0*tolp
               h0=elev(ijk)
               do while (errp.gt.tolp)
                  f=hwater(ijk)-h0
                  df=1.d0
                  do k=1,mnrec(ijk)
                     if (ht(ijk).gt.ht(mrec(k,ijk))) then
                        fact = kfint(ijk)*dt*(a(ijk)*mwrec(k,ijk))**m/mlrec(k,ijk)**n
                        f=f+fact*max(0.d0,hwater(ijk)-hwater(mrec(k,ijk)))**n
                        df=df+fact*n*max(0.d0,hwater(ijk)-hwater(mrec(k,ijk)))**(n-1.d0)
                     endif
                  enddo
                  hn=hwater(ijk)-f/df
                  errp=abs(hn-hwater(ijk))
                  hwater(ijk)=hwater(ijk)*(1.d0-omega)+hn*omega
               enddo
            endif
         endif
      enddo

   endif

   err=maxval(abs(hwater-hp))

        if (nGSStreamPowerLaw.ge.100) then
        stop 'too many Gauss-Seidl iterations in modified SPL'
        endif

        enddo

! transpose topography from lake surface to topographic surface
h = h + hwater - ht

deallocate (mrec,mwrec,mlrec,mnrec,mstack,hwater)

! stores total erosion, erosion rate and flux for output
etot=etot+ht-h
erate=(ht-h)/dt
Sedflux=ht-h

!deallocate (rhs,hn,bt,g,kf,h0)

return

end subroutine StreamPowerLaw

!----------

recursive subroutine find_stack (ij,don,ndon,nn,stack,nstack,catch)

! recursive routine to go through all nodes following donor information

implicit none
integer don(8,nn),ndon(nn),stack(nn)
double precision catch(nn)
integer k,ij,ijk,nn,nstack

  do k=1,ndon(ij)
  ijk=don(k,ij)
  nstack=nstack+1
  stack(nstack)=ijk
  catch(ijk)=catch(ij)
  call find_stack (ijk,don,ndon,nn,stack,nstack,catch)
  enddo

return
end

!-----------------------------------------------------

subroutine find_mult_rec (h,rec0,stack0,water,rec,nrec,wrec,lrec,stack,nx,ny,dx,dy,p,ibc)

! subroutine to find multiple receiver information
! in input:
! h is topography
! rec0 is single receiver information
! stack0 is stack (from bottom to top) obtained by using single receiver information
! water is the surface of the lakes (or topography where there is no lake)
! nx, ny resolution in x- and y-directions
! dx, dy grid spacing in x- and y-directions
! p is exponent to which the slope is put to share the water/sediment among receivers
! ibc boundary conditions (1111 to 0000)
! in output:
! rec: multiple receiver information
! nrec: number of receivers for each node
! wrec: weight for each receiver
! lrec: distance to each receiver
! stack: stoack order for multiple receivers (from top to bottom)

integer nx,ny,ibc
double precision h(nx*ny),wrec(8,nx*ny),lrec(8,nx*ny),dx,dy,p,water(nx*ny)
integer rec(8,nx*ny),nrec(nx*ny),stack(nx*ny),rec0(nx*ny),stack0(nx*ny)

integer :: nn,i,j,ii,jj,iii,jjj,ijk,k,ijr,nparse,nstack,ijn,i1,i2,j1,j2
double  precision :: slopemax,pp,sumweight,deltah
integer, dimension(:), allocatable :: ndon,vis,parse
integer, dimension(:,:), allocatable :: don
double precision, dimension(:), allocatable :: h0
character cbc*4
logical xcyclic,ycyclic

nn=nx*ny

! set bc

write (cbc,'(i4)') ibc
i1=1
i2=nx
j1=1
j2=ny
if (cbc(4:4).eq.'1') i1=2
if (cbc(2:2).eq.'1') i2=nx-1
if (cbc(1:1).eq.'1') j1=2
if (cbc(3:3).eq.'1') j2=ny-1
xcyclic=.FALSE.
ycyclic=.FALSE.
if (cbc(4:4).ne.'1'.and.cbc(2:2).ne.'1') xcyclic=.TRUE.
if (cbc(1:1).ne.'1'.and.cbc(3:3).ne.'1') ycyclic=.TRUE.

allocate (h0(nn))

h0=h

! fill the local minima with a nearly planar surface

deltah = 1.d-8
        do i=1,nn
        ijk=stack0(i)
        ijr=rec0(ijk)
                if (ijr.ne.0) then
                        if (h0(ijr).gt.h0(ijk)) then
                        h0(ijk)=h0(ijr)+deltah
                        endif
                endif
        enddo

water = h0

nrec=0
wrec=0.d0

! loop on all nodes
        do j=j1,j2
                do i=i1,i2
                ij = (j-1)*nx + i
                slopemax = 0.
                        do jj=-1,1
                        jjj= j + jj
                        if (jjj.lt.1.and.ycyclic) jjj=jjj+ny
                        jjj=max(jjj,1)
                        if (jjj.gt.ny.and.ycyclic) jjj=jjj-ny
                        jjj=min(jjj,ny)
                                do ii=-1,1
                                iii = i + ii
                                if (iii.lt.1.and.xcyclic) iii=iii+nx
                                iii=max(iii,1)
                                if (iii.gt.nx.and.xcyclic) iii=iii-nx
                                iii=min(iii,nx)
                                ijk = (jjj-1)*nx + iii
                                        if (h0(ij).gt.h0(ijk)) then
                                        nrec(ij)=nrec(ij)+1
                                        rec(nrec(ij),ij) = ijk
                                        lrec(nrec(ij),ij) = sqrt((ii*dx)**2 + (jj*dy)**2)
                                        wrec(nrec(ij),ij) = (h0(ij) - h0(ijk))/lrec(nrec(ij),ij)
                                        endif
                                enddo
                        enddo
                enddo
        enddo

        do ij =1,nn
        pp = p
                if (pp<0.d0) then
                slope = 0.d0
                if (nrec(ij).ne.0) slope = sum(wrec(1:nrec(ij),ij))/nrec(ij)
                pp = 0.5 + 0.6*slope
                endif
                do k=1,nrec(ij)
                wrec(k,ij) = wrec(k,ij)**pp
                enddo
        sumweight = sum(wrec(1:nrec(ij),ij))
        wrec(1:nrec(ij),ij) = wrec(1:nrec(ij),ij)/sumweight
        enddo

allocate (ndon(nn),don(8,nn))

ndon=0

        do ij=1,nn
                do k=1,nrec(ij)
                ijk = rec(k,ij)
                ndon(ijk)=ndon(ijk)+1
                don(ndon(ijk),ijk) = ij
                enddo
        enddo

allocate (vis(nn),parse(nn))

nparse=0
nstack=0
stack=0
vis=0
parse=0

! we go through the nodes
        do ij=1,nn
        ! when we find a "summit" (ie a node that has no donors)
        ! we parse it (put it in a stack called parse)
                if (ndon(ij).eq.0) then
                nparse=nparse+1
                parse(nparse)=ij
                endif
        ! we go through the parsing stack
                do while (nparse.gt.0)
                ijn=parse(nparse)
                nparse=nparse-1
                ! we add the node to the stack
                nstack=nstack+1
                stack(nstack)=ijn
                ! for each of its receivers we increment a counter called vis
                        do  ijk=1,nrec(ijn)
                        ijr=rec(ijk,ijn)
                        vis(ijr)=vis(ijr)+1
                        ! if the counter is equal to the number of donors for that node we add it to the parsing stack
                                if (vis(ijr).eq.ndon(ijr)) then
                                nparse=nparse+1
                                parse(nparse)=ijr
                                endif
                        enddo
                enddo
        enddo
if (nstack.ne.nn) stop 'error in stack'

deallocate (ndon,don,vis,parse,h0)

return

end subroutine find_mult_rec
