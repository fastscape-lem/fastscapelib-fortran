subroutine FlowRouting

  ! subroutine to compute all flow routing arrays using both single
  ! and multiple flow

  use FastScapeContext

  implicit none

  integer i, j, i1, i2, j1, j2, ii, jj, iii, jjj, ijk, ij
  double precision smax, l, slope, dx, dy
  character*4 cbc

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

  ! resets receiver and distance between node and its receiver

  do ij=1,nn
    rec(ij)=ij
    length(ij)=0.d0
  enddo

  dx = xl/(nx - 1)
  dy = yl/(ny - 1)

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
            if (l.eq.0.d0) print*,ij,ijk,l,ii,jj,dx,dy
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
  dx = xl/(nx-1)
  dy = yl/(ny-1)
  call LocalMinima (stack,rec,bc,ndon,don,h,length,nx,ny,dx,dy)

  ! computes receiver and stack information for mult-direction flow
  call find_mult_rec (h,rec,stack,hwater,mrec,mnrec,mwrec,mlrec,mstack,nx,ny,dx,dy,p,ibc,distributed_flow)

  if (.not.distributed_flow) then
    where (mnrec.ne.0) rec=mrec(1,:)
  endif

end subroutine FlowRouting

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
end subroutine find_stack

!-----------------------------------------------------

subroutine find_mult_rec (h,rec0,stack0,water,rec,nrec,wrec,lrec,stack,nx,ny,dx,dy,p,ibc,distributed_flow)

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
  double  precision :: slopemax,pp,sumweight,deltah,smax,slope,l
  integer, dimension(:), allocatable :: ndon,vis,parse
  integer, dimension(:,:), allocatable :: don
  double precision, dimension(:), allocatable :: h0
  character cbc*4
  logical xcyclic,ycyclic,distributed_flow

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

  if (distributed_flow) then

    nrec = 0
    wrec = 0.d0

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
      if (abs(pp+1.d0).lt.tiny(pp)) then
        slope = 0.d0
        if (nrec(ij).ne.0) slope = real(sum(wrec(1:nrec(ij),ij))/nrec(ij))
        pp = 0.5 + 0.6*slope
      endif
      do k=1,nrec(ij)
        wrec(k,ij) = wrec(k,ij)**pp
      enddo
      sumweight = sum(wrec(1:nrec(ij),ij))
      wrec(1:nrec(ij),ij) = wrec(1:nrec(ij),ij)/sumweight
    enddo

  else

    nrec = 0
    wrec = 0.d0

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
              slope=(h0(ij)-h0(ijk))/l
              if (slope.gt.smax) then
                smax=slope
                nrec(ij) = 1
                rec(1,ij) = ijk
                lrec(1,ij) = l
                wrec(1,ij) = 1.d0
              endif
            endif
          enddo
        enddo
      enddo
    enddo

  endif

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
