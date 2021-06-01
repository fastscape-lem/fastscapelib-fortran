#include "Error.fpp"
subroutine Marine(ierr)

  ! Marine transport component
  ! using silt and sand coupling diffusion solver
  ! developed by Xiaoping Yuan (2017-2018)

  use FastScapeContext

  implicit none

  double precision, dimension(:), allocatable :: flux,shelfdepth,ht,Fs,dh,dh1,dh2,Fmixt,mwater
  double precision, dimension(:), allocatable :: dhs, dhs1, F1, F2, zi, zo
  integer, dimension(:), allocatable :: flag,mmnrec,mmstack
  integer, dimension(:,:), allocatable :: mmrec
  double precision, dimension(:,:), allocatable :: mmwrec,mmlrec
  double precision shelfslope,ratio1,ratio2,dx,dy
  integer ij,ijr,ijk,k
  integer, intent(inout):: ierr

  allocate (flux(nn),shelfdepth(nn),ht(nn),Fs(nn),dh(nn),dh1(nn),dh2(nn),Fmixt(nn),flag(nn))
  allocate (dhs(nn),dhs1(nn),F1(nn),F2(nn),zi(nn),zo(nn))

  ! set nodes at transition between ocean and continent
  flag=0

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  ! computing flux from continental erosion
  flux=0.d0
  where (h.gt.sealevel) flux=Sedflux
  do ij=nn,1,-1
    ijk=stack(ij)
    ijr=rec(ijk)
    if (ijr.ne.ijk.and.h(ijk).gt.sealevel) then
      flux(ijr)=flux(ijr)+flux(ijk)
    endif
  enddo
  ! here the integral of erosion/deposition has been done
  ! and distributed as flux to ocean
  where (h.gt.sealevel) flux=0.d0

  ! set nodes at transition between ocean and continent
  !where (flux.gt.tiny(flux)) flag=1

  ! decompact volume of pure solid phase (silt and sand) from onshore
  ratio1=ratio/(1.d0-poro1)
  ratio2=(1.d0-ratio)/(1.d0-poro2)
  ! total volume of silt and sand after decompaction
  flux=flux*(ratio1+ratio2)

  ! modifications made by Jean for multiple flow to distribute continental flux to ocean on the shelf
  ! Dec 2018

  allocate (mmrec(8,nn),mmnrec(nn),mmwrec(8,nn),mmlrec(8,nn),mmstack(nn),mwater(nn))

  call find_mult_rec (h,rec,stack,mwater,mmrec,mmnrec,mmwrec,mmlrec,mmstack,nx,ny,dx,dy,0.d0,p_mfd_exp, &
    bounds_i1, bounds_i2, bounds_j1, bounds_j2, bounds_xcyclic, bounds_ycyclic, ierr)

  !print*,count(flux>0.and.mmnrec==0),count(flux>0),count(mmstack==0)

  ! modifications made by Jean
  ! to compute shelf depth
  shelfdepth=sealevel
  shelfslope=-1.d-4
  do ij=1,nn
    ijk=mmstack(ij)
    do k=1,mmnrec(ijk)
      ijr=mmrec(k,ijk)
      if (h(ijk).lt.sealevel) then
        shelfdepth(ijr)=min(shelfdepth(ijr),shelfdepth(ijk)+mmlrec(k,ijk)*shelfslope)
        shelfdepth(ijr)=max(shelfdepth(ijr),h(ijr))
      endif
    enddo
  enddo
  ! end modifications

  ! passes the flux across the shelf
  ! modifications made by Jean

  where (h.lt.sealevel) flux=flux+(h-shelfdepth)
  do ij=1,nn
    ijk=mmstack(ij)
    do k=1,mmnrec(ijk)
      ijr=mmrec(k,ijk)
      flux(ijr)=flux(ijr)+max(0.d0,flux(ijk)*mmwrec(k,ijk))
    enddo
  enddo
  ! modifications made by Jean

  deallocate (mmrec,mmnrec,mmwrec,mmlrec,mmstack,mwater)

  where (flux.gt.0.d0.and.h.lt.sealevel) flux=-(h-shelfdepth)
  where (flux.le.0.d0.and.h.lt.sealevel) flux=flux-(h-shelfdepth)
  where (h.ge.sealevel) flux=0.d0
  flux=max(flux,0.d0)

  ! silt fraction (after decompaction) in shelf
  Fs=0.d0
  where (flux.gt.0.d0) Fs=ratio1/(ratio1+ratio2)

  ! scales flux by time step
  flux=flux/dt

  ! stores initial height and fraction
  ht=h
  Fmixt=Fmix

  !print*,'flux',minval(flux),sum(flux)/nx/ny,maxval(flux)
  !print*,'Fmix',minval(Fmix),sum(Fmix)/nx/ny,maxval(Fmix)

  ! silt and sand coupling diffusion in ocean
  call SiltSandCouplingDiffusion (h,Fmix,flux*Fs,flux*(1.d0-Fs), &
  nx,ny,dx,dy,dt,sealevel,layer,kdsea1,kdsea2,nGSMarine,flag,bounds_ibc,ierr);FSCAPE_CHKERR(ierr)

  ! pure silt and sand during deposition/erosion
  dh1=((h-ht)*Fmix+layer*(Fmix-Fmixt))*(1.d0-poro1)
  dh2=((h-ht)*(1.d0-Fmix)+layer*(Fmixt-Fmix))*(1.d0-poro2)
  dh=dh1+dh2

  ! >>>>>>>> compaction starts added by Jean (Dec 2018)

  ! sum of pure silt and solid phase
  if (step.eq.0) then
    dhs1=dh1
    dhs=dh
  else
    dhs1=dhs1+dh1
    dhs=dhs+dh
  endif
  where (dhs1.lt.0.d0) dhs1=0.d0
  where (dhs.lt.0.d0) dhs=0.d0

  ! calculate the average silt (and sand) fraction in ocean part
  F1=0.d0;F2=0.d0
  where (h.le.sealevel.and.dhs.gt.0.d0) F1=dhs1/dhs
  F1=max(0.d0,F1);F1=min(1.d0,F1)
  where (h.le.sealevel.and.dhs.gt.0.d0) F2=1.d0-F1

  ! calculate the thickness after compaction, initial thickness of sediments
  !zi=ht-b
  !call compaction (F1,F2,poro1,poro2,zporo1,zporo2,nn,dh,zi,zo)
  ! update the elevation
  !h=b+zo

  ! >>>>>>>> compaction ends

  ! update the elevation
  !h=ht+dh
  etot=etot+ht-h
  erate=erate+(ht-h)/dt
  where (h.lt.sealevel) Sedflux=0.d0
  where (h.lt.sealevel) etot=0.d0
  where (h.lt.sealevel) erate=0.d0

  ! set the silt fraction in continent
  where (h.ge.sealevel+1.d-3) Fmix=0.d-1
  Fmix=max(0.d0,Fmix)
  Fmix=min(1.d0,Fmix)

  ! updates basement
  b=min(h,b)

  deallocate (flux,shelfdepth,ht,Fs,dh,dh1,dh2,Fmixt)

  return

end subroutine Marine

!----------------------------------------------------------------------------------

subroutine SiltSandCouplingDiffusion (h,f,Q1,Q2,nx,ny,dx,dy,dt, &
  sealevel,L,kdsea1,kdsea2,niter,flag,ibc,ierr)

  use FastScapeErrorCodes

  implicit none

  ! define the parameters
  integer i,j,ij,ipj,imj,ijp,ijm,nn,nx,ny,niter,ibc
  double precision h(nx*ny),f(nx*ny),Q1(nx*ny),Q2(nx*ny)
  double precision, dimension(:), allocatable :: hp,fp,ht,ft,hhalf,fhalf,fhalfp
  double precision, dimension(:), allocatable :: diag,sup,inf,rhs,res,tint
  integer flag(nx*ny)

  double precision dx,dy,dt,sealevel,L,kdsea1,kdsea2
  double precision K1,K2,tol,err1,err2
  double precision Ap,Bp,Cp,Dp,Ep,Mp,Np
  integer, intent(inout):: ierr

  character cbc*4

  write (cbc,'(i4)') ibc

  K1=kdsea1
  K2=kdsea2

  nn=nx*ny

  allocate (hp(nn),fp(nn),ht(nn),ft(nn),hhalf(nn),fhalf(nn),fhalfp(nn),tint(nn))

  ! initilize the elevation and silt fraction at time t and t+dt/2
  ht=h
  ft=f
  hhalf=h
  fhalf=f

  ! tolerance is in m
  tol=1.d0
  err1=2*tol
  err2=2*tol
  niter=0

  ! iteration until convergence is reached
  do while (err1.gt.tol)
    ! update the elevation and silt fraction during each iteration
    hp=h
    fp=f
    fhalfp=fhalf
    niter=niter+1
    ! calculate the elevation h in x-direction
    allocate (diag(nx),sup(nx),inf(nx),rhs(nx),res(nx))
    do j=2,ny-1
      do i=1,nx
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          if (i.eq.1) then
            if (cbc(4:4).eq.'1') then
              diag(i)=1.d0
              sup(i)=0.d0
              rhs(i)=ht(ij)
            else
              Ap=dt/2.d0*(K2+(K1-K2)*(fhalfp(ipj)+fhalfp(ij))/2.d0)/dx**2
              diag(i)=1.d0+Ap
              sup(i)=-Ap
              Cp=dt/2.d0*(K2+(K1-K2)*(ft(ijp)+ft(ij))/2.d0)*(ht(ijp)-ht(ij))/dy**2 &
              -dt/2.d0*(K2+(K1-K2)*(ft(ij)+ft(ijm))/2.d0)*(ht(ij)-ht(ijm))/dy**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(i)=Cp+ht(ij)
            endif
          elseif (i.eq.nx) then
            if (cbc(2:2).eq.'1') then
              diag(i)=1.d0
              inf(i)=0.d0
              rhs(i)=ht(ij)
            else
              Bp=-dt/2.d0*(K2+(K1-K2)*(fhalfp(ij)+fhalfp(imj))/2.d0)/dx**2
              diag(i)=1.d0-Bp
              inf(i)=Bp
              Cp=dt/2.d0*(K2+(K1-K2)*(ft(ijp)+ft(ij))/2.d0)*(ht(ijp)-ht(ij))/dy**2 &
              -dt/2.d0*(K2+(K1-K2)*(ft(ij)+ft(ijm))/2.d0)*(ht(ij)-ht(ijm))/dy**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(i)=Cp+ht(ij)
            endif
          else
            Ap=dt/2.d0*(K2+(K1-K2)*(fhalfp(ipj)+fhalfp(ij))/2.d0)/dx**2
            Bp=-dt/2.d0*(K2+(K1-K2)*(fhalfp(ij)+fhalfp(imj))/2.d0)/dx**2
            diag(i)=1.d0+Ap-Bp
            sup(i)=-Ap
            inf(i)=Bp
            Cp=dt/2.d0*(K2+(K1-K2)*(ft(ijp)+ft(ij))/2.d0)*(ht(ijp)-ht(ij))/dy**2 &
            -dt/2.d0*(K2+(K1-K2)*(ft(ij)+ft(ijm))/2.d0)*(ht(ij)-ht(ijm))/dy**2 &
            +(Q1(ij)+Q2(ij))*dt/2.d0
            rhs(i)=Cp+ht(ij)
          endif
          ! in continent
        else
          diag(i)=1.d0
          sup(i)=0.d0
          inf(i)=0.d0
          rhs(i)=ht(ij)
        endif
      enddo
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,nx)
      do i=1,nx
        ij=(j-1)*nx+i
        hhalf(ij)=res(i)
      enddo
    enddo
    tint=hhalf
    ! the corner nodes (1,1) and (1,ny)
    hhalf(1)=hhalf(2)
    hhalf((ny-1)*nx+1)=hhalf((ny-1)*nx+2)
    ! the corner nodes (nx,1) and (nx,ny)
    hhalf(nx)=hhalf(nx-1)
    hhalf(nx*ny)=hhalf(nx*ny-1)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the silt fraction F in x-direction
    allocate (diag(nx),sup(nx),inf(nx),rhs(nx),res(nx))
    do j=2,ny-1
      do i=2,nx-1
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          ! deposition
          if (hhalf(ij).ge.(1.d0+1.d-6)*ht(ij)) then
            Dp=(hhalf(ij)-ht(ij))/dt
            Ep=K1/2.d0*(hhalf(ipj)-hhalf(ij))/dx**2
            Mp=-K1/2.d0*(hhalf(ij)-hhalf(imj))/dx**2
            Np=K1/2.d0*(ft(ijp)+ft(ij))*(ht(ijp)-ht(ij))/dy**2 &
            -K1/2.d0*(ft(ij)+ft(ijm))*(ht(ij)-ht(ijm))/dy**2 &
            +Q1(ij)
            diag(i)=2.d0*L/dt+Dp-Mp-Ep
            sup(i)=-Ep
            inf(i)=-Mp
            rhs(i)=Np-Dp*ft(ij)+2.d0*L*ft(ij)/dt
            ! erosion
          else
            diag(i)=1.d0
            sup(i)=0.d0
            inf(i)=0.d0
            rhs(i)=ft(ij)
          endif
          ! in continent
        else
          diag(i)=1.d0
          sup(i)=0.d0
          inf(i)=0.d0
          rhs(i)=ft(ij)
        endif
      enddo
      ! bc on i=1
      diag(1)=1.d0
      sup(1)=-1.d0
      rhs(1)=0.d0
      ! bc on i=nx
      diag(nx)=1.d0
      inf(nx)=-1.d0
      rhs(nx)=0.d0
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,nx)
      do i=1,nx
        ij=(j-1)*nx+i
        fhalf(ij)=res(i)
      enddo
    enddo
    fhalf=max(0.d0,fhalf)
    fhalf=min(1.d0,fhalf)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the elevation h in y-direction
    allocate (diag(ny),sup(ny),inf(ny),rhs(ny),res(ny))
    do i=2,nx-1
      do j=1,ny
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          if (j.eq.1) then
            if (cbc(1:1).eq.'1') then
              diag(j)=1.d0
              sup(j)=0.d0
              rhs(j)=hhalf(ij)
            else
              Ap=dt/2.d0*(K2+(K1-K2)*(fp(ijp)+fp(ij))/2.d0)/dy**2
              diag(j)=1.d0+Ap
              sup(j)=-Ap
              Cp=dt/2.d0*(K2+(K1-K2)*(fhalf(ipj)+fhalf(ij))/2.d0)*(hhalf(ipj)-hhalf(ij))/dx**2 &
              -dt/2.d0*(K2+(K1-K2)*(fhalf(ij)+fhalf(imj))/2.d0)*(hhalf(ij)-hhalf(imj))/dx**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(j)=Cp+hhalf(ij)
            endif
          elseif (j.eq.ny) then
            if (cbc(3:3).eq.'1') then
              diag(j)=1.d0
              inf(j)=0.d0
              rhs(j)=hhalf(ij)
            else
              Bp=-dt/2.d0*(K2+(K1-K2)*(fp(ij)+fp(ijm))/2.d0)/dy**2
              diag(j)=1.d0-Bp
              inf(j)=Bp
              Cp=dt/2.d0*(K2+(K1-K2)*(fhalf(ipj)+fhalf(ij))/2.d0)*(hhalf(ipj)-hhalf(ij))/dx**2 &
              -dt/2.d0*(K2+(K1-K2)*(fhalf(ij)+fhalf(imj))/2.d0)*(hhalf(ij)-hhalf(imj))/dx**2 &
              +(Q1(ij)+Q2(ij))*dt/2.d0
              rhs(j)=Cp+hhalf(ij)
            endif
          else
            Ap=dt/2.d0*(K2+(K1-K2)*(fp(ijp)+fp(ij))/2.d0)/dy**2
            Bp=-dt/2.d0*(K2+(K1-K2)*(fp(ij)+fp(ijm))/2.d0)/dy**2
            diag(j)=1.d0+Ap-Bp
            sup(j)=-Ap
            inf(j)=Bp
            Cp=dt/2.d0*(K2+(K1-K2)*(fhalf(ipj)+fhalf(ij))/2.d0)*(hhalf(ipj)-hhalf(ij))/dx**2 &
            -dt/2.d0*(K2+(K1-K2)*(fhalf(ij)+fhalf(imj))/2.d0)*(hhalf(ij)-hhalf(imj))/dx**2 &
            +(Q1(ij)+Q2(ij))*dt/2.d0
            rhs(j)=Cp+hhalf(ij)
          endif
          ! in continent
        else
          diag(j)=1.d0
          sup(j)=0.d0
          inf(j)=0.d0
          rhs(j)=hhalf(ij)
        endif
      enddo
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,ny)
      do j=1,ny
        ij=(j-1)*nx+i
        tint(ij)=res(j)
      enddo
    enddo
    h=tint
    ! the corner nodes (1,1) and (1,ny)
    h(1)=h(2)
    h((ny-1)*nx+1)=h((ny-1)*nx+2)
    ! the corner nodes (nx,1) and (nx,ny)
    h(nx)=h(nx-1)
    h(nx*ny)=h(nx*ny-1)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the silt fraction F in y-direction
    allocate (diag(ny),sup(ny),inf(ny),rhs(ny),res(ny))
    do i=2,nx-1
      do j=2,ny-1
        ij=(j-1)*nx+i
        ipj=(j-1)*nx+i+1
        imj=(j-1)*nx+i-1
        ijp=(j)*nx+i
        ijm=(j-2)*nx+i
        ! in ocean and not at ocean-continent transition
        if (ht(ij).le.sealevel.and.flag(ij).eq.0) then
          ! deposition
          if (h(ij).ge.(1.d0+1.d-6)*hhalf(ij)) then
            Dp=(h(ij)-hhalf(ij))/dt
            Ep=K1/2.d0*(h(ijp)-h(ij))/dy**2
            Mp=-K1/2.d0*(h(ij)-h(ijm))/dy**2
            Np=K1/2.d0*(fhalf(ipj)+fhalf(ij))*(hhalf(ipj)-hhalf(ij))/dx**2 &
            -K1/2.d0*(fhalf(ij)+fhalf(imj))*(hhalf(ij)-hhalf(imj))/dx**2 &
            +Q1(ij)
            diag(j)=2.d0*L/dt+Dp-Mp-Ep
            sup(j)=-Ep
            inf(j)=-Mp
            rhs(j)=Np-Dp*fhalf(ij)+2.d0*L*fhalf(ij)/dt
            ! erosion
          else
            diag(j)=1.d0
            sup(j)=0.d0
            inf(j)=0.d0
            rhs(j)=fhalf(ij)
          endif
          ! in continent
        else
          diag(j)=1.d0
          sup(j)=0.d0
          inf(j)=0.d0
          rhs(j)=fhalf(ij)
        endif
      enddo
      ! bc on j=1
      diag(1)=1.d0
      sup(1)=-1.d0
      rhs(1)=0.d0
      ! bc on j=ny
      diag(ny)=1.d0
      inf(ny)=-1.d0
      rhs(ny)=0.d0
      ! solve a tri-diagonal system of equations
      call tridag (inf,diag,sup,rhs,res,ny)
      do j=1,ny
        ij=(j-1)*nx+i
        f(ij)=res(j)
      enddo
    enddo
    f=max(0.d0,f)
    f=min(1.d0,f)
    deallocate (diag,sup,inf,rhs,res)

    ! calculate the errors in each iteration
    err1=maxval(abs(h-hp))
    err2=maxval(abs(h-hp)/(1.d0+abs(h)))

    !print*,'niter',niter,minval(h-hp),sum(h-hp)/nn,maxval(h-hp),err1

    if (niter.gt.1000) then
      FSCAPE_RAISE_MESSAGE('Marine error: Multi-lithology diffusion not converging; decrease time step',ERR_NotConverged,ierr)
      FSCAPE_CHKERR(ierr)
    endif

    ! end of iteration
  enddo

  ! set the silt fraction for continent
  where (h.ge.sealevel+1.d-3) f=0.d-1

  deallocate (hp,fp,ht,ft,hhalf,fhalf,fhalfp)

  ! end of the subroutine
end subroutine SiltSandCouplingDiffusion

!-----------------------------------------------------

subroutine compaction (F1,F2,poro1,poro2,z1,z2,nn,dh,zi,zo)

  ! Newton iteration to calculate the thickness after compaction

  implicit none

  integer k,nn
  double precision poro1,poro2,z1,z2,fx,dfx
  double precision F1(nn),F2(nn),dh(nn),zi(nn),zo(nn)

  ! initial guess on zo
  zo=zi
  ! iteration process
  do k=1,nn
    1000 continue
    fx=zo(k)-zi(k)+F1(k)*poro1*z1*(exp(-zo(k)/z1)-exp(-zi(k)/z1)) &
    +F2(k)*poro2*z2*(exp(-zo(k)/z2)-exp(-zi(k)/z2))-dh(k)
    dfx=1.d0-F1(k)*poro1*exp(-zo(k)/z1)-F2(k)*poro2*exp(-zo(k)/z2)
    zo(k)=zo(k)-fx/dfx
    if (abs(fx/dfx).gt.1.d-6) goto 1000
  enddo

  return

end subroutine compaction
