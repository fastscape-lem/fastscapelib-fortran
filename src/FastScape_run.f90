#include "Error.fpp"
program FastScapeRUN

implicit none

integer :: nx,ny,istep,nstep,nn,nfreq
double precision, dimension(:), allocatable :: u,ux,uy,h,b,etot,erate,a,chi,catchment,sedflux,sedflux_shore,kf1,kd1
double precision, dimension(:,:), allocatable :: field
real :: time_in,time_out
double precision :: m,n,g1,g2,preci_rate,p_flow_dir_exp,kf2,kd2
double precision xl,yl,dx,dy,dt
double precision sealevel, poro1, poro2, z1, z2, ratio, L, kds1, kds2

integer i,j,ij
integer ierr

nx=201
ny=201
nn=nx*ny

call FastScape_Init(ierr);FSCAPE_CHKERR_ABORT(ierr)
call FastScape_Set_NX_NY (nx,ny,ierr);FSCAPE_CHKERR_ABORT(ierr)
call FastScape_Setup(ierr);FSCAPE_CHKERR_ABORT(ierr)

xl=200.d3
yl=200.d3
dx=xl/(nx-1)
dy=yl/(ny-1)
call FastScape_Set_XL_YL (xl,yl,ierr);FSCAPE_CHKERR_ABORT(ierr)
dt=1.d3
call FastScape_Set_DT (dt,ierr);FSCAPE_CHKERR_ABORT(ierr)
allocate (kf1(nn),kd1(nn))
kf1=1.d-5
kf2=2.d-5
!kf2=kf1
m=0.4d0
n=1.d0
p_flow_dir_exp = -2.d0
kd1=5.d-2
kd2=0.d-2
kd1 = 1.d0
kd2 = 0.d0
!kd2=kd1
g1=1.d0
g2=1.d0
preci_rate = 1.d0 ! precipitation rate
call FastScape_Set_Erosional_Parameters (kf1,kf2,m,n,kd1,kd2,g1,g2,p_flow_dir_exp,ierr);FSCAPE_CHKERR_ABORT(ierr)
!call FastScape_Set_Precipitation_Rate (preci_rate)
sealevel = 0.d0
!poro1 = 0.63d0
!poro2 = 0.49d0
poro1 = 0.0d0
poro2 = 0.0d0
ratio = 0.5d0
L = 0.5d2
kds1 = 2.d2
kds2 = 1.d2
z1 = 1.d3
z2 = 1.d3
call FastScape_Set_Marine_Parameters (sealevel, poro1, poro2, z1, z2, ratio, L, kds1, kds2,ierr);FSCAPE_CHKERR_ABORT(ierr)

call FastScape_Set_BC (1010,ierr);FSCAPE_CHKERR_ABORT(ierr)

allocate (h(nn),b(nn),u(nn),ux(nn),uy(nn),etot(nn),erate(nn),a(nn),chi(nn),catchment(nn),sedflux(nn),sedflux_shore(nn))
allocate (field(nn,2))
call random_number (h)
  do j=1,ny
    do i=1,nx
    ij=(j-1)*nx+i
      if (j.lt.ny/2) then
      h(ij)=h(ij)-200.d0
      elseif (j.gt.ny/2) then
      h(ij)=h(ij)+1000.d0
      endif
    enddo
  enddo

call FastScape_Init_H (h,ierr);FSCAPE_CHKERR_ABORT(ierr)

do j=1,ny
    do i=1,nx
        ij=(j-1)*nx+i
        u(ij)=5.d-4
        if (j.lt.ny/4) then
            u(ij)=-5.d-4*float(j-1)/(ny/4-1)
        elseif (j.gt.3*ny/4) then
            u(ij)=0.d0
        endif
    enddo
enddo
u = 0.d0
call FastScape_Set_U (u,ierr);FSCAPE_CHKERR_ABORT(ierr)

nstep=500
nfreq=100 ! frequency of output
call FastScape_View(ierr);FSCAPE_CHKERR_ABORT(ierr)
istep=nstep

call cpu_time (time_in)
  do while (istep.le.nstep)
  call FastScape_Execute_Step (ierr);FSCAPE_CHKERR_ABORT(ierr)
  call FastScape_Get_Step (istep,ierr);FSCAPE_CHKERR_ABORT(ierr)
    if (mod(istep,nfreq)==0) then
    call FastScape_Copy_H (h,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Basement (b,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Total_Erosion (etot,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Erosion_Rate (erate,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Drainage_Area (a,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Chi (chi,ierr);FSCAPE_CHKERR_ABORT(ierr)
    call FastScape_Copy_Catchment (catchment,ierr);FSCAPE_CHKERR_ABORT(ierr)
    !call FastScape_Copy_Sediment_Flux(sedflux)
    !call FastScape_Copy_Sediment_Flux_Shore(sedflux_shore)
    print*,istep
    print*,'topo',minval(h),sum(h)/nn,maxval(h)
    print*,'basement',minval(b),sum(b)/nn,maxval(b)
    print*,'etot',minval(etot),sum(etot)/nn,maxval(etot)
    print*,'erate',minval(erate),sum(erate)/nn,maxval(erate)
    print*,'a',minval(a),sum(a)/nn,maxval(a)
    print*,'chi',minval(chi),sum(chi)/nn,maxval(chi)
    print*,'catchment',minval(catchment),sum(catchment)/nn,maxval(catchment)
    call FastScape_Debug(ierr);FSCAPE_CHKERR_ABORT(ierr)
    field(:,1)=etot
    field(:,2)=erate
    call FastScape_VTK (chi, 2.d0,ierr);FSCAPE_CHKERR_ABORT(ierr) ! if value > 0, elevation plus value is written to file. If value <0, basement and sealevel + value is written to file.
    call FastScape_VTK (erate, -2.d0,ierr);FSCAPE_CHKERR_ABORT(ierr) ! if value > 0, elevation plus value is written to file. If value <0, basement and sealevel + value is written to file.
    !field(:,3)=sedflux
    !field(:,4)=sedflux_shore
    !call VTK (h,b,2,field,nx,ny,dx,dy,istep)
    endif
  enddo
call cpu_time (time_out)
print*,'Total run time',time_out-time_in

call FastScape_Destroy (ierr);FSCAPE_CHKERR_ABORT(ierr)

deallocate (u,h,b,a,etot,erate,chi,catchment,sedflux,sedflux_shore,field)
deallocate (kf1,kd1)

end program FastScapeRun
