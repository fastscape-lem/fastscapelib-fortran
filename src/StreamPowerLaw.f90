!--------------------------------------------------------------------------------------------

subroutine StreamPowerLaw ()

  ! subroutine to solve the stream power law equation following the FastScape method described
  ! in Braun and Willett, Geomorphology, 2015

  use FastScapeContext

  implicit none

  integer :: ij,ijk,ijr,k,ijr1
  double precision :: dx,dy,fact,tol,err
  double precision :: f,df,errp,h0,hn,omega,tolp,w_rcv
  double precision, dimension(:), allocatable :: ht,kfint,dh,hp
  double precision, dimension(:), allocatable :: elev
  logical, dimension(:), allocatable :: bc
  double precision, dimension(:), allocatable :: water,lake_water_volume,lake_sediment
  integer, dimension(:), allocatable :: lake_sill

  allocate (ht(nn),kfint(nn),bc(nn),dh(nn),hp(nn))
  allocate (elev(nn))
  allocate (water(nn),lake_water_volume(nn),lake_sediment(nn),lake_sill(nn))

  dx=xl/(nx-1)
  dy=yl/(ny-1)

  ! set g, dimensionless parameter for sediment transport and deposition
  ! if g1<0, skip and use g values directly from FastScapeContext (not in API!!!)
  if (g1.ge.0.d0) then
    g=g1
    if (g2.gt.0.d0) where ((h-b).gt.1.d0) g=g2
  endif

  ! set kf / kfsed
  kfint=kf
  if (kfsed.gt.0.d0) where ((h-b).gt.1.d0) kfint=kfsed

  lake_depth = hwater - h

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
  tol=1.d-4*(maxval(abs(h)) + 1.d0)
  err=2.d0*tol

  ! store the elevation at t
  ht=h

  ! Gauss-Seidel iteration
  nGSStreamPowerLaw=0

  lake_sediment=0.d0
  lake_sill=0.d0

  do while (err.gt.tol.and.nGSStreamPowerLaw.lt.99)
    nGSStreamPowerLaw=nGSStreamPowerLaw+1
    ! guess/update the elevation at t+Δt (k)
    hp=h

    ! calculate erosion/deposition at each node
    dh=ht-hp

    ! sum the erosion in stack order
    do ij=1,nn
      ijk=mstack(ij)
      ijr1=rec(ijk)
      if (ijr1.ne.ijk) then
        dh(ijk)=dh(ijk)-(ht(ijk)-hp(ijk))
        if (lake_sill(ijk).eq.ijk) then
          if (dh(ijk).le.0.d0) then
            lake_sediment(ijk)=0.d0
          else
            lake_sediment(ijk)=dh(ijk)
          endif
        endif
        dh(ijk)=dh(ijk)+(ht(ijk)-hp(ijk))
        do k=1,mnrec(ijk)
          ijr=mrec(k,ijk)
          dh(ijr)=dh(ijr)+dh(ijk)*mwrec(k,ijk)
        enddo
      else
        lake_sediment(ijk)=dh(ijk)
      endif
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
          ijr1=rec(ijk)
          if (ijr1.eq.ijk) then
            water(ijk)=ht(ijk)
            lake_sill(ijk)=ijk
            lake_water_volume(ijk)=0.d0
          else
            w_rcv=water(ijr1)
            if (elev(ijk).gt.w_rcv) then
              if (mnrec(ijk).gt.0) then
                if (h(ijk).ge.sealevel.or..not.runMarine) then
                  f = elev(ijk)
                  df = 1.d0
                  do k=1,mnrec(ijk)
                    if (ht(ijk).ge.ht(mrec(k,ijk))) then
                      fact = kfint(ijk)*dt*(a(ijk)*mwrec(k,ijk))**m/mlrec(k,ijk)
                      f = f + fact*h(mrec(k,ijk))
                      df = df + fact
                    endif
                  enddo
                  h(ijk)=f/df
                endif
              endif
              lake_sill(ijk)=ijk
              lake_water_volume(ijk)=0.d0
              if (h(ijk).lt.w_rcv) h(ijk)=w_rcv
            else
              h(ijk)=elev(ijk)
              lake_sill(ijk)=lake_sill(ijr1)
              if (lake_sill(ijk).ne.0) lake_water_volume(lake_sill(ijk)) = &
              lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))
            endif
            water(ijk)=max(w_rcv,h(ijk))
          endif
        enddo

      else

        do ij=nn,1,-1
          ijk=mstack(ij)
          ijr1=rec(ijk)
          if (ijr1.eq.ijk) then
            water(ijk)=ht(ijk)
            lake_sill(ijk)=ijk
            lake_water_volume(ijk)=0.d0
          else
            w_rcv=water(ijr1)
            if (elev(ijk).gt.w_rcv) then
              if (mnrec(ijk).gt.0) then
                if (ht(ijk).ge.sealevel.or..not.runMarine) then
                  omega=0.875d0/n
                  tolp=1.d-3
                  errp=2.d0*tolp
                  h0=elev(ijk)
                  do while (errp.gt.tolp)
                    f=h(ijk)-h0
                    df=1.d0
                    do k=1,mnrec(ijk)
                      if (ht(ijk).gt.ht(mrec(k,ijk))) then
                        fact = kfint(ijk)*dt*(a(ijk)*mwrec(k,ijk))**m/mlrec(k,ijk)**n
                        f=f+fact*max(0.d0,h(ijk)-h(mrec(k,ijk)))**n
                        df=df+fact*n*max(0.d0,h(ijk)-h(mrec(k,ijk)))**(n-1.d0)
                      endif
                    enddo
                    hn=h(ijk)-f/df
                    errp=abs(hn-h(ijk))
                    h(ijk)=h(ijk)*(1.d0-omega)+hn*omega
                  enddo
                endif
              endif
              lake_sill(ijk)=ijk
              lake_water_volume(ijk)=0.d0
              if (h(ijk).lt.w_rcv) h(ijk)=w_rcv
            else
              h(ijk)=elev(ijk)
              lake_sill(ijk)=lake_sill(ijr1)
              if (lake_sill(ijk).ne.0) lake_water_volume(lake_sill(ijk)) = &
              lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))
            endif
            water(ijk)=max(w_rcv,h(ijk))
          endif
        enddo

      endif

      err=sqrt(sum((h-hp)**2)/nn)
      if (maxval(g).lt.tiny(g)) err=0.d0
      if (nGSStreamPowerLaw.eq.99) print*,'Beware Gauss-Siedel scheme not convergent'

    enddo

    b=min(h,b)

    do ij=1,nn
      if (lake_sill(ij).ne.0) then
        if (lake_water_volume(lake_sill(ij)).gt.0.d0) h(ij)=h(ij) &
        +max(0.d0,min(lake_sediment(lake_sill(ij)),lake_water_volume(lake_sill(ij))))/ &
        lake_water_volume(lake_sill(ij))*(water(ij)-h(ij))
      endif
    enddo

    ! stores total erosion, erosion rate and flux for output
    etot=etot+ht-h
    erate=(ht-h)/dt
    Sedflux=ht-h
    !if (runMarine) where (h.lt.sealevel) Sedflux=0.d0

    return

  end subroutine StreamPowerLaw

  !--------------------------------------------------------------------------------------------

  subroutine StreamPowerLawSingleFlowDirection ()

    ! subroutine to solve the stream power law equation following the FastScape method described
    ! in Braun and Willett, Geomorphology, 2015

    use FastScapeContext

    implicit none

    integer :: ij,ijk,ijr
    double precision :: dx,dy,fact,tol,err
    double precision :: f,df,errp,h0,hn,omega,tolp,w_rcv
    double precision, dimension(:), allocatable :: ht,kfint,dh,hp
    double precision, dimension(:), allocatable :: elev
    logical, dimension(:), allocatable :: bc
    double precision, dimension(:), allocatable :: water,lake_water_volume,lake_sediment
    integer, dimension(:), allocatable :: lake_sill

    allocate (ht(nn),kfint(nn),bc(nn),dh(nn),hp(nn))
    allocate (elev(nn))
    allocate (water(nn),lake_water_volume(nn),lake_sediment(nn),lake_sill(nn))

    ! sets boundary conditions

    dx=xl/(nx-1)
    dy=yl/(ny-1)

    ! set g, dimensionless parameter for sediment transport and deposition
    ! if g1<0, skip and use g values directly from FastScapeContext (not in API!!!)
    if (g1.ge.0.d0) then
      g=g1
      if (g2.gt.0.d0) where ((h-b).gt.1.d0) g=g2
    endif

    ! set kf / kfsed
    kfint=kf
    if (kfsed.gt.0.d0) where ((h-b).gt.1.d0) kfint=kfsed

    lake_depth = hwater - h

    ! computes drainage area
    a=dx*dy*precip
    do ij=nn,1,-1
      ijk=stack(ij)
      a(rec(ijk))=a(rec(ijk))+a(ijk)
    enddo

    ! calculate the elevation / SPL, including sediment flux
    tol=1.d-4*(maxval(abs(h))+1.d0)
    err=2.d0*tol

    ! store the elevation at t
    ht=h

    ! Gauss-Seidel iteration
    nGSStreamPowerLaw=0

    lake_sediment=0.d0

    do while (err.gt.tol.and.nGSStreamPowerLaw.lt.99)
      nGSStreamPowerLaw=nGSStreamPowerLaw+1
      ! guess/update the elevation at t+Δt (k)
      hp=h

       ! calculate erosion/deposition at each node
      dh=ht-hp

      ! sum the erosion in stack order
      do ij=nn,1,-1
        ijk=stack(ij)
        ijr=rec(ijk)
        if (ijr.ne.ijk) then
          dh(ijk)=dh(ijk)-(ht(ijk)-hp(ijk))
          if (lake_sill(ijk).eq.ijk) then
            if (dh(ijk).le.0.d0) then
              lake_sediment(ijk)=0.d0
            else
              lake_sediment(ijk)=dh(ijk)
            endif
          endif
          dh(ijk)=dh(ijk)+(ht(ijk)-hp(ijk))
          dh(ijr)=dh(ijr)+dh(ijk)
        else
          lake_sediment(ijk)=dh(ijk)
        endif
      enddo

      where (bc)
      elev=ht
      elsewhere
        elev=ht+(dh-(ht-hp))*g*dx*dy/a
        endwhere

        ! apply modified stream power law using lake surface (hwater)

        if (abs(n-1.d0).lt.tiny(n)) then

          do ij=1,nn
            ijk=stack(ij)
            ijr=rec(ijk)
            if (ijr.eq.ijk) then
              water(ijk)=ht(ijk)
              lake_sill(ijk)=ijk
              lake_water_volume(ijk)=0.d0
            else
              w_rcv=water(ijr)
              if (elev(ijk).gt.w_rcv) then
                if (h(ijk).ge.sealevel.or..not.runMarine) then
                f = elev(ijk)
                df = 1.d0
                ! todo: check if we don't need those checks for single flow
!                if (ht(ijk).ge.ht(ijr)) then
                  fact = kfint(ijk)*dt*a(ijk)**m/length(ijk)
                  f = f + fact*h(ijr)
                  df = df + fact
!                endif
                h(ijk)=f/df
!                h(ijk)=min(f/df,minval(h(don(1:ndon(ijk),ijk))))
                  endif
                lake_sill(ijk)=ijk
                lake_water_volume(ijk)=0.d0
                if (h(ijk).lt.w_rcv) h(ijk)=w_rcv
              else
                h(ijk)=elev(ijk)
                lake_sill(ijk)=lake_sill(ijr)
                if (lake_sill(ijk).ne.0) lake_water_volume(lake_sill(ijk)) = &
                lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))
              endif
              water(ijk)=max(w_rcv,h(ijk))
            endif
          enddo

        else

          do ij=1,nn
            ijk=stack(ij)
            ijr=rec(ijk)
            if (ijr.eq.ijk) then
              water(ijk)=ht(ijk)
              lake_sill(ijk)=ijk
              lake_water_volume(ijk)=0.d0
            else
              w_rcv=water(ijr)
              if (elev(ijk).gt.w_rcv) then
                if (ht(ijk).ge.sealevel.or..not.runMarine) then
                omega=0.875d0/n
                tolp=1.d-3
                errp=2.d0*tolp
                h0=elev(ijk)
                do while (errp.gt.tolp)
                  f=h(ijk)-h0
                  df=1.d0
                  if (ht(ijk).gt.ht(ijr)) then
                    fact = kfint(ijk)*dt*a(ijk)**m/length(ijk)**n
                    f=f+fact*max(0.d0,h(ijk)-h(ijr))**n
                    df=df+fact*n*max(0.d0,h(ijk)-h(ijr))**(n-1.d0)
                  endif
                  hn=h(ijk)-f/df
                  errp=abs(hn-h(ijk))
                  h(ijk)=h(ijk)*(1.d0-omega)+hn*omega
                enddo
                endif
                lake_sill(ijk)=ijk
                lake_water_volume(ijk)=0.d0
                if (h(ijk).lt.w_rcv) h(ijk)=w_rcv
              else
                h(ijk)=elev(ijk)
                lake_sill(ijk)=lake_sill(ijr)
                if (lake_sill(ijk).ne.0) lake_water_volume(lake_sill(ijk)) = &
                lake_water_volume(lake_sill(ijk))+(w_rcv-h(ijk))
              endif
              water(ijk)=max(w_rcv,h(ijk))
            endif
          enddo

        endif

        err=sqrt(sum((h-hp)**2)/nn)
        if (maxval(g).lt.tiny(g)) err=0.d0
        if (nGSStreamPowerLaw.eq.99) print*,'Beware Gauss-Siedel scheme not convergent'
      enddo

      b=min(h,b)

      do ij=1,nn
        if (lake_sill(ij).ne.0) then
          if (lake_water_volume(lake_sill(ij)).gt.0.d0) h(ij)=h(ij) &
          +max(0.d0,min(lake_sediment(lake_sill(ij)),lake_water_volume(lake_sill(ij))))/ &
          lake_water_volume(lake_sill(ij))*(water(ij)-h(ij))
        endif
      enddo

      deallocate (hwater)

      ! stores total erosion, erosion rate and flux for output
      etot=etot+ht-h
      erate=(ht-h)/dt
      Sedflux=ht-h
      !if (runMarine) where (h.lt.sealevel) Sedflux=0.d0

      return

    end subroutine StreamPowerLawSingleFlowDirection

    !--------------------------------------------------------------------------------------------
