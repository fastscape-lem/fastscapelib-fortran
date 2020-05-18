module FastScapeTypes

  ! Defines derived types used in FastScape

  implicit none

  type, public :: boundaries
    integer :: ibc
    integer :: i1, i2, j1, j2
    logical :: xcyclic, ycyclic
    logical, dimension(:), allocatable :: bc
  end type boundaries

end module FastScapeTypes
