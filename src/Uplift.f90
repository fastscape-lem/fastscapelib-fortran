subroutine Uplift ()

  ! subroutine to apply an uplift step using the uplift function/array u

  use FastScapeContext

  implicit none

  h = h + u*dt
  b = b + u*dt

  return

  end subroutine Uplift
