
! Error.fpp

!
! Macros for error handling.
! Fortran preprocessor must be enabled: -fpp.
!

!
! Raise an exception by pushing the name of the exception + file name, line number to stdout
! Sets ierr to the error type specified by err_type
!
#define FSCAPE_RAISE(err_type, ierr) \
  ierr = (err_type);                 \
  print '(A,I4,A)', "*** FastScape exception *** -> " // TRIM(err_names((ierr))) // " (" // __FILE__ // ", line", __LINE__, ")";

!
! Raise an exception by pushing custom user provided message to stdout,
! along with the name of the exception + file name, line number to stdout
! Sets ierr to the error type specified by err_type
!
#define FSCAPE_RAISE_MESSAGE(msg, err_type, ierr)                    \
  print '(A,A)', "*** FastScape exception *** ", TRIM((msg)) // "!"; \
  FSCAPE_RAISE(err_type, ierr);

!
! Calls return if error code is non-zero.
!
#define FSCAPE_CHKERR(ierr) \
  if (ierr /= 0) then;      \
    return;                 \
  end if;

!
! Force stop if error code is non-zero.
! This should only be called in driver code - never within library code or the API.
!
#define FSCAPE_CHKERR_ABORT(ierr) \
  if (ierr /= 0) then;            \
    stop;                         \
  end if;
