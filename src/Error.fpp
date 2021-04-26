
! Error.fpp

!
! Macros for error handling.
! Fortran preprocessor must be enabled: -fpp.
!

!
! Initialize ierr argument in API functions
!
#define FSCAPE_INITERR(ierr, ierr_, fname)        \
  ierr_ = 0;                                      \
  if (present(ierr)) then;                        \
    ierr = ierr_;                                 \
  else;                                           \
    print '(A)', "*** Depreciation warning *** "; \
    print '(A,A,A)', "Calling ", TRIM((fname)), " without 'ierr' argument (integer) is depreciated! Please update your code!"; \
  end if;

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
! Either calls return or stop if error code is non-zero, dependening on whether the
! the ierr argument is present.
!
#define FSCAPE_CHKERR_OPT(ierr, ierr_)                                         \
  if (ierr_ /= 0) then;                                                        \
    if (present(ierr)) then;                                                   \
      ierr = ierr_;                                                            \
      return;                                                                  \
    else;                                                                      \
      print '(A)', "*** Execution being halted by calling stop!";               \
      print '(A)', "Call all API routine with 'ierr' argument to avoid this."; \
      stop;                                                                    \
    end if;                                                                    \
  end if;

!
! Force stop if error code is non-zero.
! This should only be called in driver code - never within library code or the API.
!
#define FSCAPE_CHKERR_ABORT(ierr) \
  if (ierr /= 0) then;            \
    stop;                         \
  end if;
