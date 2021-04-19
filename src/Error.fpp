
! Error.fpp

! Macros for error handling.
! Enables user to store errors and exit the subroutine in single statement.
! Fortran preprocessor must be enabled: -fpp.

#define FSCAPE_RAISE(err_type, ierr)\
ierr = (err_type);\
print '(A,A,A,I4)', "Exception: file ", __FILE__, ", line", __LINE__;\
print '(A,I4,A,A)', "Exception: code ", (err_type), " -> ", err_names((ierr));\
print '(A)', "==========";\

#define FSCAPE_RAISE1(msg, err_type, ierr)\
ierr = (err_type);\
print '(A,A,A,I4)', "Exception: file ", __FILE__, ", line", __LINE__;\
print '(A,I4,A,A)', "Exception: code ", (err_type), " -> ", err_names((ierr));\
print '(A,A)', "Exception: ", (msg);\
print '(A)', "==========";

!
! Return earlier if error code is non-zero.
!
#define FSCAPE_CHKERR(ierr)\
if (ierr /= 0) then;\
return;\
end if;

!
! Force stop if error code is non-zero.
! This should only be called in driver code - never within library code or the API.
!
#define FSCAPE_CHKERR_ABORT(ierr)\
if (ierr /= 0) then;\
stop;\
end if;
