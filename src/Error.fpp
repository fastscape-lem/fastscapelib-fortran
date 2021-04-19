! Error.fpp

! Macros for error handling.
! Enables user to store errors and exit the subroutine in single statement.
! Fortran preprocessor must be enabled: -fpp.

! Raise Error
! Store the error code and info (only if the current code is zero).
! Return from the subroutine.
#define RAISE_ERROR(ms,err,err_type)\
if (err%Code == ERR_None) then;\
err=ErrorType(Code=err_type,Message=ms);\
print *, "Exception: file ",__FILE__," line ",__LINE__;\
print *, "Exception: code ",err_type," ",err_names(err_type);\
print *, "===";\
write(*,*)'Error Raised: ',ms;end if;

! Pass Error
! Returns if there's an error.
#define HANDLE_ERROR(err)\
if (err%Code /= ERR_None) then;\
write(*,*)'Error Handled: ',err%Message;\
stop;\
end if;


#define FSCAPE_RAISE(err_type, ierr)\
ierr = (err_type);\
write (*, 666) "Exception: file", __FILE__, ", line", __LINE__;\
666 format (A,A,A,I4);\
write (*, 667) "Exception: code ", (err_type), " -> ", err_names((ierr));\
667 format (A,I4,A,A);\
write (*, 999) "==========";\
999 format (A);

#define FSCAPE_RAISE1(msg, err_type, ierr)\
ierr = (err_type);\
write (*, 666) "Exception: file", __FILE__, ", line", __LINE__;\
write (*, 667) "Exception: code ", (err_type), " -> ", err_names((ierr));\
write (*, 668) "Exception: ", (msg);\
668 format (A,A);\
write (*, 999) "==========";

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

