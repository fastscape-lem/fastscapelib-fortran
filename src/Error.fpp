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
write(*,*)'Error Raised: ',ms;end if;

! Pass Error
! Returns if there's an error.
#define HANDLE_ERROR(err)\
if (err%Code /= ERR_None) then;\
write(*,*)'Error Handled: ',err%Message;\
stop;\
end if;
