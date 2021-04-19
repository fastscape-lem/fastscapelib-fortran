
module FastScapeErrorCodes
  implicit none
  integer, parameter :: ERR_None                = 0,  &
                        ERR_Default             = 1,  &
                        ERR_FileNotFound        = 2,  &
                        ERR_FileOpenFailed      = 3,  &
                        ERR_ParameterInvalid    = 4,  &
                        ERR_ParameterOutOfRange = 5

  character(len=50), dimension(5) :: err_names = [character(len=50) :: "Default run time error", &
      "File error: File not found",        &
      "File error: File open failed",      &
      "Parameter error: Input invalid",    &
      "Parameter error: Out of range"      ]

end module FastScapeErrorCodes


