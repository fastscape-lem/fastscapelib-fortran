
module FastScapeErrorCodes
  implicit none
  integer, parameter :: ERR_None                = 0,  &
                        ERR_Default             = 1,  &
                        ERR_FileNotFound        = 2,  &
                        ERR_FileOpenFailed      = 3,  &
                        ERR_ParameterInvalid    = 4,  &
                        ERR_ParameterOutOfRange = 5,  &
                        ERR_TridiagNotSolvable  = 6,  &
                        ERR_SetupNotRun         = 7,  &
                        ERR_NotConverged        = 8

  character(len=50), dimension(8) :: err_names = [character(len=50) :: "Default run time error", &
      "File error: File not found",             &
      "File error: File open failed",           &
      "Parameter error: Input invalid",         &
      "Parameter error: Out of range" ,         &
      "Solver error: Tridiag not solvable",     &
      "Order invalid: Must call SetUp() first", &
      "Solver error: not converged"  ]

end module FastScapeErrorCodes
