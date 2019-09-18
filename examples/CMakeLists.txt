cmake_minimum_required(VERSION 3.5)

set(FASTSCAPELIB_EXAMPLES_SRC_FILES
  DippingDyke.f90
  Fan.f90
  Margin.f90
  Mountain.f90
  Strati_test.f90
  flexure_test.f90
  )

foreach(EXAMPLE_SRC_FILE ${FASTSCAPELIB_EXAMPLES_SRC_FILES})
  string(REGEX REPLACE "\\.[^.]*$" "" EXAMPLE_EXE_FILE ${EXAMPLE_SRC_FILE})

  add_executable(${EXAMPLE_EXE_FILE} ${EXAMPLE_SRC_FILE})

  target_link_libraries(${EXAMPLE_EXE_FILE} ${FASTSCAPELIB_STATIC})
endforeach(EXAMPLE_SRC_FILE)
