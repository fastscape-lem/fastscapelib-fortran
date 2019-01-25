cmake_minimum_required(VERSION 3.5)
project(fastscapelib-fortran)

enable_language(Fortran)

# Sources
# =======

set(FASTSCAPELIB_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(FASTSCAPELIB_SRC_FILES
  ${FASTSCAPELIB_SRC_DIR}/Advect.f90
  ${FASTSCAPELIB_SRC_DIR}/Diffusion.f90
  ${FASTSCAPELIB_SRC_DIR}/FastScape_api.f90
  ${FASTSCAPELIB_SRC_DIR}/FastScape_ctx.f90
  ${FASTSCAPELIB_SRC_DIR}/LocalMinima.f90
  ${FASTSCAPELIB_SRC_DIR}/Strati.f90
  ${FASTSCAPELIB_SRC_DIR}/StreamPowerLaw.f90
  ${FASTSCAPELIB_SRC_DIR}/VTK.f90
  )

# Flags
# =====

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

# set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

# Build options
# =============

OPTION(BUILD_LIBFASTSCAPE_STATIC "build fastscapelib as a static library" ON)
OPTION(BUILD_LIBFASTSCAPE_SHARED "build fastscapelib as a shared library" OFF)
OPTION(BUILD_PYTHON_MODULE "build fastscapelib library as a Python module" OFF)
OPTION(USE_FLEXURE "include flexure routines in the library" OFF)
OPTION(BUILD_EXAMPLES "build usage examples that are in the 'examples' directory" OFF)

if(BUILD_EXAMPLES)
  set(BUILD_LIBFASTSCAPE_STATIC ON)
  set(USE_FLEXURE ON)
endif()

if(USE_FLEXURE)
  set(FLEXURE_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Flexure2D_v1.0/src)
  set(FLEXURE_SRC_FILES
    ${FLEXURE_SRC_DIR}/flexure2D.f90
    ${FLEXURE_SRC_DIR}/four1.f
    ${FLEXURE_SRC_DIR}/realft.f
    ${FLEXURE_SRC_DIR}/sinft.f
    )
  list(APPEND FASTSCAPELIB_SRC_FILES ${FLEXURE_SRC_FILES})
endif()

# Fortran library (static/shared)
# ===============================

set(LIBFASTSCAPE_OBJECTS libfastscape_objects)
add_library(${LIBFASTSCAPE_OBJECTS} OBJECT ${FASTSCAPELIB_SRC_FILES})

set(LIBFASTSCAPE_NAME fastscape)

if(BUILD_LIBFASTSCAPE_STATIC)
  set(LIBFASTSCAPE_STATIC libfastscape_static)

  add_library(${LIBFASTSCAPE_STATIC} STATIC $<TARGET_OBJECTS:${LIBFASTSCAPE_OBJECTS}>)

  set_target_properties(${LIBFASTSCAPE_STATIC}
    PROPERTIES
    OUTPUT_NAME ${LIBFASTSCAPE_NAME})

  install(TARGETS ${LIBFASTSCAPE_STATIC} DESTINATION lib)
endif()

if(BUILD_LIBFASTSCAPE_SHARED)
  set(LIBFASTSCAPE_SHARED libfastscape_shared)

  # shared libraries need PIC
  set_target_properties(${LIBFASTSCAPE_OBJECTS}
    PROPERTIES
    POSITION_INDEPENDENT_CODE 1)

  add_library(${LIBFASTSCAPE_SHARED} SHARED $<TARGET_OBJECTS:${LIBFASTSCAPE_OBJECTS}>)

  set_target_properties(${LIBFASTSCAPE_SHARED}
    PROPERTIES
    OUTPUT_NAME ${LIBFASTSCAPE_NAME})

  install(TARGETS ${LIBFASTSCAPE_SHARED} DESTINATION lib)
endif()

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(BUILD_PYTHON_MODULE)
  if(NOT SKBUILD)
    set(CMAKE_MODULE_PATH
      ${PROJECT_SOURCE_DIR}/cmake
      ${CMAKE_MODULE_PATH})
  endif()

  find_package(F2PY REQUIRED)

  set(F2PY_MODULE_NAME "_fastscapelib_fortran")
  set(F2PY_MODULE_FILE
    ${CMAKE_CURRENT_BINARY_DIR}/${F2PY_MODULE_NAME}${PYTHON_EXTENSION_MODULE_SUFFIX}
    )

  add_custom_target(${F2PY_MODULE_NAME} ALL
    DEPENDS ${F2PY_MODULE_FILE}
    )

  add_custom_command(
    OUTPUT ${F2PY_MODULE_FILE}
    COMMAND ${F2PY_EXECUTABLE}
      -m ${F2PY_MODULE_NAME}
      -c
      ${FASTSCAPELIB_SRC_FILES}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

  if(SKBUILD)
    install(FILES ${F2PY_MODULE_FILE} DESTINATION fastscapelib_fortran)
  endif()
endif()
