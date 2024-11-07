cmake_minimum_required(VERSION 3.5)
project(fastscapelib-fortran C Fortran)

set(FASTSCAPELIB_NAME fastscapelib_fortran)

# Sources
# =======

set(FASTSCAPELIB_SRC_DIR src)

set(FASTSCAPELIB_SRC_FILES
  ${FASTSCAPELIB_SRC_DIR}/Advect.f90
  ${FASTSCAPELIB_SRC_DIR}/Diffusion.f90
  ${FASTSCAPELIB_SRC_DIR}/FastScape_api.f90
  ${FASTSCAPELIB_SRC_DIR}/FastScape_ctx.f90
  ${FASTSCAPELIB_SRC_DIR}/LocalMinima.f90
  ${FASTSCAPELIB_SRC_DIR}/Marine.f90
  ${FASTSCAPELIB_SRC_DIR}/Strati.f90
  ${FASTSCAPELIB_SRC_DIR}/StreamPowerLaw.f90
  ${FASTSCAPELIB_SRC_DIR}/FlowRouting.f90
  ${FASTSCAPELIB_SRC_DIR}/Uplift.f90
  ${FASTSCAPELIB_SRC_DIR}/VTK.f90
  ${FASTSCAPELIB_SRC_DIR}/TerrainDerivatives.f90
  )

set(FLEXURE_SRC_DIR Flexure2D_v1.0/src)

set(FLEXURE_SRC_FILES
    ${FLEXURE_SRC_DIR}/flexure2D.f90
    ${FLEXURE_SRC_DIR}/four1.f
    ${FLEXURE_SRC_DIR}/realft.f
    ${FLEXURE_SRC_DIR}/sinft.f
    )

# Pre-processor
# =============

if(WIN32)
  add_compile_definitions(ON_WINDOWS)
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

# Flags
# =====

if(CMAKE_Fortran_COMPILER_ID MATCHES "Flang|GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none -fall-intrinsics")
    set(ioflags "-fconvert=big-endian")
    set(f77flags "-std=legacy -ffixed-form")
    set(bounds "-fbounds-check")
    set(warnings "-Wall")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(ioflags "-convert big_endian")
    set(f77flags "-stand none -fixed")
    set(bounds "-check bounds")
    set(warnings "-warn all")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(ioflags "")
    set(f77flags "-Mfixed")
    set(bounds "-C")
    set(warnings "-Wall")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG
  "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds} ${warnings}")

# let F2PY (scikit-build) configure flags for F77 and F90 source
if(NOT SKBUILD)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect} ${ioflags}")
else()
  if(APPLE AND MAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
    # https://github.com/numpy/numpy/issues/25869 ?
    # (seems problematic only on MacOS Mx)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-error=incompatible-function-pointer-types")
  endif()
endif()

# override flags for old-school F77 files
set_source_files_properties(${FLEXURE_SRC_DIR}/four1.f PROPERTIES COMPILE_FLAGS ${f77flags})
set_source_files_properties(${FLEXURE_SRC_DIR}/realft.f PROPERTIES COMPILE_FLAGS ${f77flags})
set_source_files_properties(${FLEXURE_SRC_DIR}/sinft.f PROPERTIES COMPILE_FLAGS ${f77flags})

# Build options
# =============

OPTION(BUILD_FASTSCAPELIB_STATIC "build fastscapelib as a static library" ON)
OPTION(BUILD_FASTSCAPELIB_SHARED "build fastscapelib as a shared library" OFF)
OPTION(USE_FLEXURE "include flexure routines in the library" OFF)
OPTION(BUILD_EXAMPLES "build usage examples that are in the 'examples' directory" OFF)

if(BUILD_EXAMPLES)
  set(BUILD_FASTSCAPELIB_STATIC ON)
  set(USE_FLEXURE ON)
endif()

if(USE_FLEXURE)
  list(APPEND FASTSCAPELIB_SRC_FILES ${FLEXURE_SRC_FILES})
endif()

# Fortran library (static/shared)
# ===============================

if(BUILD_FASTSCAPELIB_STATIC OR BUILD_FASTSCAPELIB_SHARED)
  set(FASTSCAPELIB_OBJECTS libfastscape_objects)
  add_library(${FASTSCAPELIB_OBJECTS} OBJECT ${FASTSCAPELIB_SRC_FILES})
endif()

if(BUILD_FASTSCAPELIB_STATIC)
  set(FASTSCAPELIB_STATIC libfastscape_static)

  add_library(${FASTSCAPELIB_STATIC} STATIC $<TARGET_OBJECTS:${FASTSCAPELIB_OBJECTS}>)

  set_target_properties(${FASTSCAPELIB_STATIC}
    PROPERTIES
    OUTPUT_NAME ${FASTSCAPELIB_NAME})

  install(TARGETS ${FASTSCAPELIB_STATIC} DESTINATION lib)
endif()

if(BUILD_FASTSCAPELIB_SHARED)
  set(FASTSCAPELIB_SHARED libfastscape_shared)

  # shared libraries need PIC
  set_target_properties(${FASTSCAPELIB_OBJECTS}
    PROPERTIES
    POSITION_INDEPENDENT_CODE 1)

  add_library(${FASTSCAPELIB_SHARED} SHARED $<TARGET_OBJECTS:${FASTSCAPELIB_OBJECTS}>)

  set_target_properties(${FASTSCAPELIB_SHARED}
    PROPERTIES
    OUTPUT_NAME ${FASTSCAPELIB_NAME})

  install(TARGETS ${FASTSCAPELIB_SHARED} DESTINATION lib)
endif()

# Examples
# ========

if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# Python module
# =============

if(SKBUILD)
  # TODO: remove when https://github.com/scikit-build/scikit-build/pull/495 is merged
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})

  find_package(PythonLibs REQUIRED)
  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(F2PY REQUIRED)

  set(F2PY_MODULE_NAME "_fastscapelib_fortran")

  add_python_library(${F2PY_MODULE_NAME} MODULE
    SOURCES ${FASTSCAPELIB_SRC_FILES}
    INCLUDE_DIRECTORIES ${FASTSCAPELIB_SRC_DIR}
    )

  if (UNIX)
    if (APPLE)
      set_target_properties(${F2PY_MODULE_NAME} PROPERTIES
        LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
    else()
      set_target_properties(${F2PY_MODULE_NAME} PROPERTIES
        LINK_FLAGS  '-Wl,--allow-shlib-undefined')
    endif()
  endif()

  python_extension_module(${F2PY_MODULE_NAME})

  install(
    TARGETS ${F2PY_MODULE_NAME}
    LIBRARY DESTINATION "src_python/${FASTSCAPELIB_NAME}"
    )
endif()
