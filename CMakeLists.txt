#
# CMakeLists.txt file for Peridigm
#

Cmake_minimum_required(VERSION 2.8.12)

enable_testing()
project(Peridigm)

# MPI configuration
# use the FindMPI module to find a version of MPI
# expose CMAKE_CXX_COMPILER so that the user can change it to the MPI compiler
mark_as_advanced(CLEAR CMAKE_CXX_COMPILER)

mark_as_advanced(CLEAR CMAKE_CXX_FLAGS)

if(NOT CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS "-DMPICH_IGNORE_CXX_SEEK" CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif(NOT CMAKE_CXX_FLAGS)

#
# Blas and Lapack Libraries

#
find_library(Blas_LIBRARY
  NAMES blas
  PATHS ${BLAS_LIBRARY_DIRS}
)
find_library(Lapack_LIBRARY
  NAMES lapack
  PATHS ${LAPACK_LIBRARY_DIRS}
)
set(BlasLapack_Libraries ${blas} ${lapack})

link_directories(${USER_LIBRARY_DIRS})
find_library(Fortran_Interface_LIBRARY
  NAMES libinterface.so
  PATHS ${USER_LIBRARY_DIRS}
)
find_library(User_Test_LIBRARY
  NAMES libusertest.so
  PATHS ${USER_LIBRARY_DIRS}
)
find_library(User_LIBRARY
  NAMES libusermat.so
  PATHS ${USER_LIBRARY_DIRS}
)
MESSAGE("User Libraries:")
MESSAGE("   USER_LIBRARY_DIRS = ${USER_LIBRARY_DIRS}")
MESSAGE(STATUS "User_LIBRARY: [${User_LIBRARY}]")
MESSAGE("")

if(NOT User_LIBRARY)
  MESSAGE(SEND_ERROR "Did not find user lib")
endif(NOT User_LIBRARY)

set(User_Libraries ${Fortran_Interface_LIBRARY} ${User_Test_LIBRARY} ${User_LIBRARY})

#

# DAKOTA configuration
#
option(USE_DAKOTA
  "Enable tests and functionality that depend on a DAKOTA installation."
  OFF
)

#
# Trilinos configuration
#
set(TRILINOS_DIR CACHE PATH "Path to Trilinos directory")
SET(CMAKE_PREFIX_PATH ${TRILINOS_DIR} ${CMAKE_PREFIX_PATH})
FIND_PACKAGE(Trilinos REQUIRED)
MESSAGE("Trilinos installation:")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")

# MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
# MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")

# Hack for now until Trilinos exports cmake variable for binaries directory
SET(Trilinos_BINARY_PATH ${Trilinos_LIBRARY_DIRS}/../bin)

#
# Remove PyTrilinos if it's there to avoid errors, we don't use it in Peridigm
#
LIST(REMOVE_ITEM Trilinos_LIBRARIES "pytrilinos")

#
# verbose output
# MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("")

# Check for netcdf
LIST(FIND Trilinos_TPL_LIST Netcdf Netcdf_Package_Index)

IF(Netcdf_Package_Index LESS 0)
  MESSAGE("FATAL_ERROR \n\nError:  Netcdf NOT found.  Netcdf is a REQUIRED Trilinos TPL.\n\n")
ENDIF()

# Check for yaml
LIST(FIND Trilinos_TPL_LIST yaml-cpp yaml-cpp_Package_Index)

IF(yaml-cpp_Package_Index GREATER -1)
  MESSAGE("-- Trilinos was compiled with TPL_ENABLE_yaml-cpp.\n\n   Will compile Peridigm to support YAML input files.\n\n")
  ADD_DEFINITIONS(-DUSE_YAML)
  SET(HAVE_YAML TRUE)
ENDIF()

#
# Enable performance testing
#
IF(PERFORMANCE_TEST_MACHINE)
  MESSAGE("-- Performance testing is enabled, machine = ${PERFORMANCE_TEST_MACHINE}.\n")
ELSE()
  MESSAGE("-- Performance testing is NOT enabled.\n")
ENDIF()

# Setup paths for SEACAS tools needed to run some of the tests
IF(ALBANY_SEACAS)
  find_program(
    SEACAS_EPU
    NAMES epu
    PATHS ${Trilinos_BIN_DIRS} ENV PATH
  )
  find_program(
    SEACAS_EXODIFF
    NAMES exodiff
    PATHS ${Trilinos_BIN_DIRS} ENV PATH
  )
  find_program(
    SEACAS_ALGEBRA
    NAMES algebra
    PATHS ${Trilinos_BIN_DIRS} ENV PATH
  )
ENDIF()

#
# Optional dependency on Laboratory for Computational Mechanics (LCM)
#
IF(USE_LCM)
  SET(CMAKE_PREFIX_PATH ${ALBANY_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Albany REQUIRED)
  MESSAGE("Found Albany!  Here are the details: ")
  MESSAGE("   Albany_DIR = ${Albany_DIR}")
  MESSAGE("   Albany_LIBRARIES = ${Albany_LIBRARIES}")
  MESSAGE("   Albany_INCLUDE_DIRS = ${Albany_INCLUDE_DIRS}")
  MESSAGE("   Albany_LIBRARY_DIRS = ${Albany_LIBRARY_DIRS}")
  MESSAGE("End of Albany details\n")
  MESSAGE("-- LCM is enabled, compiling with -DPERIDIGM_LCM")
  MESSAGE("")
  ADD_DEFINITIONS(-DPERIDIGM_LCM)
  SET(PERIDIGM_LCM TRUE)
  include_directories(${Albany_INCLUDE_DIRS})
  link_directories(${Albany_LIBRARY_DIRS})
  set(LCM_LIBRARY ${Albany_LIBRARIES})
ELSE()
  MESSAGE("-- LCM is NOT enabled.\n")
  SET(PERIDIGM_LCM FALSE)
  SET(LCM_LIBRARY)
ENDIF()

#
# set(Eigen3_DIR CACHE PATH "Path to Eigen directory")
# SET(CMAKE_PREFIX_PATH ${EIGEN_DIR} ${CMAKE_PREFIX_PATH})
# FIND_PACKAGE(Eigen3 3.4 REQUIRED NO_MODULE)
MESSAGE("Eigen installation:")
MESSAGE("   Eigen3_DIR = ${Eigen3_DIR}")
# MESSAGE("   Eigen3_VERSION = ${Eigen3_VERSION}")  
MESSAGE("")
include_directories(${Eigen3_DIR}) 
# Enable CJL development features
#
IF(USE_CJL)
  MESSAGE("-- CJL is enabled, compiling with -DPERIDIGM_CJL.\n")
  ADD_DEFINITIONS(-DPERIDIGM_CJL)
  SET(PERIDIGM_CJL TRUE)
ELSE()
  MESSAGE("-- CJL is NOT enabled.\n")
  SET(PERIDIGM_CJL FALSE)
  SET(CJL_LIBRARY)
ENDIF()

#
# Enable Sandia internal development
#
IF(USE_SANDIA_INTERNAL)
  MESSAGE("-- SANDIA_INTERNAL is enabled, compiling with -DPERIDIGM_SANDIA_INTERNAL.\n")
  ADD_DEFINITIONS(-DPERIDIGM_SANDIA_INTERNAL)
  SET(PERIDIGM_SANDIA_INTERNAL TRUE)
ELSE()
  MESSAGE("-- SANDIA_INTERNAL is NOT enabled.\n")
  SET(PERIDIGM_SANDIA_INTERNAL FALSE)
  SET(SANDIA_INTERNAL_LIBRARY)
ENDIF()

#
# Enable partial volume calculations
#
IF(USE_PV)
  MESSAGE("-- PV is enabled, compiling with -DPERIDIGM_PV.\n")
  ADD_DEFINITIONS(-DPERIDIGM_PV)
  SET(PERIDIGM_PV TRUE)
ELSE()
  MESSAGE("-- PV is NOT enabled.\n")
  SET(PERIDIGM_PV FALSE)
ENDIF()

# Optional Installation helpers
# Note that some of this functionality depends on CMAKE > 2.8.8
SET(INSTALL_PERIDIGM FALSE)

# Set compile definitions for logging
add_compile_definitions(CMAKE_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
add_compile_definitions(CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")


IF(ENABLE_INSTALL)
  IF(${CMAKE_VERSION} VERSION_GREATER 2.8.8 OR ${CMAKE_VERSION} VERSION_EQUAL 2.8.8)
    SET(INSTALL_PERIDIGM TRUE)
    MESSAGE("")
    MESSAGE("Installation files for Peridigm will be created")
    MESSAGE("-- install location = ${CMAKE_INSTALL_PREFIX}")
    MESSAGE("")
  ELSE()
    MESSAGE("")
    MESSAGE("Installation has been requested, but cannot continue")
    MESSAGE("-- Please upgrade CMAKE to version >= 2.8.8")
    MESSAGE("")
  ENDIF()
ENDIF()

# Optional functionality for improved quadrature
if(USE_IMPROVED_QUADRATURE)
  if(IMPROVED_QUADRATURE_DIR)
    message("-- Compiling with improved quadrature:")
    message("       IMPROVED_QUADRATURE_DIR: ${IMPROVED_QUADRATURE_DIR}\n")
    add_definitions(-DPERIDIGM_IMPROVED_QUADRATURE)
    set(PERIDIGM_IMPROVED_QUADRATURE TRUE)
    set(IMPROVED_QUADRATURE_SOURCES ${IMPROVED_QUADRATURE_DIR}/scalarDiffusionQuad/nonloc.cpp ${IMPROVED_QUADRATURE_DIR}/bondBasedQuad/nonloc.cpp)
    find_package(GSL REQUIRED)
    set(IMPROVED_QUADRATURE_LIBS ${GSL_LIBRARIES})
    set(IMPROVED_QUADRATURE_INCLUDE_DIRS ${IMPROVED_QUADRATURE_DIR} ${GSL_INCLUDE_DIRS})
  else()
    message(FATAL_ERROR "\n**** Error:  Expected IMPROVED_QUADRATURE_DIR.\n")
    set(IMPROVED_QUADRATURE_SOURCES)
    set(IMPROVED_QUADRATURE_INCLUDE_DIRS)
    set(IMPROVED_QUADRATURE_LIBS)
  endif()
endif()

IF(INSTALL_PERIDIGM)
  set(INCLUDE_INSTALL_DIR include)
  set(LIB_INSTALL_DIR lib)
  include(CMakePackageConfigHelpers)
ENDIF()

# other required libraries
set(MESH_INPUT_LIBS QuickGrid)
set(PDNEIGH_LIBS PdNeigh)
set(UTILITIES_LIBS Utilities)
set(PdMathUtilitiesLib PdMathUtilities)
set(PdMaterialUtilitiesLib PdMaterialUtilities)
set(FEMUtilitiesLib FEMUtilities)
set(PdDamageUtilitiesLib PdDamageUtilities)
set(PdAdditiveUtilitiesLib PdAdditiveUtilities)

# include for all directories
# Collection of include and lib directories
set(INCLUDE_DIRS
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${IMPROVED_QUADRATURE_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${CMAKE_SOURCE_DIR}/src/core
  ${CMAKE_SOURCE_DIR}/src/io
  ${CMAKE_SOURCE_DIR}/src/io/discretization
  ${CMAKE_SOURCE_DIR}/src/io/mesh_input
  ${CMAKE_SOURCE_DIR}/src/io/mesh_input/quick_grid
  ${CMAKE_SOURCE_DIR}/src/io/mesh_output
  ${CMAKE_SOURCE_DIR}/src/io/bond_volume
  ${CMAKE_SOURCE_DIR}/src/io/utilities
  ${CMAKE_SOURCE_DIR}/src/io/pdneigh
  ${CMAKE_SOURCE_DIR}/src/contact
  ${CMAKE_SOURCE_DIR}/src/math
  ${CMAKE_SOURCE_DIR}/src/materials
  ${CMAKE_SOURCE_DIR}/src/materials/umats
  ${CMAKE_SOURCE_DIR}/src/materials/fem
  ${CMAKE_SOURCE_DIR}/src/damage
  ${CMAKE_SOURCE_DIR}/src/additive
  ${CMAKE_SOURCE_DIR}/src/api
)

include_directories(${INCLUDE_DIRS})

#
# Some convenience definitions
#
set(COMPUTE_DIR ${CMAKE_SOURCE_DIR}/src/compute)
set(CONTACT_DIR ${CMAKE_SOURCE_DIR}/src/contact)
set(CORE_DIR ${CMAKE_SOURCE_DIR}/src/core)
set(DISCRETIZATION_DIR ${CMAKE_SOURCE_DIR}/src/io/discretization)
set(IO_DIR ${CMAKE_SOURCE_DIR}/src/io)
set(MATH_DIR ${CMAKE_SOURCE_DIR}/src/math)
set(MATERIALS_DIR ${CMAKE_SOURCE_DIR}/src/materials)
set(FEM_DIR ${CMAKE_SOURCE_DIR}/src/materials/fem)
set(DAMAGE_DIR ${CMAKE_SOURCE_DIR}/src/damage)
set(ADDITIVE_DIR ${CMAKE_SOURCE_DIR}/src/additive)
set(API_DIR ${CMAKE_SOURCE_DIR}/src/api)
set(MESH_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/src/io/mesh_output)
set(QUICK_GRID_DIR ${CMAKE_SOURCE_DIR}/src/io/mesh_input/quick_grid)
set(PDNEIGH_DIR ${CMAKE_SOURCE_DIR}/src/io/pdneigh)
set(UTILITIES_DIR ${CMAKE_SOURCE_DIR}/src/io/utilities)
set(BOND_VOL_QUICK_GRID ${CMAKE_SOURCE_DIR}/src/io/bond_volume/quick_grid)
set(MESH_CONVERTER_DIR ${CMAKE_SOURCE_DIR}/src/mesh_converter)

include_directories(
  ${COMPUTE_DIR}
  ${CONTACT_DIR}
  ${CORE_DIR}
  ${MATH_DIR}
  ${MATERIALS_DIR}
  ${FEM_DIR}
  ${DAMAGE_DIR}
  ${ADDITIVE_DIR}
  ${API_DIR}
  ${QUICK_GRID_DIR}
  ${PDNEIGH_DIR}
  ${UTILITIES_DIR}
  ${BOND_VOL_QUICK_GRID}
)

set(REQUIRED_LIBS
  ${BlasLapack_Libraries}
  ${User_Libraries}
  ${Trilinos_TPL_LIBRARIES}
)

set(UT_REQUIRED_LIBS
  ${BlasLapack_Libraries}
)

add_subdirectory(src)
add_subdirectory(test/regression)
add_subdirectory(test/verification)
add_subdirectory(test/performance)

# Create symbolic links to files in directories specified.  By default it recurses directories
# Can exclude files with flag -e <filename>
string(COMPARE EQUAL ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} buildInSource)

if(NOT ${buildInSource})
  set(UT_DIRS ${UT_DIRS} -d examples)
  set(UT_DIRS ${UT_DIRS} -d test/regression)
  set(UT_DIRS ${UT_DIRS} -d test/verification)
  set(UT_DIRS ${UT_DIRS} -d test/performance)
  execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${UT_DIRS})
  set(SCRIPT_DIRS ${SCRIPT_DIRS} -d scripts)
  execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/scripts/mirror.py -s ${CMAKE_SOURCE_DIR} -b ${CMAKE_BINARY_DIR} ${SCRIPT_DIRS})
endif(NOT ${buildInSource})

# Set symbolic link to needed seacas tools in scripts dir
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/exodiff ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/epu ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/conjoin ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/exotxt ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/nem_slice ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/nem_spread ${CMAKE_BINARY_DIR}/scripts)
execute_process(COMMAND ln -sf ${Trilinos_BINARY_PATH}/decomp ${CMAKE_BINARY_DIR}/scripts)