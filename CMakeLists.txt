cmake_minimum_required(VERSION 3.7.0...${CMAKE_VERSION})

project(Omega_h VERSION 9.33.0 LANGUAGES CXX)

set(Omega_h_USE_DOLFIN_DEFAULT OFF)
option(Omega_h_USE_DOLFIN "Whether to build the DOLFIN interface" ${Omega_h_USE_DOLFIN_DEFAULT})
if (Omega_h_USE_DOLFIN)
  # The FindPETSc module will use try_compile with a C code file
  enable_language(C)
endif()

include(cmake/bob.cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

bob_begin_package()

bob_option(Omega_h_CHECK_BOUNDS "Check array bounds (makes code slow too)" OFF)
bob_option(Omega_h_THROW "Errors throw exceptions instead of abort" ${USE_XSDK_DEFAULTS})
bob_input(Omega_h_DATA "" PATH "Path to omega_h-data test files")
bob_option(Omega_h_USE_EGADS "Use EGADS from ESP for geometry" OFF)
bob_input(EGADS_PREFIX "" PATH "EGADS (or ESP) installation directory")
bob_option(Omega_h_USE_Kokkos "Use Kokkos as a backend" OFF)
bob_input(Kokkos_PREFIX "" PATH "Path to Kokkos install")
bob_option(Omega_h_USE_CUDA_AWARE_MPI "Assume MPI is CUDA-aware, make use of that" OFF)
bob_input(Omega_h_VALGRIND "" STRING "Valgrind plus arguments for testing")
bob_option(Omega_h_EXAMPLES "Compile examples" OFF)

set(Omega_h_USE_ZLIB_DEFAULT ON)
bob_add_dependency(PUBLIC NAME ZLIB TARGETS ZLIB::ZLIB)

set(Omega_h_USE_Kokkos_DEFAULT OFF)
set(KokkosCore_PREFIX_DEFAULT ${Kokkos_PREFIX})
bob_add_dependency(PUBLIC NAME Kokkos TARGETS kokkos INCLUDE_DIR_VARS Kokkos_INCLUDE_DIRS)
if (Omega_h_USE_Kokkos)
  set_target_properties(kokkos PROPERTIES INTERFACE_COMPILE_OPTIONS "")
endif()

set(libMeshb_REQUIRED_VERSION 7.24)
set(Omega_h_USE_libMeshb_DEFAULT OFF)
bob_add_dependency(PUBLIC NAME libMeshb TARGETS libMeshb::Meshb.7)

set(Omega_h_USE_Gmsh_DEFAULT OFF)
bob_add_dependency(PUBLIC NAME Gmsh
    INCLUDE_DIR_VARS Gmsh_INCLUDE_DIRS
    LIBRARY_VARS Gmsh_LIBRARIES)
if(Omega_h_USE_Gmsh)
  if(Gmsh_VERSION_STRING VERSION_LESS 4.4.1)
    message(SEND_ERROR "Gmsh minimum version is 4.4.1 but found version ${Gmsh_VERSION_STRING}")
  endif()
endif()

set(Gmodel_REQUIRED_VERSION 2.1.0)
set(Omega_h_USE_Gmodel_DEFAULT OFF)
bob_add_dependency(PRIVATE NAME Gmodel TARGETS gmodel)

set(Omega_h_USE_SEACASExodus_DEFAULT OFF)
set(SEACASExodus_PREFIX_DEFAULT ${Trilinos_PREFIX})
bob_add_dependency(PUBLIC NAME SEACASExodus
    TARGETS exodus
    INCLUDE_DIR_VARS
      SEACASExodus_INCLUDE_DIRS
      SEACASExodus_TPL_INCLUDE_DIRS
    LIBRARY_VARS SEACASExodus_TPL_LIBRARIES)
if (Omega_h_USE_SEACASExodus)
  # On Windows, we must find hdf5 targets since 
  # it is linked PUBLIC
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(Omega_h_USE_HDF5_DEFAULT ON)
    bob_add_dependency(PUBLIC NAME HDF5)
  endif()
endif()

set(Omega_h_USE_pybind11_DEFAULT OFF)
bob_public_dep(pybind11)

include(${CMAKE_CURRENT_LIST_DIR}/cmake/detect_kokkos_opts.cmake)
detect_kokkos_opts()
if (Kokkos_HAS_CUDA AND (NOT Kokkos_HAS_CUDA_LAMBDA))
  message(FATAL_ERROR
          "Please reconfigure Kokkos with -DKokkos_ENABLE_Cuda_Lambda:BOOL=ON")
endif()

bob_option(Omega_h_USE_MPI "Use MPI for parallelism" OFF)
if(Omega_h_USE_MPI)
  set(MPI_CXX_SKIP_MPICXX ON) #don't need MPI C++ bindings
  bob_add_dependency(PUBLIC NAME MPI
    INCLUDE_DIR_VARS MPI_CXX_INCLUDE_DIRS
    LIBRARY_VARS MPI_CXX_LIBRARIES)
endif()
bob_option(Omega_h_USE_OpenMP "Whether to use OpenMP" "${Kokkos_HAS_OpenMP}")
bob_option(Omega_h_USE_CUDA "Whether to use CUDA" "${Kokkos_HAS_CUDA}")

if (Omega_h_USE_CUDA)
  enable_language(CUDA)
endif()

include(cmake/osh_use_dolfin.cmake)
osh_use_dolfin()

set(Omega_h_KEY_BOOLS
    Omega_h_USE_MPI
    Omega_h_USE_Kokkos
    Omega_h_USE_OpenMP
    Omega_h_USE_CUDA
    Omega_h_USE_ZLIB
    Omega_h_USE_libMeshb
    Omega_h_USE_EGADS
    Omega_h_USE_SEACASExodus
    Omega_h_USE_DOLFIN
    Omega_h_USE_dwarf
    Omega_h_CHECK_BOUNDS
    Omega_h_THROW
    Omega_h_USE_CUDA_AWARE_MPI
    Omega_h_USE_Gmsh
   )

set(Omega_h_KEY_INTS
    Omega_h_VERSION_MAJOR
    Omega_h_VERSION_MINOR
    Omega_h_VERSION_PATCH
   )

set(Omega_h_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(Omega_h_KEY_STRINGS
    Omega_h_SEMVER
    Omega_h_COMMIT
    Omega_h_CXX_FLAGS
    Omega_h_CMAKE_ARGS
   )

bob_get_commit()
bob_get_semver()

if(Omega_h_EXAMPLES)
  find_program(GMSH_EXECUTABLE gmsh REQUIRED)
  add_subdirectory(example/field_on_square)
  add_subdirectory(example/castle)
  add_subdirectory(example/gradient)
  add_subdirectory(example/laplacian)
  if(Omega_h_USE_DOLFIN)
    set(SAVED_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    string(REPLACE " -Werror" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    add_subdirectory(example/fenics/cpp)
    set(CMAKE_CXX_FLAGS "${SAVED_CXX_FLAGS}")
  endif()
endif()
add_subdirectory(src)

if (Omega_h_USE_CUDA)
  install(PROGRAMS
      ${CMAKE_CURRENT_SOURCE_DIR}/bin/nvcc_wrapper
      DESTINATION bin)
endif()

bob_end_package()
