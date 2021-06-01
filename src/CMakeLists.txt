set(Omega_h_SOURCES
  Omega_h_adapt.cpp
  Omega_h_adj.cpp
  Omega_h_align.cpp
  Omega_h_amr.cpp
  Omega_h_amr_topology.cpp
  Omega_h_amr_transfer.cpp
  Omega_h_any.cpp
  Omega_h_approach.cpp
  Omega_h_array.cpp
  Omega_h_array_ops.cpp
  Omega_h_assoc.cpp
  Omega_h_base64.cpp
  Omega_h_bbox.cpp
  Omega_h_bcast.cpp
  Omega_h_bipart.cpp
  Omega_h_box.cpp
  Omega_h_build.cpp
  Omega_h_build_parser.cpp
  Omega_h_chartab.cpp
  Omega_h_check_collapse.cpp
  Omega_h_class.cpp
  Omega_h_cmdline.cpp
  Omega_h_coarsen.cpp
  Omega_h_coarsen_overshoot.cpp
  Omega_h_coarsen_qualities.cpp
  Omega_h_coarsen_topology.cpp
  Omega_h_collapse_rail.cpp
  Omega_h_comm.cpp
  Omega_h_compare.cpp
  Omega_h_confined.cpp
  Omega_h_conserve.cpp
  Omega_h_dist.cpp
  Omega_h_eigen.cpp
  Omega_h_expr.cpp
  Omega_h_fail.cpp
  Omega_h_fence.cpp
  Omega_h_file.cpp
  Omega_h_filesystem.cpp
  Omega_h_finite_automaton.cpp
  Omega_h_fix.cpp
  Omega_h_future.cpp
  Omega_h_ghost.cpp
  Omega_h_globals.cpp
  Omega_h_gmsh.cpp
  Omega_h_grammar.cpp
  Omega_h_graph.cpp
  Omega_h_hilbert.cpp
  Omega_h_histogram.cpp
  Omega_h_indset.cpp
  Omega_h_inertia.cpp
  Omega_h_input.cpp
  Omega_h_int128.cpp
  Omega_h_int_scan.cpp
  Omega_h_language.cpp
  Omega_h_laplace.cpp
  Omega_h_library.cpp
  Omega_h_linpart.cpp
  Omega_h_malloc.cpp
  Omega_h_map.cpp
  Omega_h_mark.cpp
  Omega_h_math_lang.cpp
  Omega_h_matrix.cpp
  Omega_h_mesh.cpp
  Omega_h_metric.cpp
  Omega_h_metric_input.cpp
  Omega_h_migrate.cpp
  Omega_h_modify.cpp
  Omega_h_owners.cpp
  Omega_h_parser.cpp
  Omega_h_parser_graph.cpp
  Omega_h_pool.cpp
  Omega_h_print.cpp
  Omega_h_profile.cpp
  Omega_h_quality.cpp
  Omega_h_reader.cpp
  Omega_h_recover.cpp
  Omega_h_refine.cpp
  Omega_h_refine_qualities.cpp
  Omega_h_refine_topology.cpp
  Omega_h_regex.cpp
  Omega_h_regress.cpp
  Omega_h_remotes.cpp
  Omega_h_reorder.cpp
  Omega_h_scatterplot.cpp
  Omega_h_shape.cpp
  Omega_h_shared_alloc.cpp
  Omega_h_simplify.cpp
  Omega_h_sort.cpp
  Omega_h_stacktrace.cpp
  Omega_h_surface.cpp
  Omega_h_swap.cpp
  Omega_h_swap2d.cpp
  Omega_h_swap2d_qualities.cpp
  Omega_h_swap2d_topology.cpp
  Omega_h_swap3d.cpp
  Omega_h_swap3d_qualities.cpp
  Omega_h_swap3d_topology.cpp
  Omega_h_tag.cpp
  Omega_h_timer.cpp
  Omega_h_transfer.cpp
  Omega_h_unmap_mesh.cpp
  Omega_h_vector.cpp
  Omega_h_vtk.cpp
  Omega_h_xml.cpp
  Omega_h_xml_lite.cpp
  Omega_h_yaml.cpp
  )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
else()
  set(Omega_h_SOURCES ${Omega_h_SOURCES}
    Omega_h_random.cpp
    Omega_h_coarsen_flip.cpp
	)
endif()

if (Omega_h_USE_CUDA)
  if (NOT Omega_h_USE_Kokkos)
    set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_for.cpp)
  endif()
else()
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_overlay.cpp)
endif()

if(Omega_h_USE_libMeshb)
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_meshb.cpp)
endif()

if(Omega_h_USE_EGADS)
  message(STATUS "Omega_h_USE_EGADS: ${Omega_h_USE_EGADS}")
  if(NOT EGADS_PREFIX)
    message(FATAL_ERROR "Omega_h_USE_EGADS=${Omega_h_USE_EGADS} but EGADS_PREFIX not set!")
  endif()
  find_path(EGADS_INCLUDE_DIR NAMES egads.h
      PATHS
      "${EGADS_PREFIX}/EngSketchPad/include"
      "${EGADS_PREFIX}/include"
      NO_DEFAULT_PATH)
  if(NOT EGADS_INCLUDE_DIR)
    message(FATAL_ERROR "could not find path to \"egads.h\"")
  endif()
  find_library(EGADS_LIBRARY NAMES egads
      PATHS
      "${EGADS_PREFIX}/EngSketchPad/lib"
      NO_DEFAULT_PATH)
  if(NOT EGADS_LIBRARY)
    message(FATAL_ERROR "could not find EGADS_LIBRARY")
  endif()
  message(STATUS "EGADS_INCLUDE_DIR: ${EGADS_INCLUDE_DIR}")
  message(STATUS "EGADS_LIBRARY: ${EGADS_LIBRARY}")
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_egads.cpp)
endif()

if(Omega_h_USE_SEACASExodus)
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_exodus.cpp)
endif()

if(Omega_h_USE_DOLFIN)
  set(Omega_h_SOURCES ${Omega_h_SOURCES} Omega_h_dolfin.cpp)
endif()

if (Omega_h_USE_CUDA AND NOT Omega_h_USE_Kokkos)
  set_source_files_properties(${Omega_h_SOURCES} PROPERTIES LANGUAGE CUDA)
endif()

add_library(omega_h ${Omega_h_SOURCES})

set_property(TARGET omega_h PROPERTY CXX_STANDARD "14")
set_property(TARGET omega_h PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET omega_h PROPERTY CXX_EXTENSIONS OFF)

bob_library_includes(omega_h)

get_filename_component(TPL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../tpl" ABSOLUTE)
target_include_directories(omega_h PRIVATE $<BUILD_INTERFACE:${TPL_PATH}>)

if (Omega_h_USE_OpenMP)
  target_compile_options(omega_h PUBLIC -fopenmp)
endif()

if (Omega_h_USE_CUDA AND NOT Omega_h_USE_Kokkos)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(omega_h PUBLIC
        "SHELL:--compiler-options -W,-Wall,-Wextra,-Werror,-Wno-noexcept-type"
        "SHELL:--Werror cross-execution-space-call,deprecated-declarations"
        --expt-extended-lambda)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(omega_h PUBLIC
        "SHELL:--compiler-options /wd6294,/wd26451,/wd26439"
        "SHELL:--Werror cross-execution-space-call,deprecated-declarations"
        --expt-extended-lambda)
  endif()
else()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(omega_h PRIVATE
        -W
        -Wall
        -Wextra
        -Werror
        -Wimplicit-fallthrough
        )
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(omega_h PRIVATE
        -W
        -Wall
        -Wextra
        -Werror
        )
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(omega_h PUBLIC
        "/wd26451" #arithmetic overflow, long(int * int), shows up in much vector access
		    "/wd26439" #noexcept advice. We've determined that noexcept is useless for performance, we don't use it
        )
  endif()
endif()

bob_link_dependency(omega_h PUBLIC Kokkos)

bob_link_dependency(omega_h PUBLIC libMeshb)

bob_link_dependency(omega_h PUBLIC Gmsh)

if(Omega_h_USE_EGADS)
  target_include_directories(omega_h PUBLIC "${EGADS_INCLUDE_DIR}")
  target_link_libraries(omega_h PUBLIC "${EGADS_LIBRARY}")
endif()

bob_link_dependency(omega_h PUBLIC SEACASExodus)

bob_link_dependency(omega_h PUBLIC ZLIB)

if (Omega_h_USE_MPI)
  target_link_libraries(omega_h PUBLIC MPI::MPI_CXX)
endif()

if (Omega_h_USE_dwarf)
  target_include_directories(omega_h PRIVATE "${LIBDWARF_INCLUDE_DIRS}")
  target_link_libraries(omega_h PUBLIC "${LIBDWARF_LIBRARIES}")
endif()

if(Omega_h_USE_DOLFIN)
  target_link_libraries(omega_h PUBLIC Boost::timer dolfin)
  target_include_directories(omega_h SYSTEM PUBLIC ${DOLFIN_INCLUDE_DIRS})
  target_include_directories(omega_h SYSTEM PUBLIC ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set_target_properties(omega_h PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

bob_export_target(omega_h)

function(osh_add_exe EXE_NAME)
  if (Omega_h_USE_CUDA AND NOT Omega_h_USE_Kokkos)
    set_source_files_properties(${EXE_NAME}.cpp PROPERTIES LANGUAGE CUDA)
  endif()
  add_executable(${EXE_NAME} ${EXE_NAME}.cpp)
  set_property(TARGET ${EXE_NAME} PROPERTY CXX_STANDARD "14")
  set_property(TARGET ${EXE_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)
  set_property(TARGET ${EXE_NAME} PROPERTY CXX_EXTENSIONS OFF)
  target_link_libraries(${EXE_NAME} PRIVATE omega_h)
endfunction(osh_add_exe)

macro(osh_add_util EXE_NAME)
  osh_add_exe(${EXE_NAME})
  bob_export_target(${EXE_NAME})
endmacro(osh_add_util)

osh_add_util(msh2osh)
osh_add_util(osh2vtk)
osh_add_util(oshdiff)
osh_add_util(vtkdiff)
osh_add_util(osh_box)
osh_add_util(osh_part)
osh_add_util(osh_scale)
osh_add_util(osh_reorder)
osh_add_util(osh_fix)
osh_add_util(osh_eval_implied)
osh_add_util(osh_calc)
if(Omega_h_USE_libMeshb)
  osh_add_util(meshb2osh)
  osh_add_util(osh2meshb)
  osh_add_util(osh_read_solb)
  osh_add_util(ugawg_stats)
  osh_add_util(ugawg_cone_cone)
  osh_add_util(ugawg_parallel_adapt)
  if (Omega_h_USE_EGADS)
    osh_add_util(ugawg_solution)
  endif()
endif()
if(Omega_h_USE_SEACASExodus)
  osh_add_util(exo2osh)
  osh_add_util(osh2exo)
  osh_add_util(osh_class_morph)
endif()
osh_add_util(osh_adapt)
osh_add_util(osh_filesystem)
osh_add_util(ascii_vtk2osh)

if(BUILD_TESTING)
  if(Omega_h_USE_MPI)
    get_filename_component(COMPILER_DIR ${CMAKE_CXX_COMPILER} PATH)
    find_program(MPIRUN mpirun PATHS ${COMPILER_DIR}
        DOC "Location of the mpirun program")
    if(NOT MPIRUN)
      message(STATUS "USE_MPI and BUILD_TESTING but MPIRUN not found!")
    endif()
  endif()

  function(test_func TEST_NAME NUM_PROCS)
    string(REPLACE " " ";" VALGRIND "${Omega_h_VALGRIND}")
    if(MPIRUN)
      set(TEST_STR ${MPIRUN} -np ${NUM_PROCS} ${VALGRIND} ${ARGN})
    else()
      if(NOT (${NUM_PROCS} EQUAL "1"))
        message(STATUS "test ${TEST_NAME} ignored because MPIRUN not found!")
        return()
      endif()
      set(TEST_STR ${VALGRIND} ${ARGN})
    endif()
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_STR})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
      set_tests_properties(${TEST_NAME} PROPERTIES ENVIRONMENT "PATH=${CMAKE_INSTALL_PREFIX}\\bin;$ENV{PATH}" )
    endif()
  endfunction(test_func)

  osh_add_exe(unit_math)
  set(TEST_EXES ${TEST_EXES} unit_math)
  test_func(run_unit_math 1 ./unit_math)
  osh_add_exe(unit_array_algs)
  set(TEST_EXES ${TEST_EXES} unit_array_algs)
  test_func(run_unit_array_algs 1 ./unit_array_algs)
  osh_add_exe(unit_mesh)
  set(TEST_EXES ${TEST_EXES} unit_mesh)
  test_func(run_unit_mesh 1 ./unit_mesh)
  osh_add_exe(unit_io)
  set(TEST_EXES ${TEST_EXES} unit_io)
  test_func(run_unit_io 1 ./unit_io)
  osh_add_exe(unit_parser)
  set(TEST_EXES ${TEST_EXES} unit_parser)
  test_func(run_unit_parser 1 ./unit_parser)
  osh_add_exe(r3d_test)
  test_func(run_test_r3d 1 ./r3d_test)
  osh_add_exe(mpi_tests)
  set(TEST_EXES ${TEST_EXES} mpi_tests)
  if(Omega_h_USE_MPI)
    test_func(run_mpi_tests 4 ./mpi_tests)
  else()
    test_func(run_mpi_tests 1 ./mpi_tests)
  endif()
  osh_add_exe(1d_test)
  test_func(serial_1d_test 1 ./1d_test)
  if(Omega_h_USE_MPI)
    test_func(parallel_1d_test 2 ./1d_test)
    test_func(run_mpi_unit_io 2 ./unit_io)
  endif()
  osh_add_exe(corner_test)
  test_func(run_corner_test 1 ./corner_test)
  osh_add_exe(coarsen_test)
  test_func(run_coarsen_test 1 ./coarsen_test)

  osh_add_exe(2d_conserve_test)
  test_func(serial_2d_conserve 1 ./2d_conserve_test)
  if(Omega_h_USE_MPI)
    test_func(parallel_2d_conserve 2 ./2d_conserve_test)
  endif()
  osh_add_exe(warp_test)
  set(TEST_EXES ${TEST_EXES} warp_test)
  test_func(warp_test_serial 1 ./warp_test)
  if(Omega_h_USE_MPI)
    test_func(warp_test_parallel 2 ./warp_test)
  endif()
  osh_add_exe(3d_conserve_test)
  if(Omega_h_DATA)
    set(TEST_EXES ${TEST_EXES} 3d_conserve_test)
    test_func(serial_3d_conserve 1 ./3d_conserve_test
              ${Omega_h_DATA}/ball_in_cube.msh)
    if(Omega_h_USE_MPI)
      test_func(parallel_3d_conserve 2 ./3d_conserve_test
                ${Omega_h_DATA}/ball_in_cube.msh)
    endif()
  endif()
  osh_add_exe(ring_test)
  if(Omega_h_DATA)
    set(TEST_EXES ${TEST_EXES} ring_test)
    test_func(run_ring_test 1 ./ring_test ${Omega_h_DATA}/ring.msh)
  endif()
  osh_add_exe(aniso_test)
  test_func(run_aniso_test 1 ./aniso_test)
  osh_add_exe(imr_cases)
  osh_add_exe(curvature_test)
  if(Omega_h_DATA)
    set(TEST_EXES ${TEST_EXES} curvature_test)
    test_func(ring_curvature 1 ./curvature_test
              ${Omega_h_DATA} ring)
    test_func(ball_curvature 1 ./curvature_test
              ${Omega_h_DATA} ball)
    if(Omega_h_USE_MPI)
      test_func(parallel_ring_curvature 2 ./curvature_test
                ${Omega_h_DATA} ring)
      test_func(parallel_ball_curvature 2 ./curvature_test
                ${Omega_h_DATA} ball)
    endif()
  endif()
  osh_add_exe(bend_test)
  osh_add_exe(test_degree)
  osh_add_exe(test_dihedral)
  if(Omega_h_USE_libMeshb)
    osh_add_exe(ugawg_linear)
    if(Omega_h_USE_EGADS)
      osh_add_exe(ugawg_cylinder)
    endif()
    osh_add_exe(ugawg_hsc)
  endif()
  osh_add_exe(slope_test)
  if(Omega_h_USE_EGADS)
    osh_add_exe(egads_test)
  endif()
  osh_add_exe(advect2d_test)
  if(Omega_h_DATA)
    set(TEST_EXES ${TEST_EXES} advect2d_test)
    test_func(advect2d_serial 1 ./advect2d_test ${Omega_h_DATA}/advect2d.osh)
    if(Omega_h_USE_MPI)
      test_func(advect2d_parallel 2 ./advect2d_test ${Omega_h_DATA}/advect2d.osh)
    endif()
  endif()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
      DEPENDS ${TEST_EXES})
  osh_add_exe(cylinder_adapt_test)
  if (Omega_h_DATA)
    test_func(cylinder_adapt_test 1 ./cylinder_adapt_test ${Omega_h_DATA}/cylinder.msh)
  endif()
  if (Omega_h_USE_DOLFIN)
    osh_add_exe(dolfin_test)
  endif()
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  else()
    osh_add_exe(random_test)
    target_include_directories(random_test PRIVATE ${TPL_PATH})
    test_func(run_random_test 1 ./random_test)
  endif()
  osh_add_exe(hypercube_test)
  osh_add_exe(amr_test)
  test_func(amr_test 1 ./amr_test)
  osh_add_exe(amr_test2)
  test_func(amr_test2 1 ./amr_test2)
  osh_add_exe(refine_scale)
  osh_add_exe(amr_mpi_test)
endif()

bob_config_header("${CMAKE_CURRENT_BINARY_DIR}/Omega_h_config.h")
set(Omega_h_HEADERS
  "${CMAKE_CURRENT_BINARY_DIR}/Omega_h_config.h"
  Omega_h_adapt.hpp
  Omega_h_adj.hpp
  Omega_h_affine.hpp
  Omega_h_align.hpp
  Omega_h_amr.hpp
  Omega_h_any.hpp
  Omega_h_array.hpp
  Omega_h_array_ops.hpp
  Omega_h_assoc.hpp
  Omega_h_base64.hpp
  Omega_h_bbox.hpp
  Omega_h_box.hpp
  Omega_h_build.hpp
  Omega_h_class.hpp
  Omega_h_cmdline.hpp
  Omega_h_comm.hpp
  Omega_h_compare.hpp
  Omega_h_dbg.hpp
  Omega_h_defines.hpp
  Omega_h_dist.hpp
  Omega_h_eigen.hpp
  Omega_h_element.hpp
  Omega_h_expr.hpp
  Omega_h_fail.hpp
  Omega_h_fence.hpp
  Omega_h_few.hpp
  Omega_h_file.hpp
  Omega_h_filesystem.hpp
  Omega_h_finite_automaton.hpp
  Omega_h_for.hpp
  Omega_h_functors.hpp
  Omega_h_future.hpp
  Omega_h_grammar.hpp
  Omega_h_graph.hpp
  Omega_h_hypercube.hpp
  Omega_h_indset.hpp
  Omega_h_indset_inline.hpp
  Omega_h_input.hpp
  Omega_h_int128.hpp
  Omega_h_int_iterator.hpp
  Omega_h_int_scan.hpp
  Omega_h_kokkos.hpp
  Omega_h_language.hpp
  Omega_h_library.hpp
  Omega_h_lie.hpp
  Omega_h_macros.h
  Omega_h_malloc.hpp
  Omega_h_map.hpp
  Omega_h_mark.hpp
  Omega_h_math_lang.hpp
  Omega_h_matrix.hpp
  Omega_h_mesh.hpp
  Omega_h_metric.hpp
  Omega_h_mpi.h
  Omega_h_owners.hpp
  Omega_h_parser.hpp
  Omega_h_print.hpp
  Omega_h_profile.hpp
  Omega_h_qr.hpp
  Omega_h_quality.hpp
  Omega_h_random.hpp
  Omega_h_rbtree.hpp
  Omega_h_reader.hpp
  Omega_h_reader_tables.hpp
  Omega_h_recover.hpp
  Omega_h_reduce.hpp
  Omega_h_remotes.hpp
  Omega_h_scalar.hpp
  Omega_h_scan.hpp
  Omega_h_scatterplot.hpp
  Omega_h_shape.hpp
  Omega_h_shared_alloc.hpp
  Omega_h_simplex.hpp
  Omega_h_sort.hpp
  Omega_h_stacktrace.hpp
  Omega_h_std_vector.hpp
  Omega_h_svd.hpp
  Omega_h_table.hpp
  Omega_h_tag.hpp
  Omega_h_template_up.hpp
  Omega_h_timer.hpp
  Omega_h_vector.hpp
  Omega_h_vtk.hpp
  Omega_h_xml_lite.hpp
  Omega_h_xml.hpp
  r3d.hpp
  )

if (NOT Omega_h_USE_CUDA)
  set(Omega_h_HEADERS ${Omega_h_HEADERS} Omega_h_overlay.hpp)
endif()

if(Omega_h_USE_DOLFIN)
  set(Omega_h_HEADERS ${Omega_h_HEADERS} Omega_h_dolfin.hpp)
endif()

install(FILES ${Omega_h_HEADERS} DESTINATION include)

if (Omega_h_USE_pybind11)
  set(PYBIND11_SOURCES
      PyOmega_h.cpp
      PyOmega_h_defines.cpp
      PyOmega_h_array.cpp
      PyOmega_h_comm.cpp
      PyOmega_h_library.cpp
      PyOmega_h_mesh.cpp
      PyOmega_h_build.cpp
      PyOmega_h_adapt.cpp
      PyOmega_h_file.cpp
      PyOmega_h_class.cpp
      )
  if (Omega_h_USE_DOLFIN)
    set(PYBIND11_SOURCES ${PYBIND11_SOURCES}
        PyOmega_h_dolfin.cpp
       )
  endif()
  pybind11_add_module(PyOmega_h
      ${PYBIND11_SOURCES})
  target_link_libraries(PyOmega_h PUBLIC omega_h)
  if (PYTHON_VERSION_MAJOR EQUAL 3)
    set(PyOmega_h_DEST "lib/python3/dist-packages")
  else()
    set(PyOmega_h_DEST "lib/python/dist-packages")
  endif()
  install(TARGETS PyOmega_h
      ARCHIVE DESTINATION "${PyOmega_h_DEST}"
      LIBRARY DESTINATION "${PyOmega_h_DEST}")
endif()

add_custom_target(test_install
    COMMAND ${CMAKE_INSTALL_PREFIX}/bin/osh_box 1 1 1 2 2 2 ${CMAKE_CURRENT_BINARY_DIR}/box.osh
    COMMAND ${CMAKE_INSTALL_PREFIX}/bin/osh_scale ${CMAKE_CURRENT_BINARY_DIR}/box.osh 100 ${CMAKE_CURRENT_BINARY_DIR}/box_100.osh
    COMMAND ${CMAKE_INSTALL_PREFIX}/bin/osh2vtk ${CMAKE_CURRENT_BINARY_DIR}/box_100.osh ${CMAKE_CURRENT_BINARY_DIR}/box_100_vtk
    COMMENT "Test installed Omega_h utilities")

bob_end_subdir()
