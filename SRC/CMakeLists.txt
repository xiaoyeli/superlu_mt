
find_package(Python3 COMPONENTS Interpreter REQUIRED)

file (GLOB sources "*.c")
add_library (superlu_mt${PLAT} ${sources})
set_target_properties (superlu_mt${PLAT} PROPERTIES POSITION_INDEPENDENT_CODE ON)

set_property(TARGET superlu_mt${PLAT} PROPERTY C_STANDARD 99)

if (HAVE_LIB_M)
  target_link_libraries (superlu_mt${PLAT} PUBLIC m)
endif ()

target_link_libraries (superlu_mt${PLAT} PUBLIC ${BLAS_LIB})
target_compile_definitions (superlu_mt${PLAT} PRIVATE Add_)

if (PLAT STREQUAL "_PTHREAD")
  target_compile_definitions (superlu_mt${PLAT} PUBLIC __PTHREAD)
  target_link_libraries (superlu_mt${PLAT} PUBLIC Threads::Threads)
elseif (PLAT STREQUAL "_OPENMP")
  target_compile_definitions (superlu_mt${PLAT} PUBLIC __OPENMP)
  target_link_libraries (superlu_mt${PLAT} PUBLIC OpenMP::OpenMP_C)
endif ()

if (LONGINT)
  target_compile_definitions (superlu_mt${PLAT} PUBLIC _LONGINT)
endif ()

target_include_directories (superlu_mt${PLAT}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${SUPERLUMT_INSTALL_INCLUDEDIR}>
)
install (TARGETS superlu_mt${PLAT} EXPORT SuperLUMTTargets DESTINATION ${CMAKE_INSTALL_LIBDIR})

file (GLOB headers "*.h")
install (FILES ${headers} DESTINATION ${SUPERLUMT_INSTALL_INCLUDEDIR})

if (enable_doubledouble)
  set(dd_generated_dir ${CMAKE_CURRENT_BINARY_DIR}/dd_generated)
  set(dd_generated_header ${dd_generated_dir}/slu_mt_dd_internal.h)
  set(dd_sources
    dgstrs.c
    dmyblas2.c
    dsp_blas2.c
    pdmemory.c
    pdgssv.c
    pdgstrf.c
    pdgstrf_bmod1D.c
    pdgstrf_bmod1D_mv2.c
    pdgstrf_bmod2D.c
    pdgstrf_bmod2D_mv2.c
    pdgstrf_column_bmod.c
    pdgstrf_column_dfs.c
    pdgstrf_copy_to_ucol.c
    pdgstrf_factor_snode.c
    pdgstrf_init.c
    pdgstrf_panel_bmod.c
    pdgstrf_panel_dfs.c
    pdgstrf_pivotL.c
    pdgstrf_snode_bmod.c
    pdgstrf_snode_dfs.c
    pdgstrf_thread.c
    pdgstrf_thread_finalize.c
    pdgstrf_thread_init.c
  )
  set(dd_generated_sources "")
  foreach(dd_source ${dd_sources})
    string(REPLACE ".c" "_dd.cpp" dd_generated_source ${dd_source})
    list(APPEND dd_generated_sources ${dd_generated_dir}/${dd_generated_source})
  endforeach()

  add_custom_command(
    OUTPUT ${dd_generated_header} ${dd_generated_sources}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${dd_generated_dir}
    COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tools/generate-superlu-dd.py
      --source-dir ${CMAKE_CURRENT_SOURCE_DIR}
      --output-dir ${dd_generated_dir}
    DEPENDS
      ${PROJECT_SOURCE_DIR}/tools/generate-superlu-dd.py
      ${CMAKE_CURRENT_SOURCE_DIR}/slu_mt_ddefs.h
      ${CMAKE_CURRENT_SOURCE_DIR}/slu_mt_dddefs.h
      ${dd_sources}
    VERBATIM
  )

  add_library (ddsuperlu_mt${PLAT}
    ${dd_generated_sources}
    ${CMAKE_CURRENT_SOURCE_DIR}/slu_mt_ddutil.cpp
  )
  add_dependencies(ddsuperlu_mt${PLAT} superlu_mt${PLAT})
  set_target_properties(ddsuperlu_mt${PLAT} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME ddSuperLU_mt${PLAT}
  )
  target_compile_features(ddsuperlu_mt${PLAT} PRIVATE cxx_std_17)
  target_compile_options(ddsuperlu_mt${PLAT} PRIVATE -fpermissive)
  target_include_directories(ddsuperlu_mt${PLAT}
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<BUILD_INTERFACE:${dd_generated_dir}>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/external/doubledouble/include>
      $<INSTALL_INTERFACE:${SUPERLUMT_INSTALL_INCLUDEDIR}>
  )
  target_link_libraries(ddsuperlu_mt${PLAT} PUBLIC superlu_mt${PLAT})
  install(TARGETS ddsuperlu_mt${PLAT} EXPORT SuperLUMTTargets DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/slu_mt_dddefs.h
    ${PROJECT_SOURCE_DIR}/external/doubledouble/include/doubledouble.h
    DESTINATION ${SUPERLUMT_INSTALL_INCLUDEDIR}
  )
endif ()
