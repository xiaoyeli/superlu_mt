include_directories(${PROJECT_SOURCE_DIR}/SRC)

# build examples for target "all" only if enable_examples is set
if(NOT enable_examples)
  set(_DEPENDENCY_ALL "EXCLUDE_FROM_ALL")
endif()

# targets to build examples
add_custom_target(examples_complex)
add_dependencies(examples_complex
		pclinsol pclinsolx pclinsolx1 pclinsolx2 pcrepeat pcspmd)
		
add_custom_target(examples_double)
add_dependencies(examples_double
		pdlinsol pdlinsolx pdlinsolx1 pdlinsolx2 pdrepeat pdspmd)
		
add_custom_target(examples_float)
add_dependencies(examples_float
		pslinsol pslinsolx pslinsolx1 pslinsolx2 psrepeat psspmd)
		
add_custom_target(examples_doublecomplex)
add_dependencies(examples_doublecomplex
		pzlinsol pzlinsolx pzlinsolx1 pzlinsolx2 pzrepeat pzspmd)

add_custom_target(examples)
add_dependencies(examples
                 examples_complex
                 examples_double
                 examples_float
                 examples_doublecomplex)


# examples for complex
add_executable(pclinsol  ${_DEPENDENCY_ALL} pclinsol.c)
target_link_libraries(pclinsol superlu_mt${PLAT})

add_executable(pclinsolx  ${_DEPENDENCY_ALL} pclinsolx.c)
target_link_libraries(pclinsolx superlu_mt${PLAT})

add_executable(pclinsolx1  ${_DEPENDENCY_ALL} pclinsolx1.c)
target_link_libraries(pclinsolx1 superlu_mt${PLAT})

add_executable(pclinsolx2  ${_DEPENDENCY_ALL} pclinsolx2.c)
target_link_libraries(pclinsolx2 superlu_mt${PLAT})

add_executable(pcrepeat  ${_DEPENDENCY_ALL} pcrepeat.c)
target_link_libraries(pcrepeat superlu_mt${PLAT})

add_executable(pcspmd  ${_DEPENDENCY_ALL} pcspmd.c)
target_link_libraries(pcspmd superlu_mt${PLAT})

# examples for double
add_executable(pdlinsol  ${_DEPENDENCY_ALL} pdlinsol.c)
target_link_libraries(pdlinsol superlu_mt${PLAT})

add_executable(pdlinsolx  ${_DEPENDENCY_ALL} pdlinsolx.c)
target_link_libraries(pdlinsolx superlu_mt${PLAT})

add_executable(pdlinsolx1  ${_DEPENDENCY_ALL} pdlinsolx1.c)
target_link_libraries(pdlinsolx1 superlu_mt${PLAT})

add_executable(pdlinsolx2  ${_DEPENDENCY_ALL} pdlinsolx2.c)
target_link_libraries(pdlinsolx2 superlu_mt${PLAT})

add_executable(pdrepeat  ${_DEPENDENCY_ALL} pdrepeat.c)
target_link_libraries(pdrepeat superlu_mt${PLAT})

add_executable(pdspmd  ${_DEPENDENCY_ALL} pdspmd.c)
target_link_libraries(pdspmd superlu_mt${PLAT})

# examples for float
add_executable(pslinsol  ${_DEPENDENCY_ALL} pslinsol.c)
target_link_libraries(pslinsol superlu_mt${PLAT})

add_executable(pslinsolx  ${_DEPENDENCY_ALL} pslinsolx.c)
target_link_libraries(pslinsolx superlu_mt${PLAT})

add_executable(pslinsolx1  ${_DEPENDENCY_ALL} pslinsolx1.c)
target_link_libraries(pslinsolx1 superlu_mt${PLAT})

add_executable(pslinsolx2  ${_DEPENDENCY_ALL} pslinsolx2.c)
target_link_libraries(pslinsolx2 superlu_mt${PLAT})

add_executable(psrepeat  ${_DEPENDENCY_ALL} psrepeat.c)
target_link_libraries(psrepeat superlu_mt${PLAT})

add_executable(psspmd  ${_DEPENDENCY_ALL} psspmd.c)
target_link_libraries(psspmd superlu_mt${PLAT})


# examples for double complex
add_executable(pzlinsol  ${_DEPENDENCY_ALL} pzlinsol.c)
target_link_libraries(pzlinsol superlu_mt${PLAT})

add_executable(pzlinsolx  ${_DEPENDENCY_ALL} pzlinsolx.c)
target_link_libraries(pzlinsolx superlu_mt${PLAT})

add_executable(pzlinsolx1  ${_DEPENDENCY_ALL} pzlinsolx1.c)
target_link_libraries(pzlinsolx1 superlu_mt${PLAT})

add_executable(pzlinsolx2  ${_DEPENDENCY_ALL} pzlinsolx2.c)
target_link_libraries(pzlinsolx2 superlu_mt${PLAT})

add_executable(pzrepeat  ${_DEPENDENCY_ALL} pzrepeat.c)
target_link_libraries(pzrepeat superlu_mt${PLAT})

add_executable(pzspmd  ${_DEPENDENCY_ALL} pzspmd.c)
target_link_libraries(pzspmd superlu_mt${PLAT})


if(MSVC AND WinGetOpt_FOUND)
  setx(NEEDS_GETOPT ${examples})

  foreach(proj ${NEEDS_GETOPT})
    target_include_directories(${proj} PRIVATE ${WinGetOpt_INCLUDE_DIR})
    target_link_libraries(${proj} ${WinGetOpt_LIBRARY})
  endforeach(proj ${NEEDS_GETOPT})
endif()
