include_directories(${PROJECT_SOURCE_DIR}/SRC)

# Libs linked to all of the tests
set(test_link_libs matgen)

add_subdirectory(MATGEN)

set(NVAL  10 19)
set(NRHS   2)
set(LWORK  0 100000000)
set(NPROCS 1 4)

function(add_superlumt_test target input)
  set(TEST_INPUT "${PROJECT_SOURCE_DIR}/EXAMPLE/${input}")
  foreach (s ${NRHS})
    foreach(l ${LWORK})
      foreach (p ${NPROCS})
          # LA tests
	  foreach (n ${NVAL})
          	  set(testName "${target}_${p}_${n}_${s}_${l}")
        	  add_test(
			NAME ${testName}_LA
          		COMMAND ${target} -t "LA" -n ${n} -s ${s} -l ${l} -p ${p} ) 
      	  endforeach()

          # SP tests
## CTest commands are not run in a shell's context, cannot use redirection,
## so, need "sh" as the command to be called with add_test()		
          set(testName "${target}_${p}_${s}_${l}")
          add_test(
		NAME ${testName}_SP
        	COMMAND sh -c ${target} -t "SP" -s ${s} -l ${l} -p ${p} < ${TEST_INPUT} )
      endforeach()
    endforeach()
  endforeach()

  # Add getopt dependency to target in case of MSVC
  if(MSVC)
    target_include_directories(${target} PRIVATE ${WinGetOpt_INCLUDE_DIR})
    target_link_libraries(${target} ${WinGetOpt_LIBRARY})
  endif()

endfunction(add_superlumt_test)

# prevent creating testing for Windows if WinGetOps is not found
if(MSVC AND NOT WinGetOpt_FOUND)
  message("Disabled tests pstest, pdtest, pctest, and pztest because optinal dependency WinGetOpt is missing")
  return()
endif()

if(enable_single)
  add_executable(pstest
		psdrive.c sp_sconvert.c psgst01.c psgst02.c psgst04.c psgst07.c
          	psgssv.c psgssvx.c sp_ienv.c )
  target_link_libraries(pstest ${test_link_libs})
  add_superlumt_test(pstest g10)
endif()


if(enable_double)
  add_executable(pdtest
		pddrive.c sp_dconvert.c pdgst01.c pdgst02.c pdgst04.c pdgst07.c
		pdgssv.c pdgssvx.c sp_ienv.c )
  target_link_libraries(pdtest ${test_link_libs})
  add_superlumt_test(pdtest g10)
endif()

if(enable_complex)
  add_executable(pctest
		pcdrive.c sp_cconvert.c pcgst01.c pcgst02.c pcgst04.c pcgst07.c
          	pcgssv.c pcgssvx.c sp_ienv.c )
  target_link_libraries(pctest ${test_link_libs})
  add_superlumt_test(pctest cg20.cua)
endif()

if(enable_complex16)
  add_executable(pztest
		pzdrive.c sp_zconvert.c pzgst01.c pzgst02.c pzgst04.c pzgst07.c
		pzgssv.c pzgssvx.c sp_ienv.c )
  target_link_libraries(pztest ${test_link_libs})
  add_superlumt_test(pztest cg20.cua)
endif()
