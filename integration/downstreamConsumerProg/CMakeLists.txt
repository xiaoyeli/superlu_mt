cmake_minimum_required(VERSION 3.22)

project(downstreamConsumerProg LANGUAGES C CXX)

enable_testing()

find_package(SuperLU_MT CONFIG REQUIRED)

foreach(candidate SuperLU_MT::superlu_mt_PTHREAD SuperLU_MT::superlu_mt_OPENMP)
  if(TARGET ${candidate})
    set(SUPERLUMT_TARGET ${candidate})
    break()
  endif()
endforeach()

foreach(candidate SuperLU_MT::ddsuperlu_mt_PTHREAD SuperLU_MT::ddsuperlu_mt_OPENMP)
  if(TARGET ${candidate})
    set(SUPERLUMT_DD_TARGET ${candidate})
    break()
  endif()
endforeach()

if(NOT SUPERLUMT_TARGET OR NOT SUPERLUMT_DD_TARGET)
  message(FATAL_ERROR "Expected both SuperLU_MT base and dd imported targets")
endif()

add_library(downstreamConsumerProg src/downstream_consumer.cpp)
target_include_directories(downstreamConsumerProg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(downstreamConsumerProg PUBLIC cxx_std_17)
target_link_libraries(downstreamConsumerProg PUBLIC ${SUPERLUMT_TARGET} ${SUPERLUMT_DD_TARGET})

add_executable(downstreamConsumerProg_test tests/downstream_consumer_test.cpp)
target_link_libraries(downstreamConsumerProg_test PRIVATE downstreamConsumerProg)
add_test(NAME downstreamConsumerProg_test COMMAND downstreamConsumerProg_test)
