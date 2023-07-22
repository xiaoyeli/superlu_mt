cmake_minimum_required (VERSION 3.13)

set (CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

project (superlu_mt C)

option (BUILD_SHARED_LIBS "shared/static" OFF)

include (GNUInstallDirs)
include (CheckLibraryExists)
check_library_exists(m sin "" HAVE_LIB_M)

find_package (BLAS REQUIRED)

set (PLAT "_PTHREAD" CACHE STRING "threading flavor _PTHREAD/_OPENMP")
if (PLAT STREQUAL "_PTHREAD")
  find_package (Threads REQUIRED)
elseif (PLAT STREQUAL "_OPENMP")
  find_package (OpenMP REQUIRED)
else ()
  message (SEND_ERROR "invalid PLAT setting")
endif ()

set (SUPERLUMT_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}/superlu_mt CACHE STRING "include dir")

add_subdirectory (SRC)
