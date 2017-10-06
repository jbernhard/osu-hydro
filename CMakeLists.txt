cmake_minimum_required(VERSION 3.4)
project(osu-hydro Fortran)

# default build type: Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# default install prefix: ~/.local
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local"
      CACHE PATH "Install path prefix, prepended onto install directories."
      FORCE)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

# By default, optimize for the system's native architecture.  Disable via the
# NATIVE option.  In addition, detect if another architecture flag is already
# set and do not override it.
option(NATIVE "compile for native architecture" ON)
if(NATIVE)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    # Intel compiler: search for -m, -x, -ax flags; if not found add -xHost.
    if(NOT CMAKE_Fortran_FLAGS MATCHES "(^| )-(m|x|ax)[^ ]+")
      string(APPEND CMAKE_Fortran_FLAGS " -xHost")
    endif()
  else()
    # All other compilers (gcc-like): search for -m flags;
    # if not found add -march=native.
    if(NOT CMAKE_Fortran_FLAGS MATCHES "(^| )-m[^ ]+")
      string(APPEND CMAKE_Fortran_FLAGS " -march=native")
    endif()
  endif()
endif()

string(APPEND CMAKE_Fortran_FLAGS " -cpp")

option(WARNINGS "enable compiler warnings" OFF)
if(WARNINGS)
  string(APPEND CMAKE_Fortran_FLAGS " -Wall -Wextra")
endif()

install(FILES ${PROJECT_NAME}.conf DESTINATION share/${PROJECT_NAME})

add_subdirectory(eos)
add_subdirectory(src)
