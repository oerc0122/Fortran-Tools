# CMake initialisation.
cmake_minimum_required(VERSION 3.1)

# Project name
project(Fortran_Tools Fortran)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${Fortran_Tools_BINARY_DIR})

if (mpi)
  find_package(MPI REQUIRED)
endif()
add_subdirectory("Timer")
add_subdirectory("Hash")
add_subdirectory("Units")
