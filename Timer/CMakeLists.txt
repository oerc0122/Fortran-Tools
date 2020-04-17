add_library(timer SHARED timer.F90)
if (mpi)
  target_link_libraries(timer ${MPI_Fortran_LIBRARIES})
endif()

add_executable(timer-test timer_test.f90)
target_link_libraries(timer-test timer m)
