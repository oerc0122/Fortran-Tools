add_library(units SHARED units.f90)
target_include_directories(units PRIVATE "${Fortran_Tools_BINARY_DIR}/Hash")
target_link_libraries(units hash m)

add_executable(units-test units_test.f90)
target_link_libraries(units-test units hash m)

add_executable(whats whats.f90)
target_link_libraries(whats units hash m)
