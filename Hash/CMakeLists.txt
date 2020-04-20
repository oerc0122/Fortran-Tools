add_library(hash SHARED hash.f90)

add_executable(test-hash test_hash.f90)
target_link_libraries(test-hash hash m)
