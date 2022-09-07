add_executable(TestLatticeExe testlattice.cpp)
target_link_libraries(TestLatticeExe PRIVATE TOOLS CORE)
target_compile_options(TestLatticeExe PRIVATE -O0 -g)
add_test(NAME TestLattice
         COMMAND TestLatticeExe)
file(COPY test_files DESTINATION ${CMAKE_CURRENT_BINARY_DIR})