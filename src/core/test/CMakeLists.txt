add_executable(TestSimExe testsim.cpp)
target_link_libraries(TestSimExe PRIVATE TOOLS CORE)
target_compile_options(TestSimExe PRIVATE -O0 -g)
add_test(NAME TestSim
         COMMAND TestSimExe)
file(COPY test_files DESTINATION ${CMAKE_CURRENT_BINARY_DIR})