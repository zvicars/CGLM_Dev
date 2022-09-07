add_executable(TestRandomExe testrandom.cpp)
target_link_libraries(TestRandomExe PRIVATE TOOLS CORE)
target_compile_options(TestRandomExe PRIVATE -O0 -g)
add_test(NAME TestRandom
         COMMAND TestRandomExe)
file(COPY test_files DESTINATION ${CMAKE_CURRENT_BINARY_DIR})