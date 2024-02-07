add_executable(TestRampingExe testramping.cpp)
target_link_libraries(TestRampingExe PRIVATE TOOLS CORE)
target_compile_options(TestRampingExe PRIVATE -O0 -g)
add_test(NAME TestRampingLG
         COMMAND TestRampingExe)
file(COPY test_files DESTINATION ${CMAKE_CURRENT_BINARY_DIR})