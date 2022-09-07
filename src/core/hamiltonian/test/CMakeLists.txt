add_executable(TestLGExe testlg.cpp)
target_link_libraries(TestLGExe PRIVATE TOOLS CORE)
target_compile_options(TestLGExe PRIVATE -O0 -g)
add_test(NAME TestLG
         COMMAND TestLGExe)
file(COPY test_files DESTINATION ${CMAKE_CURRENT_BINARY_DIR})