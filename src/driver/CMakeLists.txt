add_executable(CGLM CGLM.cpp)
target_link_libraries(CGLM PRIVATE TOOLS CORE)
target_compile_options(CGLM PRIVATE -O3 -g)