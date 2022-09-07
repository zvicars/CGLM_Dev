add_library (GENPHILIB STATIC
  ff_functions.cpp
  filehandler.cpp
)
target_include_directories(GENPHILIB PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(GENPHILIB PRIVATE -fPIC)
target_link_libraries(GENPHILIB PRIVATE TOOLS CORE)