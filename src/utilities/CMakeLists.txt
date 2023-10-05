add_executable(b2xyz b2xyz.cpp)
target_link_libraries(b2xyz PRIVATE TOOLS)
target_compile_options(b2xyz PRIVATE -O3)

add_executable(phi2xyz phi2xyz.cpp)
target_link_libraries(phi2xyz PRIVATE TOOLS)
target_compile_options(phi2xyz PRIVATE -O3)

add_subdirectory(genphilib)
add_executable(genphi genphi.cpp)
target_link_libraries(genphi PRIVATE TOOLS GENPHILIB)
target_compile_options(genphi PRIVATE -O3)

add_library(PHIFILELIB STATIC
  SimpleGroReader.cpp
)
target_link_libraries(PHIFILELIB PRIVATE TOOLS)

add_executable(gro2phifile gro2phifile.cpp)
target_link_libraries(gro2phifile PRIVATE PHIFILELIB)
target_compile_options(gro2phifile PRIVATE -O3)
