add_library (CORE STATIC "")
target_sources(CORE PRIVATE
 lattice/lattice.cpp
 lattice/lattice_1DWrap.cpp
 pv/probevolume.cpp
 pv/pv_simplerect.cpp
 bias/bias_nv.cpp
 hamiltonian/hamiltonian.cpp
 hamiltonian/hamiltonian_LG.cpp
 rng/mt19937.cpp
 #hamiltonian/hamiltonian_CGLM.cpp
 output/OutputHandler.cpp
 factory.cpp
 simulation.cpp
 inputpack.cpp
)
#target_include_directories(ANALYSIS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../extern/eigen-3.3.9/)
target_link_libraries(CORE PRIVATE TOOLS)
#if(OpenMP_CXX_FOUND)
#target_link_libraries(ANALYSIS PRIVATE GOLOSIO OpenMP::OpenMP_CXX)
#endif()
target_compile_options(CORE PRIVATE -O0)

add_subdirectory(lattice/test)
add_subdirectory(hamiltonian/test)
add_subdirectory(rng/test)
add_subdirectory(test)