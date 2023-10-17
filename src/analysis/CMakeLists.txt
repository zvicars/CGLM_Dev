add_executable(analysis Analysis.cpp)
target_link_libraries(analysis PRIVATE TOOLS INSTANT_INT)
target_sources(analysis PRIVATE
 probevolumes/ProbeVolume.cpp
 probevolumes/PV_Boolean.cpp
 probevolumes/PV_Simple.cpp
 calculations/Calculation.cpp
 calculations/Calc_Isosurface.cpp
 calculations/Calc_LeeNP.cpp
 modifiers/Modifier.cpp
 modifiers/Mod_UsePhi.cpp
 modifiers/Mod_AvgWithActivation.cpp
 factory.cpp
 AnalysisInputPack.cpp
)
target_compile_options(analysis PRIVATE -O0)