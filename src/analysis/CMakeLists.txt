add_executable(analysis Analysis.cpp)
target_link_libraries(analysis PRIVATE TOOLS INSTANT_INT GENPHILIB)
target_sources(analysis PRIVATE
 probevolumes/ProbeVolume.cpp
 probevolumes/PV_Boolean.cpp
 probevolumes/PV_Simple.cpp
 probevolumes/PV_Sphere.cpp 
 calculations/Calculation.cpp
 calculations/Calc_Isosurface.cpp
 calculations/Calc_IsosurfaceMultiphase.cpp
 calculations/Calc_LeeNP.cpp
 calculations/Calc_PhiFromLiqDens.cpp
 calculations/Calc_WriteTraj.cpp
 calculations/Calc_PoreFilling.cpp
 modifiers/Modifier.cpp
 modifiers/Mod_UsePhi.cpp
 modifiers/Mod_AvgWithActivation.cpp
 modifiers/Mod_SetInRegion.cpp
 modifiers/Mod_COMCorrectPillar.cpp
 modifiers/Mod_LoadPhi.cpp
 modifiers/Mod_Gaussian.cpp
 factory.cpp
 AnalysisInputPack.cpp
)
target_compile_options(analysis PRIVATE -O3 -g -fopenmp)