#include "factory.hpp"
#include "../tools/Assert.hpp"
#include "probevolumes/PV_Boolean.hpp"
#include "probevolumes/PV_Simple.hpp"
#include "calculations/Calc_Isosurface.hpp"
#include "calculations/Calc_IsosurfaceMultiphase.hpp"
#include "calculations/Calc_LeeNP.hpp"
#include "calculations/Calc_PhiFromLiqDens.hpp"
#include "calculations/Calc_WriteTraj.hpp"
#include "calculations/Calc_PoreFilling.hpp"
#include "modifiers/Mod_UsePhi.hpp"
#include "modifiers/Mod_AvgWithActivation.hpp"
#include "modifiers/Mod_SetInRegion.hpp"
#include "modifiers/Mod_COMCorrectPillar.hpp"
#include "modifiers/Modifier.hpp"
#include "modifiers/Mod_LoadPhi.hpp"
#include "modifiers/Mod_Gaussian.hpp"

ProbeVolume* ProbeVolume_Factory(std::string key, AnalysisInputPack& input){
  if(key == "rectilinear") return new PV_DiscreteRect(input);
  if(key == "boolean") return new PV_Boolean(input);
  FANCY_ASSERT(0, "Failed to find matching case for key: " + key);
  return 0;
}

Modifier* Modifier_Factory(std::string key, AnalysisInputPack& input){
  if(key == "average") return new Mod_AvgWithActivation(input);
  if(key == "usephi") return new Mod_UsePhi(input);
  if(key == "loadphi") return new Mod_LoadPhi(input);
  if(key == "setregion") return new Mod_SetInRegion(input);
  if(key == "compillar") return new Mod_COMCorrectPillar(input);
  if(key == "gaussian") return new Mod_Gaussian(input);
  FANCY_ASSERT(0, "Failed to find matching case for key: " + key);
  return 0;
}

Calculation* Calculation_Factory(std::string key, AnalysisInputPack& input){
  if(key == "isosurface") return new Calc_Isosurface(input);
  if(key == "lee_np") return new Calc_LeeNP(input);
  if(key == "phi_from_liq") return new Calc_PhiFromLiqDens(input);
  if(key == "write_traj") return new Calc_WriteTraj(input);
  if(key == "isosurface_multiphase") return new Calc_IsosurfaceMultiphase(input);
  if(key == "porefill") return new Calc_PoreFilling(input);
  FANCY_ASSERT(0, "Failed to find matching case for key: " + key);
  return 0;
}