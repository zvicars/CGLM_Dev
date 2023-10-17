#pragma once
#include "probevolumes/ProbeVolume.hpp"
#include "calculations/Calculation.hpp"
#include "modifiers/Modifier.hpp"

ProbeVolume* ProbeVolume_Factory(std::string key, AnalysisInputPack& input);
Modifier* Modifier_Factory(std::string key, AnalysisInputPack& input);
Calculation* Calculation_Factory(std::string key, AnalysisInputPack& input);