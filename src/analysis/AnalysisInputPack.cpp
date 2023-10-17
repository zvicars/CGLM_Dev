#include "AnalysisInputPack.hpp"
#include "probevolumes/ProbeVolume.hpp"
#include "modifiers/Modifier.hpp"
#include "calculations/Calculation.hpp"
AnalysisInputPack::~AnalysisInputPack()
{
  if(is_master_pack){ //only a pack that initialized the registries should delete them
    for(auto it = pv_registry_->begin(); it!=pv_registry_->end();) {
      delete it->second;
      it = pv_registry_->erase(it);
    }
    for(auto it = calc_registry_->begin(); it!=calc_registry_->end();) {
      delete it->second;
      it = calc_registry_->erase(it);
    }
    for(auto it = mod_registry_->begin(); it!=mod_registry_->end();) {
      delete it->second;
      it = mod_registry_->erase(it);
    }
    delete pv_registry_;
    delete calc_registry_;
    delete mod_registry_;
  }
}