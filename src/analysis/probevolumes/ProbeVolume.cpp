#include "ProbeVolume.hpp"
#include "../Analysis.hpp"
ProbeVolume::ProbeVolume(AnalysisInputPack& input) : AnalysisObject{input}{
  input.addProbeVolume(name_, this);
  return;
}