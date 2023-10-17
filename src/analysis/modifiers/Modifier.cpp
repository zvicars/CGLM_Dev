#include "Modifier.hpp"
Modifier::Modifier(AnalysisInputPack& input) : AnalysisObject{input}{
  box = input.getBox();
  return;
}