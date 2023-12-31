
#include "Mod_SetInRegion.hpp"
Mod_SetInRegion::Mod_SetInRegion(AnalysisInputPack& input):Modifier{input}{
  std::vector<real> xrange, yrange, zrange;
  input.params().readNumber("value", ParameterPack::KeyType::Required, value_);
  input.params().readVector("x_range", ParameterPack::KeyType::Required, xrange);
  input.params().readVector("y_range", ParameterPack::KeyType::Required, yrange);
  input.params().readVector("z_range", ParameterPack::KeyType::Required, zrange);
  FANCY_ASSERT(xrange.size() == 2, "xrange should have 2 entries, it has " + std::to_string(xrange.size()));
  FANCY_ASSERT(yrange.size() == 2, "yrange should have 2 entries, it has " + std::to_string(yrange.size()));
  FANCY_ASSERT(zrange.size() == 2, "zrange should have 2 entries, it has " + std::to_string(zrange.size()));
  for(int i = 0; i < 2; i++){
    x_range_[i] = xrange[i];
    y_range_[i] = yrange[i];
    z_range_[i] = zrange[i];
  }
  return;
}
//takes in a lattice, calculates some quantity, overrides the original lattice
void Mod_SetInRegion::compute(Matrix3d<real>& lattice){
  #pragma omp parallel for
  for(int i = 0; i < lattice.size_1d(); i++){
    bool flag = 0;
    auto idx3d = lattice.map1N(i);
    if(x_range_[0] < idx3d[0] && x_range_[1] > idx3d[0] &&  y_range_[0] < idx3d[1] && y_range_[1] > idx3d[1] && z_range_[0] < idx3d[2] && z_range_[1] > idx3d[2]){
        lattice.at_1d(i) = value_;
    }
  }
  return;
}
