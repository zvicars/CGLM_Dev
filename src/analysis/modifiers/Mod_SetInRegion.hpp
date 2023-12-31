//this modifier includes the phi field
#include "Modifier.hpp"
//
class Mod_SetInRegion : public Modifier{
public:
  Mod_SetInRegion(AnalysisInputPack& input);
  //takes in a lattice, calculates some quantity, overrides the original lattice
  virtual void compute(Matrix3d<real>& lattice);
protected:
  std::array<real, 2> x_range_, y_range_, z_range_;
  real value_;
};