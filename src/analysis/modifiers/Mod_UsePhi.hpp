//this modifier includes the phi field
#include "Modifier.hpp"
//
class Mod_UsePhi : public Modifier{
public:
  Mod_UsePhi(AnalysisInputPack& input);
  //takes in a lattice, calculates some quantity, overrides the original lattice
  virtual void compute(Matrix3d<real>& lattice);
protected:
  bool invert_;
  std::array<real, 2> cutoff_;
};