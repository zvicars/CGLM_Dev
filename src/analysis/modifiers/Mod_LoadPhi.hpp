//this modifier includes the phi field
#include "Modifier.hpp"
//
class Mod_LoadPhi : public Modifier{
public:
  Mod_LoadPhi(AnalysisInputPack& input);
  //takes in a lattice, calculates some quantity, overrides the original lattice
  virtual void compute(Matrix3d<real>& lattice);
protected:
  void loadFile(std::string filename);
  bool invert_=0;
  bool is_initialized_=0;
  std::array<real, 2> cutoff_;
  std::string filename_;
  Matrix3d<real> phi_internal_;
};