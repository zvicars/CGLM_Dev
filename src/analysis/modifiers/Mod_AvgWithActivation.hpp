//performs a local average and then maps it from 0->1 depending on threshold value and sigmoidal function
// A /( 1+exp(-B(x-x0)) )
#pragma once 
#include "Modifier.hpp"
class Mod_AvgWithActivation : public Modifier{
public:
  Mod_AvgWithActivation(AnalysisInputPack& input);
  //takes in a lattice, calculates some quantity, overrides the original lattice
  virtual void compute(Matrix3d<real>& lattice);
protected:
  real getLocalAverage(Vec3<std::size_t> index, Matrix3d<real>& lattice);
  bool invert_;
  real x0, A, B;
};