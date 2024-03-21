//performs a local average and then maps it from 0->1 depending on threshold value and sigmoidal function
// A /( 1+exp(-B(x-x0)) )
#pragma once 
#include "Modifier.hpp"
class Mod_COMCorrectPillar : public Modifier{
public:
  Mod_COMCorrectPillar(AnalysisInputPack& input);
  //takes in a lattice, calculates some quantity, overrides the original lattice
  virtual void compute(Matrix3d<real>& lattice);
protected:
  void offset(Matrix3d<real>& lattice, Vec3<int> offset);
  real compute_rss(Matrix3d<real>& lattice, Vec3<int> ref_pos);
  Vec<Vec3<int> > test_positions_;
  Vec<real> rss_values_;
  Vec3<int> ref_pos_;
};