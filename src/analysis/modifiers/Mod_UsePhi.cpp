
#include "Mod_UsePhi.hpp"
Mod_UsePhi::Mod_UsePhi(AnalysisInputPack& input):Modifier{input}{
  Vec<real> cutoff; //low and high cutoff
  input.params().readVector("cutoff_threshold", ParameterPack::KeyType::Required, cutoff);
  FANCY_ASSERT(cutoff.size() == 2, "This arguments needs two cutoffs, low and high ( e.g. [ -1 100 ] ).");
  cutoff_[0] = cutoff[0];
  cutoff_[1] = cutoff[1];
  //default behavior is anything outside of the bounds becomes 1 and anything inside becomes 0
  //this flag inverts it
  input.params().readFlag("invert", ParameterPack::KeyType::Optional, invert_);
  return;
}
//takes in a lattice, calculates some quantity, overrides the original lattice
void Mod_UsePhi::compute(Matrix3d<real>& lattice){
  #pragma omp parallel for
  for(int i = 0; i < lattice.size_1d(); i++){
    bool flag = 0;
    real phival = box->phi_.read_at_1d(i);
    if(phival <= cutoff_[0] || phival >= cutoff_[1]) flag = 1;
    if(invert_) flag = !flag;
    lattice.at_1d(i) = flag;
  }
  return;
}
