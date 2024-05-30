
#include "Mod_LoadPhi.hpp"
#include "../../tools/CGLMFileHelper.hpp"
Mod_LoadPhi::Mod_LoadPhi(AnalysisInputPack& input):Modifier{input}{
  Vec<real> cutoff; //low and high cutoff
  input.params().readVector("cutoff_threshold", ParameterPack::KeyType::Required, cutoff);
  FANCY_ASSERT(cutoff.size() == 2, "This arguments needs two cutoffs, low and high ( e.g. [ -1 100 ] ).");
  cutoff_[0] = cutoff[0];
  cutoff_[1] = cutoff[1];
  //default behavior is anything outside of the bounds becomes 1 and anything inside becomes 0
  //this flag inverts it
  input.params().readFlag("invert", ParameterPack::KeyType::Optional, invert_);
  input.params().readString("file", ParameterPack::KeyType::Required, filename_);
  loadFile(filename_);
  return;
}

void Mod_LoadPhi::loadFile(std::string filename){
  std::vector<real> phi_temp;
  Vec3<std::size_t> size_temp; 
  std::ifstream ifile(filename);
  FANCY_ASSERT(ifile.is_open(), "failed to open phi input file");
  binary_real_read(ifile, phi_temp, size_temp);
  bool success = ifile.fail();
  FANCY_ASSERT(!success, "failure on input phi read");
  phi_internal_.initialize(size_temp);
  for(int i = 0; i < phi_temp.size(); i++){
    phi_internal_.at_1d(i) = phi_temp[i];
  } 
  return;
}

//takes in a lattice, calculates some quantity, overrides the original lattice
void Mod_LoadPhi::compute(Matrix3d<real>& lattice){
  #pragma omp parallel for
  for(int i = 0; i < lattice.size_1d(); i++){
    bool flag = 0;
    real phival = phi_internal_.read_at_1d(i);
    if(phival <= cutoff_[0] || phival >= cutoff_[1]) flag = 1;
    if(invert_) flag = !flag;
    lattice.at_1d(i) = flag;
  }
  return;
}
