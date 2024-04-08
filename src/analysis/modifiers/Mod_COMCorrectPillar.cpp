
#include "Mod_COMCorrectPillar.hpp"
#include <cmath>
#include "../../tools/pbcfunctions.hpp"
#include <limits>
//tries to identify the offset from the set of provided offsets that's closest to having the COM of the box at a given position
Mod_COMCorrectPillar::Mod_COMCorrectPillar(AnalysisInputPack& input):Modifier{input}{
  auto position_vec = input.params().findVectors("position", ParameterPack::KeyType::Optional);
  Vec<Vec3<int> > real_pos;
  for(int i=0; i < position_vec.size(); i++){
    auto entry = *position_vec[i];
    FANCY_ASSERT(entry.size() == 3, "Invalid position entry, needs 3 values");
    Vec3<int> real_entry;
    for(int j=0; j<3; j++){
      real_entry[j] = std::stoi(entry[j]);
    }
    real_pos.push_back(real_entry);
  }
  test_positions_ = real_pos;
  Vec<int> ref_pos; 
  input.params().readVector("reference", ParameterPack::KeyType::Optional, ref_pos);
  FANCY_ASSERT(ref_pos.size() == 3, "reference position needs 3 entries");
  for(int i=0; i<3; i++){
    ref_pos_[i] = ref_pos[i];
  }
  return;
}
//takes in a lattice, calculates some quantity, overrides the original lattice
void Mod_COMCorrectPillar::compute(Matrix3d<real>& lattice){
  auto lattice_temp = lattice;
  //offset to a test point
  real min_rss_val = std::numeric_limits<real>::max(); 
  Vec3<int> best_offset = {0,0,0};
  for(auto test_pos : test_positions_){
    offset(lattice_temp, test_pos);
    real rss_val = compute_rss(lattice_temp, ref_pos_);
    if(rss_val < min_rss_val){
      best_offset = test_pos;
      min_rss_val = rss_val;
    }
    lattice_temp = lattice;
  }
  offset(lattice_temp, best_offset);
  lattice = lattice_temp;
  return;
}

void Mod_COMCorrectPillar::offset(Matrix3d<real>& lattice, Vec3<int> offset){
  auto lattice_temp = lattice;
  int numel = lattice_temp.size_1d();
  Vec3<std::size_t> box_size = lattice_temp.size();
  for(int i = 0; i < numel; i++){
    auto pos3d = lattice_temp.map1N(i);
    Vec3<std::size_t> pos3d_offset;
    for(int j = 0; j < 3; j++){
      pos3d_offset[j] = wrapIndex(pos3d[j] + offset[j], box_size[j]);
    }
    lattice.at(pos3d_offset) = lattice_temp.at(pos3d);
  }
  return;
}

real Mod_COMCorrectPillar::compute_rss(Matrix3d<real>& lattice, Vec3<int> ref_pos){
  int numel = lattice.size_1d();
  Vec3<std::size_t> box_size = lattice.size();
  Vec3<real> ref_pos_real, box_size_real;
  for(int i = 0; i < 3; i++){
    ref_pos_real[i] = (real)ref_pos[i]; 
    box_size_real[i] = box_size[i];
  }

  real rss_total = 0.0;
  for(int i = 0; i < numel; i++){
    auto pos3d = lattice.map1N(i);
    Vec3<real> real_pos3d;
    for(int j = 0; j < 3; j++){
      real_pos3d[j] = (real)pos3d[j];
    }
    real eval = lattice.at_1d(i) * getDistance(real_pos3d, ref_pos_real, box_size_real);
    rss_total += eval*eval;
  }
  rss_total = rss_total / (real)numel;
  return rss_total;
}