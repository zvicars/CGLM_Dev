
#include "Mod_Gaussian.hpp"
#include <cmath>
#include "../../tools/pbcfunctions.hpp"
#include "../../tools/smearfunctions.hpp"
#include "../../tools/stlmath.hpp"
Mod_Gaussian::Mod_Gaussian(AnalysisInputPack& input):Modifier{input}{
  Vec<real> cutoff; //low and high cutoff
  sigma_ = 1.0;
  cutoff_ = 2.0;
  span_ = 2;
  input.params().readNumber("sigma", ParameterPack::KeyType::Required, sigma_);
  cutoff_ = 2.0*sigma_;
  input.params().readNumber("cutoff", ParameterPack::KeyType::Required, cutoff_);
  span_ = ceil(cutoff_);
  return;
}
//takes in a lattice, calculates some quantity, overrides the original lattice
void Mod_Gaussian::compute(Matrix3d<real>& lattice){
  auto lattice_temp = lattice;
  #pragma omp parallel for
  for(int i = 0; i < lattice.size_1d(); i++){
    lattice_temp.at_1d(i) = getLocalCG(lattice.map1N(i), lattice);
  }
  lattice = lattice_temp;
  return;
}
real Mod_Gaussian::getLocalCG(Vec3<std::size_t> pos, Matrix3d<real>& lattice){
  //this is probably faster with a 3d FFT, but sometimes I'm lazy
  auto lattice_size = lattice.size();
  int counter = 1; 
  Vec3<int> int_pos, pos_temp, pos_temp_nowrap;
  for(int i=0; i<3; i++){
    int_pos[i] = pos[i];
  }
  double gaussian_sum=0.0;
  for(int ix = int_pos[0] - span_; ix <= int_pos[0] + span_; ix++){
    for(int iy = int_pos[1] - span_; iy <= int_pos[1] + span_; iy++){
      for(int iz = int_pos[2] - span_; iz <= int_pos[2] + span_; iz++){
        pos_temp = {wrapIndex(ix, lattice_size[0]), wrapIndex(iy, lattice_size[1]), wrapIndex(iz, lattice_size[2])};
        pos_temp_nowrap = {ix, iy, iz};
        real state_temp = lattice.at(pos_temp);
        Vec3<int> x_offset = pos_temp_nowrap-int_pos;
        gaussian_sum += state_temp * h_x(x_offset[0], 0, 1, sigma_, cutoff_)*h_x(x_offset[1], 0, 1, sigma_, cutoff_)*h_x(x_offset[2], 0, 1, sigma_, cutoff_);
      }
    }
  }
  return gaussian_sum;
}