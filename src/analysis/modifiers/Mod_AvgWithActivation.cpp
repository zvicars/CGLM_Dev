
#include "Mod_AvgWithActivation.hpp"
#include <cmath>
#include "../../tools/pbcfunctions.hpp"
Mod_AvgWithActivation::Mod_AvgWithActivation(AnalysisInputPack& input):Modifier{input}{
  Vec<real> cutoff; //low and high cutoff
  A = 1.0;
  B = 10.0; 
  x0 = 0.5;
  input.params().readNumber("A", ParameterPack::KeyType::Optional, A);
  input.params().readNumber("B", ParameterPack::KeyType::Optional, B);
  input.params().readNumber("x0", ParameterPack::KeyType::Optional, x0);
  //use 1-sigmoid instead of sigmoid
  invert_ = 0;
  input.params().readFlag("invert", ParameterPack::KeyType::Optional, invert_);
  return;
}
//takes in a lattice, calculates some quantity, overrides the original lattice
void Mod_AvgWithActivation::compute(Matrix3d<real>& lattice){
  auto lattice_temp = lattice;
  #pragma omp parallel for
  for(int i = 0; i < lattice.size_1d(); i++){
    real val = getLocalAverage(lattice.map1N(i), lattice);
    if(invert_) lattice_temp.at_1d(i) = 1.0 - A / ( 1 + std::exp(-B*(val-x0)));
    else lattice_temp.at_1d(i) = A / ( 1 + std::exp(-B*(val-x0)) );
  }
  lattice = lattice_temp;
  return;
}
real Mod_AvgWithActivation::getLocalAverage(Vec3<std::size_t> pos, Matrix3d<real>& lattice){
  real avg = lattice.at(pos);
  auto lattice_size = lattice.size();
  int counter = 1; 
  Vec3<int> int_pos, pos_temp;
  Vec3<std::size_t> u_pos_temp;
  for(int i=0; i<3; i++){
    int_pos[i] = pos[i];
  }
  for(int s=-1; s <=1; s+=2){
    for(int i=0; i<3; i++){
      pos_temp = int_pos;
      pos_temp[i] = wrapIndex(int_pos[i]+s, lattice_size[i]);
      for(int j = 0; j < 3; j++){
        u_pos_temp[j] = pos_temp[j];
      }
      avg += lattice.at(u_pos_temp);
    }
  }
  avg *= (1.0/7.0);
  return avg;
}