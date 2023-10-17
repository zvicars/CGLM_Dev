//lattice object for analyses, basically a skeleton/public version of the one from the simulation code
//modifiers can be used to set additional properties of the box
#pragma once
#include "../tools/InputParser.hpp"
#include "../tools/Stacktrace.hpp"
#include "../typedefs.hpp"
#include "../tools/Matrix.hpp"
class Box{
public:
  Box(){
    return;
  }
  void setLattice(const Vec<char>& states, const Vec3<std::size_t>& size){
    size_ = size;
    lattice_.initialize(size);
    for(int i = 0; i < states.size(); i++){
      lattice_.at_1d(i) = (states[i] == '1');
    }
    //phi is imported first in driver
    if(hasPhi_){
      FANCY_ASSERT(size_ == phi_size_, "phi file is not the same size as the simulation box");
    }
    return;
  }
  void iterateFrame(){
    frame_++;
    return;
  }
  int frame() const{
    return frame_;
  }  
  void setEpsilon(real eps){
    epsilon_ = eps; 
  }
  void setLambda(real lambda){
    lambda_ = lambda;
    real_cell_volume_ = lambda*lambda*lambda;
    real_box_volume_ = real_cell_volume_*size_[0]*size_[1]*size_[2];
    return;
  }
  void setPhi(Vec<real> phi, Vec3<std::size_t> size){
    //if loading phi with a modifier, probably only want to do that once
    if(hasPhi_ == 1) return;
    phi_size_ = size;
    phi_.initialize(size);
    for(int i = 0; i < phi.size(); i++){
      phi_.at_1d(i) = phi[i];
    }
    hasPhi_ = 1;
  }
  Vec3<std::size_t> getSize() const{
    return size_;
  }
  Vec3<int> getSizeInt() const{
    Vec3<int> output = {(int)size_[0], (int)size_[1], (int)size_[2]};
    return output;
  }
  real lambda_, epsilon_, real_cell_volume_, real_box_volume_;
  Vec3<std::size_t> size_, phi_size_;
  Matrix3d<real> lattice_;
  bool hasPhi_ = 0;
  Matrix3d<char> phi_;
  int frame_=0;
};