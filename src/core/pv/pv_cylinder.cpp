#include "pv_cylinder.hpp"
#include "../../tools/stlmath.hpp"
#include "../lattice/lattice.hpp"
PV_Cylinder::PV_Cylinder(InputPack& input) : ProbeVolume{input}{
  std::vector<int> range;
  std::vector<real> origin;
  input_->params().readNumber("height", ParameterPack::KeyType::Required, h_);
  input_->params().readVector("origin", ParameterPack::KeyType::Required, origin);
  input_->params().readNumber("radius", ParameterPack::KeyType::Required, r_);
  input_->params().readNumber("axis", ParameterPack::KeyType::Required, axis_);
  input_->params().readFlag("invert", ParameterPack::KeyType::Optional, invert_);
  FANCY_ASSERT(origin.size() == 3, "Inappropriate vector size for origin in PV_Cylinder: " + name_);
  FANCY_ASSERT( axis_ == 0 || axis_ == 1 || axis_ == 2, "Invalid axis specified in PV_Cylinder: " + name_);
  origin_ = vec2Array<real,3>(origin);
  r2_ = r_*r_;
  return;
}

real PV_Cylinder::calc_nv(const Lattice& lattice){
  stored = &lattice;
  int n = 0;
  Vec3<std::size_t> index;
  auto size = lattice.size();
  for(int i = 0; i < size[0]; i++){
    index[0] = i;
    for(int j = 0; j < size[1]; j++){
      index[1] = j;
      for(int k = 0; k < size[2]; k++){
        index[2] = k;
        if(isInside(index)) n += lattice.getState(index);
      }
    }
  }
  //biasing the real number of waters
  nv_ = n;
  return nv_;
}

real PV_Cylinder::calc_nv_compute(){
  auto& lattice = *stored;
  int n = 0;
  Vec3<std::size_t> index;
  auto size = lattice.size();
  for(int i = 0; i < size[0]; i++){
    index[0] = i;
    for(int j = 0; j < size[1]; j++){
      index[1] = j;
      for(int k = 0; k < size[2]; k++){
        index[2] = k;
        if(isInside(index)) n += lattice.getState(index);
      }
    }
  }
  return n;
}

real PV_Cylinder::calc_dnv(const Lattice& lattice){
  dnv_ = 0.0;
  auto idx = lattice.getActiveIndex();
  if(isInside(idx)){
    if(!lattice.getActiveState()) dnv_ = 1.0;
    else dnv_ = -1.0;
  }
  return dnv_;
}

ProbeVolume* PV_Cylinder::clone(){
  PV_Cylinder* ptr = new PV_Cylinder(*input_);
  return ptr;
}

std::string PV_Cylinder::printStepOutput(std::string s){
  if(s == "nv") return std::to_string(nv_);
  if(s == "nv_compute") return std::to_string(calc_nv_compute());
  FANCY_ASSERT(0, "invalid output type specified");
  return "";
}

//returns cuboid bounds of lattice sites that could possibly be in this volume
void PV_Cylinder::bounds(Vec<std::size_t>& vec_out){
  if(stored == 0){
    std::cerr << "Requested bounds in pv_cylinder: " + name_ + ", but the observation volume currently has no pointer to the box." << std::endl;
    return;
  }
  vec_out.clear();
  for(int i = 0; i < 3; i++){
    if(axis_ == i){
      vec_out[i] = std::max((int)std::floor(origin_[axis_]), 0);
      vec_out[i+3] = std::min((std::size_t)std::ceil(origin_[axis_] + h_), stored->size()[i]-1);
    }
    else{
      vec_out[i] = std::min((int)std::floor(origin_[i] - r_), 0);
      vec_out[i+3] = std::min((std::size_t)std::ceil(origin_[i] + r_), stored->size()[i]-1);
    }
  }
  return;
}

bool PV_Cylinder::isInside(Vec3<std::size_t> idx){
  //convert index to center-of-voxel coordinates
  Vec3<real> idx_real;
  for(int i=0; i<3; i++){
    idx_real[i] = idx[i] + 0.5;
  }

  if(idx[axis_] < origin_[axis_]) return 0;
  if(idx[axis_] > origin_[axis_] + h_) return 0;
  
  real nrm = 0.0;
  for(int i = 0; i < 3; i++){
    if(i == axis_) continue;
    nrm +=  (origin_[i] - idx_real[i]) * (origin_[i] - idx_real[i]);
  }
  bool retVal = nrm < r2_;
  if(invert_) retVal = !retVal;
  return retVal;
}