#include "pv_simplerect.hpp"
#include "../../tools/stlmath.hpp"
#include "../lattice/lattice.hpp"
PV_SimpleRect::PV_SimpleRect(InputPack& input) : ProbeVolume{input}{
  std::vector<int> xmin, xmax;
  input_->params().readVector("min", ParameterPack::KeyType::Required, xmin);  
  input_->params().readVector("max", ParameterPack::KeyType::Required, xmax);  
  FANCY_ASSERT(xmin.size() == 3, "Inappropriate vector size for xmin in PV_SimpleRect: " + name_);
  FANCY_ASSERT(xmax.size() == 3, "Inappropriate vector size for xmax in PV_SimpleRect: " + name_);
  min_ = vec2Array<int,3>(xmin);
  max_ = vec2Array<int,3>(xmax);
  return;
}

real PV_SimpleRect::calc_nv(const Lattice& lattice){
  stored = &lattice;
  int n = 0;
  Vec3<std::size_t> index;
  auto size = lattice.size();
  for(int i = 0; i < size[0]; i++){
    index[0] = i;
    if( i < min_[0] || i > max_[0]) continue;
    for(int j = 0; j < size[1]; j++){
      index[1] = j;
      if( j < min_[1] || j > max_[1]) continue;
      for(int k = 0; k < size[2]; k++){
        index[2] = k;
        if( k < min_[2] || k > max_[2]) continue;
        n += lattice.getState(index); //based on indexing, guaranteed to be in-bounds
      }
    }
  }
  //biasing the real number of waters
  nv_ = n;
  return nv_;
}

real PV_SimpleRect::calc_nv_compute(){
  auto& lattice = *stored;
  int n = 0;
  Vec3<std::size_t> index;
  auto size = lattice.size();
  for(int i = 0; i < size[0]; i++){
    index[0] = i;
    if( i < min_[0] || i > max_[0]) continue;
    for(int j = 0; j < size[1]; j++){
      index[1] = j;
      if( j < min_[1] || j > max_[1]) continue;
      for(int k = 0; k < size[2]; k++){
        index[2] = k;
        if( k < min_[2] || k > max_[2]) continue;
        n += lattice.getState(index); //based on indexing, guaranteed to be in-bounds
      }
    }
  }
  return n;
}

real PV_SimpleRect::calc_dnv(const Lattice& lattice){
  dnv_ = 0.0;
  auto idx = lattice.getActiveIndex();
  for(int i = 0; i < 3; i++){
    int ix = idx[i], minval = min_[i], maxval = max_[i];
    if(ix < minval || ix > maxval){
      dnv_ = 0.0;
      return 0.0;
    }
  }
  if(!lattice.getActiveState()) dnv_ = 1.0;
  else dnv_ = -1.0;
  return dnv_;
}

ProbeVolume* PV_SimpleRect::clone(){
  PV_SimpleRect* ptr = new PV_SimpleRect(*input_);
  return ptr;
}

std::string PV_SimpleRect::printStepOutput(std::string s){
  if(s == "nv") return std::to_string(nv_);
  if(s == "nv_compute") return std::to_string(calc_nv_compute());
  FANCY_ASSERT(0, "invalid output type specified");
  return "";
}

void PV_SimpleRect::bounds(Vec<std::size_t>& vec_out){
  vec_out.clear();
  #pragma unroll
  for(int i = 0; i < 3; i++){
    vec_out[i] = min_[i];
    vec_out[i+3] = max_[i];
  }
  return;
}

bool PV_SimpleRect::isInside(Vec3<std::size_t> idx){
  #pragma unroll
  for(int i = 0; i < 3; i++){
    if(idx[i] < min_[i]) return 0;
    if(idx[i] > max_[i]) return 0;
  }
  return 1;
}