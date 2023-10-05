#include "hamiltonian_CGLM.hpp"
#include "../lattice/lattice.hpp"
#include "../lattice/lattice_pbc.hpp"
#include "../../tools/stlmath.hpp"
Hamiltonian_CGLM::Hamiltonian_CGLM(InputPack& input) : Hamiltonian_LG{input}{
  input_->params().readFlag("use_aux", ParameterPack::KeyType::Optional, use_aux_);
  if(use_aux_){
    input_->params().readNumber("aux_type", ParameterPack::KeyType::Required, aux_type_);
  }
  input_->params().readString("chi_file", ParameterPack::KeyType::Required, chi_file_);
  std::ifstream ifile(chi_file_);
  FANCY_ASSERT(ifile.is_open(), "failed to open ifile");
  std::string line;
  while(getline(ifile, line)){
    line = StringTools::trimWhitespace(line);
    if(line.size() == 0) continue;
    if(line.at(0) == '#') continue;
    Vec3<int> offset; 
    real value;
    std::stringstream ss(line);
    ss >> offset[0] >> offset[1] >> offset[2] >> value;
    chi_ij_offsets_.push_back(offset);
    chi_ij_values_.push_back(value);
  }

  std::vector<std::string> ve_names;
  std::string vp_name;
  ProbeVolume* vp;
  input_->params().readVector("excluded_volumes", ParameterPack::KeyType::Optional, ve_names);
  input_->params().readString("probe_volume", ParameterPack::KeyType::Required, vp_name);
  for(auto name : ve_names){
    auto ptr = input_->find(input_->probevolumes(), name); 
    FANCY_ASSERT(ptr != 0, "failed to find probevolume: " + name + ".");
    ve_vec_.push_back(ptr->clone());
  }
  auto ptr = input_->find(input_->probevolumes(), vp_name);
  FANCY_ASSERT(ptr != 0, "failed to find probevolume: " + vp_name + ".");
  vp_ = ptr->clone();
  //should have all of the requisite data
  //build subvolume pairs
  for(int i = 0; i < ve_vec_.size(); i++){
    for(int j = i; j < ve_vec_.size(); j++){
      bool self = j==i;
      SubvolumePair svp1;
      Vec2<ProbeVolume*> pvs = {ve_vec_[i], ve_vec_[j]};
      svp1.pvs = pvs;
      svp1.self = self;
    }
  }

  return; 
}

Hamiltonian* Hamiltonian_CGLM::clone(){
  Hamiltonian_CGLM* ret_ptr = new Hamiltonian_CGLM(*input_);
  return ret_ptr;
}

real Hamiltonian_CGLM::calc_dh(const Lattice& lattice){
  real dhlg = Hamiltonian_LG::calc_dh(lattice);
  return 0.0;
}
real Hamiltonian_CGLM::calc_h(const Lattice& lattice){
  real hlg = Hamiltonian_LG::calc_h(lattice);

  return 0.0;
}

real Hamiltonian_CGLM::computeChiUU(const Lattice& lattice){
  real eval = 0.0;
  auto size = lattice.size();
  for(std::size_t i = 0; i < size[0]; i++){
    for(std::size_t j = 0; j < size[1]; j++){
      for(std::size_t k = 0; k < size[2]; k++){
        Vec3<std::size_t> ref_idx = {i,j,k};
        Vec3<int> ref_site = {i,j,k};
        if(!lattice.getState(ref_idx)) continue;
        for(int l = 0; l < chi_ij_offsets_.size(); i++){
          Vec3<int> new_site = ref_site + chi_ij_offsets_[l];
          auto new_idx = wrap3(new_site, lattice.size());
          if(!lattice.getState(new_idx)) continue;
          eval += lattice.getPhi(ref_idx)*lattice.getPhi(new_idx)*chi_ij_values_[l]; 
        }
      }
    }
  }
  return eval;
}

real computeChiPU(ProbeVolume* pv, const Lattice& lattice){
  real eval = 0.0;
  auto size = lattice.size();
  for(std::size_t i = 0; i < size[0]; i++){
    for(std::size_t j = 0; j < size[1]; j++){
      for(std::size_t k = 0; k < size[2]; k++){
        Vec3<std::size_t> ref_idx = {i,j,k};
        Vec3<int> ref_site = {i,j,k};
        if(!lattice.getState(ref_idx)) continue;
        for(int l = 0; l < chi_ij_offsets_.size(); i++){
          Vec3<int> new_site = ref_site + chi_ij_offsets_[l];
          auto new_idx = wrap3(new_site, lattice.size());
          if(!lattice.getState(new_idx)) continue;
          eval += lattice.getPhi(ref_idx)*lattice.getPhi(new_idx)*chi_ij_values_[l]; 
        }
      }
    }
  }
  return eval;
}
