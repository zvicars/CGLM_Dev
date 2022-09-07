#include "hamiltonian_CGLM.hpp"
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