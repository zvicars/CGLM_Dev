#include "hamiltonian_LG.hpp"
#include "../lattice/lattice.hpp"
#include "../bias/bias.hpp"
#include "../pv/probevolume.hpp"
Hamiltonian_LG::Hamiltonian_LG(InputPack& input) : Hamiltonian{input}{
  
  return;
}

Hamiltonian* Hamiltonian_LG::clone(){
  Hamiltonian_LG* ret_ptr = new Hamiltonian_LG(*input_);
  return ret_ptr;
}

real Hamiltonian_LG::calc_dh(const Lattice& lattice){
  real eval = -lattice.eps()*lattice.getActiveAdj() + lattice.getActiveMu() + lattice.getActivePhi();
  if( lattice.getActiveState() ) eval = -eval;
  for(Bias* bs : biases_){
    eval += bs->calc_du(lattice);
  }
  dh_ = eval;
  return eval;
}

real Hamiltonian_LG::calc_h(const Lattice& lattice){
  real eval = 0.0;
  Vec3<std::size_t> size = lattice.size(); 
  Vec3<std::size_t> index;
  for(int i = 0; i < size[0]; i++){
    index[0] = i;
    for(int j = 0; j < size[1]; j++){
      index[1] = j;
      for(int k = 0; k < size[2]; k++){
        index[2] = k;
        eval += (real)lattice.getState(index)*
        ( -0.5*(real)lattice.getAdj(index)*lattice.eps() + lattice.getMu(index) + lattice.getPhi(index) ); //in bounds due to loop
      }
    }
  }
  for(Bias* bs : biases_){
    eval += bs->calc_u(lattice);
  }
  h_ = eval;
  return eval;
}

std::string Hamiltonian_LG::printStepOutput(std::string s){
  if(s == "h") return std::to_string(h_);
  if(s.find("bias") == 0){
    std::size_t index;
    std::string s_pre = s.substr(0, s.find('.'));
    std::string s_post = s.substr(s.find('.')+1);
    int length = s_pre.find(']') - s_pre.find('[')-1;
    FANCY_ASSERT(length > 0, "failed to find index in ouput that requires an index")
    std::string s_step = s.substr(s_pre.find('[')+1, s_pre.find(']') - s_pre.find('[')-1);
    index = std::stoi(s_step);
    FANCY_ASSERT(index < biases_.size(), "invalid index specified");
    return biases_[index]->printStepOutput(s_post);
  }
  if(s.find("pv") == 0){
    std::size_t index;
    std::string s_pre = s.substr(0, s.find('.'));
    std::string s_post = s.substr(s.find('.')+1);
    int length = s_pre.find(']') - s_pre.find('[')-1;
    FANCY_ASSERT(length > 0, "failed to find index in ouput that requires an index")
    std::string s_step = s.substr(s_pre.find('[')+1, s_pre.find(']') - s_pre.find('[')-1);
    index = std::stoi(s_step);
    FANCY_ASSERT(index < pvs_.size(), "invalid index specified");
    return pvs_[index]->printStepOutput(s_post);
  }
  FANCY_ASSERT(0, "invalid output type specified");
  return "";
}