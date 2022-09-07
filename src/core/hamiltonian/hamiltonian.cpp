#include "hamiltonian.hpp"
#include "../pv/probevolume.hpp"
#include "../bias/bias.hpp"

Hamiltonian::Hamiltonian(InputPack& input) : Object{input}{
  std::vector<std::string> bias_names;
  std::vector<Bias*> bias_ptrs;
  input_->params().readVector("biases", ParameterPack::KeyType::Optional, bias_names);
  for(auto name : bias_names){
    auto ptr = input_->find(input_->biases(), name); 
    FANCY_ASSERT(ptr != 0, "Failed to find bias: " + name + ".");
    bias_ptrs.push_back(ptr);
  } 
  for(auto ptr : bias_ptrs){
    biases_.push_back(ptr->clone());
  }
  return;
}
Hamiltonian::~Hamiltonian(){
  for(auto bias : biases_){
    delete bias;
  }
  biases_.clear();
  return;
}