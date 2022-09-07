#include "lattice.hpp"
#include "../../tools/stlmath.hpp"
#include "lattice_pbc.hpp"
Lattice::Lattice(InputPack& input) : Object{input}{
  Vec<std::size_t> size;
  input_->params().readVector("size", ParameterPack::KeyType::Required, size);  
  FANCY_ASSERT(size.size() == 3, "Inappropriate vector size for size in Lattice: " + name_); 
  size_ = vec2Array<std::size_t,3>(size);
  input_->params().readNumber("epsilon", ParameterPack::KeyType::Required, eps_);
  input_->params().readNumber("mu", ParameterPack::KeyType::Required, mu_);
  input_->params().readNumber("lambda", ParameterPack::KeyType::Required, lambda_);  
  input_->params().readNumber("density", ParameterPack::KeyType::Required, density_);  
  FANCY_ASSERT(lambda_ > 0, "lambda (gridspacing in nm) must be greater than 0");
  hasStateField = input_->params().readString("state_file", ParameterPack::KeyType::Optional, sf_);
  hasMuField = input_->params().readString("mu_file", ParameterPack::KeyType::Optional, mff_);
  hasPhiField =  input_->params().readString("phi_file", ParameterPack::KeyType::Optional, pff_);
  if(!hasStateField){
    input_->params().readFlag("default_state", ParameterPack::KeyType::Optional, default_state_);
  }
  //derived classes will handle actual importation of state field file
  std::string rng_name;
  input_->params().readString("generator", ParameterPack::KeyType::Required, rng_name);
  random_generator_ = input_->find(input_->randomgenerators(), rng_name);
  FANCY_ASSERT(random_generator_ != 0, "Failed to find random generator "  + rng_name + " in lattice " + name_ + ".");
  random_generator_ = random_generator_->clone();
  return;
}

int Lattice::getAdjPBC(const Vec3<int>& index) const{
  return getAdj(wrap3(index, size_));
}

int Lattice::getAdj(const Vec3<std::size_t>& index) const{
  Vec3<std::size_t> index_temp;
  //pbc correct the index and assign value to index pbc (reference) and index temp
  #pragma unroll
  for(int i = 0; i < 3; i++){
    index_temp[i] = index[i];
  }
  int n_adj = 0;
  #pragma unroll
  for(int i = 0; i < 3; i++){
    #pragma unroll
    for(int j = -1; j <= 1; j+=2)
    {
      index_temp[i] = wrap((int)index_temp[i]+j, 0, (int)size_[i]-1);
      n_adj += getState(index_temp);
      index_temp[i] = index[i];
    }
  }
  return n_adj;
}

bool Lattice::getStatePBC(const Vec3<int>& index) const{
  return getState(wrap3(index, size_));
}
bool Lattice::flipActive(){
  return flip(getActiveIndex());
}
bool Lattice::flipPBC(const Vec3<int>& index){
  return flip(wrap3(index, size_));
}

Lattice::~Lattice(){
  delete random_generator_;  
  return;
}

void Lattice::reportInfo(string& in){
  Object::reportInfo(in);
  if(hasStateField) in = in + "state field file = " + sf_ + "\n";
  if(hasPhiField) in = in + "phi field file = " + sf_ + "\n";
  if(hasMuField) in = in + "mu field file = " + sf_ + "\n";
  in = in + "default state = " + std::to_string(default_state_) + "\n";
  in = in + "epsilon = " + std::to_string(eps_) + "\n";
  in = in + "mu = " + std::to_string(mu_) + "\n";
  in = in + "density = " + std::to_string(density_) + "\n";  
  in = in + "lambda = " + std::to_string(lambda_) + "\n";
  in = in + "lambda3 = " + std::to_string(lambda3_) + "\n";
  in = in + "size = [ " + std::to_string(size_[0]) + " " + std::to_string(size_[1]) + " " + std::to_string(size_[2]) + " ]\n";
  in = in + "random generator = " + random_generator_->getName() + "\n";
  return;
}