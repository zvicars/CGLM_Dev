#include "inputpack.hpp"
#include "pv/probevolume.hpp"
#include "lattice/lattice.hpp"
#include "simulation.hpp"
#include "rng/random.hpp"
#include "hamiltonian/hamiltonian.hpp"
#include "bias/bias.hpp"

Lattice* InputPack::findLattice(std::string name){
  auto& map = lattices_; 
  auto it = map->find(name);
  if(it != map->end()) return it->second;
  return 0; //return nullptr if search fails  
}
Hamiltonian* InputPack::findHamiltonian(std::string name){
  auto& map = hamiltonians_; 
  auto it = map->find(name);
  if(it != map->end()) return it->second;
  return 0; //return nullptr if search fails    
}
Bias* InputPack::findBias(std::string name){
   auto& map = biases_; 
  auto it = map->find(name);
  if(it != map->end()) return it->second;
  return 0; //return nullptr if search fails   
}
ProbeVolume* InputPack::findProbeVolume(std::string name){
  auto& map = probevolumes_;
  auto it = map->find(name);
  if(it != map->end()) return it->second;
  return 0; //return nullptr if search fails    
}


bool InputPack::add(std::string name, Lattice* lattice){
  lattices_->insert(std::pair<std::string, Lattice*>{name, lattice});
  return 0;
}
bool InputPack::add(std::string name, Hamiltonian* hamiltonian){
  hamiltonians_->insert(std::pair<std::string, Hamiltonian*>{name, hamiltonian});
  return 0;
}
bool InputPack::add(std::string name, Bias* bias){
  biases_->insert(std::pair<std::string, Bias*>{name, bias});
  return 0;
}
bool InputPack::add(std::string name, ProbeVolume* probevolume){
  probevolumes_->insert(std::pair<std::string, ProbeVolume*>{name, probevolume});
  return 0;
}
bool InputPack::add(std::string name, RNG* randomgenerator){
  randomgenerators_->insert(std::pair<std::string, RNG*>{name, randomgenerator});
  return 0;
}
bool InputPack::add(std::string name, Simulation* simulation){
  simulations_->insert(std::pair<std::string, Simulation*>{name, simulation});
  return 0;
}

bool InputPack::set(strmap<Lattice*>* map){
  lattices_ = map;
  return 0;
}
bool InputPack::set(strmap<Hamiltonian*>* map){
  hamiltonians_ = map;
  return 0;
}
bool InputPack::set(strmap<Bias*>* map){
  biases_ = map;
  return 0;
}
bool InputPack::set(strmap<ProbeVolume*>* map){
  probevolumes_ = map;
  return 0;
}
bool InputPack::set(strmap<RNG*>* map){
  randomgenerators_ = map;
  return 0;
}
bool InputPack::set(strmap<Simulation*>* map){
  simulations_ = map;
  return 0;
}
bool InputPack::set(const ParameterPack& params){
  params_ = params;
  return 0;
}

InputPack::~InputPack(){
  //if the object was created using the parameter pack, then it will be considered the master input pack
  if(isMasterPack){
    for(auto& ptr : *lattices_) delete ptr.second;
    delete lattices_;
    for(auto& ptr : *hamiltonians_) delete ptr.second;
    delete hamiltonians_;
    for(auto& ptr : *biases_) delete ptr.second;
    delete biases_;
    for(auto& ptr : *probevolumes_) delete ptr.second;
    delete probevolumes_;
    for(auto& ptr : *randomgenerators_) delete ptr.second;
    delete randomgenerators_;
    for(auto& ptr : *simulations_) delete ptr.second;
    delete simulations_;
  }
  return;
}