//Zachariah Vicars
//The input pack contains all of the necessary information to run a simulation in addition to
//pointers to relevant data (probevolumes, biases, etc) that are defined in the input file
//just a convenient way of lumping things together to reduce complexity for constructors
//it's also important to note that all of these objects, except for the simulations, are simply templates which
//are used to produce the full objects in the constructors
#pragma once
#include "../typedefs.hpp"
#include "../tools/InputParser.hpp"
#include <map>
#include <string>
//forward declarations of pointer types the input pack contains
class Lattice;
class Bias;
class ProbeVolume;
class Hamiltonian;
class RNG;
class Simulation;

class InputPack{
  public:
  InputPack(){return;}
  InputPack(const ParameterPack& params){
    set(params);
    bool isMasterPack = 1; //if it's a master pack, it will delete the objects pointed to
    lattices_ = new strmap<Lattice*>;
    hamiltonians_ = new strmap<Hamiltonian*>;
    biases_ = new strmap<Bias*>;
    probevolumes_ = new strmap<ProbeVolume*>;
    randomgenerators_ = new strmap<RNG*>;
    simulations_ = new strmap<Simulation*>;
    return;
  }
  ~InputPack();
  Lattice* findLattice(std::string name);
  Hamiltonian* findHamiltonian(std::string name);
  Bias* findBias(std::string name);
  ProbeVolume* findProbeVolume(std::string name);

  bool add(std::string name, Lattice* lattice);
  bool add(std::string name, Hamiltonian* hamiltonian);
  bool add(std::string name, Bias* bias);
  bool add(std::string name, ProbeVolume* probevolume);
  bool add(std::string name, RNG* randomgenerator);
  bool add(std::string name, Simulation* simulation);

  bool set(strmap<Lattice*>* map);
  bool set(strmap<Hamiltonian*>* map);
  bool set(strmap<Bias*>* map);
  bool set(strmap<ProbeVolume*>* map);
  bool set(strmap<RNG*>* map);
  bool set(strmap<Simulation*>* map);
  bool set(const ParameterPack& params);

  const ParameterPack& params(){
    return params_;
  }
  const strmap<Lattice*>* lattices(){
    return lattices_;
  }
  const strmap<Hamiltonian*>* hamiltonians(){
    return hamiltonians_;
  }
  const strmap<Bias*>* biases(){
    return biases_;
  }
  const strmap<ProbeVolume*>* probevolumes(){
    return probevolumes_;
  }
  const strmap<RNG*>* randomgenerators(){
    return randomgenerators_;
  }
  const strmap<Simulation*>* simulations(){
    return simulations_;
  }

  template <class T>
  T find(const strmap<T>* map, string name) const{
    auto it = map->find(name);
    if(it != map->end()) return it->second;
    return 0; //return nullptr if search fails    
  }

  std::vector<InputPack> buildDerivedInputPacks(std::string key){
    std::vector<InputPack> inputpacks;
    auto parameterpacks = params().findParameterPacks(key, ParameterPack::KeyType::Optional);
    inputpacks.resize(parameterpacks.size());
    for(std::size_t i = 0; i < parameterpacks.size(); i++){
      inputpacks[i].set(lattices_);
      inputpacks[i].set(hamiltonians_);
      inputpacks[i].set(biases_);
      inputpacks[i].set(probevolumes_);
      inputpacks[i].set(randomgenerators_);
      inputpacks[i].set(simulations_);
      inputpacks[i].set(*parameterpacks[i]);
    }
    return inputpacks;
  }

  private:
  bool isMasterPack = 0;
  strmap<Lattice*>* lattices_;
  strmap<Hamiltonian*>* hamiltonians_;
  strmap<Bias*>* biases_;
  strmap<ProbeVolume*>* probevolumes_;
  strmap<RNG*>* randomgenerators_;
  strmap<Simulation*>* simulations_;
  ParameterPack params_;
};