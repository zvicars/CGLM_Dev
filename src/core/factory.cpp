#include "factory.hpp"
#include "simulation.hpp"

#include "lattice/lattice.hpp"
#include "lattice/lattice_1DWrap.hpp"

#include "bias/bias.hpp"
#include "bias/bias_nv.hpp"

#include "hamiltonian/hamiltonian.hpp"
#include "hamiltonian/hamiltonian_LG.hpp"
//#include "hamiltonian/hamiltonian_CGLM.hpp"

#include "pv/probevolume.hpp"
#include "pv/pv_simplerect.hpp"
#include "rng/random.hpp"



Simulation* simulationFactory(InputPack& in){
  std::string type;
  in.params().readString("type", ParameterPack::KeyType::Required, type);
  if(type == "default") return new Simulation(in);
  return 0;
}
Lattice* latticeFactory(InputPack& in){
  std::string type;
  in.params().readString("type", ParameterPack::KeyType::Required, type);
  if(type == "1DWrap") return new Lattice_1DWrap(in);
  FANCY_ASSERT(0, "unknown type given to lattice factory, user provided: " + type + " which is not a recognized type.");
  return 0;
}
Bias* biasFactory(InputPack& in){
  std::string type;
  in.params().readString("type", ParameterPack::KeyType::Required, type);
  if(type == "nv") return new Bias_Nv(in);
  return 0;
}
Hamiltonian* hamiltonianFactory(InputPack& in){
  std::string type;
  in.params().readString("type", ParameterPack::KeyType::Required, type);
  if(type == "lattice_gas") return new Hamiltonian_LG(in);
  //if(type == "cglm") return new Hamiltonian_CGLM(in);
  return 0;
}
ProbeVolume* probevolumeFactory(InputPack& in){
  std::string type;
  in.params().readString("type", ParameterPack::KeyType::Required, type);
  if(type == "simple_rect") return new PV_SimpleRect(in);
  return 0;
}
RNG* randomFactory(InputPack& in){
  std::string type;
  in.params().readString("type", ParameterPack::KeyType::Required, type);
  if(type == "mt19937") return new RNG_mt19937(in);
  return 0;
}
