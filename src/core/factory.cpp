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
#include "pv/pv_cylinder.hpp"

#include "rng/random.hpp"



Simulation* simulationFactory(InputPack& in){
  std::string type = "default";
  in.params().readString("type", ParameterPack::KeyType::Optional, type);
  if(type == "default") return new Simulation(in);
  FANCY_ASSERT(0, "unknown type given to simulation factory, user provided: " + type + " which is not a recognized type.");
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
  FANCY_ASSERT(0, "unknown type given to bias factory, user provided: " + type + " which is not a recognized type.");
  return 0;
}
Hamiltonian* hamiltonianFactory(InputPack& in){
  std::string type;
  in.params().readString("type", ParameterPack::KeyType::Required, type);
  if(type == "lattice_gas") return new Hamiltonian_LG(in);
  //if(type == "cglm") return new Hamiltonian_CGLM(in);
  FANCY_ASSERT(0, "unknown type given to hamiltonian factory, user provided: " + type + " which is not a recognized type.");
  return 0;
}
ProbeVolume* probevolumeFactory(InputPack& in){
  std::string type;
  in.params().readString("type", ParameterPack::KeyType::Required, type);
  if(type == "simple_rect") return new PV_SimpleRect(in);
  if(type == "cylinder") return new PV_Cylinder(in);
  FANCY_ASSERT(0, "unknown type given to pv factory, user provided: " + type + " which is not a recognized type.");
  return 0;
}
RNG* randomFactory(InputPack& in){
  std::string type;
  in.params().readString("type", ParameterPack::KeyType::Required, type);
  if(type == "mt19937") return new RNG_mt19937(in);
  FANCY_ASSERT(0, "unknown type given to rng factory, user provided: " + type + " which is not a recognized type.");
  return 0;
}
