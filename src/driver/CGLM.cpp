//Zachariah Vicars
//This folder is the main cpp file for the CGLM driver
//It contains the core CGLM loop
#include "../tools/Assert.hpp"
#include "../core/inputpack.hpp"
#include "../core/factory.hpp"
#include <omp.h>
using string = std::string;
using KeyType = ParameterPack::KeyType;


void buildObjects(InputPack& master_input_pack){
  std::vector<InputPack> sims    = master_input_pack.buildDerivedInputPacks("Simulation"),
                         pvs     = master_input_pack.buildDerivedInputPacks("ProbeVolume"),
                         biases  = master_input_pack.buildDerivedInputPacks("Bias"),
                         hams    = master_input_pack.buildDerivedInputPacks("Hamiltonian"),
                         rngs     = master_input_pack.buildDerivedInputPacks("Random"),
                         lattices = master_input_pack.buildDerivedInputPacks("Lattice");

  for(auto& pv : pvs){
    auto ptr = probevolumeFactory(pv);
    master_input_pack.add(ptr->getName(), ptr);
  }
  for(auto& bias : biases){
    auto ptr = biasFactory(bias);
    master_input_pack.add(ptr->getName(), ptr);
  }
  for(auto& rng : rngs){
    auto ptr = randomFactory(rng);
    master_input_pack.add(ptr->getName(), ptr);
  }
  for(auto& ham : hams){
    auto ptr = hamiltonianFactory(ham);
    master_input_pack.add(ptr->getName(), ptr);
  }
  for(auto& lat : lattices){
    auto ptr = latticeFactory(lat);
    master_input_pack.add(ptr->getName(), ptr);
  }  
  for(auto& sim : sims){
    auto ptr = simulationFactory(sim);
    master_input_pack.add(ptr->getName(), ptr);
  }
  return;
}

void runSingleSimulation(Simulation* simulation){
  //finish object construction
  simulation->run();
  return;
}

void runSimulations(InputPack& master_input_pack){
  auto sims = master_input_pack.simulations();
  #pragma omp parallel for
  for(auto entry : *sims){
    runSingleSimulation(entry.second);
  }
  return;
}

int main(int argc, char *argv[]){
  FANCY_ASSERT(argc = 2, "Expected a single input argument containing the input file name.");
  string input_file_name = argv[1];
  InputParser input_parser;
  //parameter pack objects that will store all pointers to all objects
  ParameterPack master_parameter_pack = input_parser.parseFile(input_file_name);
  InputPack master_input_pack = InputPack(master_parameter_pack);
  //construct objects in factory in the order that makes the most sense
  buildObjects(master_input_pack);
  runSimulations(master_input_pack);
  return 0; //no issues, default return
}