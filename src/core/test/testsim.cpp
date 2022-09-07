//initial test file for full simulation
#include "../factory.hpp"

int main(){
  std::cout << "starting test" << std::endl;
  std::string input_file_name = "test_files/tsim.txt";
  InputParser input_parser;
  ParameterPack master_parameter_pack = input_parser.parseFile(input_file_name);
  InputPack master_input_pack = InputPack(master_parameter_pack);
  std::vector<InputPack> sims    = master_input_pack.buildDerivedInputPacks("Simulation"),
                         pvs     = master_input_pack.buildDerivedInputPacks("ProbeVolume"),
                         biases  = master_input_pack.buildDerivedInputPacks("Bias"),
                         hams    = master_input_pack.buildDerivedInputPacks("Hamiltonian"),
                         rngs     = master_input_pack.buildDerivedInputPacks("Random"),
                         lattices = master_input_pack.buildDerivedInputPacks("Lattice");
                         
  std::cout << "derived input packs built" << std::endl;
  for(auto& rng : rngs){
    auto ptr = randomFactory(rng);
    master_input_pack.add(ptr->getName(), ptr);
  }
  std::cout << "rng factory complete" << std::endl;
  for(auto& lat : lattices){
    auto ptr = latticeFactory(lat);
    master_input_pack.add(ptr->getName(), ptr);
  } 
  std::cout << "lattice factory complete" << std::endl;

  for(auto& pv : pvs){
    auto ptr = probevolumeFactory(pv);
    master_input_pack.add(ptr->getName(), ptr);
  } 
  std::cout << "probevolume factory complete" << std::endl;

  for(auto& bias : biases){
    auto ptr = biasFactory(bias);
    master_input_pack.add(ptr->getName(), ptr);
  } 
  std::cout << "bias factory complete" << std::endl;

  for(auto& ham : hams){
    auto ptr = hamiltonianFactory(ham);
    master_input_pack.add(ptr->getName(), ptr);
  } 
  std::cout << "hamiltonian factory complete" << std::endl;
  
  for(auto& sim : sims){
    auto ptr = simulationFactory(sim);
    master_input_pack.add(ptr->getName(), ptr);
  } 
  std::cout << "sim factory complete" << std::endl;

  auto sim_ptr = master_input_pack.find(master_input_pack.simulations(), "sim");
  FANCY_ASSERT(sim_ptr != 0, "failed to find simulation");
  sim_ptr->run();
  return 0;
}