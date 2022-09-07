//lattice testing script
//make a fake lattice, test some basic functonality
//this tests the lattice gas hamiltonian to ensure that some reasonable behaviors are observed
//(1) for a fully liquid 3x3x3 lattice, the values of h and dh are known, I compare the analytical solution to the
//computed solution
//(2) I ensure that when I flip a site, the energy is conserved (h1 + dh1 = h2) and the flip is reversible (dh1 = -dh2)
//also includes a linear bias to ensure that's being computed correctly
#include "../../factory.hpp"

int main(){
  std::cout << "starting test" << std::endl;
  std::string input_file_name = "test_files/tlg.txt";
  InputParser input_parser;
  ParameterPack master_parameter_pack = input_parser.parseFile(input_file_name);
  InputPack master_input_pack = InputPack(master_parameter_pack);
  std::vector<InputPack> rngs     = master_input_pack.buildDerivedInputPacks("Random"),
                         lattices = master_input_pack.buildDerivedInputPacks("Lattice"),
                         hams     = master_input_pack.buildDerivedInputPacks("Hamiltonian"),
                         pvs      = master_input_pack.buildDerivedInputPacks("ProbeVolume"),
                         biases    = master_input_pack.buildDerivedInputPacks("Bias");
                         
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
  //master input pack should have a single lattice
  auto lattice_ptr = master_input_pack.findLattice("test_lattice");
  std::cout << "test lattice found" << std::endl;
  //test lattice will be fully liquid
  auto ham_ptr = master_input_pack.find(master_input_pack.hamiltonians(), "ham");
  lattice_ptr->setActiveSite({0,0,0});
  auto active_idx = lattice_ptr->getActiveIndex();
  std::cout << "active site = " << active_idx[0] << "  " << active_idx[1] << "  " << active_idx[2] << "\n";
  std::cout << "active state = " << lattice_ptr->getActiveState() << "\n";
  std::cout << "active mu = " << lattice_ptr->getActiveMu() << "\n";
  std::cout << "active phi = " << lattice_ptr->getActivePhi() << "\n"; 
  std::cout << "active adj = " << lattice_ptr->getActiveAdj() << "\n";
  std::cout << std::endl;
  double h_out = ham_ptr->calc_h(*lattice_ptr); //initial call to h should also set n and u, needs to be called befor calc_dh
  double dh_out = ham_ptr->calc_dh(*lattice_ptr);
  //FANCY_ASSERT(h_out == (-3*3*3*lattice_ptr->eps() + 3*3*3*lattice_ptr->mu() + 16), "testlg has the incorrect value for h");
  bool flag1 = h_out == (-27*3*lattice_ptr->eps() + 27*lattice_ptr->getActiveMu() + 8);
  bool flag2 = dh_out == (6*lattice_ptr->eps() - lattice_ptr->getActiveMu() - 1);
  FANCY_ASSERT(flag1, "incorrect value of h_out in first step of testlg, expected = "
   + std::to_string(-27*3*lattice_ptr->eps() + 27*lattice_ptr->getActiveMu() + 8) + ", output = "
   + std::to_string(h_out) + "\n");
  FANCY_ASSERT(flag2, "incorrect value of dh_out in first step of testlg, expected = "
   + std::to_string(6*lattice_ptr->eps() - lattice_ptr->getActiveMu() - 1) + ", output = "
   + std::to_string(dh_out) + "\n");


  lattice_ptr->flipActive();
  lattice_ptr->setActiveSite({0,0,0});
  ham_ptr->flip();
  double dh_realized = ham_ptr->h() - h_out; 
  FANCY_ASSERT(dh_realized == dh_out, "incorrect value for dh_realized");
  double h_old = ham_ptr->h();
  double h_out2 = ham_ptr->calc_h(*lattice_ptr);
  double dh_out2 = ham_ptr->calc_dh(*lattice_ptr);
  FANCY_ASSERT(h_out2 ==  h_old, "hamiltonian calculation is inconsistent");
  FANCY_ASSERT(dh_out2 == -dh_out, "hamiltonian is not reversible h1 = " + std::to_string(dh_out) + ", and h2 = " + std::to_string(dh_out2));
  return 0;
}