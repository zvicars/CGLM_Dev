//lattice testing script
//make a fake lattice, test some basic functonality
//tests to make sure adj() function works properly and that periodic boundaries behave as expected
#include "../../factory.hpp"

void getCounts(Lattice* lattice_ptr, int& total_liquid_count, int& total_adjacent_count){
  total_adjacent_count = 0;
  total_liquid_count = 0;
  auto size = lattice_ptr->size();
  for(std::size_t i = 0; i < size[0]; i++){
    for(std::size_t j = 0; j < size[1]; j++){
      for(std::size_t k = 0; k < size[2]; k++){
        std::array<std::size_t, 3> idx = {i, j, k};
        total_liquid_count += (int)lattice_ptr->getState(idx);
        total_adjacent_count += (int)lattice_ptr->getAdj(idx);
      }
    }
  }

  return;
}

int main(){
  std::cout << "starting test" << std::endl;
  std::string input_file_name = "test_files/t1d.txt";
  InputParser input_parser;
  ParameterPack master_parameter_pack = input_parser.parseFile(input_file_name);
  InputPack master_input_pack = InputPack(master_parameter_pack);
  std::vector<InputPack> rngs     = master_input_pack.buildDerivedInputPacks("Random"),
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
  //master input pack should have a single lattice
  auto lattice_ptr = master_input_pack.findLattice("test_lattice");
  std::cout << "test lattice found" << std::endl;
  //test lattice will be fully liquid
  //confirm some basic facts about the lattice
  int total_adjacent_count = 0, total_liquid_count = 0;
  auto size = lattice_ptr->size();
  getCounts(lattice_ptr, total_liquid_count, total_adjacent_count);
  bool f1=0, f2=0;
  if(total_liquid_count != size[0]*size[1]*size[2]){
    std::cout << "total liquid count incorrect\n";
    f1 = 1;
  }
  if(total_adjacent_count != size[0]*size[1]*size[2]*6){
    std::cout << "total adjacent count incorrect\n";
    std::cout << "expected " << size[0]*size[1]*size[2]*6 << " got " << total_adjacent_count << "\n";;
    f2 = 1;    
  }
  if(f1 || f2 ){
    std::cout << "lattice test failed at step 1\n";
    return 1;
  }

  //modify a single site and make sure that works
  std::array<std::size_t, 3> idx = {0,0,0};
  lattice_ptr->flip(idx);
  getCounts(lattice_ptr, total_liquid_count, total_adjacent_count);
  f1=0, f2=0;
  if(total_liquid_count != size[0]*size[1]*size[2] - 1){
    std::cout << "total liquid count incorrect\n";
    f1 = 1;
  }
  if(total_adjacent_count != size[0]*size[1]*size[2]*6 - 6){
    std::cout << "total adjacent count incorrect\n";
    std::cout << "expected " << size[0]*size[1]*size[2]*6 - 6 << " got " << total_adjacent_count << "\n";;
    f2 = 1;    
  }
  if(f1 || f2 ){
    std::cout << "lattice test failed at step 2\n";
    return 1;
  } 

  //flip it back and make sure nothing changed
  lattice_ptr->flip(idx);
  getCounts(lattice_ptr, total_liquid_count, total_adjacent_count);
  f1=0, f2=0;
  if(total_liquid_count != size[0]*size[1]*size[2]){
    std::cout << "total liquid count incorrect\n";
    f1 = 1;
  }
  if(total_adjacent_count != size[0]*size[1]*size[2]*6){
    std::cout << "total adjacent count incorrect\n";
    std::cout << "expected " << size[0]*size[1]*size[2]*6 << " got " << total_adjacent_count << "\n";
    f2 = 1;    
  }
  if(f1 || f2 ){
    std::cout << "lattice test failed at step 2\n";
    return 1;
  } 

  //modify two adjacent sites and make sure it works
  std::array<std::size_t, 3> idx2 = {1,0,0};
  lattice_ptr->flip(idx);
  lattice_ptr->flip(idx2);
  getCounts(lattice_ptr, total_liquid_count, total_adjacent_count);
  //two adjacent sites, so I should have lost 12 total using the naive count
  f1=0, f2=0;
  if(total_liquid_count != size[0]*size[1]*size[2] - 2 ){
    std::cout << "total liquid count incorrect\n";
    f1 = 1;
  }
  if(total_adjacent_count != size[0]*size[1]*size[2]*6 - 12){
    std::cout << "total adjacent count incorrect\n";
    f2 = 1;    
  }
  if(f1 || f2 ){
    std::cout << "lattice test failed at step 3\n";
    return 1;
  } 
  std::cout << "all tests passed\n" << std::endl;
  std::string report_str = "";
  lattice_ptr->reportInfo(report_str);
  std::cout << report_str << std::endl;


  //test pbc functions
  std::array<int, 3> pbc_idx1 = {3,0,0};
  std::array<int, 3> pbc_idx2 = {-3,0,0};
  std::array<int, 3> pbc_idx3 = {3,3,3};
  std::array<int, 3> pbc_idx4 = {-3,-3,-3};
  std::array<std::size_t, 3> wrapped_idx = {0,0,0};
  if(lattice_ptr->getAdjPBC(pbc_idx1) != lattice_ptr->getAdj(wrapped_idx)){
    std::cout << "failed pbc test 1\n";
    return 2;
  }
  if(lattice_ptr->getAdjPBC(pbc_idx2) != lattice_ptr->getAdj(wrapped_idx)){
    std::cout << "failed pbc test 2\n";
    return 2;
  }
  if(lattice_ptr->getAdjPBC(pbc_idx3) != lattice_ptr->getAdj(wrapped_idx)){
    std::cout << "failed pbc test 3\n";
    return 2;
  }
  if(lattice_ptr->getAdjPBC(pbc_idx4) != lattice_ptr->getAdj(wrapped_idx)){
    std::cout << "failed pbc test 4\n";
    return 2;
  }


  return 0;
}