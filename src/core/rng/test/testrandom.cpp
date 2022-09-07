//lattice testing script
//make a fake lattice, test some basic functonality
//tests to ensure that random number generator is behaving reasonably
#include "../../factory.hpp"


void testRandom(RNG* r_ptr){
  FANCY_ASSERT(r_ptr != 0, "failed to find rng pointer");
  std::vector<int> histogram(10, 0);
  std::vector<double> histogram2(10,0);
  for(int i = 0; i < 1e6; i++){
    std::size_t index = r_ptr->getIndex(0, 9);
    FANCY_ASSERT(index >= 0 && index < 10, "index out of range");
    histogram[index]++;
  }
  for(int i = 0; i < 10; i++){
    histogram2[i] = (double)histogram[i]/(double)1e6;
    FANCY_ASSERT(std::fabs(histogram2[i] - 0.10) < 0.001, "non-uniform distribution detected for rng getIndex()");
  }
  
  for(auto& val : histogram) val = 0;
  for(int i = 0; i < 1e6; i++){
    double val = r_ptr->getReal(0.0, 10.0);
    if(val == 10.0) val = 0.0;
    FANCY_ASSERT(val >= 0.0 && val < 10.0, "value out of range");
    std::size_t index = floor(val);
    histogram[val]++;
  }
  for(int i = 0; i < 10; i++){
    histogram2[i] = (double)histogram[i]/(double)1e6;
    FANCY_ASSERT(std::fabs(histogram2[i] - 0.10) < 0.001, "non-uniform distribution detected for rng getIndex()");
  } 
  return;  
}

int main(){
  std::cout << "starting test" << std::endl;
  std::string input_file_name = "test_files/trand.txt";
  InputParser input_parser;
  ParameterPack master_parameter_pack = input_parser.parseFile(input_file_name);
  InputPack master_input_pack = InputPack(master_parameter_pack);
  std::vector<InputPack> rngs     = master_input_pack.buildDerivedInputPacks("Random");
  std::cout << "derived input packs built" << std::endl;
  for(auto& rng : rngs){
    auto ptr = randomFactory(rng);
    master_input_pack.add(ptr->getName(), ptr);
  }
  std::cout << "rng factory complete" << std::endl;
  auto r_ptr = master_input_pack.find(master_input_pack.randomgenerators(), "rng");
  auto rng_map = master_input_pack.randomgenerators();
  for(auto& entry : *rng_map){
    auto rng_ptr = entry.second;
    auto cloned_rng_ptr = rng_ptr->clone();
    testRandom(cloned_rng_ptr);
    delete cloned_rng_ptr;
  }

  return 0;
}