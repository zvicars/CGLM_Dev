//ramping test script
//compares output of ramping header to reference data that was manually checked
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "../ramped_parameter.hpp"
#define TOL 0.0001

//just contains a raw list of output values
std::vector<real> loadRefVals(std::string filename){
  std::ifstream input(filename);
  FANCY_ASSERT(input.is_open(), "Failed to open input file for test.");
  std::vector<real> retVec;
  std::string line;
  while(getline(input, line)){
    std::stringstream ss(line);
    real val;
    while(ss >> val){
      retVec.push_back(val);
    }
  }
  return retVec;
}

int main(){
  std::vector<real> xy, x, y;
  std::vector<real> ref_values = loadRefVals("test_files/ref.txt");
  xy = {  2,   5,   10,  20,  40,  80, 160, 
        100, 150,  160, 180, 140, 100, 0};
  x = {2,   5,   10,  20,  40,  80, 160};
  y = {100, 150,  160, 180, 140, 100, 0};
  RampedParameter R1(xy);
  RampedParameter R2(x, y);
  std::size_t ref_idx = 0;
  for(real t = 0.0; t <= 200.0; t+= 0.1){
    real r1_ret = R1.compute(t);
    real r2_ret = R2.compute(t);
    if(r1_ret != r2_ret){
      std::cout << "different constructors give different outputs!" << std::endl;
      throw 1; 
    }
    if(std::fabs(r1_ret - ref_values[ref_idx]) > TOL){
      std::cout << "disagreement with reference for t =" << t << std::endl;
    }
    ref_idx++;
  }
  std::cout << "All is well!" << std::endl;
  return 0;
}