#include "../tools/CGLMFileHelper.hpp"
#include "../tools/Matrix.hpp"
#include <limits>
#include <iostream>
bool readFrameBinary(Vec3<std::size_t>& size, Vec<bool>& states, std::ifstream& ifile){
  binary_bool_read(ifile, states, size);
  return ifile.fail();
}
bool writeFrameXYZ(Vec3<std::size_t>& size, Vec<bool>& states, std::ofstream& ofile){
  std::size_t natoms = states.size();
  ofile << natoms << std::endl;
  Matrix3d<char> newMat;
  //easiest way to access the appropriate 1D->3D mapping function
  newMat.initialize(size);
  ofile << "generated from b2xyz" << std::endl;
  for(int i = 0; i < states.size(); i++){
    auto idx = newMat.map1N(i); 
    if(states[i] == 1){
      ofile << "LIQ  " << idx[0] << "  " << idx[1] << "  " << idx[2] << "\n"; 
    }
    else{
      ofile << "VAP  " << idx[0] << "  " << idx[1] << "  " << idx[2] << "\n"; 
    }
  }
  return 0;
}
int main(int argc, char** argv){
  int begin = 0 , end = std::numeric_limits<int>::max(), frequency = 1;
  std::string file_in = "traj.binary", file_out = "traj.xyz";
  for(int i = 1; i < argc; i++){
    std::string arg = argv[i];
    if(arg == "-b"){
      begin = std::stoi(argv[i+1]);
      i++;
      continue;
    }
    if(arg == "-e"){
      end = std::stoi(argv[i+1]);
      i++;
      continue;
    }
    if(arg == "-s"){
      frequency = std::stoi(argv[i+1]);
      i++;
      continue;
    }
    if(arg == "-f"){
      file_in = argv[i+1];
      i++;
      continue;
    }
    if(arg == "-o"){
      file_out = argv[i+1];
      i++;
      continue;
    }    
  }
  Vec3<std::size_t> size;
  std::vector<bool> states;
  int counter = 0; 
  std::ifstream ifile(file_in, std::ios::binary);
  std::ofstream ofile(file_out);
  if(!ofile.is_open()){
    std::cout << "failed to open output file" << std::endl;
    throw 1;
  }
  while(true){
    bool fail = readFrameBinary(size,  states, ifile);
    if(fail) break;
    if(counter >= begin && counter <= end && (counter-begin)%frequency == 0){
      writeFrameXYZ(size, states, ofile);
    }
    counter++;
  }
  return 0;
}