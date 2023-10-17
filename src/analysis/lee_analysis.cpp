#include "lee_analysis.hpp"
#include <map>

bool readFrameBinary(Vec3<std::size_t>& size, Vec<real>& states, std::ifstream& ifile){
  binary_real_read(ifile, states, size);
  return ifile.fail();
}
bool readFrameBinary(Vec3<std::size_t>& size, Vec<char>& states, std::ifstream& ifile){
  binary_char_read(ifile, states, size);
  return ifile.fail();
}

int main(int argc, char** argv){
  int begin = 0 , end = std::numeric_limits<int>::max(), frequency = 1;
  //remove isolated vapors/waters
  bool denoise=0;
  std::string file_in = "traj.traj", file_out = "traj.xyz", file_phi = "in.phi";
  for(int i = 1; i < argc; i++){
    std::string arg = argv[i];
    //first frame
    if(arg == "-b"){
      begin = std::stoi(argv[i+1]);
      i++;
      continue;
    }
    //last frame
    if(arg == "-e"){
      end = std::stoi(argv[i+1]);
      i++;
      continue;
    }
    //frequency with which to consider trajectory frames
    if(arg == "-s"){
      frequency = std::stoi(argv[i+1]);
      i++;
      continue;
    }
    //trajectory file
    if(arg == "-f"){
      file_in = argv[i+1];
      i++;
      continue;
    }
    //file to output data to
    if(arg == "-o"){
      file_out = argv[i+1];
      i++;
      continue;
    }
    //file containing the external potentials
    if(arg == "-phi"){
      file_phi = argv[i+1];
      i++;
      continue;
    }    
    if(arg == "-denoise"){
      denoise = 1;
      continue;
    }
  }
  //initialize necessary variables
  //box-sizes in number of lattice sites in x,y,z dimensions
  Vec3<std::size_t> size, size_phi;
  std::vector<char> states;
  std::vector<real> phi;
  int counter = 0; 
  //get file handle for trajectory
  std::ifstream ifile(file_in, std::ios::binary);
  //get file handle for phi file
  std::ifstream phi_file(file_phi, std::ios::binary);
  //read phi file
  readFrameBinary(size_phi,  phi, ifile);

  std::ofstream ofile(file_out);
  if(!ofile.is_open()){
    std::cout << "failed to open output file" << std::endl;
    throw 1;
  }
  while(true){
    bool fail = readFrameBinary(size,  states, ifile);
    if(fail) break;
    if(counter >= begin && counter <= end && (counter-begin)%frequency == 0){
      if(denoise) denoise_states(states);
    }
    counter++;
  }
  ofile.close();
  return 0;
}
