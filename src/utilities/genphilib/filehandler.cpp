#include "filehandler.hpp"
#include <fstream>
#include <sstream>
#include "../../tools/StringTools.hpp"
#include "../../tools/CGLMFileHelper.hpp"
void loadAtoms(std::string filename, Vec<AtomFF>& atoms, Vec3<real> size){
  atoms.clear();
  std::ifstream ifile(filename);
  FANCY_ASSERT(ifile.is_open(), "failed to open input file");
  std::string line;
  while(getline(ifile, line)){
    StringTools::trimWhitespace(line);
    if(line.size() == 0) continue;
    if(line.at(0) == '#') continue;
    std::stringstream ss(line);
    Vec<real> params;
    real cutoff;
    Vec3<real> pos;
    int funct;
    real temp;
    ss >> pos[0] >> pos[1] >> pos[2];
    ss >> cutoff;
    ss >> funct;
    while(ss >> temp){
      params.push_back(temp);
    }
    atoms.push_back(AtomFF(pos, cutoff, funct, params, size));
  }
  return;
}

void writePhiField(std::string output_file, Matrix3d<real>& phi){
  auto size = phi.size();
  std::ofstream ofile(output_file);
  FANCY_ASSERT(ofile.is_open(), "failed to open output file");
  std::vector<real> phi_1d(phi.size_1d());
  for(int i = 0; i < phi.size_1d(); i++){
    phi_1d[i] = phi.at_1d(i);
  }
  binary_real_write(ofile, phi_1d, size);
  return;
}