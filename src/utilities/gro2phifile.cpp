//Zachariah Vicars
//This folder is the main cpp file for the CGLM driver
//It contains the core CGLM loop
#include <omp.h>
#include <cmath>
#include "../tools/Assert.hpp"
#include "../tools/InputParser.hpp"
#include "../tools/Matrix.hpp"
#include "../tools/pbcfunctions.hpp"
#include "SimpleGroReader.hpp"

using string = std::string;
using KeyType = ParameterPack::KeyType;

struct Atomdef{
  std::string name, atomname, resname;
  std::vector<double> params;
  double cutoff;
  int function;
};

bool isInAtomdef(const SimpleAtom& atom, const Atomdef& atom_def){
  if(atom.name == atom_def.atomname == 1 && atom_def.resname == atom.resname) return 1;
  return 0;
}

int main(int argc, char *argv[]){
  FANCY_ASSERT(argc == 2, "Expected a single input argument containing the input file name.");
  string input_file_name = argv[1], output_file;
  string grofile;
  InputParser input_parser;
  //parameter pack objects that will store all pointers to all objects
  ParameterPack mpp = input_parser.parseFile(input_file_name);
  mpp.readString("OutputFile", ParameterPack::KeyType::Required, output_file);
  //construct objects in factory in the order that makes the most sense
  auto atomdefs_input = mpp.findParameterPacks("Atomdef", ParameterPack::KeyType::Required);
  std::vector<Atomdef> atomdefs;
  for(auto atomdef : atomdefs_input){
    std::string name, atomname, resname;
    std::vector<double> params;
    int funct = -1;
    double cutoff = 0;
    double r0;
    atomdef->readString("name",ParameterPack::KeyType::Optional, name);
    //any atom name
    atomdef->readString("atomname", ParameterPack::KeyType::Required, atomname);
    atomdef->readString("resname", ParameterPack::KeyType::Required, resname);
    atomdef->readVector("params", ParameterPack::KeyType::Required, params);
    atomdef->readNumber("function", ParameterPack::KeyType::Required, funct);
    atomdef->readNumber("cutoff", ParameterPack::KeyType::Required, cutoff);
    if(funct == 2) atomdef->readNumber("r0", ParameterPack::KeyType::Optional, r0);
    Atomdef new_atomdef; new_atomdef.atomname = atomname;
    new_atomdef.resname = resname;
    new_atomdef.params = params;
    new_atomdef.function = funct;
    new_atomdef.cutoff = cutoff;
    atomdefs.push_back(new_atomdef);
  }
  std::cout << "Found " << atomdefs.size() << " atomdefs." << std::endl;
  double groscale=1.0;
  mpp.readString("GroFile", ParameterPack::KeyType::Required, grofile);
  mpp.readNumber("GroScale", ParameterPack::KeyType::Optional, groscale);
  auto grofileAtoms = simpleReadGro(grofile);

  std::ofstream ofile(output_file);
  //loop over atomdefs, build a cell-grid
  for(const auto& atom : grofileAtoms){
    int atom_def_index = -1;
    for(int i = 0; i < atomdefs.size(); i++){
      const auto& atomdef  = atomdefs[i];
      if(isInAtomdef(atom, atomdef)){
        atom_def_index = i;
        break;
      }
    }
    if(atom_def_index == -1) continue;
    ofile << groscale*atom.x[0] << "  " << groscale*atom.x[1] << "  " << groscale*atom.x[2] << "  ";
    ofile << atomdefs[atom_def_index].cutoff << "  " << atomdefs[atom_def_index].function << "  ";
    for(auto ad_val : atomdefs[atom_def_index].params){
      ofile << ad_val << "  ";
    }
    ofile << "\n";
  }
  ofile.close();
  return 0; //no issues, default return
}