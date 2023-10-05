#include "SimpleGroReader.hpp"
#include "../tools/StringTools.hpp"
#include <fstream>
#include <iostream>
std::vector<SimpleAtom> simpleReadGro(std::string filename){
  std::vector<SimpleAtom> ret_atom_list;
  std::ifstream ifile(filename);
  if(!ifile.is_open()){
    std::cout << "Failed to open gro file " << filename << "std::endl";
    throw 1;
  }
  std::string line;
  //line 1 is a comment
  std::getline(ifile, line);
  //line 2 contains total number of atoms
  std::getline(ifile, line);
  int natoms = std::stoi(line);
  for(int i = 0; i < natoms; i++){
    SimpleAtom a;
    std::getline(ifile, line); //these should all be normal grofile lines in fixed-column format
    int resnr = std::stoi(line.substr(0, 5));
    std::string resname = line.substr(5, 5);
    resname = StringTools::trimWhitespace(resname);
    std::string atomname = line.substr(10, 5);
    atomname = StringTools::trimWhitespace(atomname);
    int atomnumber = std::stoi(line.substr(15, 5));
    std::array<double,3> position, velocity;
    position[0] = std::stod(line.substr(20, 8));
    position[1] = std::stod(line.substr(28, 8));
    position[2] = std::stod(line.substr(36, 8));
    if(line.length() >= 68){
      velocity[0] = std::stod(line.substr(44, 8));
      velocity[1] = std::stod(line.substr(52, 8));
      velocity[2] = std::stod(line.substr(60, 8));
    }
    else{
      velocity[0] = 0.0;
      velocity[1] = 0.0;
      velocity[2] = 0.0;
    }
    a.name = atomname;
    a.x = position;
    a.v = velocity;
    a.resname = resname;
    ret_atom_list.push_back(a);
  }
  return ret_atom_list;
}