//script to generate a phi potential field with a given box size and a given set of input parameters
//input file contains a list of atom positions and lennard jones parameters
#include "genphilib/filehandler.hpp"
#include "../tools/stlmath.hpp"
#include "../tools/pbcfunctions.hpp"
#include <fstream>
real computePhiForSingleCell(const Vec3<std::size_t>& index, const Vec3<std::size_t>& box_size, real spacing, const AtomFF& atom_in){
  real eval=0.0;
  Vec3<real> point;
  Vec3<real> size;
  for(int i=0; i<3; i++){
    point[i] = ((real)index[i]+0.5)*spacing; //place test point in middle of lattice cell
    size[i] = box_size[i]*spacing;
  }
  eval = atom_in.funct_ptr(point, atom_in.x, size, atom_in.params);
  return eval;
}

real computePhiForSingleCellIntegrated(const Vec3<std::size_t>& index, const Vec3<std::size_t>& box_size, real spacing, const AtomFF& atom_in){
  Vec3<real> point;
  Vec3<real> size;
  for(int i=0; i<3; i++){
    point[i] = ((real)index[i] + 0.5)*spacing; //place test point in middle of lattice cell
    size[i] = box_size[i]*spacing;
  }
  real eval=0.0;
  real counter = 0.0;
  for(real i = 0.0; i <= 1.0; i+=0.2){
  for(real j = 0.0; j <= 1.0; j+=0.2){
  for(real k = 0.0; k <= 1.0; k+=0.2){
        Vec3<real> new_point = point + spacing*Vec3<real>{i,j,k};
        real tempEval = atom_in.funct_ptr(new_point, atom_in.x, size, atom_in.params);
        real weight = 1;//exp(-tempEval);
        eval += weight*tempEval; 
        counter+=weight;
      }
    }
  }
  eval = eval / counter;
  return eval;
}

Vec3<int> pos2idx(Vec3<real> x, Vec3<real> size, real spacing){
  placeInsideBox(x, size);
  Vec3<int> eval;
  for(int i = 0; i < 3; i++){
    eval[i] = floor(x[i]/spacing);
  }
  return eval;
}

void computePhiField(Matrix3d<real>& phi, const Vec<AtomFF>& atoms, Vec3<std::size_t> box_size, real spacing, bool integrate){
  phi.initialize(box_size);
  phi.fill(0.0);
  Vec3<real> size;
  for(int i = 0; i < 3; i++){
    size[i] = box_size[i]*spacing;
  }
  //get the range of cells to consider
  for(auto& atom : atoms){
    int range = atom.cutoff / spacing  + 1;
    Vec3<int> at_idx = pos2idx(atom.x, size, spacing);
    Vec3<int> min, max;
    for(int i = 0; i < 3; i++){
      min[i] = at_idx[i] - range;
      max[i] = at_idx[i] + range;
      if(max[i] - min[i] > 0.5*box_size[i]){
        min[i] = 0;
        max[i] = box_size[i]-1;
      }
    }
    for(int i = min[0]; i <= max[0]; i++){
      std::size_t i_temp = wrapIndex(i, box_size[0]);
      for(int j = min[1]; j <= max[1]; j++){
        std::size_t j_temp = wrapIndex(j, box_size[1]);
        for(int k = min[2]; k <= max[2]; k++){
          std::size_t k_temp = wrapIndex(k, box_size[2]);
          Vec3<std::size_t> real_idx = {i_temp, j_temp, k_temp};
          double eval;
          if(integrate) eval = computePhiForSingleCellIntegrated(real_idx, box_size, spacing, atom);
          else eval = computePhiForSingleCell(real_idx, box_size, spacing, atom);
          if(eval > 100) eval = 100;
          phi.at(real_idx) += eval;
        }
      }
    }
  }
  return;
}



int main (int argc, char** argv){
  real spacing;
  std::size_t box_x, box_y, box_z;
  std::string input_file, output_file;
  bool integrate = 0;
  for(int i = 1; i < argc; i++){
    std::string arg = argv[i];
    if(arg == "-box"){
      FANCY_ASSERT(argc >= i+4, "not enough input arguments for -box, need -box x y z");
      box_x = std::stoi(argv[i+1]);
      box_y = std::stoi(argv[i+2]);
      box_z = std::stoi(argv[i+3]);
      i+=3;
      continue;
    }
    if(arg == "-gs"){
      FANCY_ASSERT(argc >= i+2, "not enough input arguments for -box, need -gs <float>");
      spacing = std::stod(argv[i+1]);
      i++;
      continue;
    }
    if(arg == "-i"){
      FANCY_ASSERT(argc >= i+2, "not enough input arguments for -i, need -i <input_file>");
      input_file = argv[i+1];
      i++;
      continue;   
    }
    if(arg == "-o"){
      FANCY_ASSERT(argc >= i+2, "not enough input arguments for -i, need -o <output_file>");
      output_file = argv[i+1];
      i++;
      continue;
    }
    if(arg == "-int"){
      FANCY_ASSERT(argc >= i+2, "not enough input arguments for -int, need -o yes/no/y/n");
      std::string arg = argv[i+1];
      if(arg.find('y') == 0) integrate = 1;
      else integrate = 0;
    }
  }  

  Matrix3d<real> phi;
  Vec<AtomFF> atoms;
  loadAtoms(input_file, atoms);
  computePhiField(phi, atoms, {box_x, box_y, box_z}, spacing, integrate);
  writePhiField(output_file, phi);
  return 0;
}