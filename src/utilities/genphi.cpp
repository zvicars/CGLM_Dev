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
  eval = atom_in.ff->compute(point);
  return eval;
}

real computePhiForSingleCellIntegrated(const Vec3<std::size_t>& index, const Vec3<std::size_t>& box_size, real spacing, const AtomFF& atom_in, int npoints){
  Vec3<real> point;
  Vec3<real> size;
  real point_spacing = 1.0;
  Vec3<real> starting_point = {0.0, 0.0, 0.0};
  if(npoints == 1) starting_point = {0.5, 0.5, 0.5};
  if(npoints > 1){
    point_spacing = 1.0 / (real)(npoints - 1);
  }
  for(int i=0; i<3; i++){
    point[i] = (real)index[i]*spacing; //place test point at the origin of the lattice cell
    size[i] = box_size[i]*spacing;
  }
  real eval=0.0;
  real counter = 0.0;
  for(int i = 0; i < npoints; i++){
  for(int j = 0; j < npoints; j++){
  for(int k = 0; k < npoints; k++){
        Vec3<real> ref_point = point + (starting_point*spacing); //place point in middle of cell if npoints = 1 
        Vec3<real> rel_pos = {(real)i * point_spacing, (real)j * point_spacing, (real)k * point_spacing};
        Vec3<real> new_point = ref_point + spacing*rel_pos;
        real tempEval = atom_in.ff->compute(new_point);
        real weight = 1;
        //real weight = exp(-tempEval);
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

void computePhiField(Matrix3d<real>& phi, const Vec<AtomFF>& atoms, Vec3<std::size_t> box_size, real spacing, bool integrate, int integrate_npoints){
  phi.initialize(box_size);
  phi.fill(0.0);
  Vec3<real> size;
  for(int i = 0; i < 3; i++){
    size[i] = box_size[i]*spacing;
  }
  //get the range of cells to consider
  for(auto& atom : atoms){
    Vec<real> bounding_box = atom.ff->getBoundingBox();
    int range = atom.cutoff / spacing  + 1;
    Vec3<int> min, max;
    for(int i = 0; i < 3; i++){
      min[i] = std::round(bounding_box[i] / spacing) - 1;
      max[i] = std::round(bounding_box[i+3] / spacing) + 1;
      if(max[i] - min[i] > 0.5*box_size[i]){
        min[i] = 0;
        max[i] = box_size[i]-1;
      }
    }
    for(int i = min[0]; i <= max[0]; i++){
      for(int j = min[1]; j <= max[1]; j++){
        for(int k = min[2]; k <= max[2]; k++){
          Vec3<std::size_t> wrapped_idx = {(std::size_t)wrapIndex(i, box_size[0]), (std::size_t)wrapIndex(j, box_size[1]), (std::size_t)wrapIndex(k, box_size[2])};
          double eval;
          if(integrate) eval = computePhiForSingleCellIntegrated(wrapped_idx, box_size, spacing, atom, integrate_npoints);
          else eval = computePhiForSingleCell(wrapped_idx, box_size, spacing, atom);
          phi.at(wrapped_idx) += eval;
          if(phi.at(wrapped_idx) > 100) phi.at(wrapped_idx) = 100;
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
  int integrate_npoints = 1;
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
      FANCY_ASSERT(argc >= i+2, "not enough input arguments for -int, need -int <npoints>");
      std::string arg = argv[i+1];
      integrate = 1;
      integrate_npoints = std::stoi(arg);
      continue;
    }
  }  

  //size computation, needed when defining atoms (kludgy, but WIP)
  Vec3<std::size_t> box_size = {box_x, box_y, box_z};
  Vec3<real> size;
  for(int i=0; i<3; i++){
    size[i] = (real)box_size[i]*spacing;
  }

  Matrix3d<real> phi;
  Vec<AtomFF> atoms;
  loadAtoms(input_file, atoms, size);
  computePhiField(phi, atoms, {box_x, box_y, box_z}, spacing, integrate, integrate_npoints);
  writePhiField(output_file, phi);
  return 0;
}
