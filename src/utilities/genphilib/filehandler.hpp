#pragma once
#include "../../tools/Matrix.hpp"
#include "../../tools/typedefs.hpp"
#include "../../tools/Assert.hpp"
#include "ff_functions.hpp"

struct AtomFF{
  Vec3<real> x;
  real cutoff;
  int funct;
  real (*funct_ptr)(Vec3<real>, Vec3<real>, Vec3<real>, Vec<real>); //x, box_size, params
  Vec<real> params;
  AtomFF(Vec3<real> x_in, double cutoff, int funct_in, Vec<real> params_in){
    funct = funct_in;
    x = x_in;
    params = params_in;
    if(funct == 1){
      funct_ptr = *LJ_6_12;
    }
    else throw 1;
    return;
  }
};

void loadAtoms(std::string filename, Vec<AtomFF>& atoms);
void writePhiField(std::string output_file, Matrix3d<real>& phi);