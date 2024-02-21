#pragma once
#include "../../tools/Matrix.hpp"
#include "../../tools/typedefs.hpp"
#include "../../tools/Assert.hpp"
#include "ff_functions.hpp"

struct AtomFF{
  Vec3<real> x;
  real cutoff;
  int funct;
  real (*funct_ptr)(Vec3<real>, Vec3<real>, Vec3<real>, real, Vec<real>); //x, box_size, params
  Vec<real> params;
  AtomFF(Vec3<real> x_in, double cutoff_in, int funct_in, Vec<real> params_in){
    funct = funct_in;
    cutoff = cutoff_in;
    x = x_in;
    params = params_in;
    //standard lj particles
    if(funct == 1){
      funct_ptr = *LJ_6_12_default;
    }
    //big particles
    else if(funct == 2){
      funct_ptr = *LJ_6_12_offset;
    }
    //box-shaped particles
    else if(funct == 3){
      funct_ptr = *LJ_6_12_offset_box;
    }
    //cylinder-shaped particles
    else if(funct == 4){
      funct_ptr= *LJ_6_12_offset_cylinder;
    }
    //lj 3-9 big particles
    else if(funct == 5){
      funct_ptr = *LJ_3_9_offset;
    }    
    //lj 3-9 box-shaped particles
    else if(funct == 6){
      funct_ptr = *LJ_3_9_offset_box;
    }
    //lj 3-9 cylinder-shaped particles
    else if(funct == 7){
      funct_ptr= *LJ_3_9_offset_cylinder;
    }
    else throw 1;
    return;
  }
};

void loadAtoms(std::string filename, Vec<AtomFF>& atoms);
void writePhiField(std::string output_file, Matrix3d<real>& phi);