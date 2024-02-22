#pragma once
#include "../../tools/Matrix.hpp"
#include "../../tools/typedefs.hpp"
#include "../../tools/Assert.hpp"
#include "ff_functions.hpp"

struct AtomFF{
  Vec3<real> x; //ref x
  real cutoff; //cutoff
  int funct; //function type
  FF_function* ff=0; //pointer to FF_function object
  Vec<real> params; //parameter list
  Vec3<real> size; //box size 
  AtomFF(Vec3<real> x_in, real cutoff_in, int funct_in, Vec<real> params_in, Vec3<real> size_in){
    funct = funct_in;
    cutoff = cutoff_in;
    x = x_in;
    params = params_in;
    size = size_in;
    setPtr();
    //standard lj particles
    return;
  }
  AtomFF(const AtomFF& t){
    funct = t.funct;
    cutoff = t.cutoff;
    x = t.x;
    params = t.params;
    size = t.size;
    setPtr();
  }
  ~AtomFF(){
    delete ff;
  }
  void setPtr(){
    if(funct == 1){
      ff = new LJ_6_12_default(x, size, cutoff, params);
    }
    //big particles
    else if(funct == 2){
      ff = new LJ_6_12_offset(x, size, cutoff, params);
    }
    //box-shaped particles
    else if(funct == 3){
      ff = new LJ_6_12_offset_box(x, size, cutoff, params);
    }
    //cylinder-shaped particles
    else if(funct == 4){
      ff = new LJ_6_12_offset_cylinder(x, size, cutoff, params);
    }
    //lj 3-9 big particles
    else if(funct == 5){
      ff = new LJ_3_9_offset(x, size, cutoff, params);
    }    
    //lj 3-9 box-shaped particles
    else if(funct == 6){
      ff = new LJ_3_9_offset_box(x, size, cutoff, params);
    }
    //lj 3-9 cylinder-shaped particles
    else if(funct == 7){
      ff = new LJ_6_12_offset_cylinder(x, size, cutoff, params);
    }
    else throw 1;
  }
};

void loadAtoms(std::string filename, Vec<AtomFF>& atoms, Vec3<real> size);
void writePhiField(std::string output_file, Matrix3d<real>& phi);