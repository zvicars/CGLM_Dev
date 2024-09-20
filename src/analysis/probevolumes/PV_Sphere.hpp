#pragma once
#include "ProbeVolume.hpp"
class PV_Sphere : public ProbeVolume{
public: 
  PV_Sphere(AnalysisInputPack& input);
  PV_Sphere(Vec3<real> pos, real radius, Vec3<real> box_size){
    type_ = "pv_sphere";
    pos_ = pos;
    radius_ = radius;
    box_size_ = box_size; 
    bSetBoxsize = 1;
    return;
  }
  virtual double compute(Vec3<double> position);
  virtual double compute_periodic(Vec3<double> position);
  virtual double compute_nonperiodic(Vec3<double> position); 
  void setBoxSize(Vec3<real> box_size){
    box_size_ = box_size;
    bSetBoxsize = 1;
    return;
  }
  Vec<real> getBounds(){
    Vec<real> output(6);
    for(int i = 0; i < 3; i++){
      output[i] = pos_[i] - radius_;
      output[i+3] = pos_[i] + radius_;
    }
    return output;
  }
  Vec3<real> getPos(){
    return pos_;
  }
private:
  std::array<real,3> pos_;
  real radius_;
  Vec3<real> box_size_;
  bool bSetBoxsize = 0;
};