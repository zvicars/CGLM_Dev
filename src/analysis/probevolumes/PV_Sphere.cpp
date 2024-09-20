#include "PV_Sphere.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/pbcfunctions.hpp"
PV_Sphere::PV_Sphere(AnalysisInputPack& input):ProbeVolume{input}
{
  Vec<real> posvec;
  input.params().readVector("pos", KeyType::Required, posvec);
  FANCY_ASSERT(posvec.size() == 3, "pos argument should have x, y, and z position");
  for(int i = 0; i < 3; i++){
    pos_[i] = posvec[i];
  }
  input.params().readNumber("radius", KeyType::Required, radius_);
  return;
}


//pbc correction should happen outside of function call
double PV_Sphere::compute(Vec3<double> position){
  if(bSetBoxsize){
    return compute_periodic(position);
  }
  return compute_nonperiodic(position);
}

double PV_Sphere::compute_nonperiodic(Vec3<double> position){
  Vec3<real> rvec = position - pos_;
  real r = norm2(rvec);
  if(r <= radius_) return 1.0;
  return 0.0;
}

double PV_Sphere::compute_periodic(Vec3<double> position){
  real r = getDistance(position, pos_, box_size_);
  if(r <= radius_) return 1.0;
  return 0.0;
}
