#include "ff_functions.hpp"
#include "../../tools/SimpleForcefields.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/pbcfunctions.hpp"
#include <iostream>

//https://stackoverflow.com/questions/5254838/calculating-distance-between-a-point-and-a-rectangular-box-nearest-point
real nearestDistanceToBox(Vec3<real> p, Vec3<real> xref, Vec3<real> min, Vec3<real> max, Vec3<real> size){
  Vec3<real> dx;
  getNearestImage3D(p, xref, size);
  for(int i = 0; i < 3; i++){
    dx[i] = std::max(min[i]-p[i], std::max(0.0, p[i]-max[i]));
  }
  return norm2(dx);
}

Vec<real> getCylindricalBounds(Vec3<real> base, real radius, real height, int axis, real cutoff){
  Vec<real> output(6);
  for(int i = 0; i < 3; i++){
    if(i == axis) continue;
    output[i] = base[i] - radius - cutoff;
  }
  for(int i = 3; i < 6; i++){
    if(i == axis+3) continue;
    output[i] = base[i] + radius + cutoff;
  }
  output[axis] = base[axis] - cutoff;
  output[axis+3] = base[axis] + height + cutoff; 
  return output;
}

Vec<real> getBoxBounds(Vec3<real> box_min, Vec3<real> box_max, real cutoff){
  Vec<real> output(6);
  for(int i = 0; i < 3; i++){
    output[i] = box_min[i]-cutoff;
    output[i+3] = box_max[i]+cutoff;
  }
  return output;
}

Vec<real> getSphereBounds(Vec3<real> xref, real radius){
  Vec<real> output(6);
  for(int i = 0; i < 3; i++){
    output[i] = xref[i] - radius;
    output[i+3] = xref[i] + radius;
  }
  return output;
}

LJ_6_12_default::LJ_6_12_default(Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params) : FF_function(xref, size, cutoff, params){
  epsilon_ = params[0];
  sigma_ = params[1];
  return;
}

real LJ_6_12_default::compute(Vec3<real> x) const{
  getNearestImage3D(x, xref_, size_);
  auto r = norm2(x-xref_);
  if(r > cutoff_) return 0.0;
  return LJ_6_12(r, epsilon_, sigma_);
}

Vec<real> LJ_6_12_default::getBoundingBox() const{
  return getSphereBounds(xref_, cutoff_);
}

LJ_6_12_offset::LJ_6_12_offset(Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params) : FF_function(xref, size, cutoff, params){
  epsilon_ = params[0];
  sigma_ = params[1];
  offset_ = params[2];
  return;
}

real LJ_6_12_offset::compute(Vec3<real> x) const{
  //models hard sphere with lj component
  //params are epsilon, sigma, dr
  getNearestImage3D(x, xref_, size_);
  auto r = norm2(x-xref_) - offset_;
  if(r <= 0 ) return 100.0;
  if(r > cutoff_) return 0.0;
  return LJ_6_12(r, epsilon_, sigma_);
}

Vec<real> LJ_6_12_offset::getBoundingBox() const{
  return getSphereBounds(xref_, cutoff_ + offset_);
}

//probably not the most efficient algorithm, but it works...
LJ_6_12_offset_cylinder::LJ_6_12_offset_cylinder(Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params) : FF_function(xref, size, cutoff, params){
  epsilon_ = params[0];
  sigma_ = params[1];
  offset_ = params[2];
  height_ = params[3];
  axis_ = std::round(params[4]);
  return;
}

real LJ_6_12_offset_cylinder::compute(Vec3<real> x) const{
  //params are epsilon, sigma, r, h, and axis
  real eval = 0.0;
  real r2_cyl = offset_*offset_;
  real r2 = 0.0;
  for(int i = 0; i < 3; i++){
    if(i == axis_) continue;
    r2 += (x[i] - xref_[i])*(x[i] - xref_[i]);
  }
  //first, check if point is within cylinder
  if(r2 <= r2_cyl && x[axis_] - xref_[axis_] <= height_ && x[axis_] - xref_[axis_] >= 0) return 100.0;
  //if not, compute closest distance
  real r = std::sqrt(r2);
  real rd = r - offset_; //radial distance to radial extent of cylinder
  //if it's within the cylinder height, closest distance will be distance to exterior of cylinder
  if(x[axis_] - xref_[axis_] <= height_ && x[axis_] - xref_[axis_] >= 0){
    return LJ_6_12(rd, epsilon_, sigma_); //can't be < 0 because of previous condition
  }
  real vd2t = x[axis_] - xref_[axis_] - height_; //vertical distance to top circle
  real vd2b = x[axis_] - xref_[axis_]; //vertical distance to bottom circle
  real top_distance=0.0, bottom_distance=0.0;
  if(rd < 0){ //if it's within the circle radius, it's the distance along the axis
    top_distance = vd2t;
    bottom_distance = vd2b;
  }
  else{ //if it's outside the radius, it's the hypotenuse of the r and y component
    top_distance = std::sqrt( (vd2t*vd2t) + (rd*rd) );
    bottom_distance = std::sqrt( (vd2b*vd2b) + (rd*rd) );
  }
  real r_final = std::min(top_distance, bottom_distance);
  if(r_final > cutoff_) return 0.0; 
  return LJ_6_12(r_final, epsilon_, sigma_); //return the minimum of these values
}

Vec<real> LJ_6_12_offset_cylinder::getBoundingBox() const{
  return getCylindricalBounds(xref_, offset_, height_, axis_, cutoff_);
}

LJ_6_12_offset_box::LJ_6_12_offset_box(Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params) : FF_function(xref, size, cutoff, params){
  epsilon_ = params[0];
  sigma_ = params[1];
  offset_ = params[2];
  Vec3<real> dx;
  for(std::size_t i = 0; i < 3; i++){
    dx[i] = params[i+2];
  }
  xmin_ = xref_ - dx;
  xmax_ = xref_ + dx;
  return;
}

real LJ_6_12_offset_box::compute(Vec3<real> x) const{
  //models cuboidal volume with lj component computed based on nearest distance to surface
  //params are epsilon, sigma, dx, dy, dz
  auto r = nearestDistanceToBox(x, xref_, xmin_, xmax_, size_);
  if(r <= 0.0) return 100.0;
  if(r > cutoff_) return 0.0;
  return LJ_6_12(r, epsilon_, sigma_);
}

Vec<real> LJ_6_12_offset_box::getBoundingBox() const{
  return getBoxBounds(xmin_, xmax_, cutoff_);
}


/////3-9
////
////
LJ_3_9_offset::LJ_3_9_offset(Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params) : FF_function(xref, size, cutoff, params){
  epsilon_ = params[0];
  sigma_ = params[1];
  offset_ = params[2];
  return;
}

real LJ_3_9_offset::compute(Vec3<real> x) const{
  //models hard sphere with lj component
  //params are epsilon, sigma, dr
  getNearestImage3D(x, xref_, size_);
  auto r = norm2(x-xref_) - offset_;
  if(r <= 0 ) return 100.0;
  if(r > cutoff_) return 0.0;
  return LJ_3_9(r, epsilon_, sigma_);
}

Vec<real> LJ_3_9_offset::getBoundingBox() const{
  return getSphereBounds(xref_, cutoff_ + offset_);
}

//probably not the most efficient algorithm, but it works...
LJ_3_9_offset_cylinder::LJ_3_9_offset_cylinder(Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params) : FF_function(xref, size, cutoff, params){
  epsilon_ = params[0];
  sigma_ = params[1];
  offset_ = params[2];
  height_ = params[3];
  axis_ = std::round(params[4]);
  return;
}

real LJ_3_9_offset_cylinder::compute(Vec3<real> x) const{
  //params are epsilon, sigma, r, h, and axis
  real eval = 0.0;
  real r2_cyl = offset_*offset_;
  real r2 = 0.0;
  for(int i = 0; i < 3; i++){
    if(i == axis_) continue;
    r2 += (x[i] - xref_[i])*(x[i] - xref_[i]);
  }
  //first, check if point is within cylinder
  if(r2 <= r2_cyl && x[axis_] - xref_[axis_] <= height_ && x[axis_] - xref_[axis_] >= 0) return 100.0;
  //if not, compute closest distance
  real r = std::sqrt(r2);
  real rd = r - offset_; //radial distance to radial extent of cylinder
  //if it's within the cylinder height, closest distance will be distance to exterior of cylinder
  if(x[axis_] - xref_[axis_] <= height_ && x[axis_] - xref_[axis_] >= 0){
    return LJ_3_9(rd, epsilon_, sigma_); //can't be < 0 because of previous condition
  }
  real vd2t = x[axis_] - xref_[axis_] - height_; //vertical distance to top circle
  real vd2b = x[axis_] - xref_[axis_]; //vertical distance to bottom circle
  real top_distance=0.0, bottom_distance=0.0;
  if(rd < 0){ //if it's within the circle radius, it's the distance along the axis
    top_distance = vd2t;
    bottom_distance = vd2b;
  }
  else{ //if it's outside the radius, it's the hypotenuse of the r and y component
    top_distance = std::sqrt( (vd2t*vd2t) + (rd*rd) );
    bottom_distance = std::sqrt( (vd2b*vd2b) + (rd*rd) );
  }
  real r_final = std::min(top_distance, bottom_distance);
  if(r_final > cutoff_) return 0.0; 
  return LJ_3_9(r_final, epsilon_, sigma_); //return the minimum of these values
}

Vec<real> LJ_3_9_offset_cylinder::getBoundingBox() const{
  return getCylindricalBounds(xref_, offset_, height_, axis_, cutoff_);
}

LJ_3_9_offset_box::LJ_3_9_offset_box(Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params) : FF_function(xref, size, cutoff, params){
  epsilon_ = params[0];
  sigma_ = params[1];
  offset_ = params[2];
  Vec3<real> dx;
  for(std::size_t i = 0; i < 3; i++){
    dx[i] = params[i+2];
  }
  xmin_ = xref_ - dx;
  xmax_ = xref_ + dx;
  return;
}

real LJ_3_9_offset_box::compute(Vec3<real> x) const{
  //models cuboidal volume with lj component computed based on nearest distance to surface
  //params are epsilon, sigma, dx, dy, dz
  auto r = nearestDistanceToBox(x, xref_, xmin_, xmax_, size_);
  if(r <= 0.0) return 100.0;
  if(r > cutoff_) return 0.0;
  return LJ_3_9(r, epsilon_, sigma_);
}

Vec<real> LJ_3_9_offset_box::getBoundingBox() const{
  return getBoxBounds(xmin_, xmax_, cutoff_);
}
