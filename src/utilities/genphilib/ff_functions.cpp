#include "ff_functions.hpp"
#include "../../tools/SimpleForcefields.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/pbcfunctions.hpp"
#include <iostream>
real LJ_6_12_default(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params){
  getNearestImage3D(x, xref, size);
  auto r = norm2(x-xref);
  if(r > cutoff) return 0.0;
  return LJ_6_12(r, params[0], params[1]);
}
real LJ_6_12_shiftscale(Vec3<real> x1, Vec3<real> x2, Vec3<real> size, real cutoff, Vec<real> params){
  return 0.0;
}
real LJ_6_12_offset(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params){
  //models hard sphere with lj component
  //params are epsilon, sigma, dr
  getNearestImage3D(x, xref, size);
  auto r = norm2(x-xref) - params[2];
  if(r <= 0 ) return 100.0;
  if(r > cutoff) return 0.0;
  return LJ_6_12(r, params[0], params[1]);
}

//probably not the most efficient algorithm, but it works...
real LJ_6_12_offset_cylinder(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params){
  //params are epsilon, sigma, r, h, and axis
  real eval = 0.0;
  int axis = std::round(params[4]);
  real epsilon = params[0], sigma = params[1], r_cyl = params[2], r2_cyl = params[2]*params[2], h = params[3];
  real r2 = 0.0;
  for(int i = 0; i < 3; i++){
    if(i == axis) continue;
    r2 += (x[i] - xref[i])*(x[i] - xref[i]);
  }
  //first, check if point is within cylinder
  if(r2 <= r2_cyl && x[axis] - xref[axis] <= h && x[axis] - xref[axis] >= 0) return 100.0;
  //if not, compute closest distance
  real r = std::sqrt(r2);
  real rd = r - r_cyl; //radial distance to radial extent of cylinder
  //if it's within the cylinder height, closest distance will be distance to exterior of cylinder
  if(x[axis] - xref[axis] <= h && x[axis] - xref[axis] >= 0){
    return LJ_6_12(rd, epsilon, sigma); //can't be < 0 because of previous condition
  }
  real vd2t = x[axis] - xref[axis] - h; //vertical distance to top circle
  real vd2b = x[axis] - xref[axis]; //vertical distance to bottom circle
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
  if(r_final > cutoff) return 0.0; 
  return LJ_6_12(r_final, epsilon, sigma); //return the minimum of these values
}

//https://stackoverflow.com/questions/5254838/calculating-distance-between-a-point-and-a-rectangular-box-nearest-point
real nearestDistanceToBox(Vec3<real> p, Vec3<real> xref, Vec3<real> min, Vec3<real> max, Vec3<real> size){
  Vec3<real> dx;
  getNearestImage3D(p, xref, size);
  for(int i = 0; i < 3; i++){
    dx[i] = std::max(min[i]-p[i], std::max(0.0, p[i]-max[i]));
  }
  return norm2(dx);
}

real LJ_6_12_offset_box(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params){
  //models cuboidal volume with lj component computed based on nearest distance to surface
  //params are epsilon, sigma, dx, dy, dz
  Vec3<real> dx;
  for(std::size_t i = 0; i < 3; i++){
    dx[i] = params[i+2];
  }
  auto xmin = xref - dx;
  auto xmax = xref + dx;
  //for(int i = 0; i < 3; i++){
  //  if(xmin[i] < 0 ) xmin[i] = 0;
  //  if(xmax[i] > size[i]) xmax[i] = size[i];
  //}
  //get shortest distance to 
  auto r = nearestDistanceToBox(x, xref, xmin, xmax, size);
  if(r <= 0.0) return 100.0;
  if(r > cutoff) return 0.0;
  return LJ_6_12(r, params[0], params[1]);
}

real LJ_3_9_offset_box(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params){
  //models cuboidal volume with lj component computed based on nearest distance to surface
  //params are epsilon, sigma, dx, dy, dz
  Vec3<real> dx;
  real epsilon = params[0], sigma = params[1];
  for(std::size_t i = 0; i < 3; i++){
    dx[i] = params[i+2];
  }
  auto xmin = xref - dx;
  auto xmax = xref + dx;
  auto r = nearestDistanceToBox(x, xref, xmin, xmax, size);
  if(r <= 0.0) return 100.0;
  if(r > cutoff) return 0.0;
  return LJ_3_9(r, epsilon, sigma);
}

real LJ_3_9_offset_cylinder(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params){
  //params are epsilon, sigma, r, h, and axis
  real eval = 0.0;
  int axis = std::round(params[4]);
  real epsilon = params[0], sigma = params[1], r_cyl = params[2], r2_cyl = params[2]*params[2], h = params[3];
  real r2 = 0.0;
  for(int i = 0; i < 3; i++){
    if(i == axis) continue;
    r2 += (x[i] - xref[i])*(x[i] - xref[i]);
  }
  //first, check if point is within cylinder
  if(r2 <= r2_cyl && x[axis] - xref[axis] <= h && x[axis] - xref[axis] >= 0) return 100.0;
  //if not, compute closest distance
  real r = std::sqrt(r2);
  real rd = r - r_cyl; //radial distance to radial extent of cylinder
  //if it's within the cylinder height, closest distance will be distance to exterior of cylinder
  if(x[axis] - xref[axis] <= h && x[axis] - xref[axis] >= 0){
    return LJ_3_9(rd, epsilon, sigma); //can't be < 0 because of previous condition
  }
  real vd2t = x[axis] - xref[axis] - h; //vertical distance to top circle
  real vd2b = x[axis] - xref[axis]; //vertical distance to bottom circle
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
  if(r_final > cutoff) return 0.0; 
  return LJ_3_9(r_final, epsilon, sigma); //return the minimum of these values
}

real LJ_3_9_offset(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params){
  //models hard sphere with lj component
  //params are epsilon, sigma, dr
  getNearestImage3D(x, xref, size);
  auto r = norm2(x-xref) - params[2];
  if(r <= 0 ) return 100.0;
  if(r > cutoff) return 0.0;
  return LJ_3_9(r, params[0], params[1]);
}