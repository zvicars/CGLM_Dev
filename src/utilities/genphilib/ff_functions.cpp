#include "ff_functions.hpp"
#include "../../tools/SimpleForcefields.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/pbcfunctions.hpp"
#include <iostream>
real LJ_6_12(Vec3<real> x, Vec3<real> xref, Vec3<real> size, Vec<real> params){
  getNearestImage3D(x, xref, size);
  auto r = norm2(x-xref);
  return LJ_6_12(r, params[0], params[1]);
}
real LJ_6_12_shiftscale(Vec3<real> x1, Vec3<real> x2, Vec3<real> size, Vec<real> params){
  return 0.0;
}
real LJ_6_12_offset(Vec3<real> x, Vec3<real> xref, Vec3<real> size, Vec<real> params){
  //models hard sphere with lj component
  //params are epsilon, sigma, dr
  getNearestImage3D(x, xref, size);
  auto r = norm2(x-xref) - params[2];
  if(r <= 0 ) return 100.0;
  return LJ_6_12(r, params[0], params[1]);
}

//https://stackoverflow.com/questions/5254838/calculating-distance-between-a-point-and-a-rectangular-box-nearest-point
real nearestDistanceToBox(Vec3<real> p, Vec3<real> min, Vec3<real> max, Vec3<real> size){
  Vec3<real> pos_distances = {0.0, 0.0, 0.0};
  for(int i = 0; i < 3; i++){
    //find nearest boundary
    if(p[i] >= min[i] && p[i] <= max[i]) continue; //it's inside, distance is 0
    real mindist = min[i] - getNearestImage1D(p[i], min[i], size[i]);
    real maxdist = getNearestImage1D(p[i], max[i], size[i]) - max[i];
    if(mindist > 0.0 || maxdist > 0.0){
      if(mindist < 0.0) pos_distances[i] = maxdist;
      else pos_distances[i] = mindist;
    }
  }
  real sum = 0.0;
  for(int i = 0; i < 3; i++){
    //std::cout << "xdist = " << pos_distances[0] << " ydist = " << pos_distances[1] << " zdist = " << pos_distances[2] << std::endl;
    sum += pos_distances[i]*pos_distances[i];
  }
  sum = std::sqrt(sum);
  return sum;
}

real LJ_Box(Vec3<real> x, Vec3<real> xref, Vec3<real> size, Vec<real> params){
  //models cuboidal volume with lj component computed based on nearest distance to surface
  //params are epsilon, sigma, dx, dy, dz
  Vec3<real> dx;
  dx[0] = params[2];
  dx[1] = params[3];
  dx[2] = params[4];
  auto xmin = xref - dx;
  auto xmax = xref + dx;
  for(int i = 0; i < 3; i++){
    if(xmin[i] < 0 ) xmin[i] = 0;
    if(xmax[i] > size[i]) xmax[i] = size[i];
  }
  //get shortest distance to 
  auto r = nearestDistanceToBox(x, xmin, xmax, size);
  //std::cout << x[0] << "  " << x[1] << "  " << x[2] << std::endl;
  //std::cout << dx[0] << "  " << dx[1] << "  " << dx[2] << std::endl;
  //std::cout << xmin[0] << "  " << xmin[1] << "  " << xmin[2] << std::endl;
  //std::cout << xmax[0] << "  " << xmax[1] << "  " << xmax[2] << std::endl;
  //std::cin.get();
  if(r == 0.0 ) return 100.0;
  return LJ_6_12(r, params[0], params[1]);
}