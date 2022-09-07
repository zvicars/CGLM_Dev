#include "ff_functions.hpp"
#include "../../tools/SimpleForcefields.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/pbcfunctions.hpp"
real LJ_6_12(Vec3<real> x, Vec3<real> xref, Vec3<real> size, Vec<real> params){
  getNearestImage3D(x, xref, size);
  auto r = norm2(x-xref);
  return LJ_6_12(r, params[0], params[1]);
}
real LJ_6_12_shiftscale(Vec3<real> x1, Vec3<real> x2, Vec3<real> size, Vec<real> params){
  return 0.0;
}