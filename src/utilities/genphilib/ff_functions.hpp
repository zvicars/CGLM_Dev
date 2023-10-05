#pragma once
#include "../../tools/typedefs.hpp"
//real (*funct_ptr)(Vec3<real>, Vec3<real>, Vec<real>)
//params: epsilon, sigma, cutoff
real LJ_6_12(Vec3<real> x, Vec3<real> x_ref, Vec3<real> size, Vec<real> params);
real LJ_6_12_shiftscale(Vec3<real> x, Vec3<real> x_ref, Vec3<real> size, Vec<real> params);
real LJ_6_12_offset(Vec3<real> x, Vec3<real> xref, Vec3<real> size, Vec<real> params);
real LJ_Box(Vec3<real> x, Vec3<real> xref, Vec3<real> size, Vec<real> params);