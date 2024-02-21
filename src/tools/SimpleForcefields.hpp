//some simple forcefields, LJ, Couloumb interactions, etc.
//Zachariah Vicars, 7-7-2022
#pragma once
#include <cmath>
#include "typedefs.hpp"
static inline real LJ_6_12(real r, real epsilon, real sigma){
  real sigma6 = std::pow(sigma/r, 6.0);
  return -4.0*epsilon*(sigma6 - sigma6*sigma6);
}

static inline real LJ_3_9(real r, real epsilon, real sigma){
  //2.598076211 = 3sqrt(3)/2 http://www.sklogwiki.org/SklogWiki/index.php/9-3_Lennard-Jones_potential
  real sigma3 = std::pow(sigma/r, 3.0);
  return 2.598076211*epsilon*(sigma3*sigma3*sigma3 - sigma3);
}