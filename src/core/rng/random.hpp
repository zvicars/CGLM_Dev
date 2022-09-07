//random number generator parent class
//performs a simple role, generates random numbers for selecting lattice sites
#pragma once
#include "../object.hpp"
#include <random>
class RNG : public Object{
  public:
  RNG(InputPack& input):Object{input}{
    return;
  }
  virtual size_t getIndex(size_t min, size_t max) = 0;
  virtual real getReal(real min, real max) = 0;
  virtual RNG* clone() = 0;
};

class RNG_mt19937 : public RNG{
  public:
  RNG_mt19937(InputPack& input);
  virtual size_t getIndex(size_t min, size_t max);
  virtual real getReal(real min, real max);
  virtual RNG* clone(){
    return new RNG_mt19937(*input_);
  }
  private:
  //returns an unsigned int, will need to cast to size_t
  std::mt19937 generator;
};
