#include "random.hpp"
#include <ctime>
RNG_mt19937::RNG_mt19937(InputPack& input):RNG{input}{
  int seed=0;
  bool seed_found = input.params().readNumber("seed", ParameterPack::KeyType::Optional, seed);
  if(!seed_found) seed = time(NULL);
  generator.seed(seed);
  return;
}

size_t RNG_mt19937::getIndex(size_t min, size_t max){
  std::uniform_int_distribution<std::size_t> dist(min, max);
  return dist(generator);
}

real RNG_mt19937::getReal(real min, real max){
  std::uniform_real_distribution<real> dist(min, max);
  return dist(generator);
}