#pragma once
#include <array>
//solution adapted from https://stackoverflow.com/questions/707370/clean-efficient-algorithm-for-wrapping-integers-in-c
//user MartinStettner

static inline int wrap(int i, int lb, int ub){
  int range = ub - lb + 1;
  i = ((i-lb) % range);
  if (i<0)
    return ub + 1 + i;
  else
    return lb + i;
  return i;
}

//used to convert int index into pbc-wrapped size_t index
static inline std::array<std::size_t,3> wrap3(std::array<int,3> in, std::array<std::size_t, 3> ub){
  std::array<std::size_t, 3> out;
  #pragma unroll
  for(int i = 0; i < 3; i++){
    out[i] = wrap(in[i], 0, ub[i]-1);
  }
  return out;
}