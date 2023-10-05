#pragma once
#include <array>
#include <vector>
#include <string>
struct SimpleAtom{
  std::array<double,3> x;
  std::array<double,3> v;
  std::string name, resname; 
  bool hasPosition = 0;
  bool hasVelocity = 0;
};

std::vector<SimpleAtom> simpleReadGro(std::string filename);