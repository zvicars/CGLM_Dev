#pragma once
#include "../../tools/Assert.hpp"
#include "../../tools/StringTools.hpp"
#include "../../typedefs.hpp"
#include "../AnalysisInputPack.hpp"
#include "../Analysis.hpp"

class Modifier : public AnalysisObject{
public:
  Modifier(AnalysisInputPack& input);
  //takes in a lattice, calculates some quantity, overrides the original lattice
  virtual void compute(Matrix3d<real>& lattice) = 0;
protected:
  //needs read access to objects inside box, even if it modifies a different copy of the lattice
  const Box* box;
};