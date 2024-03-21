#pragma once
#include "Calculation.hpp"
class Calc_WriteTraj : public Calculation{
public:
  Calc_WriteTraj(AnalysisInputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
  virtual void output();
  virtual void update();
protected:
  virtual void printOutput();
  std::ofstream ofile;
};
