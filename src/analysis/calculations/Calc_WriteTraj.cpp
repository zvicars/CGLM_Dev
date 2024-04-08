#include "Calc_WriteTraj.hpp"
#include "../../tools/CGLMFileHelper.hpp"
Calc_WriteTraj::Calc_WriteTraj(AnalysisInputPack& input):Calculation{input}
{
  std::string filepath = name_ + "_out" + ".traj";
  ofile_.open(filepath, std::ios::binary);
  FANCY_ASSERT(ofile_.is_open(), "Failed to open output file.");
  return;
}
void Calc_WriteTraj::update(){
  Calculation::update();
  return;
}
void Calc_WriteTraj::calculate(){
  if(!doCalculate()) return;
  return;
}

void Calc_WriteTraj::output(){
  if(doOutput()){
    printOutput();
  } 
  return;
}

std::string Calc_WriteTraj::printConsoleReport(){
  return "";
}

void Calc_WriteTraj::printOutput(){
  int numel = internalLattice.size_1d();
  auto size3d = internalLattice.size();
  Vec<char> values(numel);
  for(int i = 0; i < values.size(); i++){
    auto iVal = internalLattice.at_1d(i);
    if(iVal < 0.5) values[i] = '0';
    else values[i] = '1';
  }
  binary_char_write(ofile_, values, size3d);
};
void Calc_WriteTraj::finalOutput(){
  ofile_.close();
  return;
};