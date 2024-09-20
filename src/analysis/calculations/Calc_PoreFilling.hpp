#pragma once
#include "Calculation.hpp"
#include "../probevolumes/PV_Sphere.hpp"
class Calc_PoreFilling : public Calculation{
public:
  Calc_PoreFilling(AnalysisInputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
  virtual void output();
  virtual void update();
protected:
  real computeN(PV_Sphere& probevolume);
  void loadVolumeFile();
  virtual void printOutput();
  std::ofstream ofile_;
  std::vector<PV_Sphere> probeVolumes_;
  Vec<real> n_frame_; //vector of water counts for each volume, will be compared to absolute count post-analysis

  //input
  std::string volfile_;
  bool bFileLoaded_=0;
  //output
  std::ofstream nOutput_;
  Vec<real> nVectorStep_, nVectorAverage_;
  real nVectorCount_; 
  bool b_nOutputOpened_=0;

};
