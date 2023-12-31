#include "Calculation.hpp"
class Calc_PhiFromLiqDens : public Calculation{
public:
  Calc_PhiFromLiqDens(AnalysisInputPack& input);
  virtual void calculate();
  virtual void finalOutput();
  virtual void output();
  virtual void update();
protected:
  real closestLiqSite(Vec3<std::size_t> dPos);
  int cutoff_;
  real real_cut_;
  real sigma_;
  real eps_;
  bool combine_phi_;
  Matrix3d<real> avg_;
  bool isInitialized = 0;
  int frame_counter_;
};