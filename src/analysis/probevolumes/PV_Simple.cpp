#include "PV_Simple.hpp"
PV_DiscreteRect::PV_DiscreteRect(AnalysisInputPack& input):ProbeVolume{input}
{
  Vec<double> x_range, y_range, z_range;
  input.params().readVector("x_range", KeyType::Required, x_range);
  input.params().readVector("y_range", KeyType::Required, y_range);
  input.params().readVector("z_range", KeyType::Required, z_range);
  FANCY_ASSERT(x_range.size() == 2, "Invalid x-range given for Simple Rectangular PV.");
  FANCY_ASSERT(y_range.size() == 2, "Invalid y-range given for Simple Rectangular PV.");
  FANCY_ASSERT(z_range.size() == 2, "Invalid z-range given for Simple Rectangular PV.");
  for(int i = 0; i < 2; i++){
    axis_ranges_[i*3] = x_range[i];
    axis_ranges_[i*3+1] = y_range[i];
    axis_ranges_[i*3+2] = z_range[i];
  }
  return;
}

double PV_DiscreteRect::compute(Vec3<double> position){
  for(int i = 0; i < 3; i++){
    if( position[i] < axis_ranges_[i] || position[i] > axis_ranges_[i+3] ) return 0.0; 
  }
  return 1.0;
}