#pragma once
#include "Calculation.hpp"
#include "Calc_Isosurface.hpp"
#include "../../tools/Matrix.hpp"
#include "../../tools/smearfunctions.hpp"
class Calc_IsosurfaceMultiphase : public Calculation{
public:
  Calc_IsosurfaceMultiphase(AnalysisInputPack& input);
  virtual void calculate();
  virtual void finalOutput();
  virtual void output();
  virtual void update();  
protected:
  Vec3<std::size_t> box_size_;
  int axis_, axis1_, axis2_;
  real min_val_=0; 
  Vec<Calc_Isosurface*> isosurface_calculations_;
  //minimum of ray-triangle intersection from +axis
  //pretty complicated way of estimating the area
  //may wait on this
  //Vec< Matrix<real, 2> > areas;
  bool initialized_=0;
  int num_groups_ = 0, frame_counter_=0;
  std::vector<double> areas_, distance_rmax_, distance_sigmas_, num_neighbor_thresholds_, num_neighbor_sigmas_;
  std::vector<Mesh> meshes_;
  std::vector<std::string> filepaths_;
  std::vector<std::string> hex_colors_;
  std::vector<double> times_;
  std::vector<std::vector<double> > avecs_;
  bool cons_areas_ = 0;
  bool outputPLY_ = 0;
};