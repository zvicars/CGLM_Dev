#pragma once
#include "Calculation.hpp"
#include "../../external/cubes/MarchingCubesInterface.hpp"
class Calc_Isosurface : public Calculation{
public:
  Calc_Isosurface(AnalysisInputPack& input);
  virtual void calculate();
  virtual std::string printConsoleReport();
  virtual void finalOutput();
  virtual void output();
  virtual void update();
  Mesh& getMesh(){
    return mesh_;
  }
protected:
  virtual void printOutput();
  //number of actual frames computed
  int frame_counter_=0;
  bool initialized_ = 0;
  std::string method_;
  //total surface area of the isosurface
  std::vector<double> areas_;
  //voxel grid input information
  Vec3<int> npoints_;
  double area_, isovalue_;
  //voxel grid output
  Mesh mesh_;
  VoxelGrid frame_;
  VoxelGrid average_;
  Vec3<int> size_;
  ProbeVolume* pv_;
  bool haspv_=0, computeCurvature_=0;
};
