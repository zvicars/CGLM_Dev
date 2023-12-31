//depends on two isosurface calculations, one using the nanoparticles and one using the waters
#pragma once
#include "Calculation.hpp"
#include "Calc_Isosurface.hpp"
#include "../../tools/Matrix.hpp"
class Calc_LeeNP : public Calculation{
public:
  Calc_LeeNP(AnalysisInputPack& input);
  virtual void calculate();
  virtual void finalOutput();
  virtual void output();
  virtual void update();  
protected:
  real findNearestIntersection(Vec3<real> p0, Vec3<real> n0, const Mesh& mesh) const;
  bool rayTriangleIntersection(Vec3<real> p0, Vec3<real> n0, Vec3<Vec3<real> > triangle, Vec3<real>& intersection_point) const;
  void writeImage(int width, int height, Matrix<unsigned char, 2>& r, Matrix<unsigned char, 2>& g, Matrix<unsigned char, 2>& b, std::string filename);
  int axis_, axis1_, axis2_;
  real min_val_=0; 
  Vec<Calc_Isosurface*> isosurface_calculations_;
  //minimum of ray-triangle intersection from +axis
  Vec< Matrix<real, 2> > distances_;
  //pretty complicated way of estimating the area
  //may wait on this
  //Vec< Matrix<real, 2> > areas;
  bool initialized=0;
};