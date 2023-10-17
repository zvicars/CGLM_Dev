#include "Calc_Isosurface.hpp"
Calc_Isosurface::Calc_Isosurface(AnalysisInputPack& input):Calculation{input}
{
  frame_counter_ = 0;
  input.params().readNumber("isovalue", KeyType::Required, isovalue_);
  FANCY_ASSERT(isovalue_ > 0, "Invalid isovalue given for isosurface calculation.");
  method_ = "golosio";
  input.params().readString("method", KeyType::Optional, method_);
  FANCY_ASSERT(method_ == "golosio", "Invalid method chosen for instantaneous interface calculation, valid options are \'golosio\' and \'rchandra\'.");
  std::string pvname;
  haspv_ = input.params().readString("probevolume", KeyType::Optional, pvname);
  if(haspv_){
    pv_ = input.findProbeVolume(pvname);
    FANCY_ASSERT(pv_ != 0, "probe volume not found");
  }
  input.params().readFlag("compute_curvature", ParameterPack::KeyType::Optional, computeCurvature_);
  return;
}
void Calc_Isosurface::update(){
  Calculation::update();
  size_ = box->getSizeInt();
  if(!initialized_){
    average_.initialize(size_, isovalue_);
    frame_.initialize(size_, isovalue_);
    initialized_ = 1;
  }
  else{
    frame_.clear();
  }
  return;
}
void Calc_Isosurface::calculate(){
  if(!doCalculate()) return;
  //for lattice model, actually have guarantee of no data races
  int nsites = box->lattice_.size_1d();
  #pragma omp parallel for
  for(int i = 0; i < nsites; i++ ){
    auto dPos = internalLattice.map1N(i);
    Vec3<real> realPos;
    for(int j = 0; j < 3; j++){
      realPos[j] = (real)dPos[j];
    }
    if(haspv_ && pv_->compute(realPos) == 0.0) continue;
    real lattice_value = internalLattice.read_at_1d(i);
    frame_.add_discrete(dPos, lattice_value);
  }
  marchingCubes(method_, frame_, mesh_);
  average_.sumInPlace(frame_);
  frame_counter_++;
  return;
}

void Calc_Isosurface::output(){
  if(doOutput()){
    printOutput();
  }  
}

std::string Calc_Isosurface::printConsoleReport(){
  return "";
}

void Calc_Isosurface::printOutput(){
  std::string filepath = name_ + "_interface"+ ".stla";
  std::ofstream ofile(filepath, std::ofstream::app);
  FANCY_ASSERT(ofile.is_open(), "Failed to open output file for instantaneous interface step calculation.");
  std::string output;
  printSTL(mesh_, output);
  ofile << output;
  ofile.close();
};
void Calc_Isosurface::finalOutput(){
  //make the final mesh no matter what
  average_.scalarMult(1.0/(double)frame_counter_);
  marchingCubes(method_, average_, mesh_);
  std::string avgstlname = name_ + "_average.stl";
  std::string output;
  printSTL(mesh_, output);
  std::ofstream ofile2(avgstlname);
  ofile2 << output;
  ofile2.close();
  final_output_flag_ = 1;
  if(computeCurvature_){
    std::string filepath = name_ + "_average.ply";
    std::ofstream ofile(filepath);
    FANCY_ASSERT(ofile.is_open(), "Failed to open output file for instantaneous interface average.");
    std::string output;
    std::vector<std::array<double, 2> > curv;
    curv = computeMeshCurvature(mesh_, 3);
    std::vector<double> ac(curv.size()), gc(curv.size());
    for(int i = 0; i < curv.size(); i++){
      ac[i] = 0.5*(curv[i][0]+curv[i][1]);
      gc[i] = curv[i][0]*curv[i][1];
    }
    std::string ac_info = printPLYWithCurvature(mesh_, ac); 
    std::string gc_info = printPLYWithCurvature(mesh_, gc);   
    ofile << ac_info;
    ofile.close(); 
    filepath = name_ + "_gaussian.ply";
    ofile.open(filepath);
    FANCY_ASSERT(ofile.is_open(), "Failed to open output file for instantaneous interface average.");
    ofile << gc_info;
    ofile.close();
  }
  return;
};