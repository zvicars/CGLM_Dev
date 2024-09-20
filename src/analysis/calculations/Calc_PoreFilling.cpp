#include "Calc_PoreFilling.hpp"
#include "../../tools/StringTools.hpp"
#include "../../tools/pbcfunctions.hpp"
#include "../../tools/stlmath.hpp"
#include <cmath>
Calc_PoreFilling::Calc_PoreFilling(AnalysisInputPack& input):Calculation{input}
{
  input.params().readString("volumes", ParameterPack::KeyType::Required, volfile_);
  return;
}
void Calc_PoreFilling::update(){
  Calculation::update();
  if(!bFileLoaded_){
    loadVolumeFile();
    bFileLoaded_ = 1;
    nVectorAverage_.resize(probeVolumes_.size(), 0.0);
  }
  return;
}

void Calc_PoreFilling::loadVolumeFile(){
  std::ifstream ifile(volfile_);
  FANCY_ASSERT(ifile.is_open(), "Failed to open input file " + volfile_ + "in Calculation " + name_);
  std::string line;
  Vec3<real> rBoxSize = {(real)box->size_[0], (real)box->size_[1], (real)box->size_[2]};
  while(getline(ifile, line)){
    std::string trLine = StringTools::trimWhitespace(line);
    if(trLine.at(0) == '#' || trLine.at(0) == '%' || trLine.at(0) == ';') continue;
    std::stringstream ss(trLine);
    double x, y, z, r;
    ss >> x >> y >> z >> r;
    Vec3<real> pos = {x,y,z};
    FANCY_ASSERT(!ss.fail(), "stringstream failed in calculation " + name_ + ", invalid line likely read from input volume file");
    probeVolumes_.emplace_back(PV_Sphere(pos, r, rBoxSize));
  }
  return;
}

real Calc_PoreFilling::computeN(PV_Sphere& pv){
  Vec<real> bounds = pv.getBounds();
  Vec2<int> xb, yb, zb;
  xb[0] = std::floor(bounds[0]); xb[1] = std::ceil(bounds[3]);
  yb[0] = std::floor(bounds[1]); yb[1] = std::ceil(bounds[4]);    
  zb[0] = std::floor(bounds[2]); zb[1] = std::ceil(bounds[5]);
  real count = 0.0;
  Vec3<int> box_size_int = box->getSizeInt();
  for(int i = xb[0]; i <= xb[1]; i++){
  for(int j = yb[0]; j <= yb[1]; j++){
  for(int k = zb[0]; k <= zb[1]; k++){
        Vec3<real> site_vec = {(real)i + 0.5, (real)j + 0.5, (real)k + 0.5};
        std::size_t index1 = wrapIndex(i, box_size_int[0]);
        std::size_t index2 = wrapIndex(j, box_size_int[1]);
        std::size_t index3 = wrapIndex(k, box_size_int[2]);
        Vec3<std::size_t> wrappedIdx = {index1, index2, index3};
        real val = internalLattice.at(wrappedIdx) * pv.compute(site_vec);
        count += val;
      }
    }
  }
  return count;
}

void Calc_PoreFilling::calculate(){
  if(!doCalculate()) return;
  std::size_t nVols = probeVolumes_.size();
  std::vector<real> calcEvals(nVols);
  for(std::size_t volIdx = 0; volIdx < nVols; volIdx++){
    auto& pv =  probeVolumes_[volIdx];
    calcEvals[volIdx] = computeN(pv);
  }
  nVectorStep_ = calcEvals;
  nVectorAverage_ = nVectorAverage_ + calcEvals;
  nVectorCount_ += 1.0;
  return;
}

void Calc_PoreFilling::output(){
  if(doOutput()){
    printOutput();
  } 
  return;
}

std::string Calc_PoreFilling::printConsoleReport(){
  return "";
}

void Calc_PoreFilling::printOutput(){
  if(!b_nOutputOpened_){
    nOutput_.open(name_ + "_" + type_ + ".txt");
    nOutput_ <<  "#num volumes = " << probeVolumes_.size() << "\n";
    nOutput_ << "# frame  n_i" << std::endl;
    b_nOutputOpened_ = 1;
  }
  nOutput_ << box->frame_ << "  "; 
  for(int i = 0; i < probeVolumes_.size(); i++){
    nOutput_ << nVectorStep_[i] << "   ";
  }
  nOutput_ << "\n";
  return;
};
void Calc_PoreFilling::finalOutput(){

  nOutput_ << "Av: "; 
  for(int i = 0; i < probeVolumes_.size(); i++){
    nOutput_ << nVectorAverage_[i] / nVectorCount_ << "   ";
  }  
  nOutput_.close();
  return;
};