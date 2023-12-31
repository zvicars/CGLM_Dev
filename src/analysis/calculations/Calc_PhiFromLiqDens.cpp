#include "Calc_PhiFromLiqDens.hpp"
#include "../../tools/pbcfunctions.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/SimpleForcefields.hpp"
#include "../../utilities/genphilib/filehandler.hpp"
Calc_PhiFromLiqDens::Calc_PhiFromLiqDens(AnalysisInputPack& input) : Calculation{input}{
  cutoff_ = 6; //check all sites within 6 LU of the reference point by default (0.184*6 = 1.104 nm)
  input.params().readNumber("cutoff", ParameterPack::KeyType::Optional, cutoff_);
  real_cut_ = (real)cutoff_; 
  sigma_ = 1.630; //0.3 nm / 0.184 nm for default args
  input.params().readNumber("well_dist", ParameterPack::KeyType::Optional, sigma_); //effectively a Lennard Jones sigma, but the way the function computes is strange
  eps_ = 0.2;
  input.params().readNumber("well_depth", ParameterPack::KeyType::Optional, eps_); //effectively a LJ eps, but only using r associated with closest liquid-like site, no prefactor
  combine_phi_ = 0;
  input.params().readFlag("combine_phi", ParameterPack::KeyType::Optional, combine_phi_); //merge it with the supplied phi file
  return;
}

void Calc_PhiFromLiqDens::calculate(){
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
    real lattice_value = internalLattice.read_at_1d(i);
  }
  avg_.add_inPlace(internalLattice);
  frame_counter_++;
  return;
}

void Calc_PhiFromLiqDens::finalOutput(){
  //loop over all cells within cutoff range
  avg_.div_inPlace((real)frame_counter_);
  Matrix3d<real> newLattice(avg_.size());
  newLattice.fill(0.0);
  int nsites = box->lattice_.size_1d();
  #pragma omp parallel for
  for(int i = 0; i < nsites; i++ ){
    if(avg_.at_1d(i) > 0.5){
        newLattice.at_1d(i) = 100.0;
        continue;
    }
    auto dPos = avg_.map1N(i);
    real r = closestLiqSite(dPos);
    if( r > real_cut_ ){
        newLattice.at_1d(i) = 0.0;
        if(combine_phi_) newLattice.at_1d(i) += box->phi_.read_at_1d(i);
        continue;
    }
    real potential = LJ_6_12(r, eps_, sigma_); 
    if(potential > 10 ) potential = 100;
    newLattice.at_1d(i) = potential;
    if(combine_phi_) newLattice.at_1d(i) += box->phi_.read_at_1d(i);
  }
  std::string filepath = name_ + ".phi";
  writePhiField(filepath, newLattice);
  return;
}

real Calc_PhiFromLiqDens::closestLiqSite(Vec3<std::size_t> pos){
    real r = 1e10;
    Vec3<real> realPosRef;
    for(int dim = 0; dim < 3; dim++){
        realPosRef[dim] = pos[dim];
    }
    auto size = avg_.size();
    for(int i = (int)pos[0] - cutoff_; i <= (int)pos[0] + cutoff_; i++){
        for(int j = (int)pos[1] - cutoff_; j <= (int)pos[1] + cutoff_; j++){
            for(int k = (int)pos[2] - cutoff_; k <= (int)pos[2] + cutoff_; k++){
                Vec3<real> realPosInner  = {(real) i,(real) j,(real) k};
                Vec3<std::size_t> wrappedPos = {(std::size_t) wrapIndex(i, (int)size[0]), (std::size_t) wrapIndex(j, (int)size[1]) , (std::size_t) wrapIndex(k, (int)size[2]) };
                if(avg_.at(wrappedPos) < 0.5) continue;
                real r_inner = norm2(realPosRef - realPosInner);
                if(r_inner < r) r = r_inner;
            }
        }
    }
    return r;
}

void Calc_PhiFromLiqDens::output(){
  return;
}

void Calc_PhiFromLiqDens::update(){
  Calculation::update();
  if(!isInitialized){
    avg_.initialize(box->size_);
    avg_.fill(0.0);
    frame_counter_=0;
    isInitialized=1;
  }
  return;  
}