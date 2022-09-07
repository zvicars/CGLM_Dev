#pragma once
#include "bias.hpp"
#include "../pv/probevolume.hpp"
#include "bias_functions.hpp"
class Bias_Nv : public Bias{
  public:
  Bias_Nv(InputPack& input);
  ~Bias_Nv(){
    delete pv_;
    delete bf_;
  };
  virtual real calc_u(const Lattice& lattice){
    u_ = (*bf_)(pv_->calc_nv(lattice));
    return u_;
  }
  //calc_du function call needs calc_u to have been called at some point before
  virtual real calc_du(const Lattice& lattice){
    du_ = (*bf_)(pv_->nv(), pv_->calc_dnv(lattice));
    return du_;
  }
  virtual void flip(){
    Bias::flip();
    pv_->flip();
    u_ = u_ + du_;
    du_ = 0;
    return;
  }
  virtual Bias* clone(){
    return new Bias_Nv(*input_);
  }
  virtual std::string printStepOutput(std::string s);
  private:
  ProbeVolume* pv_;
  BiasFunction* bf_;
};