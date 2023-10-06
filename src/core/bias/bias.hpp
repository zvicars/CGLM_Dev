#pragma once
#include "../object.hpp"
class Bias : public Object{
  public:
  Bias(InputPack& input):Object{input}{
    return;
  }
  virtual real calc_u(const Lattice& lattice)=0;
  virtual real calc_du(const Lattice& lattice)=0;
  real u(){
    return u_;
  }
  real du(){
    return du_;
  }
  virtual void flip(){
    u_ += du_;
    du_ = 0;
    return;
  }
  virtual void sweepUpdate(real time){
    return;
  }
  virtual Bias* clone()=0;
  protected:
  //value for u so there's some updated value of u that updates with du
  real u_;
  //value for du to apply if the site flip is successful
  real du_;
};
