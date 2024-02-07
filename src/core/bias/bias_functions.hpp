#pragma once
#include "../../tools/InputParser.hpp"
#include "../../typedefs.hpp"
#include "ramped_parameter.h"
class BiasFunction{
  public:
    virtual real operator()(real n) = 0;
    virtual real operator()(real n, real dn) = 0;
    virtual void sweepUpdate(real time){
      return;
    }
};
class HarmonicPotential : public BiasFunction{
  public:
    HarmonicPotential(const ParameterPack& input){
      std::vector<double> kr, xr;
      dokappaRamp_ = input.readVector("ramp_kappa", ParameterPack::KeyType::Optional, kr);
      if(dokappaRamp_){
        FANCY_ASSERT(kr.size()%2 == 0, "kappa ramping needs to have an even number of elements in the vector");
      }
      else{
        input.readNumber("kappa", ParameterPack::KeyType::Required, kappa_);
        kr = { 0, 1, kappa_, kappa_ };
      }
      doxstarRamp = input.readVector("ramp_xstar", ParameterPack::KeyType::Optional, xr);
      std::cout << doxstarRamp << std::endl; 
      if(doxstarRamp){
        FANCY_ASSERT(xr.size()%2 == 0, "xstar ramping needs to have an even number of elements in the vector");
      }
      else{
        input.readNumber("xstar", ParameterPack::KeyType::Required, xstar_);
        xr = { 0, 1, xstar_, xstar_};
      }
      xstar_ramp_.initialize(xr);
      kappa_ramp_.initialize(kr);
      sweepUpdate(0.0);
      return;
    }
    HarmonicPotential(real kappa, real xstar){
      xstar_ =  xstar;
      kappa_ = kappa;
      halfkappa_ = kappa*0.5;
    }

    virtual void sweepUpdate(real time){
      xstar_ = xstar_ramp_.compute(time);
      kappa_ = kappa_ramp_.compute(time);
      halfkappa_ = 0.5*kappa_;
      return;
    }

    //u
    virtual real operator()(real x){
      return halfkappa_*(x-xstar_)*(x-xstar_);
    }
    //u2-u1 for a change of dn
    virtual real operator()(real x, real dx){
      return kappa_*dx*(x - xstar_ + 0.5*dx);
    }

  private:
    real xstar_, kappa_, halfkappa_;
    RampedParameter xstar_ramp_, kappa_ramp_;
    bool dokappaRamp_=0, doxstarRamp=0;
};

class LinearPotential : public BiasFunction{
  public:
    LinearPotential(const ParameterPack& input){
      input.readNumber("phi", ParameterPack::KeyType::Required, phi_);
      return;
    }
    LinearPotential(real phi){
      phi_ = phi;
    }
    //u
    virtual real operator()(real x){
      return phi_*x;
    }
    //u2-u1 for a change of dn
    virtual real operator()(real x, real dx){
      return phi_*dx;
    }
  private:
    real phi_;
};

class HarmonicPlusLinear : public BiasFunction{
  public:
    HarmonicPlusLinear(const ParameterPack& input) : p1{input}, p2{input}{
      return;
    }
    HarmonicPlusLinear(real phi, real kappa, real xstar) : p1{phi}, p2{kappa,xstar}{
      return;
    }
    virtual real operator()(real x){
      return p1(x) + p2(x);
    }
    //u2-u1 for a change of dn
    virtual real operator()(real x, real dx){
      return p1(x, dx) + p2(x,dx);
    }    
  protected:
  LinearPotential p1;
  HarmonicPotential p2;
};

static BiasFunction* bias_function_factory(const ParameterPack& params){
  std::string type;
  params.readString("type", ParameterPack::KeyType::Required, type);
  if(type == "harmonic") return new HarmonicPotential(params);
  if(type == "linear") return new LinearPotential(params);
  if(type == "linearharmonic") return new HarmonicPlusLinear(params);
  throw 999; //just being lazy here
  return 0;
}