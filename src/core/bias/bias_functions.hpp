#pragma once
#include "../../tools/InputParser.hpp"
#include "../../typedefs.hpp"
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
        FANCY_ASSERT(kr.size() == 4, "kappa ramping needs [ t_start(sweeps) t_end(sweeps) val_start val_end ]");
      }
      else{
        input.readNumber("kappa", ParameterPack::KeyType::Required, kappa_);
        kr = { 0, 1, kappa_, kappa_ };
      }
      kt0_ = kr[0]; ktf_ = kr[1]; k0_ = kr[2]; kf_ = kr[3];
      kslope_ = (kr[3]-kr[2])/(kr[1]-kr[0]);

      doxstarRamp = input.readVector("ramp_xstar", ParameterPack::KeyType::Optional, xr);
      if(doxstarRamp){
        FANCY_ASSERT(xr.size() == 4, "kappa ramping needs [ t_start(sweeps) t_end(sweeps) val_start val_end ]");
      }
      else{
        input.readNumber("xstar", ParameterPack::KeyType::Required, xstar_);
        xr = { 0, 1, xstar_, xstar_};
      }
      xt0_ = xr[0]; xtf_ = xr[1]; x0_ = xr[2]; xf_ = xr[3];
      xslope_ = (xr[3]-xr[2])/(xr[1]-xr[0]);
      halfkappa_ = kappa_*0.5;
      return;
    }
    HarmonicPotential(real kappa, real xstar){
      xstar_ =  xstar;
      kappa_ = kappa;
      halfkappa_ = kappa*0.5;
    }

    virtual void sweepUpdate(real time){
      time_ = time;
      if(time_ < xt0_){
        xstar_ = x0_;
      }
      else if (time_ > xtf_){
        xstar_ = xf_;
      }
      else{
        xstar_ = xslope_*(time - xt0_) + x0_;
      }
      if(time_ < kt0_){
        kappa_ = k0_;
      }
      else if (time_ > ktf_){
        kappa_ = kf_;
      }
      else{
        kappa_ = kslope_*(time - kt0_) + k0_;
      }
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
    real xstar_, halfkappa_, kappa_, time_;
    real xslope_, kslope_;
    real  x0_, xf_, xt0_, xtf_;
    real k0_, kf_, kt0_, ktf_;
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