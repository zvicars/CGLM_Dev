#pragma once
#include "../../tools/InputParser.hpp"
#include "../../typedefs.hpp"
class BiasFunction{
  public:
    virtual real operator()(real n) = 0;
    virtual real operator()(real n, real dn) = 0;
};
class HarmonicPotential : public BiasFunction{
  public:
    HarmonicPotential(const ParameterPack& input){
      input.readNumber("xstar", ParameterPack::KeyType::Required, xstar_);
      input.readNumber("kappa", ParameterPack::KeyType::Required, kappa_);
      halfkappa_ = kappa_*0.5;
      return;
    }
    HarmonicPotential(real kappa, real xstar){
      xstar_ =  xstar;
      kappa_ = kappa;
      halfkappa_ = kappa*0.5;
    }
    //u
    virtual real operator()(real x){
      return halfkappa_*(x-xstar_)*(x-xstar_);
    }
    //u2-u1 for a change of dn
    virtual real operator()(real x, real dx){
      return halfkappa_*dx*(2.0*(x - xstar_) - dx);
    }
  private:
    real xstar_, halfkappa_, kappa_;
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