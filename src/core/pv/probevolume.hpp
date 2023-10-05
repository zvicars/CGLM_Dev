#pragma once
#include "../object.hpp"
class ProbeVolume : public Object{
  public:
  ProbeVolume(InputPack& input);
  virtual ~ProbeVolume() = default;
  virtual real calc_nv(const Lattice& lattice)=0;
  virtual real calc_dnv(const Lattice& lattice)=0;
  real nv(){
    return nv_;
  }
  real dnv(){
    return dnv_;
  }
  void flip(){
    nv_ += dnv_;
    dnv_ = 0.0;
  }
  virtual void bounds(Vec<std::size_t>&) = 0;
  virtual bool isInside(Vec3<std::size_t> idx) = 0;
  virtual ProbeVolume* clone() = 0;
  protected:
  real nv_, dnv_;
  const Lattice* stored=0;
};