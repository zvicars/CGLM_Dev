#pragma once
#include "probevolume.hpp"
class PV_SimpleRect : public ProbeVolume{
  public:
    PV_SimpleRect(InputPack& input);
    ~PV_SimpleRect() = default;
    virtual ProbeVolume* clone();
    virtual real calc_nv(const Lattice& lattice);
    virtual real calc_nv_compute();
    virtual real calc_dnv(const Lattice& lattice);
    virtual std::string printStepOutput(std::string s);
  protected:
    Vec3<int> min_, max_;
};