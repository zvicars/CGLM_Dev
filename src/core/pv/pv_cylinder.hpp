#pragma once
#include "probevolume.hpp"
class PV_Cylinder : public ProbeVolume{
  public:
    PV_Cylinder(InputPack& input);
    ~PV_Cylinder() = default;
    virtual ProbeVolume* clone();
    virtual real calc_nv(const Lattice& lattice);
    virtual real calc_nv_compute();
    virtual real calc_dnv(const Lattice& lattice);
    virtual std::string printStepOutput(std::string s);
    virtual void bounds(Vec<std::size_t>&);
    virtual bool isInside(Vec3<std::size_t> idx);
  protected:
    Vec3<real> origin_;
    int axis_;
    real h_, r_, r2_;
    bool invert_=0;
};