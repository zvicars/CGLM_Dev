//header for the cglm hamiltonian
//this will include a contribution from the lattice gas in tandem with 
//the chi_ij matrix that accounts for excluded volumes
//in the event that there are no excluded volumes, it will report the relevant
//covariance terms for getting the appropriate histograms from probe volumes
//on output steps
#pragma once
#include "hamiltonian_LG.hpp"
class SubvolumePair{
public:
  double compute_chi(const Lattice& lattice);
  double compute_chi_init(const Lattice& lattice);
  double compute_dchi(const Lattice& lattice);
  Vec2<ProbeVolume*> pvs;
  real chi_12;
  bool self; //does i == j?
private:
  double chi, dchi, chi_uu, dchi_uu;
  Vec<real> values_;
  Vec<Vec3<int> > offsets_;
};

class Hamiltonian_CGLM : public Hamiltonian_LG{
  public:
  Hamiltonian_CGLM(InputPack& input);
  ~Hamiltonian_CGLM(){
    for(std::size_t i = 0; i < ve_vec_.size(); i++){
      delete ve_vec_[i];
    }
    delete vp_;
  }
  virtual real calc_dh(const Lattice& lattice);
  virtual real calc_h(const Lattice& lattice);
  virtual Hamiltonian* clone();
  protected:
  real computeChiUU(const Lattice& lattice);
  real computeChiPU(ProbeVolume* pv, const Lattice& lattice);
  bool use_aux_;
  int aux_type_;
  std::string chi_file_;
  Vec<real> chi_ij_values_;
  Vec<Vec3<int> > chi_ij_offsets_;
  Vec<ProbeVolume*> ve_vec_;
  Vec<real> chi_pu_;
  ProbeVolume* vp_;

  std::vector<std::string> ve_names;
};