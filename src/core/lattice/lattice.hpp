//definition of lattice, contains all information that a hamiltonian will need to compute relevant energies
//all energies in kT and all lengthscales in lattice units
#pragma once
#include "../object.hpp"
#include "../rng/random.hpp"
#include "../../tools/Matrix.hpp"
class Lattice : public Object{
  public:
  Lattice(InputPack& input);
  virtual ~Lattice();
  //non-const functions
  virtual void chooseActiveSite() = 0; //do rng to find new site
  virtual void setActiveSite(const Vec3<int>& idx) = 0; //pick one automatically, using int input to allow for pbc wrap
  virtual Lattice* clone() = 0;
  virtual bool flip(const Vec3<std::size_t>& idx) = 0;
  bool flipPBC(const Vec3<int>& index);
  bool flipActive();
  //system info getters
  Vec3<std::size_t> size() const { return size_;}
  real eps() const { return eps_; }
  real lambda() const { return lambda_;}
  real lambda3() const { return lambda3_;}
  real density() const { return density_;}
  virtual void TestPrint()=0;
  //site info getters
  //any get active site must be in-bounds
  Vec3<std::size_t> getActiveIndex() const { return active_index_; }
  virtual bool getState(const Vec3<std::size_t>& idx) const = 0;
  bool getStatePBC(const Vec3<int>& index) const;
  virtual bool getActiveState() const { return active_state_; }
  virtual real getMu(const Vec3<std::size_t>& idx) const = 0;
  real getMuPBC(const Vec3<int>& index) const;
  virtual real getActiveMu() const { return active_mu_; }
  virtual real getPhi(const Vec3<std::size_t>& idx) const = 0;
  real getPhiPBC(const Vec3<int>& index) const;
  virtual real getActivePhi() const { return active_phi_; }
  int getAdj(const Vec3<std::size_t>& index) const; //num of adjacent sites
  int getAdjPBC(const Vec3<int>& index) const; //num of adjacent sites
  virtual int getActiveAdj() const { return active_adj_; }
  virtual int sweepSize(){
    return size_[0]*size_[1]*size_[2];
  }
  //output functions
  virtual void reportInfo(string& in);
  protected:
  std::string sf_, mff_, pff_; //filenames
  bool hasPhiField=0, hasMuField=0, hasStateField=0; //various flags
  bool default_state_ = 0;
  real mu_, eps_, density_; //default mu/eps
  real lambda_, lambda3_; //gridspacings
  Vec3<std::size_t> size_;
  Vec3<bool> pbc_;
  RNG* random_generator_ = 0;
  //active quantities, these update whenever an active site is chosen
  Vec3<std::size_t> active_index_;
  char active_state_;
  int active_adj_;
  double active_phi_, active_mu_;
};