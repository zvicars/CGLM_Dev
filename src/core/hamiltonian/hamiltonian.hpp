//Zachariah Vicars
//Hamiltonian header is an abstract class that can compute the total Hamiltonian for a lattice
#pragma once
#include "../object.hpp"
#include "../pv/probevolume.hpp"
#include "../bias/bias.hpp"
class Hamiltonian : public Object{
  public:
  Hamiltonian(InputPack& input);
  virtual ~Hamiltonian();
  virtual real calc_dh(const Lattice& lattice) = 0;
  virtual real calc_h(const Lattice& lattice) = 0; 
  //general update function, do any time-dependent stuff
  virtual void sweepUpdate(real time){
    for(auto bias : biases_){
      bias->sweepUpdate(time);
    }
  }
  real h(){
    return h_;
  }
  real dh(){
    return dh_;
  }
  virtual Hamiltonian* clone() = 0;
  virtual void flip(){
    h_ += dh_;
    dh_ = 0;
    for(auto bias : biases_){
      bias->flip();
    }
  }
  protected:
  std::vector<Bias*> biases_;
  double h_;
  double dh_;
};