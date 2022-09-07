#pragma once
#include "hamiltonian.hpp"
class Hamiltonian_LG : public Hamiltonian{
  public:
  Hamiltonian_LG(InputPack& input);
  ~Hamiltonian_LG() = default;
  virtual real calc_dh(const Lattice& lattice);
  virtual real calc_h(const Lattice& lattice); 
  virtual Hamiltonian* clone();
  //output control
  virtual std::string printStepOutput(std::string s);
  protected:
};