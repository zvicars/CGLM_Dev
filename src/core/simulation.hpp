#pragma once
#include "object.hpp"
#include "output/OutputHandler.hpp"
#include <vector>
class Simulation : public Object{
  public:
    Simulation(InputPack& input);
    ~Simulation();
    bool initialize();
    bool step();
    void sweep();
    void run();
    protected:
    Lattice* lattice_;
    RNG* rng_;
    Hamiltonian* ham_;
    std::size_t current_sweep_, max_sweeps_;
    real accept_ratio_;
    OutputHandler output_;
  friend OutputHandler;
};