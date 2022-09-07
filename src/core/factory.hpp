#pragma once
#include "inputpack.hpp"
#include "object.hpp"
#include "bias/bias.hpp"
#include "hamiltonian/hamiltonian.hpp"
#include "lattice/lattice.hpp"
#include "pv/probevolume.hpp"
#include "rng/random.hpp"
#include "simulation.hpp"
Simulation* simulationFactory(InputPack& in);
Lattice* latticeFactory(InputPack& in);
Bias* biasFactory(InputPack& in);
Hamiltonian* hamiltonianFactory(InputPack& in);
ProbeVolume* probevolumeFactory(InputPack& in);
RNG* randomFactory(InputPack& in);