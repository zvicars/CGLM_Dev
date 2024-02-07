#include "simulation.hpp"
#include "hamiltonian/hamiltonian.hpp"
#include "lattice/lattice.hpp"
#include "rng/random.hpp"

Simulation::Simulation(InputPack& input) : Object{input}{
  //get all probe volumes
  input_ = &input;
  std::string lattice_name, rng_name, ham_name;
  //read required parameters
  input_->params().readString("lattice", ParameterPack::KeyType::Required, lattice_name);
  input_->params().readString("rng", ParameterPack::KeyType::Required, rng_name);
  input_->params().readString("hamiltonian", ParameterPack::KeyType::Required, ham_name);
  print_freq_ = 1;
  input_->params().readNumber("print_freq", ParameterPack::KeyType::Optional, print_freq_);
  if(print_freq_ < 1) print_freq_ = std::numeric_limits<int>::max();
  FANCY_ASSERT(print_freq_ > 0, "print frequency needs to be a positive integer");
  //search input pack for constructed versions of parameters
  Lattice* lattice_ptr;
  RNG* rng_ptr;
  Hamiltonian* hamiltonian_ptr;
  lattice_ptr = input_->find(input_->lattices(), lattice_name);
  FANCY_ASSERT(lattice_ptr != 0, "Failed to find lattice: " + lattice_name + ".");
  hamiltonian_ptr = input_->find(input_->hamiltonians(), ham_name);
  FANCY_ASSERT(hamiltonian_ptr != 0, "Failed to find hamiltonian: " + ham_name + ".");
  rng_ptr = input_->find(input_->randomgenerators(), rng_name);
  FANCY_ASSERT(rng_ptr != 0, "Failed to find rng: " + rng_name + ".");

  //copy construct all of the pointers into internal members
  lattice_ = lattice_ptr->clone();
  ham_ = hamiltonian_ptr->clone();
  rng_ = rng_ptr->clone();
  //load output handler
  auto p_out = input_->params().findParameterPacks("Output", ParameterPack::KeyType::Optional);
  if(p_out.size() > 0){
    FANCY_ASSERT(p_out.size() == 1, "only one output type supported per simulation");
    output_.initialize(*p_out[0]);
  }
  input_->params().readNumber("sweeps", ParameterPack::KeyType::Required, max_sweeps_);
  current_sweep_ = 0;
  total_step_counter_ = 0;
  return;
};

//build all relevant objects
bool Simulation::initialize(){
  ham_->calc_h(*lattice_); //should set all of the initial values for probe volumes and the like as well
  return 0;
}
bool Simulation::step(){
  lattice_->chooseActiveSite();
  double dh = ham_->calc_dh(*lattice_);
  bool accept = rng_->getReal(0.0, 1.0) < exp(-dh);
  if(accept){
    lattice_->flipActive();
    ham_->flip();
    return 1; 
  }
  return 0;
};

void Simulation::sweep(){
  std::size_t nflips_max = lattice_->sweepSize(), nflips = 0, nsucc = 0;
  while(nflips < nflips_max){
    //attempt flip
    nsucc += step();
    //iterate
    nflips++;
    total_step_counter_++;
  }
  accept_ratio_ = (real)nsucc/(real)nflips;
  return;
}

void Simulation::run(){
  initialize();
  total_step_counter_ = 0;
  for(current_sweep_ = 0; current_sweep_ < max_sweeps_; current_sweep_++){
    //pretty much just updates ramped biases, once per sweep because it's kinda wasteful to do it every flip
    ham_->sweepUpdate( current_sweep_ );
    //per-step updates
    output_.output_traj(this, current_sweep_);
    output_.output_ts(this, current_sweep_);
    if((current_sweep_+1)%print_freq_ == 0){
      std::cout << current_sweep_+1 << " ";
      std::cout << std::flush;
    }

    //sweep
    sweep();
  }
  output_.finalize_output(this, current_sweep_);
  //output_.close_outputs();
  return;
}

Simulation::~Simulation(){
  delete lattice_;
  delete rng_;
  delete ham_;
  return;
}