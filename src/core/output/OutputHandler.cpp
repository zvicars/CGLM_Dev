#include "OutputHandler.hpp"
#include "../simulation.hpp"
#include "../hamiltonian/hamiltonian.hpp"
#include "../bias/bias.hpp"
#include "../pv/probevolume.hpp"
#include "../rng/random.hpp"
#include "../lattice/lattice.hpp"
#include "../../tools/CGLMFileHelper.hpp"
void OutputHandler::output_ts(Simulation* sim, std::size_t step){
  if(step < begin_ts_) return;
  if(step > end_ts_) return;
  std::size_t steps_from_begin = step - begin_ts_;
  if(steps_from_begin%freq_ts_ != 0) return;
  std::vector<std::string> output_list(output_commands_.size());
  for(std::size_t i = 0; i < output_commands_.size(); i++){
    std::string s = output_commands_[i];
    if(s.find("h") == 0){
      std::string s_post = s.substr(s.find('.')+1);
      output_list[i] = sim->ham_->printStepOutput(s_post);
      continue;
    }
    if(s == "ratio"){
      output_list[i] = std::to_string(sim->accept_ratio_);
    }    
  }
  *ofile_handle_ts_ << sim->current_sweep_ << "  ";
  for(auto entry : output_list){
    *ofile_handle_ts_ << entry << "  ";
  }
  *ofile_handle_ts_ << "\n";

  return;
};

void OutputHandler::output_traj(Simulation* sim, std::size_t step){
  if(step < begin_traj_) return;
  if(step > end_traj_) return;
  std::size_t steps_from_begin = step - begin_traj_;
  if(steps_from_begin%freq_ts_ != 0) return;
  auto lattice = sim->lattice_;
  auto size = lattice->size();
  std::vector<bool> state_out(size[0]*size[1]*size[2]);
  std::size_t idx = 0; 
  for(std::size_t i = 0; i < size[0]; i++)
  for(std::size_t j = 0; j < size[1]; j++)
  for(std::size_t k = 0; k < size[2]; k++){
    state_out[idx] = lattice->getState({i,j,k});
    idx++; 
  }
  binary_bool_write(*ofile_handle_traj_, state_out, size);
  return;
};