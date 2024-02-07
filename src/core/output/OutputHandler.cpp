#include "OutputHandler.hpp"
#include "../simulation.hpp"
#include "../hamiltonian/hamiltonian.hpp"
#include "../bias/bias.hpp"
#include "../pv/probevolume.hpp"
#include "../rng/random.hpp"
#include "../lattice/lattice.hpp"
#include "../../tools/CGLMFileHelper.hpp"
void OutputHandler::output_ts(Simulation* sim, std::uint64_t step){
  if(step < begin_ts_) return;
  if(step > end_ts_) return;
  std::size_t steps_from_begin = step - begin_ts_;
  if(steps_from_begin%freq_ts_ != 0) return;
  //std::cout << "Outputting timeseries..." << std::endl;
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
  *ofile_handle_ts_ << step << "  " << sim->current_sweep_ << "  ";
  for(auto entry : output_list){
    *ofile_handle_ts_ << entry << "  ";
  }
  *ofile_handle_ts_ << "\n";

  return;
};

void OutputHandler::output_traj(Simulation* sim, std::uint64_t step){
  if(step < begin_traj_) return;
  if(step > end_traj_) return;
  std::size_t steps_from_begin = step - begin_traj_;
  if(steps_from_begin%freq_traj_ != 0) return;
  //std::cout << "Outputting trajectory..." << std::endl;
  //sim->lattice_->TestPrint();
  auto size = sim->lattice_->size();
  int size_1d = size[0]*size[1]*size[2];
  std::vector<char> state_out;
  state_out.reserve(size_1d);
  std::size_t idx = 0; 
  Matrix3d<char> temp(size);
  temp.fill(0.0);
  for(int i = 0; i < size_1d; i++){
    auto pos = temp.map1N(i);
    if(sim->lattice_->getState(pos)) state_out.push_back('1');
    else state_out.push_back('0');
  }
  binary_char_write(*ofile_handle_traj_, state_out, size);
  return;
};


void OutputHandler::finalize_output(Simulation* sim, std::uint64_t step){
  if(ofile_final_.empty()) return;
  std::ofstream final(ofile_final_);
  auto size = sim->lattice_->size();
  int size_1d = size[0]*size[1]*size[2];
  std::vector<char> state_out;
  state_out.reserve(size_1d);
  std::size_t idx = 0; 
  Matrix3d<char> temp(size);
  temp.fill(0.0);
  for(int i = 0; i < size_1d; i++){
    auto pos = temp.map1N(i);
    if(sim->lattice_->getState(pos)) state_out.push_back('1');
    else state_out.push_back('0');
  }
  binary_char_write(final, state_out, size); 
  final.close();
}