#pragma once
#include <map>
#include <string>
#include <vector>
#include <limits>
#include <stdint.h>
#include "../inputpack.hpp"
//output handler class will be hand-written to allow outputting information from a simulation at a given timestep
//each object will have a function called sendOutput(), which sends a set of key-value pairs that the output handler
//can parse
//how to manage nested objects? 
//simulation <sim>
//-hamiltonian <h>
//--bias <b>
//---probevolume <v>
//sim.h.h
//sim.h.dh
//sim.h.bias[0].u
//sim.h.bias[0].pv.n
//needs to be navigable somehow
//register value-type and object-type
//if object-type is specified, gets the next 
class OutputHandler{
  public:
    OutputHandler() = default;
    ~OutputHandler(){
      ofile_handle_ts_->close();
      ofile_handle_traj_->close();
      delete ofile_handle_ts_;
      delete ofile_handle_traj_;
    }
    void initialize(const ParameterPack& p){
      p.readVector("args", ParameterPack::KeyType::Required, output_commands_);
      p.readNumber("begin_ts", ParameterPack::KeyType::Optional, begin_ts_);
      p.readNumber("end_ts", ParameterPack::KeyType::Optional, end_ts_);
      p.readNumber("freq_ts", ParameterPack::KeyType::Optional, freq_ts_);
      p.readString("output_file", ParameterPack::KeyType::Required, ofile_ts_);
      ofile_handle_ts_ = new std::ofstream(ofile_ts_);
      *ofile_handle_ts_ << "#sweeps  ";
      for(auto cmd : output_commands_){
        *ofile_handle_ts_ << cmd << "  ";
      }
      *ofile_handle_ts_ << std::endl;
      p.readNumber("begin_traj", ParameterPack::KeyType::Optional, begin_traj_);
      p.readNumber("end_traj", ParameterPack::KeyType::Optional, end_traj_);
      p.readNumber("freq_traj", ParameterPack::KeyType::Optional, freq_traj_);
      p.readString("output_traj", ParameterPack::KeyType::Required, ofile_traj_);
      p.readString("final_config", ParameterPack::KeyType::Optional, ofile_final_);
      ofile_handle_traj_ = new std::ofstream(ofile_traj_, std::ios::binary);      
      return;
    }
    void output_ts(Simulation* sim, std::uint64_t step);
    void output_traj(Simulation* sim, std::uint64_t step);
    void finalize_output(Simulation* sim, std::uint64_t step);
    void close_outputs(){
      ofile_handle_ts_->close();
      ofile_handle_traj_->close();
    }
  private:
    std::vector<std::string> output_commands_;
    std::string ofile_ts_, ofile_traj_, ofile_final_;
    uint64_t begin_ts_ = 0, end_ts_ = std::numeric_limits<uint64_t>::max(), freq_ts_ = 1;
    uint64_t begin_traj_ = 0, end_traj_ = std::numeric_limits<uint64_t>::max(), freq_traj_ = 1;
    std::ofstream* ofile_handle_ts_, * ofile_handle_traj_;
};