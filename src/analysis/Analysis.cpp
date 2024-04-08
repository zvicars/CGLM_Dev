
#include "../tools/Assert.hpp"
#include "../tools/InputParser.hpp"
#include "../tools/CGLMFileHelper.hpp"
#include "factory.hpp"
#include "Box.hpp"

bool readFrameBinary(Vec3<std::size_t>& size, Vec<real>& states, std::ifstream& ifile){
  binary_real_read(ifile, states, size);
  return ifile.fail();
}

bool readFrameBinary(Vec3<std::size_t>& size, Vec<char>& states, std::ifstream& ifile){
  binary_char_read(ifile, states, size);
  return ifile.fail();
}

bool nextFrame(Box& box_ptr, std::ifstream& box_file){
  std::vector<char> states_temp;
  Vec3<std::size_t> size_temp;
  bool fail = readFrameBinary(size_temp, states_temp, box_file);
  for(auto size : size_temp){
    if(size > 1000000){
      std::cout << "unreasonable box size" << std::endl;
      return 0;
    }
  }
  if(fail) return 0;
  box_ptr.setLattice(states_temp, size_temp);
  return 1; 
}

bool readPhi(Box& box_ptr, std::ifstream& phi_file){
  std::vector<real> phi_temp;
  Vec3<std::size_t> size_temp; 
  bool success = readFrameBinary(size_temp, phi_temp, phi_file);
  box_ptr.setPhi(phi_temp, size_temp);
  return success; 
}

int main(int argc, char **argv)
{
  FANCY_ASSERT(argc == 2, "Analysis code only accepts a single input that specifies the analysis input file.");
  std::string op_input_file_ = argv[1];
  std::string trajectory_file_, phi_file_;
  InputParser input_parser;
  ParameterPack master_pack = input_parser.parseFile(op_input_file_);
  using KeyType = ParameterPack::KeyType;
  bool trajectory_found = master_pack.readString("trajectory", KeyType::Required, trajectory_file_);
  bool phi_found = master_pack.readString("phi", KeyType::Optional, phi_file_);
  Box b1;
  std::ifstream box_stream(trajectory_file_);
  FANCY_ASSERT(box_stream.is_open(), "Failed to open trajectory input file");
  if(phi_found){
    std::ifstream phistream(phi_file_);
    FANCY_ASSERT(phistream.is_open(), "Failed to open phi input file.");
    readPhi(b1, phistream);
    phistream.close();
  }
  AnalysisInputPack master_input_pack = AnalysisInputPack(&master_pack, &b1);

  std::vector<AnalysisInputPack> mod_packs = master_input_pack.buildDerivedInputPacks("Modifier");
  for(std::size_t i = 0; i < mod_packs.size(); i++){
    std::string type, name;
    mod_packs[i].params().readString("type", ParameterPack::KeyType::Required, type);
    mod_packs[i].params().readString("name", ParameterPack::KeyType::Required, name);
    auto mod_ptr = Modifier_Factory(type, mod_packs[i]);
    master_input_pack.addModifier(name, mod_ptr);
  }
  auto mod_reg_ = master_input_pack.ModifierMap();
  std::vector<AnalysisInputPack> pv_packs = master_input_pack.buildDerivedInputPacks("ProbeVolume");
  for(std::size_t i = 0; i < pv_packs.size(); i++){
    std::string type, name;
    pv_packs[i].params().readString("type", ParameterPack::KeyType::Required, type);
    pv_packs[i].params().readString("name", ParameterPack::KeyType::Required, name);
    auto pv_ptr = ProbeVolume_Factory(type, pv_packs[i]);
    master_input_pack.addProbeVolume(name, pv_ptr);
  } 
  std::vector<AnalysisInputPack> calc_packs = master_input_pack.buildDerivedInputPacks("Calculation");
  for(std::size_t i = 0; i < calc_packs.size(); i++){
    std::string type, name;
    calc_packs[i].params().readString("type", ParameterPack::KeyType::Required, type);
    calc_packs[i].params().readString("name", ParameterPack::KeyType::Required, name);
    auto calc_ptr = Calculation_Factory(type, calc_packs[i]);
    master_input_pack.addCalculation(name, calc_ptr);
  } 

  int step_iterator = 0;
  
  auto pv_reg_ = master_input_pack.ProbeVolumeMap();
  auto calc_reg_ = master_input_pack.CalculationMap();
  while(nextFrame(b1, box_stream)){
    if(box_stream.fail()) break;
    std::cout << step_iterator << "\n";
    //run all calculations
    int finished_counter = 0;
    for(auto& i : calc_reg_){
      i.second->update();
    }    
    for(auto& i : calc_reg_){
      i.second->calculate();
      finished_counter += (int)i.second->isFinished(); //count the number of calculations that are completely finished
    }
    //perform all outputs
     for(auto& i : calc_reg_){
      i.second->printConsoleReport();
      i.second->output();
    } 
    step_iterator++;
    b1.iterateFrame(); 
    for(auto& i : calc_reg_){
      i.second->finish();
    }
    if(finished_counter >= calc_reg_.size()){
      break;
    }
  }
  for(auto& i : calc_reg_){
    i.second->finalOutput();
  }
  box_stream.close();
  return 0;
}