#include "Calculation.hpp"
#include "../modifiers/Modifier.hpp"
#include <limits>
Calculation::Calculation(AnalysisInputPack& input):AnalysisObject{input}
{
    output_freq_ = -1; //will never output
    input.params().readNumber("output_frequency", KeyType::Optional, output_freq_);
    if(output_freq_ > 0) calc_freq_ = output_freq_;
    else calc_freq_ = 1;
    input.params().readNumber("calculation_frequency", KeyType::Optional, calc_freq_);
    equilibration_ = 0;
    FANCY_ASSERT( output_freq_%calc_freq_ == 0 || output_freq_ < 0, "Can only output data on a calculation step! Output frequency should be a multiple of calculation frequency.");
    //equilibration in ps
    input.params().readNumber("equilibration", KeyType::Optional, equilibration_);
    end_ = std::numeric_limits<int>::max();
    input.params().readNumber("end", KeyType::Optional, end_);
    //probably advantageous to make modifiers site-wise, functions
    Vec<std::string> modifier_names;
    input.params().readVector("modifiers", KeyType::Optional, modifier_names);
    for(auto name : modifier_names){
      auto mod_ptr = input.findModifier(name);
      FANCY_ASSERT(mod_ptr != NULL, "could not find the specified modifier name.");
      modifiers_.push_back(mod_ptr);
    }
    input.addCalculation(name_, this);
    box = input.getBox();
    update_flag_ = 0;
    calculate_flag_ = 0;
    return;
}

bool Calculation::doCalculate(){
  if(hasCalculated()) return 0;
  if(calc_freq_ == -1) return 0;
  if(current_time_ < equilibration_ || current_time_ > end_) return 0;
  if(current_time_%calc_freq_ != 0) return 0;
  calculate_flag_ = 1;
  return 1;
}

//another version that doesn't check to see if a calculation has already been called, prevents child class from calculating
//when parent class hasn't calculated
bool Calculation::doCalculateNoFlag() const{
  if(current_time_ < equilibration_ || current_time_ > end_) return 0;
  if(current_time_%calc_freq_ != 0) return 0;
  return 1;
}

bool Calculation::doOutput(){
  if(output_freq_ == -1) return 0;
  if(current_time_ < equilibration_ || current_time_ > end_) return 0;
  if(current_time_%output_freq_ != 0) return 0;
  return 1;
}

void Calculation::update(){
    if(update_flag_ == 1) return;
    current_time_ = box->frame();
    update_flag_ = 1;
    //rest of this is pretty expensive, so really only need to do it when I'm actually calculating something
    if(!doCalculateNoFlag()) return;
    //apply modifiers to the lattice for this calculation
    internalLattice.clear();
    internalLattice.initialize(box->size_);
    for(int i = 0; i < internalLattice.size_1d(); i++){
      internalLattice.at_1d(i) = box->lattice_.read_at_1d(i);
    }
    for(auto modifier : modifiers_){
      modifier->compute(internalLattice);
    }
    return;
}