#include "bias_nv.hpp"
Bias_Nv::Bias_Nv(InputPack& input) : Bias{input}{
  std::string pv_name, function_type;
  input_->params().readString("probevolume", ParameterPack::KeyType::Required, pv_name);
  auto pv_temp = input_->findProbeVolume(pv_name);
  FANCY_ASSERT(pv_temp != 0, "Failed to find specified probevolume: " + pv_name + " in bias: " + name_);
  //clone the probevolume template in the input pack to make the internal probe volume for the bias
  pv_ = pv_temp->clone();
  auto function_params = *(input_->params().findParameterPack("function", ParameterPack::KeyType::Required));
  bf_ = bias_function_factory(function_params);
  return;
}

std::string Bias_Nv::printStepOutput(std::string s){
  if(s == "u") return std::to_string(u_);
  if(s.find("pv") == 0){
    std::string s_post = s.substr(s.find('.')+1);
    return pv_->printStepOutput(s_post);    
  }
  FANCY_ASSERT(0, "invalid output type specified");
  return "";
}