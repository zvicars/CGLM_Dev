
#pragma once
#include "AnalysisInputPack.hpp"
#include <map>

class AnalysisObject{
  public:
  using string = std::string;
  AnalysisObject(AnalysisInputPack& input){
    input.params().readString("name", ParameterPack::KeyType::Required, name_);
    input.params().readString("type", ParameterPack::KeyType::Required, type_);
    input_ = &input;
    return;
  }
  string getName(){
    return name_;
  }
  string getType(){
    return type_;
  }
  protected:
  string name_="", type_="";
  //all objects will have modify-access on the original input pack
  AnalysisInputPack* input_=0;
  bool isInitialized=0;
};
