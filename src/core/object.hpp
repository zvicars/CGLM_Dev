#pragma once
#include "inputpack.hpp"
#include <map>

class Object{
  public:
  using string = std::string;
  Object(InputPack& input){
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
  virtual void reportInfo(string& in){
    in = "";
    in = in + "name = " + name_ + "\ntype = " + type_ + "\ninitialized = " + std::to_string(isInitialized) + "\n";
    return;
  }
  virtual std::string printStepOutput(std::string){
    FANCY_ASSERT(0, "object type does not have any outputs");
    return "";
  }
  protected:
  string name_="", type_="";
  //all objects will have modify-access on the original input pack
  InputPack* input_=0;
  bool isInitialized=0;
};