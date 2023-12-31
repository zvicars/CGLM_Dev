#include "lattice_1DWrap.hpp"
#include "lattice_pbc.hpp"
#include "../../tools/CGLMFileHelper.hpp"
Lattice_1DWrap::Lattice_1DWrap(InputPack& input) : Lattice{input}{
  state_.initialize(size_);
  state_.fill((char)default_state_);
  if(hasStateField){
    loadFileIntoDataArray(sf_, state_);
  }
  if(hasPhiField){
    phi_field_.initialize(size_);
    phi_field_.fill(0.0);
    loadFileIntoDataArray(pff_, phi_field_);
  }
  if(hasMuField){
    mu_field_.initialize(size_);
    mu_field_.fill(0.0);
    loadFileIntoDataArray(mff_, mu_field_);
  }
  size_1d_ = state_.size_1d();
  fixed_ = 0;
  input_->params().readFlag("nonfixed_list", ParameterPack::KeyType::Optional, fixed_);
  if(fixed_){
    input_->params().readNumber("nonfixed_list_threshold", ParameterPack::KeyType::Required, nfl_cutoff_); 
    setNonfixedList();
  }
  return;
}

void Lattice_1DWrap::setNonfixedList(){
  fixed_ = 1;
  nfl_.reserve(size_1d_);
  for(int i=0; i < size_1d_; i++){
    if( getPhi(i) + getMu(i) < -nfl_cutoff_){
      state_.at_1d(i) = 1;
    }
    else if( getPhi(i) + getMu(i) > nfl_cutoff_){
      state_.at_1d(i) = 0;
    }    
    else{
      nfl_.push_back(i);
    }
  }
  nfl_size_ = nfl_.size();
}
/*
Lattice_1DWrap::Lattice_1DWrap(InputPack& input, bool placeholder):Lattice{input}{
  state_.initialize(size_);
  size_1d_ = state_.size_1d();
  input_->params().readFlag("nonfixed_list", ParameterPack::KeyType::Required, fixed_);
  return;
}
*/

Lattice_1DWrap::Lattice_1DWrap(InputPack& input, bool placeholder):Lattice{input}{
  state_.initialize(size_);
  state_.fill((char)default_state_);
  if(hasStateField){
    loadFileIntoDataArray(sf_, state_);
  }
  if(hasPhiField){
    phi_field_.initialize(size_);
    phi_field_.fill(0.0);
    loadFileIntoDataArray(pff_, phi_field_);
  }
  if(hasMuField){
    mu_field_.initialize(size_);
    mu_field_.fill(0.0);
    loadFileIntoDataArray(mff_, mu_field_);
  }
  size_1d_ = state_.size_1d();
  fixed_ = 0;
  input_->params().readFlag("nonfixed_list", ParameterPack::KeyType::Optional, fixed_);
  if(fixed_){
    input_->params().readNumber("nonfixed_list_threshold", ParameterPack::KeyType::Required, nfl_cutoff_); 
    setNonfixedList();
  }
  return;
}

void Lattice_1DWrap::chooseActiveSite(){
  if(!fixed_) active_index_1d_ = random_generator_->getIndex(0, size_1d_-1);
  else active_index_1d_ = nfl_[random_generator_->getIndex(0, nfl_size_-1)];
  active_index_ = state_.map1N(active_index_1d_);
  active_state_ = getState(active_index_1d_);
  active_mu_ = getMu(active_index_1d_);
  active_phi_ = getPhi(active_index_1d_);
  active_adj_ = getAdj(active_index_);
  return;
}
void Lattice_1DWrap::setActiveSite(const Vec3<int>& idx){
  auto idx2 = wrap3(idx, size_);
  active_index_1d_ = state_.mapN1(idx2);
  active_index_ = idx2;
  active_state_ = getState(active_index_1d_);
  active_mu_ = getMu(active_index_1d_);
  active_phi_ = getPhi(active_index_1d_);
  active_adj_ = getAdj(active_index_);
  return;
}
bool Lattice_1DWrap::getState(const Vec3<std::size_t>& idx) const{
  return state_.read_at(idx);
}
bool Lattice_1DWrap::getState(std::size_t idx) const{
  return state_.read_at_1d(idx);
}
real Lattice_1DWrap::getMu(const Vec3<std::size_t>& idx) const{
  if(hasMuField) return mu_field_.read_at(idx);
  return mu_;
}
real Lattice_1DWrap::getMu(std::size_t idx) const{
  if(hasMuField) return mu_field_.read_at_1d(idx);
  return mu_;
}
real Lattice_1DWrap::getPhi(const Vec3<std::size_t>& idx) const{
  if(hasPhiField) return phi_field_.read_at(idx);
  return 0.0;
}
real Lattice_1DWrap::getPhi(std::size_t idx) const{
  if(hasPhiField) return phi_field_.read_at_1d(idx);
  return 0.0;
}
Lattice* Lattice_1DWrap::clone(){
  //being a bit lazy here and just using the stored input pack to deal
  //with the initialation of the Lattice() part, hence the superfluous 0
  //to distinguish the constructors
  auto ptr = new Lattice_1DWrap(*input_, 0);
  //ptr->setMuField(mu_field_);
  //ptr->setStates(state_);
  //ptr->setPhiField(phi_field_);
  //if(fixed_) ptr->setNonfixedList();
  return ptr;
}

bool Lattice_1DWrap::flip(const Vec3<std::size_t>& idx){
  auto& state = state_.at(idx);
  if(state == 1) state = 0; 
  else state = 1;
  return state;
}

void Lattice_1DWrap::TestPrint(){
  std::ofstream ofile("test_out.xyz");
  ofile << "xyz file data" << std::endl;
  ofile << size_1d_ << std::endl;
  for(int i = 0; i < size_1d_; i++){
    auto pos = state_.map1N(i);
    ofile << pos[0] << "  " << pos[1] << "  " << pos[2] << "  " << (state_.at(pos)==1) << std::endl;
  }
  ofile.close();
}

void Lattice_1DWrap::reportInfo(string& in){
  Lattice::reportInfo(in);
  return;
}

bool Lattice_1DWrap::loadFileIntoDataArray(string filename, Matrix3d<real>& field){
  std::ifstream ifile(filename, std::ios::binary);
  FANCY_ASSERT(ifile.is_open(), "failed to open input file, " + filename);
  readFileIntoArray(ifile, size_, field);
  ifile.close();
  return 0;
}
bool Lattice_1DWrap::loadFileIntoDataArray(string filename, Matrix3d<char>& field){
  std::ifstream ifile(filename, std::ios::binary);
  FANCY_ASSERT(ifile.is_open(), "failed to open input file, " + filename);
  readFileIntoArray(ifile, size_, field);
  ifile.close();
  return 0;
}