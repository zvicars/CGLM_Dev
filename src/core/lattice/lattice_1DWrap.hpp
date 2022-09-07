#pragma once
#include "lattice.hpp"

class Lattice_1DWrap : public Lattice{
  public:
    Lattice_1DWrap(InputPack& input);
    //initialize empty for cloning to prevent unecessary file io
    Lattice_1DWrap(InputPack& input, bool placeholder);
    ~Lattice_1DWrap() = default;
    Lattice* clone();
    virtual bool flip(const Vec3<std::size_t>& idx);
    virtual bool getState(const Vec3<std::size_t>& idx) const;
    virtual bool getState(std::size_t idx) const;
    virtual real getMu(const Vec3<std::size_t>& idx) const;
    virtual real getMu(std::size_t idx) const;
    virtual real getPhi(const Vec3<std::size_t>& idx) const;
    virtual real getPhi(std::size_t idx) const;
    virtual void setStates(Matrix3d<char>& states){
      state_ = states;
      return;
    }
    virtual void setMuField(Matrix3d<real>& mus){
      mu_field_ = mus;
      return;
    }
    virtual void setPhiField(Matrix3d<real>& phis){
      phi_field_ = phis;
      return;
    }
    virtual void chooseActiveSite();
    virtual void setActiveSite(const Vec3<int>& idx);
    virtual void reportInfo(string& in);
  protected:
    virtual bool loadFileIntoDataArray(string, Matrix3d<real>&);
    virtual bool loadFileIntoDataArray(string, Matrix3d<char>&);
    //lattice implementations typically just change the way the underlying matrices are stored
    //mostly just for evaluating performance
    Matrix3d<char> state_;
    Matrix3d<real> phi_field_, mu_field_;  
    std::size_t active_index_1d_, size_1d_;
};