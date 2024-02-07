#pragma once
#include <vector>
#include <limits>
#include "../../typedefs.hpp"
#include "../../tools/Assert.hpp"
class RampedParameter{
public:
  RampedParameter() = default; 
  RampedParameter(const std::vector<real>& xy){
    initialize(xy);
    return;
  }
  RampedParameter(const std::vector<real>& x, const std::vector<real>& y){
    initialize(x,y);
    return;
  }
  void initialize(const std::vector<real>& xy){
    FANCY_ASSERT(xy.size() % 2 == 0, "Ramped parameters need to have a list of time points followed by a list of values\
    and thus require an even number of elements.");
    FANCY_ASSERT(xy.size() > 0, "RampedParameters need at least 2 entries.");
    std::size_t ref_index = (xy.size() / 2);
    x_.insert(x_.end(), xy.begin(), xy.begin() + ref_index);
    y_.insert(y_.end(), xy.begin() + ref_index, xy.end());
    FANCY_ASSERT(x_.size() == y_.size(), "Sanity check failed, x_ and y_ are different sizes.");
    checkIncreasing();
    setMaxMin();
    isInitialized = 1;
    return;    
  }
  void initialize(const std::vector<real>& x, const std::vector<real>& y){
    x_ = x;
    y_ = y;
    FANCY_ASSERT(x_.size() == y_.size(), "Sanity check failed, x_ and y_ are different sizes.");
    checkIncreasing();
    setMaxMin();
    isInitialized = 1;    
  }
  real compute(real x){
    FANCY_ASSERT(isInitialized, "Ramped parameter is not initialized!");
    real ret_val;
    //will assume that this vector is small ( < 10 elements), so no fancy containers
    //find the last element that's less than t_min
    std::size_t index = 0;
    while(x > x_[index] && index < x_.size()){
        index++;
    }
    if(index == 0){
        return min_y;
    }
    if(index == x_.size()){
        return max_y;
    }
    //if it's in bounds, the value should be linearly interpolated from y_[index-1] and y_[index]
    double x0 = x_[index-1], xf = x_[index], y0 = y_[index-1], yf = y_[index];
    double slope = (yf-y0)/(xf-x0);
    ret_val = slope*(x-x0) + y0;
    return ret_val;
  }
  
protected:
  void checkIncreasing(){
    real previous = std::numeric_limits<real>::lowest();
    for(int i = 0; i < x_.size(); i++){
        real val = x_[i];
        FANCY_ASSERT(previous < val, "x-points are not in increasing order");
        previous = val;
    }
    return;
  }
  void setMaxMin(){
    min_x = x_.front();
    max_x = x_.back();
    min_y = y_.front();
    max_y = y_.back();
    return;
  }
  bool isInitialized=0;
  std::vector<real> x_, y_;
  real min_x, min_y, max_x, max_y;
};