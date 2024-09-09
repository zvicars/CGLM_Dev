#pragma once
#include "../../tools/typedefs.hpp"
#include "../../tools/Matrix.hpp"

class FF_function{
public:
    FF_function(Vec3<real> x_ref, Vec3<real> size, real cutoff, Vec<real> params){
        cutoff_ = cutoff;
        xref_ = x_ref;
        size_ = size;
        return;
    };
    virtual real compute(Vec3<real> x) const = 0;
    virtual Vec<real> getBoundingBox() const = 0;
protected:
  real cutoff_;
  Vec3<real> xref_, size_;
};

class LJ_6_12_default : public FF_function{
public:
    LJ_6_12_default(Vec3<real> x_ref, Vec3<real> size, real cutoff, Vec<real> params);
    real compute(Vec3<real> x) const;
    virtual Vec<real> getBoundingBox() const;
private:
  real sigma_, epsilon_; 
};

class LJ_6_12_offset : public FF_function{
public:
    LJ_6_12_offset(Vec3<real> x_ref, Vec3<real> size, real cutoff, Vec<real> params);
    real compute(Vec3<real> x) const;
    virtual Vec<real> getBoundingBox() const;
private:
  real sigma_, epsilon_, offset_; 
};

class LJ_6_12_offset_box : public FF_function{
public:
    LJ_6_12_offset_box(Vec3<real> x_ref, Vec3<real> size, real cutoff, Vec<real> params);
    real compute(Vec3<real> x) const;
    virtual Vec<real> getBoundingBox() const;
private:
  real sigma_, epsilon_, offset_;
  Vec3<real> xmin_, xmax_;
};


class LJ_6_12_offset_cylinder : public FF_function{
public:
    LJ_6_12_offset_cylinder(Vec3<real> x_ref, Vec3<real> size, real cutoff, Vec<real> params);
    real compute(Vec3<real> x) const;
    virtual Vec<real> getBoundingBox() const;
private:
  real sigma_, epsilon_, offset_, height_;
  int axis_;
};

class LJ_3_9_offset : public FF_function{
public:
    LJ_3_9_offset(Vec3<real> x_ref, Vec3<real> size, real cutoff, Vec<real> params);
    real compute(Vec3<real> x) const;
    virtual Vec<real> getBoundingBox() const;
private:
  real sigma_, epsilon_, offset_; 
};

class LJ_3_9_offset_box : public FF_function{
public:
    LJ_3_9_offset_box(Vec3<real> x_ref, Vec3<real> size, real cutoff, Vec<real> params);
    real compute(Vec3<real> x) const;
    virtual Vec<real> getBoundingBox() const;
private:
  real sigma_, epsilon_;
  Vec3<real> xmin_, xmax_;
};


class LJ_3_9_offset_cylinder : public FF_function{
public:
    LJ_3_9_offset_cylinder(Vec3<real> x_ref, Vec3<real> size, real cutoff, Vec<real> params);
    real compute(Vec3<real> x) const;
    virtual Vec<real> getBoundingBox() const;
private:
  real sigma_, epsilon_, offset_, height_;
  int axis_;
};

class LJ_3_9_offset_ellipse : public FF_function{
public:
    LJ_3_9_offset_ellipse(Vec3<real> x_ref, Vec3<real> size, real cutoff, Vec<real> params);
    real compute(Vec3<real> x) const;
    virtual Vec<real> getBoundingBox() const;
private:
  real sigma_, epsilon_, offset_;
  //LAMMPS outputs a quaternion which should be convertible to a rotation matrix that orients a point in world-space to the object-space
  Vec3<Vec3<real>> rotationMatrix_;
  Vec3<real> shape_;
  int axis_;
};


//real (*funct_ptr)(Vec3<real>, Vec3<real>, Vec<real>)
//params: epsilon, sigma, cutoff
//real LJ_6_12_default(Vec3<real> x, Vec3<real> x_ref, Vec3<real> size, real cutoff, Vec<real> params);
//real LJ_6_12_shiftscale(Vec3<real> x, Vec3<real> x_ref, Vec3<real> size, real cutoff, Vec<real> params);


//real LJ_6_12_offset(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params);
//real LJ_6_12_offset_box(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params);
//real LJ_6_12_offset_cylinder(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params);

//real LJ_3_9_offset(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params);
//real LJ_3_9_offset_box(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params);
//real LJ_3_9_offset_cylinder(Vec3<real> x, Vec3<real> xref, Vec3<real> size, real cutoff, Vec<real> params);