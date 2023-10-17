#include "VoxelGrid.hpp"
#include <omp.h>
VoxelGrid::VoxelGrid(Vec3<int> size, double isovalue) //must be orthorhombic box
{
  initialize(size, isovalue);
  return;
}
void VoxelGrid::initialize(Vec3<int> size, double isovalue){
  int dx = size[0];
  int dy = size[1];
  int dz = size[2];
  resize_grid(dx, dy, dz);
  isovalue_ = isovalue;
  return;
}
int VoxelGrid::resize_grid(int dim_x, int dim_y, int dim_z)
{
    grid_density_.clear();
    grid_density_.resize(dim_x,std::vector<std::vector<double> >(dim_y,std::vector<double>(dim_z,0.0)));
    sz[0] = dim_x; sz[1]= dim_y; sz[2] = dim_z;
    return 1;
}
void VoxelGrid::pbcidx(int &x, int& y, int& z)
{
    if(x >= sz[0]) x = x%sz[0];
    else if(x < 0) x = sz[0] - (-x)%sz[0];
    if(y >= sz[1]) y = y%sz[1];
    else if(y < 0) y = sz[1] - (-y)%sz[1];
    if(z >= sz[2]) z = z%sz[2];
    else if(z < 0) z = sz[2] - (-z)%sz[2];
    return;
}

void VoxelGrid::clear(){
    #pragma omp for collapse(3)
    for(int i = 0; i < sz[0]; i++){
        for(int j = 0; j < sz[1]; j++){
            for(int k = 0; k < sz[2]; k++){
                grid_density_[i][j][k] = 0.0; 
            }
        }
    }
    return;
}

void VoxelGrid::scalarMult(double rhs) // compound assignment (does not need to be a member,
{                           // but often is, to modify the private members)
    #pragma omp parallel for collapse(3)
    for(int i = 0; i < sz[0]; i++){
        for(int j = 0; j < sz[1]; j++){
            for(int k = 0; k < sz[2]; k++){
                grid_density_[i][j][k] *= rhs;
            }
        }
    }
    return;
}

void VoxelGrid::sumInPlace(const VoxelGrid& vg){
    #pragma omp parallel for collapse(3)
    for(int i = 0; i < sz[0]; i++){
        for(int j = 0; j < sz[1]; j++){
            for(int k = 0; k < sz[2]; k++){
                grid_density_[i][j][k] += vg.getGridVal(i,j,k);
            }
        }
    }
    return;
}

double VoxelGrid::getTot() const{
    double sum = 0.0;
    for(int i = 0; i < sz[0]; i++){
        for(int j = 0; j < sz[1]; j++){
            for(int k = 0; k < sz[2]; k++)
            {
                sum += grid_density_[i][j][k];
            }
        }
    }
    return sum;
}