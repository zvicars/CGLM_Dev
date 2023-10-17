//Interface to several marching cubes algorithms
//Core functionality is to take a voxel grid and convert it to a series of triangles that can be analyzed/outputted
//Author: Zachariah Vicars
//Marching cubes algorithms come from various sources
#pragma once
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <iostream>
class VoxelGrid
{
	using real = double;
public:
	template <class T>
	using Vec3 = std::array<T, 3>;
	VoxelGrid(){
		return;
	}
	VoxelGrid(Vec3<int> size, double isovalue);
	void initialize(Vec3<int> size, double isovalue);
	int resize_grid(int dim_x, int dim_y, int dim_z);
	Vec3<int> getSize() const{
		return sz;
	}
	double getGridVal(int i, int j, int k) const{
		return grid_density_[i][j][k];
	}
	double getIsovalue() const{
		return isovalue_;
	}
	void sumInPlace(const VoxelGrid& vg);
	double getTot() const;
	Vec3<double> getWeightedCOM(){
		Vec3<double> sum = {0.0, 0.0, 0.0};
		for(int i = 0; i < sz[0]; i++){
			for(int j = 0; j < sz[1]; j++){
				for(int k = 0; k < sz[2]; k++)
				{
					sum[0] += grid_density_[i][j][k] * (i+0.5);
					sum[1] += grid_density_[i][j][k] * (j+0.5);
					sum[2] += grid_density_[i][j][k] * (k+0.5);
				}
			}
		}
		return sum;
	}	
	void scalarMult(double rhs); // compound assignment (does not need to be a member,
	void pbcidx(int &x, int& y, int& z);
	void add_discrete(Vec3<std::size_t> x_in, real amount){
		grid_density_[x_in[0]][x_in[1]][x_in[2]] += amount;
	};
	void clear();
private:
	std::vector<std::vector<std::vector<double > > > grid_density_;
	Vec3<int> sz;
	double isovalue_;
};