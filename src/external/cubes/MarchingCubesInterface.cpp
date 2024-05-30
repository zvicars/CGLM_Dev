#include "MarchingCubesInterface.hpp"
#include <sstream>
#include "golosio/TriangulateGlobal.hpp"

void golosio_mc(const VoxelGrid& v, Mesh& output_mesh){
	int nx = v.getSize()[0];
	int ny = v.getSize()[1];
	int nz = v.getSize()[2];
	double box_vec[3] = {nx, ny, nz};
	std::vector<float> volume_data(nx*ny*nz, 0);
	std::vector<int> triangle_data;
	std::vector<float> vertex_data; 
	int iterator = 0;
	double cell_volume = 1.0;
	double density = 1.0;
	for(int k = 0; k < nz; k++)
	for(int j = 0; j < ny; j++)
	for(int i = 0; i < nx; i++)
	{
		volume_data[iterator] = v.getGridVal(i, j, k)/( cell_volume * density );
		iterator++;
	}
	int nvtx; 
	int ntri; 
	int n[3] = {nx, ny, nz};

	float s[3] = {1.0, 1.0, 1.0}; 
	float thresh = (float)v.getIsovalue(); 
	TriangulateGlobal g1;
	g1.Triangulate(&volume_data[0], vertex_data, triangle_data, n, s, &nvtx, &ntri, thresh, 0);
	std::vector<Vec3<double> > vertices_(nvtx);
	std::vector<Vec3<double> > gradients_(nvtx);
	std::vector<Triangle> triangles_(ntri);
	for(int i = 0; i < nvtx; i++){
			for(int j = 0; j < 3; j++){
					vertices_[i][j] = vertex_data[6*i + j] + 0.5*n[j];
					gradients_[i][j] = vertex_data[6*i + j + 3];   
			}
	}
	for(int i = 0; i < ntri; i++){
		Vec3<int> triangle_indices_;
		for(int j = 0; j < 3; j++){
				triangle_indices_[j] = triangle_data[3*i+j];   
		}
		triangles_[i].indices = triangle_indices_;
	}
	output_mesh.normals = gradients_;
	output_mesh.vertices = vertices_;
	output_mesh.triangles = triangles_;
	output_mesh.nvtx = nvtx;
	output_mesh.ntri = ntri;
	return;
}

void marchingCubes(std::string type, const VoxelGrid& v, Mesh& output_mesh)
{
	if(type == "golosio"){
		golosio_mc(v, output_mesh);
	}
	else{
		std::cout << "Invalid marching cubes type selected." << std::endl;
		throw 0;
	}
	return;
}

void printSTL(const Mesh& mesh, std::string& frame)
{
	std::stringstream ofile;
	ofile << "solid " << "frame" << "\n";
	for(int i = 0; i < mesh.ntri; i++)
	{
        int vidx1 = mesh.triangles[i].indices[0];
        int vidx2 = mesh.triangles[i].indices[1];
        int vidx3 = mesh.triangles[i].indices[2];
        //just using the first vertex normal as the face normal, lazy but shouldn't affect much
        ofile << "facet normal " << mesh.normals[vidx1][0] << "  " << mesh.normals[vidx1][1] << "  " << mesh.normals[vidx1][2] << "\n"; 
		ofile << "    outer loop\n";
		ofile << "vertex " << mesh.vertices[vidx1][0] << " " << mesh.vertices[vidx1][1] << " " << mesh.vertices[vidx1][2] << "\n";
 		ofile << "vertex " << mesh.vertices[vidx2][0] << " " << mesh.vertices[vidx2][1] << " " << mesh.vertices[vidx2][2] << "\n";
		ofile << "vertex " << mesh.vertices[vidx3][0] << " " << mesh.vertices[vidx3][1] << " " << mesh.vertices[vidx3][2] << "\n";       
		ofile << "    endloop\n";
		ofile << "endfacet\n";
	}
	ofile << "endsolid " << frame << std::endl;
	frame = ofile.str();
  return;
}

std::string printPLYWithRGB(const Mesh& mesh, const std::vector<RGB>& cv){
  std::stringstream ss;
  ss << "ply\nformat ascii 1.0\n";
  ss << "element vertex " << mesh.nvtx << "\n";
  ss << "property float x\n";
  ss << "property float y\n";
  ss << "property float z\n";
  ss << "property uchar red\n";
  ss << "property uchar green\n";
  ss << "property uchar blue\n";
  ss << "element face " << mesh.ntri << "\n";
  ss << "property list uchar int vertex_index\n";
  ss << "end_header\n";
  for(int i = 0; i < mesh.nvtx; i++){
    ss << mesh.vertices[i][0] << " " << mesh.vertices[i][1] << " " << mesh.vertices[i][2] << " ";
    std::array<unsigned char,3> channels;
    channels[0] = cv[i].r;
    channels[1] = cv[i].g;
    channels[2] = cv[i].b;
    for(int j = 0; j < 3; j++){
      ss << (int)channels[j] << " ";
    }
    ss << "\n";
  }
  for(int i = 0; i < mesh.ntri; i++){
    ss << "3 " << mesh.triangles[i].indices[0] << " " << mesh.triangles[i].indices[1] << " " << mesh.triangles[i].indices[2] << "\n";
  }

  return ss.str();
}