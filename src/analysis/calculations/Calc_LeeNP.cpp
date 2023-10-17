#include "Calc_LeeNP.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/GenerateBitmap.hpp"
#include <cstring>
Calc_LeeNP::Calc_LeeNP(AnalysisInputPack& input) : Calculation{input}{
  axis_ = 2;
  input.params().readNumber("axis", ParameterPack::KeyType::Optional, axis_);
  axis1_ = (axis_ + 1) % 3;
  axis2_ = (axis_ + 2) % 3;
  min_val_ = 0.0;
  input.params().readNumber("min_val", ParameterPack::KeyType::Optional, axis_);
  FANCY_ASSERT(axis_ >= 0 && axis_ <= 3, "axis needs to be 0, 1, or 2 for x, y, or z");
  std::vector<std::string> calculations;
  input.params().readVector("calculations", ParameterPack::KeyType::Optional, calculations);
  for(auto calc : calculations){
    Calculation* calc_object = input.findCalculation(calc);
    Calc_Isosurface* cast_calc = dynamic_cast<Calc_Isosurface*>(calc_object);
    FANCY_ASSERT(cast_calc != 0, "Isosurface calculation not found");
    isosurface_calculations_.push_back(cast_calc);
  }
  distances_.resize(isosurface_calculations_.size());
  return;
}
//returns the nearest distance to the mesh
real Calc_LeeNP::findNearestIntersection(Vec3<real> p0, Vec3<real> n0, const Mesh& mesh) const{
  real eval = std::numeric_limits<real>::max();
  int closest_triangle=-1;
  //std::cout << "num_triangles = " << mesh.triangles.size() << std::endl;
  for(int i = 0; i < mesh.triangles.size(); i++){
    Vec3<real> intersection_point;
    Vec3<Vec3<real> > points;
    auto indices = mesh.triangles[i].indices;
    for(int j = 0; j < indices.size(); j++){
      points[j] = mesh.vertices[indices[j]];
    }
    bool intersects = rayTriangleIntersection(p0, n0, points, intersection_point);
    if(!intersects) continue;
    //if(intersection_point[axis_] < min_val_) continue;
    //std::cout << "int point = " << intersection_point[0] << "  " << intersection_point[1] << "  " << intersection_point[2] << std::endl;
    real d = norm2(p0 - intersection_point);
    if(d < eval){
      eval = d;
      closest_triangle = i;
    }
  }
  return eval;
}

//Muller-Trombore stolen from Wikipedia
bool Calc_LeeNP::rayTriangleIntersection(Vec3<real> p0, Vec3<real> n0, Vec3<Vec3<real> > triangle, Vec3<real>& intersection_point) const{
  const real EPSILON = 0.0000001;
  Vec3<real> vertex0 = triangle[0];
  Vec3<real> vertex1 = triangle[1];  
  Vec3<real> vertex2 = triangle[2];
  Vec3<real> edge1, edge2, h, s, q;
  real a, f, u, v;
  edge1 = vertex1 - vertex0;
  edge2 = vertex2 - vertex0;
  h = cross(n0, edge2);
  a = dot(edge1, h);
  if (a > -EPSILON && a < EPSILON)
      return 0;    // This ray is parallel to this triangle.
  f = 1.0 / a;
  s = p0 - vertex0;
  u = f * dot(s,h);
  if (u < 0.0 || u > 1.0)
    return 0;
  q = cross(s, edge1);
  v = f * dot(n0, q);
  if (v < 0.0 || u + v > 1.0) return 0;
  // At this stage we can compute t to find out where the intersection point is on the line.
  real t = f * dot(edge2, q);
  if (t > EPSILON) // ray intersection
  {
    intersection_point = p0 + (t * n0);
    return 1;
  }
  else // This means that there is a line intersection but not a ray intersection.
    return 0;  
}

void Calc_LeeNP::calculate(){
  if(!doCalculate()) return;  
  for(int i = 0; i < isosurface_calculations_.size(); i++){
    if(!isosurface_calculations_[i]->hasCalculated()) isosurface_calculations_[i]->calculate();
    std::array<std::size_t, 2> idx = {box->size_[axis1_], box->size_[axis2_]};
    distances_[i].initialize(idx);
  }
  std::size_t size = internalLattice.size()[axis_];
  std::size_t size1 = internalLattice.size()[axis1_];
  std::size_t size2 = internalLattice.size()[axis2_];
  real test_pos = size+10;
  //#pragma omp parallel for unroll(2)
  for(std::size_t i = 0; i < size1; i++){
    for(std::size_t j = 0; j < size2; j++){
      std::array<std::size_t, 2> idx = {i,j};
      //determine test ray
      Vec3<real> p0; 
      p0[axis_] = test_pos; 
      p0[axis1_] = (real)i + 0.5; 
      p0[axis2_] = (real)j + 0.5;
      Vec3<real> dir; dir[axis_] = -1.0; dir[axis1_] = 0.0; dir[axis2_] = 0.0;
      for(int k = 0; k < isosurface_calculations_.size(); k++){
        auto& mesh = isosurface_calculations_[k]->getMesh();
        real distance = findNearestIntersection(p0, dir, mesh);
        distances_[k].at(idx) = distance;
      }
    }
  }
  return;
}


void Calc_LeeNP::writeImage(int width, int height, Matrix<unsigned char, 2>& r, 
                            Matrix<unsigned char, 2>& g, Matrix<unsigned char, 2>& b, std::string filename){
    unsigned char image[height][width][BYTES_PER_PIXEL];
    // assigning value to string s 
    const int length = filename.length(); 
    // declaring character array (+1 for null terminator) 
    char* char_array = new char[length + 1]; 
    strcpy(char_array, filename.c_str()); 
    std::size_t i, j;
    for (i = 0; i < height; i++) {
        for (j = 0; j < width; j++) {
          std::array<std::size_t, 2> idx = {j,i};
          image[i][j][2] = r.at(idx);   //red
          image[i][j][1] = g.at(idx);   //green
          image[i][j][0] = b.at(idx);   //blue
        }
    }
    generateBitmapImage((unsigned char*) image, height, width, char_array);
    return;
}


void Calc_LeeNP::finalOutput(){
  for(int i = 0; i < isosurface_calculations_.size(); i++){
    if(!isosurface_calculations_[i]->hasPerformedFinalOutput()) isosurface_calculations_[i]->finalOutput();
    std::array<std::size_t, 2> idx = {box->size_[axis1_], box->size_[axis2_]};
    distances_[i].initialize(idx);
  }
  std::size_t size = internalLattice.size()[axis_];
  std::size_t size1 = internalLattice.size()[axis1_];
  std::size_t size2 = internalLattice.size()[axis2_];
  real test_pos = size+10;
  #pragma omp parallel for unroll(2)
  for(std::size_t i = 0; i < size1; i++){
    for(std::size_t j = 0; j < size2; j++){
      std::array<std::size_t, 2> idx = {i,j};
      //determine test ray
      Vec3<real> p0; 
      p0[axis_] = test_pos; 
      p0[axis1_] = (real)i + 0.5; 
      p0[axis2_] = (real)j + 0.5;
      Vec3<real> dir; dir[axis_] = -1.0; dir[axis1_] = 0.0; dir[axis2_] = 0.0;
      for(int k = 0; k < isosurface_calculations_.size(); k++){
        auto& mesh = isosurface_calculations_[k]->getMesh();
        real distance = findNearestIntersection(p0, dir, mesh);
        distances_[k].at(idx) = distance;
      }
    }
  }

 std::string filepath = name_ + "_average.png";
  Matrix<unsigned char, 2> r, g, b;
  int axis1 = (axis_ + 1) % 3;
  int axis2 = (axis_ + 2) % 3;
  r.initialize({internalLattice.size()[axis1], internalLattice.size()[axis2]});
  g.initialize({internalLattice.size()[axis1], internalLattice.size()[axis2]});
  b.initialize({internalLattice.size()[axis1], internalLattice.size()[axis2]});

  Vec<Vec3<unsigned char> > color_map;
  color_map.push_back({255, 0, 0 }); //0
  color_map.push_back({0, 255, 0 }); //1
  color_map.push_back({0, 0, 255 }); //2
  color_map.push_back({255, 255, 0 }); //3 
  color_map.push_back({255, 0, 255 }); //4
  color_map.push_back({0, 255, 255 }); //5 
  color_map.push_back({255, 255, 255}); //6

  for(std::size_t i = 0; i < internalLattice.size()[axis1]; i++){
    for(std::size_t j = 0; j < internalLattice.size()[axis2]; j++){
      std::array<std::size_t, 2> idx = {i,j};
      real closestDistance = 1000;
      int closest_entry = 6; //default to white
      for(int k = 0; k < distances_.size(); k++){
        if(distances_.at(k).at(idx) < closestDistance){
          closestDistance = distances_.at(k).at(idx);
          closest_entry = k;
        }
      }
      r.at(idx) = color_map[closest_entry][0];
      g.at(idx) = color_map[closest_entry][1];
      b.at(idx) = color_map[closest_entry][2];
    }
  }
  writeImage(internalLattice.size()[axis1], internalLattice.size()[axis2], r, g, b, filepath);

  return;
}


void Calc_LeeNP::output(){
  if(!doOutput()) return;
  std::string filepath = name_ + "_" + std::to_string(box->frame_) + ".png";
  Matrix<unsigned char, 2> r, g, b;
  int axis1 = (axis_ + 1) % 3;
  int axis2 = (axis_ + 2) % 3;
  r.initialize({internalLattice.size()[axis1], internalLattice.size()[axis2]});
  g.initialize({internalLattice.size()[axis1], internalLattice.size()[axis2]});
  b.initialize({internalLattice.size()[axis1], internalLattice.size()[axis2]});

  Vec<Vec3<unsigned char> > color_map;
  color_map.push_back({255, 0, 0 }); //0
  color_map.push_back({0, 255, 0 }); //1
  color_map.push_back({0, 0, 255 }); //2
  color_map.push_back({255, 255, 0 }); //3 
  color_map.push_back({255, 0, 255 }); //4
  color_map.push_back({0, 255, 255 }); //5 
  color_map.push_back({255, 255, 255}); //6

  for(std::size_t i = 0; i < internalLattice.size()[axis1]; i++){
    for(std::size_t j = 0; j < internalLattice.size()[axis2]; j++){
      std::array<std::size_t, 2> idx = {i,j};
      real closestDistance = 1000;
      int closest_entry = 6; //default to white
      for(int k = 0; k < distances_.size(); k++){
        if(distances_.at(k).at(idx) < closestDistance){
          closestDistance = distances_.at(k).at(idx);
          closest_entry = k;
        }
      }
      r.at(idx) = color_map[closest_entry][0];
      g.at(idx) = color_map[closest_entry][1];
      b.at(idx) = color_map[closest_entry][2];
    }
  }
  writeImage(internalLattice.size()[axis1], internalLattice.size()[axis2], r, g, b, filepath);
  return;
}

void Calc_LeeNP::update(){
  Calculation::update();
  return;  
}