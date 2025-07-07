#ifndef Geometry_h
#define Geometry_h 1

#include <string>
#include <vector>
#include <random>
#include "Global.hh"
#include "Triangle.hh"

#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags

class Geometry
{
public:
  Geometry(std::mt19937& MT);
  virtual ~Geometry();
  void LoadCAD(std::string name);
  int NTriangle(){return triangle.size();}
  Triangle GetTriangle(int n){return triangle[n];}
  bool InSolid(const Position& pos);
  double Round(double p0);
  bool IntersectsAABB(const Position& origin, const Position& end);
  bool InAABB(const Position& pos);
  std::mt19937& mt;
  std::vector<Triangle> triangle;
  std::uniform_real_distribution<double> unirand;
  std::normal_distribution<double> gausrand;
  Position box_max, box_min;
};

#endif
