#ifndef Geometry_h
#define Geometry_h 1

#include <string>
#include <iostream>
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
  Geometry();
  virtual ~Geometry();

public:
  std::mt19937 mt;
  std::uniform_real_distribution<double> unirand;
  std::normal_distribution<double> gausrand;
  std::vector<Triangle> triangle;

  /*
  Assimp::Importer importer;
  const aiScene* scene;
  aiMesh* aim;
  */
  
  void LoadCAD(std::string name);
  int NTriangle(){return triangle.size();}
  Triangle GetTriangle(int n){return triangle[n];}
  bool InSolid(double pos[3]);
};

#endif
