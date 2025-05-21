#ifndef Charged_h
#define Charged_h 1

#include "Global.hh"
#include "Material.hh"
#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <float.h>

class Charged
{
public:
  Charged(std::mt19937 MT, std::vector<Material*> &MAT);
  virtual ~Charged();
  void Simulate(const Position& pos, const Direction& vec);
  void Generate(Position& pos, Direction& vec, int n);
  int GetNPhotons();
  
private:
  double Distance(const Position& p1, const Position& p2);
  Position Interpolate(const Position& p1, const Position& p2, double t);
  Direction Isotropic();
  Direction CherenkovDir(double index, const Direction& vec);
  std::mt19937 mt;
  std::uniform_real_distribution<double> unirand;
  std::vector<Material*> mat;
  Position crs,ptn,far;
  Direction dir;
  std::vector<Position> cross;
  std::vector<std::pair<Position, Direction>> particle;
  double p_sci,p_che,dist,s;
  int n_sci,n_che;
};

#endif
