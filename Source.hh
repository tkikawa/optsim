#ifndef Source_h
#define Source_h 1

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include "Global.hh"
#include "Geometry.hh"


class Source : public Geometry
{
public:
  Source(std::mt19937 MT, Config config);
  virtual ~Source();
  void Generate(Position& pos, Direction& vec);
  Position PointInSource();
  Direction ParticleDir();
  bool ChargedMode();
  double GetBeta();

private:
  void Compare(double &A_max, double &A_min, double A);
  void Normalize(double &vx, double &vy, double &xz);
  double CalcBeta(double m, double T);
  std::string sourcemode;
  std::string particlemode;
  std::string directionmode;
  bool charged;
  double x_min, x_max, y_min, y_max, z_min, z_max, x_center, y_center, z_center;
  double r_min, r_max, phi_min, phi_max, r;
  double topsurf, insurf, outsurf, xsurf, ysurf, zsurf, totsurf;
  std::string sourcefile;
  double v_x, v_y, v_z;
  double cost,sint,cosp,sinp;
  double rz,vr,ang;
  double mass,energy,beta;
};

#endif
