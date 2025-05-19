#ifndef Source_h
#define Source_h 1

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <math.h>
#include "Global.hh"
#include "Geometry.hh"
#include "Material.hh"


class Source : public Geometry
{
public:
  Source(std::mt19937 MT, Config config, std::vector<Material*> &MAT);
  virtual ~Source();
  void Generate(double *pos, double *vec);
  void PointInSource(double *pos);
  void Direction(double *vec);

private:
  void Compare(double &A_max, double &A_min, double A);
  void Normalize(double &vx, double &vy, double &xz);
  std::vector<Material*> mat;
  std::string sourcemode;
  std::string directionmode;
  double x_min, x_max, y_min, y_max, z_min, z_max;
  double r_min, r_max, phi_min, phi_max;
  double topsurf, insurf, outsurf, totsurf;
  std::string sourcefile;
  double v_x, v_y, v_z, phi;
  double cost,sint,cosp,sinp;
  double v[3],rz,vr,ang;
};

#endif
