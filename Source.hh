#ifndef Source_h
#define Source_h 1

#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include "Global.hh"
#include "Geometry.hh"

class Source : public Geometry
{
public:
  Source(std::mt19937& MT, Config config);
  virtual ~Source();
  void Generate(Position& pos, Direction& vec, Direction& pol);
  Position PointInSource();
  Direction ParticleDir();
  Direction Polarize(const Direction& v);
  bool ChargedMode();
  double GetBeta();
  bool IsElectron();

private:
  void Compare(double &A_max, double &A_min, double A);
  double CalcBeta(double m, double T);
  std::string sourcemode;
  std::string particlemode;
  std::string directionmode;
  bool charged;
  double x_min, x_max, y_min, y_max, z_min, z_max, x_center, y_center, z_center;
  double r_min, r_max, phi_min, phi_max, r;
  double topsurf, insurf, outsurf, xsurf, ysurf, zsurf, totsurf;
  std::string sourcefile;
  double cost,sint,cosp,sinp, rz, beta,polar;
  bool is_electron;
};

#endif
