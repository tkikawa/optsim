#ifndef Charged_h
#define Charged_h 1

#include "TRandom3.h"
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
  void Generate(Position& pos, Position& vec, Position& pol, double& t, int n);
  int GetNPhotons();
  void SetBeta(double BETA);
  void SetScinti(std::string sci_type, double YIELD, double DELAY);
  void SetCherenkov(double WLMIN, double WLMAX);
    
private:
  double Distance(const Position& p1, const Position& p2);
  Position Interpolate(const Position& p1, const Position& p2, double t);
  void ScintiDir(Direction& vsci, Direction& psci);
  void CherenkovDir(double index, const Direction& vec, Direction& vche, Direction& pche);
  int NumPhoton(double d, double prob, bool is_scinti);
  void CalcSciPar();
  double CalcCheProb(double n);
  double ScintiDelay();
  std::mt19937 mt;
  std::uniform_real_distribution<double> unirand;
  std::vector<Material*> mat;
  Position crs,ptn,far;
  Direction dir, plr;
  double time;
  std::vector<Position> cross;
  std::vector<std::tuple<Position, Direction, Direction, double>> particle;
  double cheprob,dedx,width_pp,p_sci,p_che,dist,s;
  int n_sci,n_che;
  double wlmin, wlmax, yield, A, Z, I, density, beta, scinti_lifetime;
  bool act_sci,act_che;
  TRandom3 rng;
};

#endif
