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
  double Generate(Position& pos, Direction& vec, int n);
  int GetNPhotons();
  void SetBeta(double BETA);
  void SetScinti(std::string sci_type, double YIELD);
  void SetCherenkov(double WLMIN, double WLMAX);
    
private:
  double Distance(const Position& p1, const Position& p2);
  Position Interpolate(const Position& p1, const Position& p2, double t);
  Direction Isotropic();
  Direction CherenkovDir(double index, const Direction& vec);
  int NumPhoton(double d, double prob, bool is_scinti);
  void CalcSciPar();
  double CalcCheProb(double n);
  std::mt19937 mt;
  std::uniform_real_distribution<double> unirand;
  std::vector<Material*> mat;
  Position crs,ptn,far;
  Direction dir;
  double time;
  std::vector<Position> cross;
  std::vector<std::tuple<Position, Position, double>> particle;
  double cheprob,dedx,width_pp,p_sci,p_che,dist,s;
  int n_sci,n_che;
  double wlmin, wlmax, yield,A,Z,I,density,beta;
  bool act_sci,act_che;
  TRandom3 rng;
};

#endif
