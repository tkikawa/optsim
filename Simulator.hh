#ifndef Simulator_h
#define Simulator_h 1

#include "Global.hh"
#include "Source.hh"
#include "Material.hh"
#include "Charged.hh"
#include <string>
#include <iostream>
#include <vector>
#include <random>
#include <TFile.h>
#include <TTree.h>

#include "TCanvas.h"
#include "TView.h"
#include "TPolyLine3D.h"

class Simulator
{
public:
  Simulator(std::mt19937 MT, Config config, std::string OUTPUT, Source *SRC, std::vector<Material*> &MAT);
  virtual ~Simulator();
  void Run();
  void Display();
  
private:
  int PointMaterial(const Position& p);
  void Initialize();
  int FType(int bt);
  void Summary();
  bool Fresnel(const Direction& v, Direction& newv, Direction norm, double idx_in, double idx_out);
  void Specular(const Direction& v, Direction& newv, const Direction& norm);
  void Lambert(const Direction& v, Direction& newv, Direction norm);
  void Rayleigh(const Direction& v, Direction& newv);
  void Draw(double vtx[3][3], int type);
  void Track(const Position& p0, const Position& p1, bool charged = false);
  void Compare(double &A_max, double &A_min, double A);
  void Normalize(Direction& v);
  void Isotropic(Direction& v);
  std::mt19937 mt;
  std::uniform_real_distribution<double> unirand;
  std::string output;
  Source *src;
  Charged *chg;
  std::vector<Material*> mat;
  int nevt;
  int nph;
  double index0;
  Position pos,vec,cross,cand,newpos;
  Direction newvec,normal;
  double index,newindex,attlen,newattlen,scatlen,newscatlen;
  int matid, newmatid, btype, mn, newmn;
  double pl,apl,spl;
  bool displaymode;
  int count[5];
  double x_min, x_max, y_min, y_max, z_min, z_max, r_max;
  TFile *file;
  TTree *tree;
  double ipos[3],fpos[3],ivec[3],fvec[3],time,length;//Branch variables for TTree *tree
  int imat, fmat, ftype, nref, npas;                 //Branch variables for TTree *tree
};

#endif
