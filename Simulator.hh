#ifndef Simulator_h
#define Simulator_h 1

#include "Source.hh"
#include "Material.hh"
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
  int PointMaterial(double p[3]);
  void Initialize();
  int FType(int bt);
  void Summary();
  bool Fresnel(double v[3], double *newv, double norm[3], double idx_in, double idx_out);
  void Specular(double v[3], double *newv, double norm[3]);
  void Lambert(double v[3], double *newv, double norm[3]);
  void Rayleigh(double v[3], double *newv);
  void Draw(double vtx[3][3], int type);
  void Track(double p0[3], double p1[3]);
  void Compare(double &A_max, double &A_min, double A);
  void Normalize(double *v);
  void Isotropic(double *v);
  std::mt19937 mt;
  std::uniform_real_distribution<double> unirand;
  std::string output;
  Source *src;  
  std::vector<Material*> mat;
  int nevt;
  double index0;
  double pos[3],vec[3],cross[3],cand[3],newpos[3],newvec[3],normal[3],index,newindex,attlen,newattlen,scatlen,newscatlen;
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
