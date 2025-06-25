#ifndef Simulator_h
#define Simulator_h 1

#include "Global.hh"
#include "Source.hh"
#include "Material.hh"
#include "Charged.hh"
#include "Gui.hh"
#include <string>
#include <iostream>
#include <vector>
#include <random>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TView.h>
#include <TPolyLine3D.h>

class Simulator
{
public:
  Simulator(std::mt19937 MT, Config config, std::string OUTPUT, Source *SRC, std::vector<Material*> &MAT);
  virtual ~Simulator();
  void Run();
  void SetDisplay(Gui *GUI);
  Material* GetMaterial(int i){return mat[i];};
  int GetNMat(){return mat.size();};
  void SetNevt(int NEVT){nevt=NEVT;};
  Source* GetSource(){return src;};
  std::array<int, 5> GetResult(){return count;};
  void check_defined(const std::string& input);
  
private:
  int PointMaterial(const Position& p);
  void Initialize();
  int FType(int bt);
  void Summary();
  bool Fresnel(const Direction& v, const Direction& p, Direction& newv, Direction& newp, Direction norm, double idx_in, double idx_out);
  void Specular(const Direction& v, const Direction& p, Direction& newv, Direction& newp, Direction norm);
  void Lambert(const Direction& v, Direction& newv, Direction& newp, Direction norm);
  void Rayleigh(const Direction& v, const Direction& p, Direction& newv, Direction& newp);
  void Mie(const Direction& v, const Direction& p, Direction& newv, Direction& newp);
  std::mt19937 mt;
  std::uniform_real_distribution<double> unirand;
  std::string output;
  Source *src;
  Charged *chg;
  Gui *gui=nullptr;
  std::vector<Material*> mat;
  int nevt, nph;
  double index0, yield, delay, wlmin, wlmax,gmie;
  std::string sci_type;
  Position pos,vec,pol,cross,cand,newpos;
  Direction newvec,newpol,normal;
  double index,newindex,attlen,newattlen,scatlen,newscatlen;
  int matid, newmatid, btype, mn, newmn;
  double pl,apl,spl;
  bool displaymode,act_scinti,act_cherenkov,usemie;
  std::array<int, 5> count;
  double x_min, x_max, y_min, y_max, z_min, z_max, r_max;
  TFile *file;
  TTree *tree;
  TTree *charged;
  double ipos[3],fpos[3],ivec[3],fvec[3],ipol[3],fpol[3],time,length; //Branch variables for TTree *tree
  int imat, fmat, ftype, nref, npas, id;              //Branch variables for TTree *tree
};

#endif
