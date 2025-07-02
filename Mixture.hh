#ifndef Mixture_h
#define Mixture_h 1

#include <string>
#include <vector>
#include "Material.hh"

class Mixture : public Material
{
public:
  Mixture(std::mt19937& MT, int ID, std::string MATFILE, std::vector<std::pair<double, std::string>> MIX, double INDEX, double ATTLEN, double SCATLEN);
  virtual ~Mixture();
  int Type();
  int TType(){return type;}
  bool Is_Scinti(){return is_scinti;}

private:
  int type;
  bool is_scinti;
  std::vector<std::pair<double, int>> mix;
};

#endif
