#ifndef Material_h
#define Material_h 1

#include <string>
#include <iostream>
#include <vector>
#include "Geometry.hh"

class Material : public Geometry
{
public:
  Material(std::mt19937 MT, int ID, std::string MATFILE, std::string TYPE, double INDEX, double ATTLEN, double SCATLEN);
  virtual ~Material();
  virtual int Type(){return type;}
  virtual int TType(){return type;}
  int ID(){return id;}
  double Index(){return index;}
  double AttLen(){return attlen;}
  double ScatLen(){return scatlen;}
  virtual bool Is_Scinti(){return is_scinti;}

private:
  int id, type;
  double index, attlen, scatlen;
  bool is_scinti;
};

#endif
