#ifndef Material_h
#define Material_h 1

#include <string>
#include <iostream>
#include <vector>
#include "Geometry.hh"

class Material : public Geometry
{
public:
  Material(std::mt19937 MT, int ID, std::string MATFILE, std::string TYPE, double INDEX, double ATTLEN, double SCATLEN, double SCIPROB, double CHEPROB);
  virtual ~Material();
  int Type(){return type;}
  int ID(){return id;}
  double Index(){return index;}
  double AttLen(){return attlen;}
  double ScatLen(){return scatlen;}
  double SciProb(){return sciprob;}
  double CheProb(){return cheprob;}

private:
  int id, type;
  double index, attlen, scatlen, sciprob, cheprob;
};

#endif
