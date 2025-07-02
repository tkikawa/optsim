#ifndef Material_h
#define Material_h 1

#include <string>
#include <vector>
#include "Geometry.hh"

class Material : public Geometry
{
public:
  Material(std::mt19937& MT, int ID, std::string MATFILE, std::string TYPE, double INDEX, double ATTLEN, double SCATLEN);
  virtual ~Material();
  virtual int Type() const {return type;}
  virtual int TType() const {return type;}
  int ID() const {return id;}
  double Index() const {return index;}
  double AttLen() const {return attlen;}
  double ScatLen() const {return scatlen;}
  virtual bool Is_Scinti() const {return is_scinti;}

private:
  int id, type;
  double index, attlen, scatlen;
  bool is_scinti;
};

#endif
