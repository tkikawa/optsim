#ifndef Triangle_h
#define Triangle_h 1

#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <float.h>
#include "Global.hh"

class Triangle
{
public:
  Triangle(double vtx[3][3]);
  virtual ~Triangle();
  double X(int n){return vertex[n][0];}
  double Y(int n){return vertex[n][1];}
  double Z(int n){return vertex[n][2];}
  double A(){return a;}
  double B(){return b;}
  double C(){return c;}
  double D(){return d;}
  double GetArea(){return area;}
  bool Collision(double s[3], double t[3], double *p);
  void GetSurfPoint(double *p, double r1, double r2);
  void GetNormal(double *norm);

private:
  double vertex[3][3],normal[3];
  double a,b,c,d;// Plane equation: ax+by+cz+d=0
  double area;// Surface area
};

#endif
