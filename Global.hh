#ifndef Global_h
#define Global_h 1

#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <random>

static const long double pi = 3.1415926535897932384626L; //Circular constant
static const long double conv = 3.1415926535897932384626L/180.L; //Deg to rad conversion factor
static const long double c_0 = 299.792458L; //Light speed [mm/ns]
static const double world=1e4; //World size in +/- X, Y, Z (mm) [region in which simulation is done.]
static const double nano=1e-6;
static const double micro=1e-3;
static const double cadunit=1; //Unit of CAD files [1 for mm, 10 for cm, 1000 for m]
static const double nreflimit=1e5; //Limit of the number of reflections
static const double hsize=1000;//Horizontal size of event display
static const double vsize=1000;//Vertical size of event display
static const double n_cosmic=2;//Cosmic ray angular distribution ~ cos^n_cosmic(theta)

using Config = std::map<std::string, std::map<std::string, std::string>>;
using Position = std::array<double, 3>;
using Direction = std::array<double, 3>;

Position CustomSource(std::mt19937 mt);
Direction CustomDirection(std::mt19937 mt);
void ReadInFile(const char *inpath, Config &vars);//Read variables from input card file into map
void Normalize(Direction& v);
void Isotropic(std::mt19937 mt, Direction& v);
void RandomPolarization(std::mt19937 mt, const Direction& v, Direction& p);
double Dot(const Direction& a, const Direction& b);
void Cross(const Direction& a, const Direction& b, Direction& out);
double Distance(const Position& p1, const Position& p2);
void Compare(double &A_max, double &A_min, double A);

#endif
