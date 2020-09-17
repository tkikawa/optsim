#ifndef Global_h
#define Global_h 1

#include <iostream>
#include <string>
#include <map>
#include <fstream>

static const long double pi = 3.1415926535897932384626L; //circular constant
static const long double conv = 3.1415926535897932384626L/180.L; //deg to rad conversion factor
static const long double c_0 = 299.792458L; ///< light speed [mm/ns]
static const double world=1e4;// world size in +/- X, Y, Z (mm) [region in which simulation is done.]
static const double micro=1e-6;
static const double cadunit=1;// Unit of CAD files [1 for mm, 10 for cm, 1000 for m]

typedef std::map<std::string, std::map<std::string, std::string> > Config;

void ReadInFile(const char *inpath, Config &vars);//read variables from input card file into map
void VectNorm(double *vec);

#endif
