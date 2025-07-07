#include "Global.hh"

Position CustomSource(std::mt19937 mt){//User's custom function for primary particle position
  std::uniform_real_distribution<double> uniform;
  std::normal_distribution<double> gauss;
  std::exponential_distribution<double> expo;
  Position p;
  p[0]=-1+2*uniform(mt);
  p[1]=-1+2*uniform(mt);
  p[2]=-1+2*uniform(mt);
  return p;
}

Direction CustomDirection(std::mt19937 mt){//User's custom function for primary particle direction
  std::uniform_real_distribution<double> uniform;
  std::normal_distribution<double> gauss;
  std::exponential_distribution<double> expo;  
  Direction v;
  double ang=uniform(mt)*2*pi;
  v[0]=cos(ang);
  v[1]=sin(ang);
  v[2]=0;
  return v;
}

void ReadInFile(const char *inpath, Config &vars){//Read out the input card file.
  std::ifstream infile(inpath);
  if(!infile){
    std::cerr<<"Error: Input card file "<<inpath<<" is not found."<<std::endl;
    exit(1);
  }
  char c;
  std::string rest,section,key;
  while (infile && (infile >> std::ws) && (c = infile.peek())){
    if (c == '[' && infile.ignore()){
      if (infile.peek() == '/'){
	section = "";
      }
      else{
	std::getline(infile, section, ']');
      }
      std::getline(infile,rest);
    }
    else if (c == '#')
      std::getline(infile,rest);
    else if (section != ""){
      infile >> key;
      std::getline(infile,rest);
      if (infile){
	std::string::size_type l = rest.find('#');
	if (l == std::string::npos)
	  vars[section][key] = rest;
	else
	  vars[section][key] = rest.substr(0,l);
      }
    }
    else
      std::getline(infile,rest);
  }
}
void Normalize(Direction& v){//Normalize the vector norm.
  double norm=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  for(int i=0;i<3;i++){
    v[i]/=norm;
  }
}
double Dot(const Direction& a, const Direction& b){
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
void Cross(const Direction& a, const Direction& b, Direction& out){
  out[0] = a[1]*b[2] - a[2]*b[1];
  out[1] = a[2]*b[0] - a[0]*b[2];
  out[2] = a[0]*b[1] - a[1]*b[0];
}
void Isotropic(std::mt19937 mt, Direction& v){//Randomely and isotropically determine the direction
  std::uniform_real_distribution<double> uniform;
  v[2]=1-2*uniform(mt);
  double ang=uniform(mt)*2*pi;
  double vr=sqrt(1-v[2]*v[2]);
  v[0]=vr*cos(ang);
  v[1]=vr*sin(ang);
  Normalize(v);
}
void RandomPolarization(std::mt19937 mt, const Direction& v, Direction& p) {
  std::uniform_real_distribution<double> uniform;
  // Construct an arbitrary vector perpendicular to v
  Direction base;
  if (fabs(v[2]) < 0.99) {
    base = {-v[1], v[0], 0.0};  // Not parallel to z-axis
  } else {
    base = {0.0, -v[2], v[1]};  // For nearly z-axis directions
  }
  Normalize(base);

  // Orthogonal vector to complete the basis
  Direction ortho;
  Cross(v, base, ortho);
  Normalize(ortho);

  // Random angle on the plane orthogonal to v
  double phi = 2.0 * pi * uniform(mt);
  for (int i = 0; i < 3; ++i) {
    p[i] = cos(phi) * base[i] + sin(phi) * ortho[i];
  }
  Normalize(p);
}
double Distance(const Position& p1, const Position& p2)
{
  double dx = p1[0] - p2[0];
  double dy = p1[1] - p2[1];
  double dz = p1[2] - p2[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}
void Compare(double &A_max, double &A_min, double A){//Update the minimum and maximum points of a coordinate
  if(A_min > A) A_min = A;
  if(A_max < A) A_max = A;
}
