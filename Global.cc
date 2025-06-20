#include "Global.hh"

Position CustomSource(std::mt19937 mt){//User's custom function for primary particle position
  std::uniform_real_distribution<double> uniform;
  std::normal_distribution<double> gauss;
  std::exponential_distribution<double> exp;
  Position p;
  p[0]=0;
  p[1]=0;
  p[2]=0;
  return p;
}

Direction CustomDirection(std::mt19937 mt){//User's custom function for primary particle direction
  std::uniform_real_distribution<double> uniform;
  std::normal_distribution<double> gauss;
  std::exponential_distribution<double> exp;
  Direction v;
  v[0]=0;
  v[1]=0;
  v[2]=1;
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
