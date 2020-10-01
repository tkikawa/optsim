#include "Global.hh"

void ReadInFile(const char *inpath, Config &vars){
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
