#include "Mixture.hh"

Mixture::Mixture(std::mt19937& MT, int ID, std::string MATFILE, std::vector<std::pair<double, std::string>> MIX, double INDEX=1, double ATTLEN=0, double SCATLEN=0)
  : Material(MT, ID, MATFILE, "mixture", INDEX, ATTLEN, SCATLEN),
    is_scinti(false)
{
  double total=0;
  for(unsigned int i=0;i<MIX.size();i++){
    total+=MIX[i].first;
  }
  int tp;
  type=7;
  std::cout<<"The mixture consists of ";
  for(unsigned int i=0;i<MIX.size();i++){
    if(MIX[i].second=="medium")         tp = 0;
    else if(MIX[i].second=="converter") tp = 1;
    else if(MIX[i].second=="mirror")    tp = 2;
    else if(MIX[i].second=="diffuser")  tp = 3;
    else if(MIX[i].second=="absorber")  tp = 4;
    else if(MIX[i].second=="detector")  tp = 5;
    else if(MIX[i].second=="scintillator"){
      tp = 0;
      is_scinti=true;
    }
    else{
      std::cerr<<std::endl;
      std::cerr<<"Error: Material type: "<<MIX[i].second<<" for ID="<<ID<<" is invalid."<<std::endl;
      std::cerr<<"Material type must be any one of medium, mirror, diffuser, absober or detector."<<std::endl;
      exit(1);
    }
    std::cout<<MIX[i].first/total*100<<"% of "<<MIX[i].second;
    if(i<MIX.size()-1){
      std::cout<<", ";
    }
    mix.push_back(std::make_pair(MIX[i].first/total,tp));
    if(tp<2)type=6;//Mixture includes medium or converter
  }
  std::cout<<"."<<std::endl;
}
Mixture::~Mixture()
{
}
int Mixture::Type(){
  double r = unirand(mt);
  double cumulative = 0.0;
  for (unsigned int i=0;i<mix.size();i++) {
    cumulative += mix[i].first;
    if (r < cumulative) {
      return mix[i].second;
    }
  }
  return mix.back().second;
}
