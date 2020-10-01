#include "Material.hh"

Material::Material(std::mt19937 MT, int ID, std::string MATFILE, std::string TYPE, double INDEX=1, double ATTLEN=0)
{
  mt = MT;
  id = ID;
  if(id<=0){
    std::cerr<<"Error: Matrial ID is invalid."<<std::endl;
    std::cerr<<"Material ID must be positive integer."<<std::endl;
    exit(1);
  }
  index = INDEX;
  attlen = ATTLEN;
  if(TYPE=="medium")        type = 0;
  else if(TYPE=="mirror")   type = 1;
  else if(TYPE=="diffuser") type = 2;
  else if(TYPE=="absorber") type = 3;
  else if(TYPE=="detector") type = 4;
  else{
    std::cerr<<"Error: Material type: "<<TYPE<<" for ID="<<ID<<" is invalid."<<std::endl;
    std::cerr<<"Material type must be any one of medium, mirror, diffuser, absober or detector."<<std::endl;
    exit(1);
  }
  if(type==0&&index<=0){
    std::cerr<<"Error: Refractive index: "<<INDEX<<" for ID="<<ID<<" is invalid."<<std::endl;
    std::cerr<<"Refractive index of medium must be > 0."<<std::endl;
    exit(1);
  }
  if(type==0&&attlen<0){
    std::cerr<<"Error: Attenuation length: "<<ATTLEN<<" for ID="<<ID<<" is invalid."<<std::endl;
    std::cerr<<"Attenuation lenght of medium must be > 0 (or = 0) to activate (inactivate) the absorption in the medium."<<std::endl;
    exit(1);
  }

  LoadCAD(MATFILE);
  
  std::cout<<"Material ID="<<ID<<", File="<<MATFILE<<", Type="<<TYPE;
  if(type==0){
    if(ATTLEN>0)std::cout<<" (Index="<<INDEX<<", Att.length="<<ATTLEN<<"mm)";
    else std::cout<<" (Index="<<INDEX<<", no attenuation)";
  }
  std::cout<<" was added."<<std::endl;
}
Material::~Material()
{
}
