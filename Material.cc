#include "Material.hh"

Material::Material(std::mt19937 MT, int ID, std::string MATFILE, std::string TYPE, double INDEX=1, double ATTLEN=0, double SCATLEN=0)
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
  scatlen = SCATLEN;
  if(TYPE=="medium")         type = 0;
  else if(TYPE=="converter") type = 1;
  else if(TYPE=="mirror")    type = 2;
  else if(TYPE=="diffuser")  type = 3;
  else if(TYPE=="absorber")  type = 4;
  else if(TYPE=="detector")  type = 5;
  else{
    std::cerr<<"Error: Material type: "<<TYPE<<" for ID="<<ID<<" is invalid."<<std::endl;
    std::cerr<<"Material type must be any one of medium, mirror, diffuser, absober or detector."<<std::endl;
    exit(1);
  }
  if(type==0 || type==1){
    if(index<=0){
      std::cerr<<"Error: Refractive index: "<<INDEX<<" for ID="<<ID<<" is invalid."<<std::endl;
      std::cerr<<"Refractive index of medium must be > 0."<<std::endl;
      exit(1);
    }
    if(attlen<0){
      std::cerr<<"Error: Attenuation length: "<<ATTLEN<<" for ID="<<ID<<" is invalid."<<std::endl;
      std::cerr<<"Attenuation lenght of medium must be > 0 (or = 0) to activate (inactivate) the absorption in the medium."<<std::endl;
      exit(1);
    }
    if(scatlen<0){
      std::cerr<<"Error: Scattering length: "<<SCATLEN<<" for ID="<<ID<<" is invalid."<<std::endl;
      std::cerr<<"Scattering lenght of medium must be > 0 (or = 0) to activate (inactivate) the scattering in the medium."<<std::endl;
    exit(1);
    }
  }
  LoadCAD(MATFILE);
  
  std::cout<<"Material ID="<<ID<<", File="<<MATFILE<<", Type="<<TYPE;
  if(type==0 || type==1){
    std::cout<<" (Index="<<INDEX<<", ";
    if(ATTLEN>0)std::cout<<"Att.length="<<ATTLEN<<"mm, ";
    else std::cout<<"no attenuation, ";
    if(SCATLEN>0)std::cout<<"Scat.length="<<SCATLEN<<"mm)";
    else std::cout<<"no scattering)";
  }
  std::cout<<" was added."<<std::endl;
}
Material::~Material()
{
}
