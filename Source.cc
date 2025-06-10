#include "Source.hh"

Source::Source(std::mt19937 MT, Config config)
  : Geometry(MT), mass(0), energy(0), beta(1), charge(0)
{
  if(config["Source"].size()==0){
    std::cerr<<"Error: Souce is not defined in the input card file."<<std::endl;
    exit(1);
  }  
  sourcemode = config["Source"].begin()->first; // only first source in input card file is read in
  std::istringstream sourceconf(config["Source"].begin()->second);
  std::cout<<"Source type is "<<sourcemode<<"."<<std::endl;
  if (sourcemode == "boxvolume"){
    sourceconf >> x_min >> x_max >> y_min >> y_max >> z_min >> z_max;
  }
  else if (sourcemode == "boxsurface"){
    sourceconf >> x_min >> x_max >> y_min >> y_max >> z_min >> z_max;
    xsurf=(y_max - y_min)*(z_max - z_min);
    ysurf=(x_max - x_min)*(z_max - z_min);
    zsurf=(x_max - x_min)*(y_max - y_min);
    totsurf=xsurf+ysurf+zsurf;
  }
  else if (sourcemode == "cylvolume"){
    sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max;
    phi_max*=conv;
    phi_min*=conv;
  }
  else if (sourcemode == "cylsurface"){
    sourceconf >> r_min >> r_max >> phi_min >> phi_max >> z_min >> z_max;
    phi_max*=conv;
    phi_min*=conv;
    topsurf=(r_max*r_max - r_min*r_min)*(phi_max - phi_min)/2;
    insurf=r_min*(phi_max - phi_min)*(z_max - z_min);
    outsurf=r_max*(phi_max - phi_min)*(z_max - z_min);
    totsurf=topsurf*2+insurf+outsurf;
  }
  else if (sourcemode == "sphvolume"){
    sourceconf >> x_center >> y_center >> z_center >> r;
  }
  else if (sourcemode == "sphsurface"){
    sourceconf >> x_center >> y_center >> z_center >> r;
  }
  else if (sourcemode == "CADvolume"){
    sourceconf >> sourcefile;
    LoadCAD(sourcefile);
    x_min = y_min = z_min = world;
    x_max = y_max = z_max = -world;
    for(unsigned int i=0; i < triangle.size(); i++){
      for(unsigned int j=0; j < 3; j++){
	Compare(x_max, x_min, triangle[i].X(j));
	Compare(y_max, y_min, triangle[i].Y(j));
	Compare(z_max, z_min, triangle[i].Z(j));
      }
    }
  }
  else if (sourcemode == "CADsurface"){
    sourceconf >> sourcefile;
    LoadCAD(sourcefile);
    totsurf = 0;
    for(unsigned int i=0; i < triangle.size(); i++){
      totsurf += triangle[i].GetArea();
    }
  }
  else if(sourcemode != "custom"){
    std::cerr<<"Error: Invalid source type."<<std::endl;
    std::cerr<<"Source type must be any one of CADvolume, CADsurface, boxvolume, boxsurface, cylvolume, cylsurface, sphvolume, sphsurface or custom."<<std::endl;
    exit(1);
  }

  if(config["Particle"].size()==0){
    std::cerr<<"Error: Particle is not defined in the input card file."<<std::endl;
    exit(1);
  }  
  particlemode = config["Particle"].begin()->first; // only first particle in input card file is read in
  std::istringstream particleconf(config["Particle"].begin()->second);
  std::cout<<"Particle type is "<<particlemode<<"."<<std::endl;  
  if (particlemode == "photon"){
    charged = false;
  }
  else if (particlemode == "electron"){
    charged = true;
    mass=0.511;
    charge=1;
    particleconf >> energy;
  }
  else if (particlemode == "muon"){
    charged = true;
    mass=105.658;
    charge=1;
    particleconf >> energy;
  }
  else if (particlemode == "pion"){
    charged = true;
    mass=139.570;
    charge=1;
    particleconf >> energy;
  }
  else if (particlemode == "proton"){
    charged = true;
    mass=938.272;
    charge=1;
    particleconf >> energy;
  }
  else if (particlemode == "alpha"){
    charged = true;
    mass=3727.38;
    charge=2;
    particleconf >> energy;
  }  
  else{
    std::cerr<<"Error: Invalid particle type."<<std::endl;
    std::cerr<<"Particle type must be either photon, electron, muon, pion, proton or alpha."<<std::endl;
    exit(1);
  }
  if(charged){
    std::cout<<"Particle kinetic energy is "<<energy<<" MeV."<<std::endl;  
    beta = CalcBeta(mass,energy);
  }
  if(config["Direction"].size()==0){
    std::cerr<<"Error: Direction is not defined in the input card file."<<std::endl;
    exit(1);
  }
  directionmode = config["Direction"].begin()->first; // only first direction in input card file is read in
  std::istringstream directionconf(config["Direction"].begin()->second);
  std::cout<<"Direction type is "<<directionmode<<"."<<std::endl;
  if (directionmode == "isotropic"){
    v_x=0; v_y=0; v_z=1; phi=180;
    phi*=conv;
  }
  else if (directionmode == "lambert"){
    directionconf >> v_x >> v_y >> v_z;
  }
  else if (directionmode == "cosmic"){
    directionconf >> v_x >> v_y >> v_z;
  }  
  else if (directionmode == "flat"){
    directionconf >> v_x >> v_y >> v_z >> phi;
    phi*=conv;
  }
  else if (directionmode == "gauss"){
    directionconf >> v_x >> v_y >> v_z >> phi;
    phi*=conv;
  }
  else if(directionmode != "custom"){
    std::cerr<<"Error: Invalid direction type."<<std::endl;
    std::cerr<<"Direction type must be any one of isotropic, lambert, cosmic, flat, gauss, custom"<<std::endl;
    exit(1);
  }
  
  Normalize(v_x, v_y, v_z);  

  cost=v_z;
  sint=sqrt(1-cost*cost);
  if(sint!=0){
    cosp=v_x/sint;
  }
  else{
    cosp=1;
  }
  if(v_y >= 0){
    sinp=sqrt(1-cosp*cosp);
  }
  else{
    sinp=-sqrt(1-cosp*cosp);
  }

  rz=cos(phi);
}
Source::~Source()
{
}
void Source::Generate(Position& pos, Direction& vec)
{
  pos=PointInSource();
  vec=ParticleDir();
}
Position Source::PointInSource(){//Randomely determine a point in the source volume or surface.
  
  double phi_r;
  Position pos;

  if (sourcemode == "boxvolume"){
    pos[0] = x_min + unirand(mt)*(x_max - x_min);
    pos[1] = y_min + unirand(mt)*(y_max - y_min);
    pos[2] = z_min + unirand(mt)*(z_max - z_min);
  }
  else if (sourcemode == "boxsurface"){
    double randsurf=unirand(mt)*totsurf;
    if(randsurf<xsurf){
      pos[1] = y_min + unirand(mt)*(y_max - y_min);
      pos[2] = z_min + unirand(mt)*(z_max - z_min);
      if(unirand(mt)<0.5) pos[0]=x_min;
      else                pos[0]=x_max;
    }
    else if(randsurf<xsurf+ysurf){
      pos[0] = x_min + unirand(mt)*(x_max - x_min);
      pos[2] = z_min + unirand(mt)*(z_max - z_min);
      if(unirand(mt)<0.5) pos[1]=y_min;
      else                pos[1]=y_max;
    }
    else{
      pos[0] = x_min + unirand(mt)*(x_max - x_min);
      pos[1] = y_min + unirand(mt)*(y_max - y_min);
      if(unirand(mt)<0.5) pos[2]=z_min;
      else                pos[2]=z_max;
    }
  }
  else if (sourcemode == "cylvolume"){
    r = sqrt(unirand(mt)*(r_max*r_max - r_min*r_min) + r_min*r_min); // weighting because of the volume element and a r^2 probability outwards
    phi_r = phi_min + unirand(mt)*(phi_max - phi_min);
    pos[0] = r*cos(phi_r);
    pos[1] = r*sin(phi_r);
    pos[2] = z_min + unirand(mt)*(z_max - z_min);
  }
  else if (sourcemode == "cylsurface"){
    double randsurf=unirand(mt)*totsurf;

    if(randsurf<outsurf){
      phi_r = phi_min + unirand(mt)*(phi_max - phi_min);
      pos[0] = r_max*cos(phi_r);
      pos[1] = r_max*sin(phi_r);
      pos[2] = z_min + unirand(mt)*(z_max - z_min);
    }
    else if(randsurf<outsurf+insurf){
      phi_r = phi_min + unirand(mt)*(phi_max - phi_min);
      pos[0] = r_min*cos(phi_r);
      pos[1] = r_min*sin(phi_r);
      pos[2] = z_min + unirand(mt)*(z_max - z_min);
    }
    else if(randsurf<outsurf+insurf+topsurf){
      r = sqrt(unirand(mt)*(r_max*r_max - r_min*r_min) + r_min*r_min);
      phi_r = phi_min + unirand(mt)*(phi_max - phi_min);
      pos[0] = r*cos(phi_r);
      pos[1] = r*sin(phi_r);
      pos[2] = z_max;
    }
    else{
      r = sqrt(unirand(mt)*(r_max*r_max - r_min*r_min) + r_min*r_min);
      phi_r = phi_min + unirand(mt)*(phi_max - phi_min);
      pos[0] = r*cos(phi_r);
      pos[1] = r*sin(phi_r);
      pos[2] = z_min;
    }
  }
  else if (sourcemode == "sphvolume"){
    double theta = acos(1.0 - 2.0 * unirand(mt));
    double phi = 2.0 * pi * unirand(mt);
    double rr = r * cbrt(unirand(mt));
    pos[0] = x_center + rr * sin(theta) * cos(phi);
    pos[1] = y_center + rr * sin(theta) * sin(phi);
    pos[2] = z_center + rr * cos(theta);    
  }  
  else if (sourcemode == "sphsurface"){
    double theta = acos(1.0 - 2.0 * unirand(mt));
    double phi = 2.0 * pi * unirand(mt);
    pos[0] = x_center + r * sin(theta) * cos(phi);
    pos[1] = y_center + r * sin(theta) * sin(phi);
    pos[2] = z_center + r * cos(theta);
  }  
  else if (sourcemode == "CADsurface"){
    double randsurf = unirand(mt)*totsurf;
    double tmpsurf = 0;
    for(unsigned int i=0; i < triangle.size(); i++){
      tmpsurf += triangle[i].GetArea();
      if(randsurf < tmpsurf){
	double r1=unirand(mt);
	double r2=unirand(mt);
	pos = triangle[i].GetSurfPoint(r1,r2);
	break;
      }
    }
  }
  else if (sourcemode == "CADvolume"){
    while(1){
      pos[0] = x_min+unirand(mt)*(x_max-x_min);
      pos[1] = y_min+unirand(mt)*(y_max-y_min);
      pos[2] = z_min+unirand(mt)*(z_max-z_min);
      if(InSolid(pos))break;
    }
  }
  else{// sourcemode == "custom"
    pos=CustomSource(mt);
  }
  return pos;
}
Direction Source::ParticleDir(){//Randomely determine the initial direction of the primary particle
  Direction vec,v;
  if (directionmode == "isotropic" || directionmode == "flat"){
    v[2]=1+unirand(mt)*(rz-1);
    vr=sqrt(1-v[2]*v[2]);
  }
  else if(directionmode == "gauss"){
    v[2]=-2;
    while(v[2]<-1){
      v[2]=1+fabs(gausrand(mt))*(rz-1);
    }
    vr=sqrt(1-v[2]*v[2]);
  }
  else if(directionmode == "lambert"){
    vr=sqrt(unirand(mt));
    v[2]=sqrt(1 - vr*vr);
  }
  else if(directionmode == "cosmic"){
    v[2]=pow(unirand(mt), 1.0 / (n_cosmic + 1));
    vr=sqrt(1-v[2]*v[2]);
  }
  else{//directionmode == "custom"
    vec=CustomDirection(mt);
    Normalize(vec[0], vec[1], vec[2]);
    return vec;
  }
  ang=unirand(mt)*2*pi;
  v[0]=vr*cos(ang);
  v[1]=vr*sin(ang);
  vec[0]=-sinp*v[0] +cost*cosp*v[1] +sint*cosp*v[2];
  vec[1]=cosp*v[0]  +cost*sinp*v[1] +sint*sinp*v[2];
  vec[2]=           -sint*v[1]      +cost*v[2];  
  return vec;
}
bool Source::ChargedMode(){
  return charged;
}
double Source::GetBeta(){
  return beta;
}
void Source::Compare(double &A_max, double &A_min, double A){//Update the minimum and maximum points of a coordinate.
  if(A_min > A) A_min = A;
  if(A_max < A) A_max = A;
}
void Source::Normalize(double &vx, double &vy, double &vz){//Normalize the vector norm.
  double norm=vx*vx+vy*vy+vz*vz;
  vx/=norm;
  vy/=norm;
  vz/=norm;
}
double Source::CalcBeta(double m, double T){
  double E = T + m;
  double gamma = E / m;
  double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
  return beta;
}
