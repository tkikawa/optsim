#include "Simulator.hh"

Simulator::Simulator(std::mt19937 MT, Config config, std::string OUTPUT, Source *SRC, std::vector<Material*> &MAT)
  : mt(MT),
    output(OUTPUT),
    src(SRC),
    mat(MAT),
    nevt(10000),
    index0(1),
    yield(0),
    delay(0),
    wlmin(-99999),
    wlmax(-99999),
    gmie(-99999),
    sci_type(""),
    displaymode(false),
    act_scinti(false),
    act_cherenkov(false),
    usemie(false)
{
  std::istringstream(config["Global"]["Number"]) >> nevt;
  if(nevt<=0){
    std::cout<<"Number of primary particles is "<<nevt<<std::endl;
    std::cout<<"It must be more than 0."<<std::endl;
  }
  std::istringstream(config["Global"]["Index"]) >> index0;
  if(index0<1){
    std::cout<<"Refractive index of surroundings is"<<index0<<std::endl;
    std::cout<<"It must be 1 or more."<<std::endl;
  }
  std::istringstream(config["Global"]["Scintillation"]) >> sci_type >> yield >> delay;
  if(!sci_type.empty()){
    check_defined(sci_type);
    if(yield==0){
      std::cout<<"Scintillation is activated with "<<sci_type <<" scintillator with default light yield."<<std::endl;
    }
    else{
      std::cout<<"Scintillation is activated with "<<sci_type <<" scintillator with "<<yield<<" photons/MeV."<<std::endl;
    }
    act_scinti=true;
  }
  std::istringstream(config["Global"]["Cherenkov"]) >> wlmin >> wlmax;
  if(wlmin>100&&wlmin<wlmax&&wlmax<700){
    std::cout<<"Cherenkov radiation is activated for wavelendth range of "<<wlmin<<" nm to "<<wlmax<<" nm."<<std::endl;
    act_cherenkov=true;
  }
  else if(wlmin!=-99999||wlmax!=-99999){
    std::cout<<"Setting of wavelength range for Cherenkov radiation is strange."<<std::endl;
    exit(1);
  }
  std::istringstream(config["Global"]["Mie"]) >> gmie;
  if(gmie>-1&&gmie<1){
    std::cout<<"Mie scattering is activated with asymmetry parameter = "<<gmie<<std::endl;
    usemie=true;
  }
  else if(gmie!=-99999){
    std::cout<<"Asymmetry parameter in Mie scattering is "<<gmie<<std::endl;
    std::cout<<"It must be between -1 and 1."<<gmie<<std::endl;
    exit(1);
  }
  chg = new Charged(mt,mat);
  if(act_scinti||act_cherenkov)chg->SetBeta(src->GetBeta());
  if(act_scinti)chg->SetScinti(sci_type,yield,delay);
  if(act_cherenkov)chg->SetCherenkov(wlmin,wlmax);
}
Simulator::~Simulator()
{
}
void Simulator::Run(){//Run simulation
  file = new TFile(output.c_str(), "recreate");//Output ROOT file
  tree = new TTree("photon","photon");
  tree->Branch("ipos[3]",&ipos,"ipos[3]/D");//Initial position of optical photon (mm). [0], [1], [2] for x, y, z.
  tree->Branch("fpos[3]",&fpos,"fpos[3]/D");//Final position of optical photon (mm). [0], [1], [2] for x, y, z.
  tree->Branch("ivec[3]",&ivec,"ivec[3]/D");//Initial direction of optical photon (mm). [0], [1], [2] for x, y, z.
  tree->Branch("fvec[3]",&fvec,"fvec[3]/D");//Final direction of optical photon (mm). [0], [1], [2] for x, y, z.
  tree->Branch("time",&time,"time/D");//Time from generation to end of optical photon (ns).
  tree->Branch("length",&length,"length/D");//Total path length traveled by optical photon (mm).
  tree->Branch("imat",&imat,"imat/I");//Material ID in which optical photon was generated.
  tree->Branch("fmat",&fmat,"fmat/I");//Material ID in which optical photon ended.
  tree->Branch("ftype",&ftype,"ftype/I");//ID which identifies how optical photon ended. (See README)
  tree->Branch("nref",&nref,"nref/I");//Number of reflections of optical photon in boundaries.
  tree->Branch("npas",&npas,"npas/I");//Number of transmissions of optical photon in boundaries.
  if(src->ChargedMode()){
    tree->Branch("id",&id,"id/I");//ID of corresponding charged particle.

    charged = new TTree("charged","charged");
    charged->Branch("ipos[3]",&pos,"ipos[3]/D");//Initial position of charged particle (mm). [0], [1], [2] for x, y, z.
    charged->Branch("ivec[3]",&vec,"ivec[3]/D");//Initial direction of charged particle (mm). [0], [1], [2] for x, y, z.
    charged->Branch("nph",&nph,"nph/I");//Number of produced photons
    charged->Branch("id",&id,"id/I");//ID for charged particles.
  }
  
  for(int j=0;j<5;j++)count[j]=0;
  if(!displaymode){
    if(src->ChargedMode()){
      std::cout<<"Generate "<<nevt<<" charged particles."<<std::endl;
    }
    else{
      std::cout<<"Generate "<<nevt<<" optical photons."<<std::endl;
    }
  }
  bool progbar=false;
  if(nevt>=100&&!displaymode)progbar=true;
  if(progbar){
    std::cout<<"0%       25%       50%       75%       100%"<<std::endl;
    std::cout<<"+---------+---------+---------+---------+"<<std::endl;;
  }
  double ndev40=nevt/40.;
  int nprog=0;
  for(int i=0;i<nevt;i++){
    if(progbar){
      if(i>=nprog*ndev40){
	std::cout<<"#"<<std::flush;
	nprog++;
      }
      else if(i==nevt-1){
	std::cout<<"#"<<std::endl;
      }
    }
    if(src->ChargedMode()){
      src->Generate(pos, vec);//Determine the initial position and direction of charged particle
      if(displaymode){
	for(int j=0;j<3;j++)cross[j]=pos[j]+vec[j]*world;
	gui->DrawTrack(pos,cross,1);
      }
      chg->Simulate(pos, vec);//Simualtion for charged particles.
      nph=chg->GetNPhotons();
      id=i;
      charged->Fill();
    }
    else nph=1;
    for(int p=0;p<nph;p++){
      if(src->ChargedMode()){
	time=chg->Generate(pos, vec, p);//Determine the initial position, direction and time
	mn=PointMaterial(pos);//Check the material of initial position
      }
      else{
	time=0;
	while(1){
	  src->Generate(pos, vec);//Determine the initial position and direction
	  mn=PointMaterial(pos);//Check the material of initial position
	  if(mn==-1)break;//Check if initial position of out of material volumes (i.e. air).
	  else if(mat[mn]->TType()<2||mat[mn]->TType()==6)break;//Check if the material type is medium, converter or mixture including medium or converter
	}
      }
      if(mn!=-1){//Initial position is in the defined material volumes.
	matid=mat[mn]->ID();
	index=mat[mn]->Index();
	attlen=mat[mn]->AttLen();
	scatlen=mat[mn]->ScatLen();
      }
      else{//Initial position is out of the defined material volumes and is in surrounding (air).
	matid=0;
	index=index0;
	attlen=0;
	scatlen=0;
      }
      Initialize();
      while(1){
	btype = -2;
	for(int j=0;j<3;j++)cross[j]=pos[j]+vec[j]*world;//far point on the straight-line of the track
	for(unsigned int m=0;m<mat.size();m++){//loop for Materials
	  if(matid!=0&&matid<mat[m]->ID()&&btype>-1){
	    goto ENDLOOP;
	  }
	  else if(matid!=0&&matid<mat[m]->ID()&&btype==-1){
	    if(mat[m]->InSolid(cross)){
	      newmatid=mat[m]->ID();
	      newindex=mat[m]->Index();
	      newattlen=mat[m]->AttLen();
	      newscatlen=mat[m]->ScatLen();
	      btype=mat[m]->Type();
	      newmn=m;
	      goto ENDLOOP;
	    }
	    continue;
	  }
	  for(int t=0;t<mat[m]->NTriangle();t++){//loop for Triangles
	    if(mat[m]->GetTriangle(t).Collision(pos,cross,cand)){//Check if a line between two points collides with the triangle or not
	      for(int j=0;j<3;j++){
		newpos[j]=cand[j];
		if(matid==0){
		  cross[j]=cand[j];
		}
		else if(matid==mat[m]->ID()){
		  cross[j]=cand[j]+vec[j]*micro;
		}
		else{
		  cross[j]=cand[j]-vec[j]*micro;
		}
	      }
	      normal=mat[m]->GetTriangle(t).GetNormal();
	      if(matid==mat[m]->ID()){
		newmatid=0;
		newindex=index0;
		newattlen=0;
		newscatlen=0;
		btype=-1;
		newmn=-1;
	      }
	      else{
		newmatid=mat[m]->ID();
		newindex=mat[m]->Index();
		newattlen=mat[m]->AttLen();
		newscatlen=mat[m]->ScatLen();
		btype=mat[m]->Type();
		newmn=m;
	      }
	    }//end of if
	  }//end of for
	}//end of for
 
      ENDLOOP:
	if(btype!=-2){
	  pl=sqrt(pow(newpos[0]-pos[0],2)+pow(newpos[1]-pos[1],2)+pow(newpos[2]-pos[2],2));
	  if(attlen>0||scatlen>0){//When attlen==0 (scatlen==0), absorption (scattering) in the medium does not occur.
	    if(attlen>0)apl=-attlen*log(unirand(mt));
	    else apl=DBL_MAX;
	    if(scatlen>0)spl=-scatlen*log(unirand(mt));
	    else spl=DBL_MAX;
	    if(apl<pl || spl<pl){
	      if(apl<spl){//Absorption in the medium
		pl=apl;
		btype=-3;
		newmatid=matid;
	      }
	      else{//Rayleigh or Mie scattering in the medium
		pl=spl;
		btype=10;
		if(!usemie) Rayleigh(vec,newvec);
		else        Mie(vec,newvec);
	      }
	      for(int j=0;j<3;j++){
		newpos[j]=pos[j]+vec[j]*pl;
	      }
	    }
	  }
	  length+=pl;
	  time+=pl/c_0*index;
	  if(displaymode){
	    gui->DrawTrack(pos,newpos);
	  }
	}
	else{
	  if(displaymode){
	    for(int j=0;j<3;j++){
	      cand[j]=pos[j]+vec[j]*world*2;
	    }
	    gui->DrawTrack(pos,cand);
	  }
	}
	if((btype >= -1 && btype <4) || btype==10){
	  if(btype <= 1){
	    if(Fresnel(vec,newvec,normal,index,newindex)){//Determine if the optical photon is reflected or transmitted with refraction following the Fresnel equation
	      mn=newmn;
	      matid=newmatid;
	      index=newindex;
	      attlen=newattlen;
	      scatlen=newscatlen;
	      npas++;
	      if(btype == 1){//Isotropic scattering for converter
		Isotropic(newvec);
	      }
	    }
	    else{
	      nref++;
	    }
	  }
	  else if(btype == 2){
	    Specular(vec,newvec,normal);//Specular reflection
	    nref++;
	  }
	  else if(btype == 3){
	    Lambert(vec,newvec,normal);//Diffusion following the Lambert's cosine law
	    nref++;
	  }
	  Normalize(newvec);
	  for(int j=0;j<3;j++){
	    pos[j]=newpos[j]+newvec[j]*nano;
	    vec[j]=newvec[j];
	  }
	}
	else break;//Absorbed, detected or go out of world volume.
	if(nref>nreflimit){
	  btype=-4;
	  break;
	}
      }//end of while
      for(int j=0;j<3;j++){
	fpos[j]=newpos[j];
	fvec[j]=newvec[j];
      }
      fmat=newmatid;
      ftype=FType(btype);//Conversion of type ID.
      tree->Fill();
      count[ftype]++;
    }//end of for
  }//end of for
  if(!displaymode){
    Summary();
  }
  file->Write();
  file->Close();  
}
void Simulator::SetDisplay(Gui *GUI){
  gui = GUI;
  gui->SetSimulator(this);
  displaymode=true;
}
int Simulator::PointMaterial(const Position& p){//Check the material of initial position
  for(unsigned int m=0;m<mat.size();m++){
    if(mat[m]->InSolid(p)){
      return m;
    }
  }
  return -1;
}
void Simulator::Initialize(){//Initialize the variables.
  for(int j=0;j<3;j++){
    ipos[j]=pos[j];
    ivec[j]=vec[j];
    newpos[j]=pos[j];
    newvec[j]=vec[j];
  }
  imat=matid;
  nref=0;
  npas=0;
  length=0;
}
int Simulator::FType(int bt){//Conver the temporary type ID to final type ID
  if(bt==-2)return 0;     //Go out of world volume.
  else if(bt==-4)return 1;//Exceed the limit of reflections
  else if(bt==-3)return 2;//Absorbed in normal medium.
  else if(bt==4)return 3; //Absorbed by absorber
  else if(bt==5)return 4; //Detected by detector
  else{
    std::cerr<<"Error in FType"<<std::endl;
    exit(1);
  }
}
void Simulator::Summary(){//Display the summary of simulation
  std::cout<<"Go out of world volume     : "<<count[0]<<std::endl;
  std::cout<<"Exceed limit of reflections: "<<count[1]<<std::endl;
  std::cout<<"Absorved in normal medium  : "<<count[2]<<std::endl;
  std::cout<<"Absorbed by absorber       : "<<count[3]<<std::endl;
  std::cout<<"Detected by detector       : "<<count[4]<<std::endl;
}
bool Simulator::Fresnel(const Direction& v, Direction& newv, Direction norm, double idx_in, double idx_out){//Determine if the optical photon is reflected or transmitted with refraction following the Fresnel equation
  double idx_ratio=idx_out/idx_in;
  double cos_in=v[0]*norm[0]+v[1]*norm[1]+v[2]*norm[2];
  double sin_in=sqrt(1-cos_in*cos_in);
  if(cos_in>0){
    for(int i=0;i<3;i++){
      norm[i]*=-1;
    }
  }
  else{
    cos_in*=-1;
  }
  if(sin_in>idx_ratio){//Total reflection
    Specular(v,newv,norm);
    return false;
  }
  else{
    double sin_out=sin_in/idx_ratio;//Snell's law
    double cos_out=sqrt(1-sin_out*sin_out);
    double ref_prob=0.5*(pow(cos_in-idx_ratio*cos_out,2)/pow(cos_in+idx_ratio*cos_out,2)+pow(idx_ratio*cos_in-cos_out,2)/pow(idx_ratio*cos_in+cos_out,2));//Fresnel equations
    if(unirand(mt)<ref_prob){//Reflection
      Specular(v,newv,norm);
      return false;
    }
    else{//Transmission with refraction
      double para[3];
      for(int i=0;i<3;i++){
	para[i]=v[i]+norm[i]*cos_in;
	if(sin_in!=0){
	  newv[i]=-norm[i]*cos_out+para[i]/sin_in*sin_out;
	}
	else{
	  newv[i]=v[i];
	}
      }
      return true;
    }
  } 
}
void Simulator::Specular(const Direction& v, Direction& newv, const Direction& norm){//Specular reflection
  double cp=v[0]*norm[0]+v[1]*norm[1]+v[2]*norm[2];
  for(int i=0;i<3;i++){
    newv[i]=v[i]-2*norm[i]*cp;
  }
}
void Simulator::Lambert(const Direction& v, Direction& newv, Direction norm){//Diffusion following the Lambert's cosine law
  double cp=v[0]*norm[0]+v[1]*norm[1]+v[2]*norm[2];
  if(cp>0){
    for(int i=0;i<3;i++){
      norm[i]*=-1;
    }
  }

  double cost,sint,cosp,sinp;
  double ranv[3],vr,ang;
  
  cost=norm[2];
  sint=sqrt(1-cost*cost);
  if(sint!=0){
    cosp=norm[0]/sint;
  }
  else{
    cosp=1;
  }
  if(norm[1]>0){
    sinp=sqrt(1-cosp*cosp);
  }
  else{
    sinp=-sqrt(1-cosp*cosp);
  }

  ranv[2]=sqrt(unirand(mt));//sqrt is give because of the cosine law.
  ang=unirand(mt)*2*pi;
  vr=sqrt(1-ranv[2]*ranv[2]);
  ranv[0]=vr*cos(ang);
  ranv[1]=vr*sin(ang);
  newv[0]=-sinp*ranv[0] +cost*cosp*ranv[1] +sint*cosp*ranv[2];
  newv[1]=cosp*ranv[0]  +cost*sinp*ranv[1] +sint*sinp*ranv[2];
  newv[2]=              -sint*ranv[1]      +cost*ranv[2];
}
void Simulator::Rayleigh(const Direction& v, Direction& newv){//Rayleigh scattering in the medium
  double costheta;

  // Sampling using cos\theta from P(cos\theta) \propto 1 + cos^2\theta using rejection sampling
  while (true) {
    double x = 2.0 * unirand(mt) - 1.0; // x in [-1, 1]
    double y = unirand(mt) * 2.0;       // y in [0, 2]
    if (y <= 1 + x * x) {
      costheta = x;
      break;
    }
  }
  double sintheta = sqrt(1.0 - costheta * costheta);
  double phi = 2 * pi * unirand(mt);
  double cosphi = cos(phi);
  double sinphi = sin(phi);

  // Define local coordinate system
  double ux = v[0], uy = v[1], uz = v[2];
  double vx[3], wx[3];
  if (fabs(uz) < 0.99) {
    vx[0] = -uy;
    vx[1] = ux;
    vx[2] = 0;
  } else {
    vx[0] = 0;
    vx[1] = -uz;
    vx[2] = uy;
  }
  // vx normalization
  double norm = sqrt(vx[0]*vx[0] + vx[1]*vx[1] + vx[2]*vx[2]);
  vx[0] /= norm; vx[1] /= norm; vx[2] /= norm;

  // wx = v × vx
  wx[0] = uy * vx[2] - uz * vx[1];
  wx[1] = uz * vx[0] - ux * vx[2];
  wx[2] = ux * vx[1] - uy * vx[0];

  // newv = sin\theta cos\phi * vx + sin\theta sin\phi * wx + cos\theta * v
  for (int i = 0; i < 3; ++i) {
    newv[i] = sintheta * cosphi * vx[i] + sintheta * sinphi * wx[i] + costheta * v[i];
  }

  Normalize(newv);
}
void Simulator::Mie(const Direction& v, Direction& newv){//Mie scattering in the medium
  // Sampling cos\theta from Henyey-Greenstein distribution
  double xi = unirand(mt);
  double costheta;

  if (fabs(gmie) < 1e-6) {
    costheta = 2.0 * xi - 1.0; // Isotropic scattering case
  } else {
    double sq = (1.0 - gmie * gmie) / (1.0 - gmie + 2.0 * gmie * xi);
    costheta = (1.0 + gmie * gmie - sq * sq) / (2.0 * gmie);
  }

  double sintheta = sqrt(std::max(0.0, 1.0 - costheta * costheta));
  double phi = 2.0 * pi * unirand(mt);
  double cosphi = cos(phi);
  double sinphi = sin(phi);

  // Define local coordinate system
  double ux = v[0], uy = v[1], uz = v[2];
  double vx[3], wx[3];
  if (fabs(uz) < 0.99) {
    vx[0] = -uy;
    vx[1] = ux;
    vx[2] = 0;
  } else {
    vx[0] = 0;
    vx[1] = -uz;
    vx[2] = uy;
  }
  // vx normalization
  double norm = sqrt(vx[0]*vx[0] + vx[1]*vx[1] + vx[2]*vx[2]);
  vx[0] /= norm; vx[1] /= norm; vx[2] /= norm;

  // wx = v × vx
  wx[0] = uy * vx[2] - uz * vx[1];
  wx[1] = uz * vx[0] - ux * vx[2];
  wx[2] = ux * vx[1] - uy * vx[0];

  // newv = sin\theta cos\phi * vx + sin\theta sin\phi * wx + cos\theta * v
  for (int i = 0; i < 3; ++i) {
    newv[i] = sintheta * cosphi * vx[i] + sintheta * sinphi * wx[i] + costheta * v[i];
  }

  Normalize(newv);
}
void Simulator::Normalize(Direction& v){//Normalize the vector norm.
  double norm=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
  for(int i=0;i<3;i++){
    v[i]/=norm;
  }
}
void Simulator::Isotropic(Direction& v){//Randomely and isotropically determine the direction
  v[2]=1-2*unirand(mt);
  double ang=unirand(mt)*2*pi;
  double vr=sqrt(1-v[2]*v[2]);
  v[0]=vr*cos(ang);
  v[1]=vr*sin(ang);
}
void Simulator::check_defined(const std::string& input){
  const std::string candidates[6] = {
     "organic", "NaI", "CsI", "BGO", "LYSO", "LaBr3"
  };
  
  for (int i = 0; i < 6; ++i) {
    if (input == candidates[i]) {
      return;
    }
  }
  std::cerr<<"Error: scintillator type must be selected from organic, NaI, CsI, BGO, LYSO or LaBr3."<<std::endl;
  exit(1);
}
