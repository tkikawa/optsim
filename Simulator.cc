#include "Simulator.hh"

Simulator::Simulator(std::mt19937 MT, Config config, std::string OUTPUT, Source *SRC, std::vector<Material*> &MAT)
  : mt(MT),
    output(OUTPUT),
    src(SRC),
    mat(MAT),
    nevt(10000),
    index0(1),
    displaymode(false)
{
  std::istringstream(config["Global"]["Number"]) >> nevt;
  std::istringstream(config["Global"]["Index"]) >> index0;
  chg = new Charged(mt,mat);
}
Simulator::~Simulator()
{
}
void Simulator::Run(){//Run simulation
  file = new TFile(output.c_str(), "recreate");
  tree = new TTree("tree","tree");
  tree->Branch("ipos[3]",&ipos,"ipos[3]/D");
  tree->Branch("fpos[3]",&fpos,"fpos[3]/D");
  tree->Branch("ivec[3]",&ivec,"ivec[3]/D");
  tree->Branch("fvec[3]",&fvec,"fvec[3]/D");
  tree->Branch("time",&time,"time/D");
  tree->Branch("length",&length,"length/D");
  tree->Branch("imat",&imat,"imat/I");
  tree->Branch("fmat",&fmat,"fmat/I");
  tree->Branch("ftype",&ftype,"ftype/I");
  tree->Branch("nref",&nref,"nref/I");
  tree->Branch("npas",&npas,"npas/I");
  for(int j=0;j<5;j++)count[j]=0;
  if(src->ChargedMode()){
    std::cout<<"Generate "<<nevt<<" charged particles."<<std::endl;
  }
  else{
    std::cout<<"Generate "<<nevt<<" optical photons."<<std::endl;
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
	Track(pos,cross,1);
      }
      chg->Simulate(pos, vec);//Simualtion for charged particles.
      nph=chg->GetNPhotons();
    }
    else nph=1;
    for(int p=0;p<nph;p++){
      if(src->ChargedMode()){
	chg->Generate(pos, vec, p);//Determine the initial position and direction
	mn=PointMaterial(pos);//Check the material of initial position
      }
      else{
	while(1){
	  src->Generate(pos, vec);//Determine the initial position and direction
	  mn=PointMaterial(pos);//Check the material of initial position
	  if(mn==-1)break;//Check if initial position of out of material volumes (i.e. air).
	  else if(mat[mn]->Type()<=1)break;//Check if the material type is medium or converter
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
	      }
	      else{//Rayleigh scattering in the medium
		pl=spl;
		btype=4;
		Rayleigh(vec,newvec);
	      }
	      for(int j=0;j<3;j++){
		newpos[j]=pos[j]+vec[j]*pl;
	      }
	    }
	  }
	  length+=pl;
	  time+=pl/c_0*index;
	  if(displaymode){
	    Track(pos,newpos);
	  }
	}
	else{
	  if(displaymode){
	    for(int j=0;j<3;j++){
	      cand[j]=pos[j]+vec[j]*world*2;
	    }
	    Track(pos,cand);
	  }
	}
	if(btype >= -1 && btype <4){
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
  Summary();  
  file->Write();
  file->Close();  
}
void Simulator::Display(){//Event display mode
  std::map<std::pair<Position, Position>, std::vector<Direction>> edge_to_norm;
  displaymode=true;
  TCanvas *c1 = new TCanvas("Event display","Event display",hsize,vsize);
  TView *view = TView::CreateView(1);
  x_min = y_min = z_min = world;
  x_max = y_max = z_max = -world;
  for(unsigned int m=0;m<mat.size();m++){//loop for Materials
    edge_to_norm.clear();
    for(int t=0;t<mat[m]->NTriangle();t++){//loop for Triangles
      for(int j=0;j<3;j++){
	Position p1 = Round({ mat[m]->GetTriangle(t).X(j),       mat[m]->GetTriangle(t).Y(j),       mat[m]->GetTriangle(t).Z(j) });
	Position p2 = Round({ mat[m]->GetTriangle(t).X((j+1)%3), mat[m]->GetTriangle(t).Y((j+1)%3), mat[m]->GetTriangle(t).Z((j+1)%3) });
	if(ComparePosition(p1,p2))std::swap(p1, p2);
	edge_to_norm[{p1, p2}].push_back(mat[m]->GetTriangle(t).GetNormal());
	Compare(x_max, x_min, p1[0]);
	Compare(y_max, y_min, p1[1]);
	Compare(z_max, z_min, p1[2]);
      }//end of for      
    }//end of for
    for (auto it = edge_to_norm.begin(); it != edge_to_norm.end(); ++it) {
      if (it->second.size() != 2) continue;
      if (!Parallel(it->second[0], it->second[1])){
	Draw(it->first.first,it->first.second,mat[m]->Type());
      }
    }
  }//end of for
  r_max = (x_max-x_min > y_max-y_min) ? x_max-x_min : y_max-y_min;
  r_max = (z_max-z_min > r_max) ? z_max-z_min : r_max;
  view->SetRange((x_max-x_min-r_max)/2, (y_max-y_min-r_max)/2, (z_max-z_min-r_max)/2,
		 (x_max-x_min+r_max)/2, (y_max-y_min+r_max)/2, (z_max-z_min+r_max)/2);

  std::cout<<"******** Display mode ********"<<std::endl;
  std::cout<<"Red:    normal medium material"<<std::endl;
  std::cout<<"Magenta:converter material"<<std::endl;
  std::cout<<"Yellow: mirror material"<<std::endl;
  std::cout<<"Green:  diffuser material"<<std::endl;
  std::cout<<"Cyan:   absorber material"<<std::endl;
  std::cout<<"Blue:   detector material"<<std::endl;
  std::string msg;
  while(1){
    std::cout<<"Waiting for input."<<std::endl;
    std::cout<<"\"quit\" or \"exit\": Terminate the display mode."<<std::endl;
    std::cout<<"\"gui\"           : Switch to GUI control mode."<<std::endl;
    std::cout<<"\"save\"          : Save the event display."<<std::endl;
    std::cout<<"Positive integer: Generate photons and display tracks."<<std::endl;
    c1->Update();
    std::cin >> msg;
    if(msg=="quit" || msg=="exit")break;
    else if(msg=="gui"){
      std::cout<<"Switched to GUI control mode."<<std::endl;
      std::cout<<"Event display can be rotated or cotrolled with mouse."<<std::endl;
      std::cout<<"Double click the event display to switch back to CUI control mode."<<std::endl;
      c1->WaitPrimitive();
      std::cout<<"Switched to CUI control mode."<<std::endl;
    }
    else if(msg=="save"){
      std::cout<<"File name :"<<std::flush;
      std::cin >> msg;
      c1->Print(msg.c_str());
    }
    else if(std::all_of(msg.cbegin(), msg.cend(), isdigit)){
      if(stoi(msg)>0){
	nevt=stoi(msg);
	Run();
      }
    }
  }
  c1->Close();
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
  time=0;
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
  if(unirand(mt)<2./3.){//Angle-independent term
    Isotropic(newv);
  }
  else{//cos^2(theta)-dependent term
  double cost,sint,cosp,sinp;
  double ranv[3],vr,ang;  
  cost=v[2];
  sint=sqrt(1-cost*cost);
  if(sint!=0){
    cosp=v[0]/sint;
  }
  else{
    cosp=1;
  }
  if(v[1]>0){
    sinp=sqrt(1-cosp*cosp);
  }
  else{
    sinp=-sqrt(1-cosp*cosp);
  }

  if(unirand(mt)<0.5){
    ranv[2]=sqrt(sqrt(unirand(mt)));
  }
  else{
    ranv[2]=-sqrt(sqrt(unirand(mt)));
  }
  ang=unirand(mt)*2*pi;
  vr=sqrt(1-ranv[2]*ranv[2]);
  ranv[0]=vr*cos(ang);
  ranv[1]=vr*sin(ang);
  newv[0]=-sinp*ranv[0] +cost*cosp*ranv[1] +sint*cosp*ranv[2];
  newv[1]=cosp*ranv[0]  +cost*sinp*ranv[1] +sint*sinp*ranv[2];
  newv[2]=              -sint*ranv[1]      +cost*ranv[2];    
  }
}
void Simulator::Draw(const Position& p0, const Position& p1, int type){//Draw materials in the event display
  TPolyLine3D *l = new TPolyLine3D(2);
  l->SetPoint(0,p0[0],p0[1],p0[2]);
  l->SetPoint(1,p1[0],p1[1],p1[2]);
  l->SetPoint(0,p0[0],p0[1],p0[2]);
  if(type==0)     l->SetLineColor(2);
  else if(type==1)l->SetLineColor(6);
  else if(type==2)l->SetLineColor(5);
  else if(type==3)l->SetLineColor(8);
  else if(type==4)l->SetLineColor(7);
  else if(type==5)l->SetLineColor(4);
  l->Draw();
}
void Simulator::Track(const Position& p0, const Position& p1, bool charged){//Draw track of optical photon in the event display
  TPolyLine3D *l = new TPolyLine3D(2);
  l->SetPoint(0,p0[0],p0[1],p0[2]);
  l->SetPoint(1,p1[0],p1[1],p1[2]);
  if(charged){
    l->SetLineColor(2);
    l->SetLineWidth(3);
  }
  else{
    l->SetLineColor(1);
    l->SetLineWidth(2);
  }
  l->Draw();
}
void Simulator::Compare(double &A_max, double &A_min, double A){//Update the minimum and maximum points of a coordinate
  if(A_min > A) A_min = A;
  if(A_max < A) A_max = A;
}
bool Simulator::ComparePosition(const Position& p0, const Position& p1){
    for (int i = 0; i < 3; ++i) {
        if (p0[i] != p1[i]) return p0[i] < p1[i];
    }
    return false;
}
bool Simulator::Parallel(const Direction& n0, const Direction& n1){
  return fabs(n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2]) > 0.999;
}
Position Simulator::Round(const Position& p0){
  double eps = 1e-5;
  return {
	  std::round(p0[0] / eps) * eps,
	  std::round(p0[1] / eps) * eps,
	  std::round(p0[2] / eps) * eps
  };
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
